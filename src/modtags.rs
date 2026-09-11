//! Rewriting of base modification tags (`MM`/`ML`/`MN`) for trimmed reads.
//!
//! `samtools fastq -T MM,ML` writes SAM tags as tab-separated fields after the
//! read name. Those tags encode modified base positions as counts along the
//! read, so trimming or splitting a read invalidates them. This module
//! recomputes them for a subsequence.
//!
//! FASTQ records are always in the original read orientation, so the canonical
//! base is counted left to right and the `+`/`-` strand character in `MM` is
//! metadata that is carried through unchanged. (In BAM the reverse-strand case
//! is driven by the record's reverse flag, which has no FASTQ equivalent.)

use std::fmt;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Number of reads seen that carried an `MM` tag. Used to warn about
/// `--update-mods` on input that has no base modification tags at all.
static READS_WITH_MODS: AtomicUsize = AtomicUsize::new(0);

pub fn reads_with_mods() -> usize {
    READS_WITH_MODS.load(Ordering::Relaxed)
}

/// Reasons a read's base modification tags could not be rewritten. Every one of
/// these means the input tags are broken or no longer describe the sequence
/// they are attached to, so they are reported as errors rather than silently
/// passed on.
#[derive(Debug, PartialEq, Eq)]
pub enum ModTagError {
    MalformedMm(String),
    MalformedMl(String),
    /// `MM` and `ML` disagree on how many modification calls there are.
    CountMismatch {
        mm: usize,
        ml: usize,
    },
    /// `MM` skips past more canonical bases than the sequence contains, which
    /// means the tags describe a different (usually untrimmed) sequence.
    BeyondSequence {
        base: char,
    },
    MalformedMn(String),
    /// `MN` states a sequence length that the record does not have.
    StaleMn {
        mn: usize,
        len: usize,
    },
    /// `ML` holds probabilities for calls that no `MM` tag describes.
    MlWithoutMm,
}

impl fmt::Display for ModTagError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ModTagError::MalformedMm(s) => write!(f, "malformed MM tag: {s}"),
            ModTagError::MalformedMl(s) => write!(f, "malformed ML tag: {s}"),
            ModTagError::CountMismatch { mm, ml } => write!(
                f,
                "MM lists {mm} modification call(s) but ML holds {ml} probability value(s)"
            ),
            ModTagError::BeyondSequence { base } => write!(
                f,
                "MM tag refers to more '{base}' bases than the sequence contains; \
                 the tags do not describe this read"
            ),
            ModTagError::MalformedMn(s) => write!(f, "malformed MN tag: {s}"),
            ModTagError::MlWithoutMm => write!(
                f,
                "ML tag holds modification probabilities but there is no MM tag to place them"
            ),
            ModTagError::StaleMn { mn, len } => write!(
                f,
                "MN tag says the modification tags describe a {mn} base sequence, \
                 but the read is {len} bases"
            ),
        }
    }
}

/// One `MM` sublist, e.g. `C+mh?,5,12`.
struct ModList<'a> {
    /// Everything before the deltas: canonical base, strand, codes, optional flag.
    prefix: &'a str,
    /// Canonical base to count, with `U` normalised to `T`.
    base: u8,
    /// Number of modification codes, i.e. `ML` values per modified position.
    stride: usize,
    deltas: Vec<usize>,
}

fn parse_mm(mm: &str) -> Result<Vec<ModList<'_>>, ModTagError> {
    let mut lists = Vec::new();
    for chunk in mm.split(';') {
        if chunk.is_empty() {
            continue;
        }
        let (head, deltas_str) = match chunk.find(',') {
            Some(i) => (&chunk[..i], &chunk[i + 1..]),
            None => (chunk, ""),
        };
        let bytes = head.as_bytes();
        if bytes.len() < 3 {
            return Err(ModTagError::MalformedMm(format!("'{chunk}' is too short")));
        }
        if bytes[1] != b'+' && bytes[1] != b'-' {
            return Err(ModTagError::MalformedMm(format!(
                "'{head}' has no '+' or '-' strand"
            )));
        }
        let base = if bytes[0] == b'U' { b'T' } else { bytes[0] };
        if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
            return Err(ModTagError::MalformedMm(format!(
                "'{}' is not a canonical base",
                bytes[0] as char
            )));
        }
        // The codes run to the end of the head, minus the optional '.'/'?' flag.
        let mut codes = &head[2..];
        if codes.ends_with('.') || codes.ends_with('?') {
            codes = &codes[..codes.len() - 1];
        }
        // A numeric ChEBI id denotes a single modification; letters are one
        // modification each.
        let stride = if codes.starts_with(|c: char| c.is_ascii_digit()) {
            if !codes.bytes().all(|c| c.is_ascii_digit()) {
                return Err(ModTagError::MalformedMm(format!(
                    "'{codes}' is neither a ChEBI id nor a list of modification codes"
                )));
            }
            1
        } else {
            if !codes.bytes().all(|c| c.is_ascii_alphabetic()) || codes.is_empty() {
                return Err(ModTagError::MalformedMm(format!(
                    "'{head}' has no modification code"
                )));
            }
            codes.len()
        };

        let mut deltas = Vec::new();
        if !deltas_str.is_empty() {
            for d in deltas_str.split(',') {
                deltas.push(d.parse::<usize>().map_err(|_| {
                    ModTagError::MalformedMm(format!("'{d}' is not a valid skip count"))
                })?);
            }
        }
        lists.push(ModList {
            prefix: head,
            base,
            stride,
            deltas,
        });
    }
    Ok(lists)
}

fn parse_ml(ml: &str) -> Result<Vec<u8>, ModTagError> {
    if ml.is_empty() {
        return Ok(Vec::new());
    }
    ml.split(',')
        .map(|v| {
            v.parse::<u8>()
                .map_err(|_| ModTagError::MalformedMl(format!("'{v}' is not a value in 0-255")))
        })
        .collect()
}

/// Walks to the modified base that follows `skip` further occurrences of
/// `base`, counting the occurrences that fall inside `[start, end)` on the way.
/// The canonical base `N` stands for any base.
fn next_mod_pos(
    seq: &[u8],
    base: u8,
    from: usize,
    skip: usize,
    start: usize,
    end: usize,
    window_occ: &mut usize,
) -> Option<usize> {
    let mut remaining = skip;
    for (i, &b) in seq.iter().enumerate().skip(from) {
        // An `N` list is counted over every base, not just literal Ns, matching
        // htslib (`freq[15] = l_qseq; // all bases count as N for base mods`).
        if base == b'N' || b.to_ascii_uppercase() == base {
            if start <= i && i < end {
                *window_occ += 1;
            }
            if remaining == 0 {
                return Some(i);
            }
            remaining -= 1;
        }
    }
    None
}

/// Recomputes `MM` and `ML` for the subsequence `seq[start..end)`.
///
/// `ML` is optional: `MM` on its own is valid and carries positions without
/// probabilities. Lists that lose all their calls are kept as empty lists so
/// that "this modification was called and none was found here" survives.
fn subset(
    mm: &str,
    ml: Option<&[u8]>,
    seq: &[u8],
    start: usize,
    end: usize,
) -> Result<(String, Vec<u8>), ModTagError> {
    let lists = parse_mm(mm)?;
    let calls: usize = lists.iter().map(|l| l.stride * l.deltas.len()).sum();
    if let Some(ml) = ml {
        if calls != ml.len() {
            return Err(ModTagError::CountMismatch {
                mm: calls,
                ml: ml.len(),
            });
        }
    }

    let mut new_mm = String::with_capacity(mm.len());
    let mut new_ml = Vec::with_capacity(ml.map_or(0, |m| m.len()));
    let mut ml_cursor = 0;

    for list in &lists {
        new_mm.push_str(list.prefix);
        let mut pos = 0;
        // Occurrences of the canonical base seen so far inside [start, end).
        let mut window_occ = 0;
        // Value of window_occ at the last modified base that was kept.
        let mut last_kept_occ = 0;
        for delta in &list.deltas {
            let p = next_mod_pos(seq, list.base, pos, *delta, start, end, &mut window_occ).ok_or(
                ModTagError::BeyondSequence {
                    base: list.base as char,
                },
            )?;
            pos = p + 1;
            if start <= p && p < end {
                // window_occ counts p itself, so the new skip count is the
                // number of occurrences between the previous kept call and p.
                new_mm.push(',');
                new_mm.push_str(&(window_occ - 1 - last_kept_occ).to_string());
                last_kept_occ = window_occ;
                if let Some(ml) = ml {
                    new_ml.extend_from_slice(&ml[ml_cursor..ml_cursor + list.stride]);
                }
            }
            ml_cursor += list.stride;
        }
        new_mm.push(';');
    }
    Ok((new_mm, new_ml))
}

/// Strips a case-insensitive tag prefix, e.g. `MM:Z:` also matches `Mm:Z:`.
fn strip_tag<'a>(field: &'a str, prefix: &str) -> Option<&'a str> {
    let (head, rest) = field.split_at_checked(prefix.len())?;
    head.eq_ignore_ascii_case(prefix).then_some(rest)
}

/// Rewrites the `MM`, `ML` and `MN` tags among the tab-separated `fields` so
/// that they describe `seq[start..end)`. Other tags are passed through
/// untouched, including ones that trimming also invalidates (`qs`, `ns`, `ts`,
/// `du`, ...).
pub fn rewrite_tags(
    fields: &[&str],
    seq: &[u8],
    start: usize,
    end: usize,
) -> Result<Vec<String>, ModTagError> {
    let mut out: Vec<String> = fields.iter().map(|s| s.to_string()).collect();

    if let Some(i) = fields.iter().position(|f| strip_tag(f, "MN:i:").is_some()) {
        let mn: usize = strip_tag(fields[i], "MN:i:")
            .expect("checked above")
            .parse()
            .map_err(|_| ModTagError::MalformedMn(format!("'{}' is not a length", fields[i])))?;
        if mn != seq.len() {
            return Err(ModTagError::StaleMn { mn, len: seq.len() });
        }
        out[i] = format!("MN:i:{}", end - start);
    }

    let Some(mm_i) = fields.iter().position(|f| strip_tag(f, "MM:Z:").is_some()) else {
        // ML without MM would otherwise be written through with stale
        // probabilities, which is the failure this module exists to prevent.
        if fields.iter().any(|f| strip_tag(f, "ML:B:").is_some()) {
            return Err(ModTagError::MlWithoutMm);
        }
        return Ok(out);
    };
    READS_WITH_MODS.fetch_add(1, Ordering::Relaxed);
    let mm = strip_tag(fields[mm_i], "MM:Z:").expect("checked above");

    // ML is optional: MM on its own is a valid way to report positions.
    let ml_i = fields.iter().position(|f| strip_tag(f, "ML:B:").is_some());
    let ml = match ml_i {
        Some(i) => {
            let body = strip_tag(fields[i], "ML:B:").expect("checked above");
            let values = body.strip_prefix('C').ok_or_else(|| {
                ModTagError::MalformedMl(format!("'{}' is not a uint8 array", fields[i]))
            })?;
            let values = match values.strip_prefix(',') {
                Some(v) => v,
                None if values.is_empty() => values,
                None => {
                    return Err(ModTagError::MalformedMl(format!(
                        "'{}' is missing the comma after the array type",
                        fields[i]
                    )))
                }
            };
            Some(parse_ml(values)?)
        }
        None => None,
    };

    let (new_mm, new_ml) = subset(mm, ml.as_deref(), seq, start, end)?;
    out[mm_i] = format!("{}{}", &fields[mm_i][..5], new_mm);
    if let Some(i) = ml_i {
        let mut field = String::from(&fields[i][..6]);
        for v in &new_ml {
            field.push(',');
            field.push_str(&v.to_string());
        }
        out[i] = field;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    // C bases at 1, 4, 8, 12, 16; G bases at 2, 9, 10, 17.
    const SEQ: &[u8] = b"ACGTCAAACGGTCAAACGT";

    /// `MM:Z:C+m?,1,1;` places calls on the C at 4 and the C at 12.
    fn subset_ok(mm: &str, ml: &[u8], start: usize, end: usize) -> (String, Vec<u8>) {
        subset(mm, Some(ml), SEQ, start, end).expect("should rewrite")
    }

    #[test]
    fn untrimmed_range_is_unchanged() {
        let (mm, ml) = subset_ok("C+m?,1,1;", &[200, 10], 0, SEQ.len());
        assert_eq!(("C+m?,1,1;".to_string(), vec![200, 10]), (mm, ml));
    }

    #[test]
    fn headcrop_shifts_the_first_skip_count() {
        // Dropping the C at position 1 leaves one fewer C before the first call.
        let (mm, ml) = subset_ok("C+m?,1,1;", &[200, 10], 2, SEQ.len());
        assert_eq!(("C+m?,0,1;".to_string(), vec![200, 10]), (mm, ml));
    }

    #[test]
    fn tailcrop_drops_calls_and_their_probabilities() {
        let (mm, ml) = subset_ok("C+m?,1,1;", &[200, 10], 0, 10);
        assert_eq!(("C+m?,1;".to_string(), vec![200]), (mm, ml));
    }

    #[test]
    fn segment_keeps_only_the_calls_it_contains() {
        let (mm, ml) = subset_ok("C+m?,1,1;", &[200, 10], 10, SEQ.len());
        assert_eq!(("C+m?,0;".to_string(), vec![10]), (mm, ml));
    }

    #[test]
    fn segment_without_calls_keeps_an_empty_list() {
        // The list is kept so that "called here, nothing found" is preserved.
        let (mm, ml) = subset_ok("C+m?,1,1;", &[200, 10], 5, 12);
        assert_eq!(("C+m?;".to_string(), Vec::<u8>::new()), (mm, ml));
    }

    #[test]
    fn multiple_codes_keep_their_probabilities_interleaved() {
        // C+mh has two codes, so ML holds two values per modified position.
        let (mm, ml) = subset("C+mh?,1,1;", Some(&[200, 5, 10, 3]), SEQ, 10, SEQ.len()).unwrap();
        assert_eq!(("C+mh?,0;".to_string(), vec![10, 3]), (mm, ml));
    }

    #[test]
    fn multiple_lists_are_rewritten_independently() {
        let (mm, ml) =
            subset("C+m?,1,1;G-m,0;", Some(&[200, 10, 250]), SEQ, 10, SEQ.len()).unwrap();
        // The single G call is on the G at position 2, which is trimmed away.
        assert_eq!(("C+m?,0;G-m;".to_string(), vec![10]), (mm, ml));
    }

    #[test]
    fn chebi_codes_are_a_single_modification() {
        let (mm, ml) = subset("C+76792?,1,1;", Some(&[200, 10]), SEQ, 0, SEQ.len()).unwrap();
        assert_eq!(("C+76792?,1,1;".to_string(), vec![200, 10]), (mm, ml));
    }

    #[test]
    fn mm_without_ml_is_valid() {
        let (mm, ml) = subset("C+m?,1,1;", None, SEQ, 10, SEQ.len()).unwrap();
        assert_eq!(("C+m?,0;".to_string(), Vec::<u8>::new()), (mm, ml));
    }

    #[test]
    fn u_is_treated_as_t() {
        // T bases at 3, 11, 18; the first call is on the T at 11.
        let (mm, _) = subset("U+m?,1;", Some(&[200]), SEQ, 0, SEQ.len()).unwrap();
        assert_eq!("U+m?,1;", mm);
    }

    #[test]
    fn mismatched_mm_and_ml_counts_are_an_error() {
        assert_eq!(
            Err(ModTagError::CountMismatch { mm: 2, ml: 1 }),
            subset("C+m?,1,1;", Some(&[200]), SEQ, 0, SEQ.len())
        );
    }

    #[test]
    fn tags_describing_a_longer_sequence_are_an_error() {
        // Four skips need five Cs after the first call; the read has five in total.
        assert_eq!(
            Err(ModTagError::BeyondSequence { base: 'C' }),
            subset("C+m?,1,4;", Some(&[200, 10]), SEQ, 0, SEQ.len())
        );
    }

    #[test]
    fn malformed_mm_is_an_error() {
        for mm in ["C,1,1;", "X+m,1;", "C+,1;", "C+m,x;", "C;"] {
            assert!(
                matches!(
                    subset(mm, None, SEQ, 0, SEQ.len()),
                    Err(ModTagError::MalformedMm(_))
                ),
                "expected {mm} to be rejected"
            );
        }
    }

    #[test]
    fn rewrite_tags_updates_mm_ml_and_mn_and_leaves_others_alone() {
        let fields = ["MM:Z:C+m?,1,1;", "ML:B:C,200,10", "MN:i:19", "qs:f:20.1"];
        let out = rewrite_tags(&fields, SEQ, 10, SEQ.len()).unwrap();
        assert_eq!(
            vec!["MM:Z:C+m?,0;", "ML:B:C,10", "MN:i:9", "qs:f:20.1"],
            out
        );
    }

    #[test]
    fn rewrite_tags_accepts_the_old_lowercase_tag_names() {
        let fields = ["Mm:Z:C+m?,1,1;", "Ml:B:C,200,10"];
        let out = rewrite_tags(&fields, SEQ, 10, SEQ.len()).unwrap();
        assert_eq!(vec!["Mm:Z:C+m?,0;", "Ml:B:C,10"], out);
    }

    #[test]
    fn rewrite_tags_empties_ml_when_no_call_survives() {
        let fields = ["MM:Z:C+m?,1,1;", "ML:B:C,200,10"];
        let out = rewrite_tags(&fields, SEQ, 5, 12).unwrap();
        assert_eq!(vec!["MM:Z:C+m?;", "ML:B:C"], out);
    }

    #[test]
    fn rewrite_tags_passes_through_reads_without_mod_tags() {
        let fields = ["qs:f:20.1", "RG:Z:abc"];
        let out = rewrite_tags(&fields, SEQ, 2, 10).unwrap();
        assert_eq!(vec!["qs:f:20.1", "RG:Z:abc"], out);
    }

    #[test]
    fn stale_mn_is_an_error() {
        // MN says the tags describe a 25 base read, but this one is 19.
        let fields = ["MM:Z:C+m?,1,1;", "ML:B:C,200,10", "MN:i:25"];
        assert_eq!(
            Err(ModTagError::StaleMn { mn: 25, len: 19 }),
            rewrite_tags(&fields, SEQ, 0, SEQ.len())
        );
    }

    #[test]
    fn malformed_ml_is_an_error() {
        let fields = ["MM:Z:C+m?,1;", "ML:B:C,300"];
        assert!(matches!(
            rewrite_tags(&fields, SEQ, 0, SEQ.len()),
            Err(ModTagError::MalformedMl(_))
        ));
    }

    #[test]
    fn n_counts_every_base_not_just_literal_n() {
        // htslib counts all bases towards an N list, so `N+n?,3` lands on the
        // base at index 3 regardless of what the sequence holds.
        let (mm, ml) = subset("N+n?,3;", Some(&[200]), SEQ, 0, SEQ.len()).unwrap();
        assert_eq!(("N+n?,3;".to_string(), vec![200]), (mm, ml));

        // Cropping three bases from the front leaves the call at the new start.
        let (mm, ml) = subset("N+n?,3;", Some(&[200]), SEQ, 3, SEQ.len()).unwrap();
        assert_eq!(("N+n?,0;".to_string(), vec![200]), (mm, ml));
    }

    #[test]
    fn n_lists_do_not_fail_on_reads_without_literal_n() {
        // Counting only literal Ns would make this a fatal error on a read that
        // htslib reads without complaint.
        assert!(subset("N+n?,5,5;", Some(&[200, 10]), SEQ, 0, SEQ.len()).is_ok());
    }

    #[test]
    fn ml_without_mm_is_an_error() {
        assert_eq!(
            Err(ModTagError::MlWithoutMm),
            rewrite_tags(&["ML:B:C,200", "MN:i:19"], SEQ, 2, 10)
        );
    }

    #[test]
    fn malformed_mn_is_reported_as_an_mn_error() {
        assert!(matches!(
            rewrite_tags(&["MM:Z:C+m?,1;", "MN:i:x"], SEQ, 0, SEQ.len()),
            Err(ModTagError::MalformedMn(_))
        ));
    }

    #[test]
    fn ml_without_the_separating_comma_is_an_error() {
        assert!(matches!(
            rewrite_tags(&["MM:Z:C+m?,1;", "ML:B:C200"], SEQ, 0, SEQ.len()),
            Err(ModTagError::MalformedMl(_))
        ));
    }
}
