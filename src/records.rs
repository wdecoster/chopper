use std::borrow::Cow;
use std::fmt;
use std::io::{BufWriter, Write};

use bio::io::fastq;
use noodles_sam::alignment::{
    record::Flags,
    record_buf::{QualityScores, Sequence},
    RecordBuf,
};

use crate::formats::{text_to_tags, Format, InputRecord};
use crate::modtags::{rewrite_tags, ModTagError};

/// A kept read or segment, already in the output format.
pub enum WritableRecord {
    Fastq(String),
    Bam(RecordBuf),
}

/// Why a read could not be written out.
#[derive(Debug)]
pub enum RecordError {
    ModTags(ModTagError),
    /// The read cannot be represented as a BAM record.
    Bam(String),
}

impl From<ModTagError> for RecordError {
    fn from(e: ModTagError) -> Self {
        RecordError::ModTags(e)
    }
}

impl fmt::Display for RecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RecordError::ModTags(e) => e.fmt(f),
            RecordError::Bam(e) => f.write_str(e),
        }
    }
}

impl WritableRecord {
    /// Creates a `WritableRecord` from an input read, restricting it
    /// to the subsequence defined by `start..end`.
    ///
    /// # Arguments
    /// * `record` - Original input read.
    /// * `start`  - Start index (inclusive).
    /// * `end`    - End index (exclusive).
    /// * `update_mods` - Recompute `MM`/`ML`/`MN` tags for the trimmed range.
    /// * `format` - Output format.
    pub fn new(
        record: &InputRecord,
        start: usize,
        end: usize,
        total_segments: usize,
        segment_idx: usize,
        update_mods: bool,
        format: Format,
    ) -> Result<Self, RecordError> {
        match format {
            Format::Fastq => Ok(WritableRecord::Fastq(record_to_string(
                &record.fastq,
                start,
                end,
                total_segments,
                segment_idx,
                update_mods,
            )?)),
            Format::Bam => Ok(WritableRecord::Bam(record_to_bam(
                record,
                start,
                end,
                total_segments,
                segment_idx,
                update_mods,
            )?)),
        }
    }

    /// Writes a FASTQ record to the provided buffer for stdout output.
    ///
    /// Uses `write_all` rather than `write`: a record at or above the buffer's
    /// capacity bypasses the buffer and reaches the writer in a single `write`
    /// call, which is free to consume only part of it. FASTQ records for long
    /// reads routinely exceed the 8 KiB default, so a short write there would
    /// truncate a read mid-sequence.
    pub fn write_on_buffer<W: Write>(&self, buf: &mut BufWriter<W>) -> Result<(), std::io::Error> {
        match self {
            WritableRecord::Fastq(record) => buf.write_all(record.as_bytes()),
            WritableRecord::Bam(_) => unreachable!("BAM records go through the BAM writer"),
        }
    }
}

/// Room reserved for the `_segment_N` suffix so the header needs one allocation.
const SEGMENT_SUFFIX_CAPACITY: usize = 12;

/// The header line as it was in the input, without the leading `@`.
fn full_header(record: &fastq::Record) -> Cow<'_, str> {
    // The fastq parser splits the header on the first space, which may fall
    // inside a tag value, so SAM tags written as tab-separated fields
    // (`samtools fastq -T MM,ML,...`) can end up spread across the id and the
    // description. Reassemble the original header before picking it apart on
    // tabs, so that every tag is seen wherever the space happened to land.
    match record.desc() {
        Some(d) => Cow::Owned(format!("{} {}", record.id(), d)),
        None => Cow::Borrowed(record.id()),
    }
}

/// Splits a header into the read name, the free text after it, and the
/// tab-separated SAM tags.
fn split_header(full: &str) -> (&str, Option<&str>, Option<&str>) {
    // Tags follow the read name as tab-separated fields. Without a tab the
    // whole header is a plain name, optionally followed by a description.
    let (name_field, tags) = match full.split_once('\t') {
        Some((name_field, tags)) => (name_field, Some(tags)),
        None => (full, None),
    };
    // The segment suffix belongs on the read name, not on a description or tag.
    let (name, name_rest) = match name_field.split_once(' ') {
        Some((name, rest)) => (name, Some(rest)),
        None => (name_field, None),
    };
    (name, name_rest, tags)
}

/// Converts a `fastq::record` into a valid FASTQ string within the range `[start..end]`.
fn record_to_string(
    record: &fastq::Record,
    start: usize,
    end: usize,
    total_segments: usize,
    segment_idx: usize,
    update_mods: bool,
) -> Result<String, ModTagError> {
    let full = full_header(record);
    let (name, name_rest, tags) = split_header(&full);

    let mut header = String::with_capacity(1 + full.len() + SEGMENT_SUFFIX_CAPACITY);
    header.push('@');
    header.push_str(name);
    if total_segments > 1 {
        // Add suffix for multiple segments
        header.push_str("_segment_");
        header.push_str(&(segment_idx + 1).to_string());
    }
    if let Some(rest) = name_rest {
        header.push(' ');
        header.push_str(rest);
    }
    if let Some(t) = tags {
        if update_mods {
            let fields: Vec<&str> = t.split('\t').collect();
            for field in rewrite_tags(&fields, record.seq(), start, end)? {
                header.push('\t');
                header.push_str(&field);
            }
        } else {
            header.push('\t');
            header.push_str(t);
        }
    }

    // Apply the trimming to both sequence and quality data
    let seq_slice = &record.seq()[start..end];
    let qual_slice = &record.qual()[start..end];

    Ok(format!(
        "{}\n{}\n+\n{}\n",
        header,
        unsafe { std::str::from_utf8_unchecked(seq_slice) },
        unsafe { std::str::from_utf8_unchecked(qual_slice) }
    ))
}

/// Builds a BAM record for the subsequence `[start..end]` of a read.
///
/// For BAM input this starts from the original record, so its flags and every
/// tag keep their exact type and order; only the name, sequence and qualities
/// are replaced, and under `update_mods` the modification tags. For FASTQ
/// input the tab-separated SAM tags in the header become typed BAM tags, and
/// any free text after the read name is kept in a `CO:Z:` comment tag.
fn record_to_bam(
    record: &InputRecord,
    start: usize,
    end: usize,
    total_segments: usize,
    segment_idx: usize,
    update_mods: bool,
) -> Result<RecordBuf, RecordError> {
    let fastq = &record.fastq;
    let full = full_header(fastq);
    let (name, name_rest, tags) = split_header(&full);

    let fields: Vec<&str> = tags.map(|t| t.split('\t').collect()).unwrap_or_default();
    let fields: Vec<Cow<str>> = if update_mods {
        rewrite_tags(&fields, fastq.seq(), start, end)?
            .into_iter()
            .map(Cow::Owned)
            .collect()
    } else {
        fields.into_iter().map(Cow::Borrowed).collect()
    };

    let mut out = match &record.bam {
        Some(original) => {
            // The FASTQ view of a BAM record only carries the tags that
            // `--update-mods` rewrites; put the rewritten values in place.
            let mut out = original.clone();
            let rewritten = text_to_tags(&fields).map_err(RecordError::Bam)?;
            for (tag, value) in rewritten.iter() {
                out.data_mut().insert(tag, value.clone());
            }
            out
        }
        None => {
            let mut fields = fields;
            if let Some(text) = name_rest {
                if fields.iter().any(|f| f.starts_with("CO:")) {
                    return Err(RecordError::Bam(format!(
                        "the header has both free text ('{text}') and a CO tag, \
                         and BAM can only hold one comment"
                    )));
                }
                fields.push(Cow::Owned(format!("CO:Z:{text}")));
            }
            let mut out = RecordBuf::builder().set_flags(Flags::UNMAPPED).build();
            *out.data_mut() = text_to_tags(&fields).map_err(RecordError::Bam)?;
            out
        }
    };

    let mut new_name = String::with_capacity(name.len() + SEGMENT_SUFFIX_CAPACITY);
    new_name.push_str(name);
    if total_segments > 1 {
        new_name.push_str("_segment_");
        new_name.push_str(&(segment_idx + 1).to_string());
    }
    *out.name_mut() = Some(new_name.into());
    *out.sequence_mut() = Sequence::from(fastq.seq()[start..end].to_vec());
    let qual = fastq.qual()[start..end]
        .iter()
        .map(|&q| q.checked_sub(33))
        .collect::<Option<Vec<u8>>>()
        .ok_or_else(|| RecordError::Bam("quality string has characters below '!'".into()))?;
    *out.quality_scores_mut() = QualityScores::from(qual);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use bio::io::fastq;

    use crate::records::record_to_string;

    #[test]
    fn test_completed_record_to_string() {
        let record = fastq::Record::with_attrs("10-bases", None, b"AAAAAAAAAA", b"IIIIIIIIII");

        let start = 0;
        let end = 10;
        let total_segments = 1;
        let segment_idx = 0;

        let expected = String::from("@10-bases\nAAAAAAAAAA\n+\nIIIIIIIIII\n");

        let actual =
            record_to_string(&record, start, end, total_segments, segment_idx, false).unwrap();

        assert_eq!(expected, actual);
    }

    #[test]
    fn test_record_to_string_one_segment() {
        let record = fastq::Record::with_attrs("10-bases", None, b"TTAAAAAATT", b"KKIIIIIIKK");

        let start = 2;
        let end = 8;
        let total_segments = 1;
        let segment_idx = 0;

        let expected = String::from("@10-bases\nAAAAAA\n+\nIIIIII\n");

        let actual =
            record_to_string(&record, start, end, total_segments, segment_idx, false).unwrap();

        assert_eq!(expected, actual);
    }

    #[test]
    #[should_panic]
    fn test_record_to_string_with_no_valid_segment() {
        let record = fastq::Record::with_attrs("10-bases", None, b"TTAAAAAATT", b"KKIIIIIIKK");

        let start = 8;
        let end = 2;
        let total_segments = 1;
        let segment_idx = 0;

        let _ = record_to_string(&record, start, end, total_segments, segment_idx, false).unwrap();
    }

    #[test]
    fn test_record_to_string_multiple_segments() {
        let record = fastq::Record::with_attrs("10-bases", None, b"TTAAAAAATT", b"KKIIIIIIKK");

        let start = 2;
        let end = 8;
        let total_segments = 2;
        let segment_idx = 1;

        let expected = String::from("@10-bases_segment_2\nAAAAAA\n+\nIIIIII\n");

        let actual =
            record_to_string(&record, start, end, total_segments, segment_idx, false).unwrap();

        assert_eq!(expected, actual);
    }

    #[test]
    fn test_record_to_string_multiple_segments_with_desc() {
        let record = fastq::Record::with_attrs(
            "10-bases",
            Some("description"),
            b"TTAAAAAATT",
            b"KKIIIIIIKK",
        );

        let start = 2;
        let end = 8;
        let total_segments = 2;
        let segment_idx = 1;

        let expected = String::from("@10-bases_segment_2 description\nAAAAAA\n+\nIIIIII\n");

        let actual =
            record_to_string(&record, start, end, total_segments, segment_idx, false).unwrap();

        assert_eq!(expected, actual);
    }

    #[test]
    fn test_record_to_string_keeps_sam_tags_on_single_segment() {
        // `samtools fastq -T MM,ML` writes tags as tab-separated header fields,
        // which the parser leaves inside the id.
        let record = fastq::Record::with_attrs(
            "10-bases\tMM:Z:C+m?,1;\tML:B:C,200",
            None,
            b"TTAAAAAATT",
            b"KKIIIIIIKK",
        );

        let expected = String::from("@10-bases\tMM:Z:C+m?,1;\tML:B:C,200\nAAAAAA\n+\nIIIIII\n");

        let actual = record_to_string(&record, 2, 8, 1, 0, false).unwrap();

        assert_eq!(expected, actual);
    }

    #[test]
    fn test_record_to_string_segment_suffix_goes_on_name_not_tags() {
        let record = fastq::Record::with_attrs(
            "10-bases\tMM:Z:C+m?,1;\tML:B:C,200",
            None,
            b"TTAAAAAATT",
            b"KKIIIIIIKK",
        );

        let expected =
            String::from("@10-bases_segment_2\tMM:Z:C+m?,1;\tML:B:C,200\nAAAAAA\n+\nIIIIII\n");

        let actual = record_to_string(&record, 2, 8, 2, 1, false).unwrap();

        assert_eq!(expected, actual);
    }

    #[test]
    fn test_record_to_string_segment_suffix_with_tags_and_desc() {
        let record = fastq::Record::with_attrs(
            "10-bases\tMM:Z:C+m?,1;",
            Some("description"),
            b"TTAAAAAATT",
            b"KKIIIIIIKK",
        );

        let expected =
            String::from("@10-bases_segment_2\tMM:Z:C+m?,1; description\nAAAAAA\n+\nIIIIII\n");

        let actual = record_to_string(&record, 2, 8, 2, 1, false).unwrap();

        assert_eq!(expected, actual);
    }

    #[test]
    fn test_record_to_string_updates_mod_tags_when_requested() {
        // C bases at 1, 4, 8; MM places a call on the C at 4.
        let record = fastq::Record::with_attrs(
            "read1\tMM:Z:C+m?,1;\tML:B:C,200\tqs:f:20.1",
            None,
            b"ACGTCAAAC",
            b"IIIIIIIII",
        );

        // Cropping the first two bases removes the C at 1.
        let expected =
            String::from("@read1\tMM:Z:C+m?,0;\tML:B:C,200\tqs:f:20.1\nGTCAAAC\n+\nIIIIIII\n");

        let actual = record_to_string(&record, 2, 9, 1, 0, true).unwrap();

        assert_eq!(expected, actual);
    }

    #[test]
    fn test_record_to_string_leaves_mod_tags_alone_by_default() {
        let record = fastq::Record::with_attrs(
            "read1\tMM:Z:C+m?,1;\tML:B:C,200",
            None,
            b"ACGTCAAAC",
            b"IIIIIIIII",
        );

        let expected = String::from("@read1\tMM:Z:C+m?,1;\tML:B:C,200\nGTCAAAC\n+\nIIIIIII\n");

        let actual = record_to_string(&record, 2, 9, 1, 0, false).unwrap();

        assert_eq!(expected, actual);
    }

    #[test]
    fn test_record_to_string_reports_broken_mod_tags() {
        // MM skips more C bases than the read holds.
        let record = fastq::Record::with_attrs(
            "read1\tMM:Z:C+m?,9;\tML:B:C,200",
            None,
            b"ACGTCAAAC",
            b"IIIIIIIII",
        );

        assert!(record_to_string(&record, 0, 9, 1, 0, true).is_err());
    }

    #[test]
    fn test_record_to_string_finds_tags_after_a_tag_value_with_a_space() {
        // The parser splits on the first space, so a tag value containing one
        // pushes every later tag into the description. They must still be found.
        let record = fastq::Record::with_attrs(
            "read1\tCO:Z:some",
            Some("comment\tMM:Z:C+m?,1;\tML:B:C,200"),
            b"ACGTCAAAC",
            b"IIIIIIIII",
        );

        let expected = String::from(
            "@read1\tCO:Z:some comment\tMM:Z:C+m?,0;\tML:B:C,200\nGTCAAAC\n+\nIIIIIII\n",
        );

        let actual = record_to_string(&record, 2, 9, 1, 0, true).unwrap();

        assert_eq!(expected, actual);
    }

    #[test]
    fn test_record_to_string_suffixes_the_name_not_the_description() {
        let record = fastq::Record::with_attrs(
            "read1",
            Some("a free text description"),
            b"ACGTCAAAC",
            b"IIIIIIIII",
        );

        let expected =
            String::from("@read1_segment_2 a free text description\nACGTCAAAC\n+\nIIIIIIIII\n");

        let actual = record_to_string(&record, 0, 9, 2, 1, false).unwrap();

        assert_eq!(expected, actual);
    }
}
