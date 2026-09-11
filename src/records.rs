use std::borrow::Cow;
use std::io::{BufWriter, Write};

use bio::io::fastq;

use crate::modtags::{rewrite_tags, ModTagError};

pub struct WritableRecord {
    record: String,
}

impl WritableRecord {
    /// Creates a `WritableRecord` from a FASTQ `Record`, restricting it
    /// to the subsequence defined by `start..end`.
    ///
    /// # Arguments
    /// * `record` - Original FASTQ record.
    /// * `start`  - Start index (inclusive).
    /// * `end`    - End index (exclusive).
    /// * `update_mods` - Recompute `MM`/`ML`/`MN` tags for the trimmed range.
    pub fn new(
        record: &fastq::Record,
        start: usize,
        end: usize,
        total_segments: usize,
        segment_idx: usize,
        update_mods: bool,
    ) -> Result<Self, ModTagError> {
        let record =
            record_to_string(record, start, end, total_segments, segment_idx, update_mods)?;

        Ok(WritableRecord { record })
    }

    /// Writes the record to the provided buffer for stdout output.
    pub fn write_on_buffer<W: Write>(
        &self,
        buf: &mut BufWriter<W>,
    ) -> Result<usize, std::io::Error> {
        buf.write(self.record.as_bytes())
    }
}

/// Room reserved for the `_segment_N` suffix so the header needs one allocation.
const SEGMENT_SUFFIX_CAPACITY: usize = 12;

/// Converts a `fastq::record` into a valid FASTQ string within the range `[start..end]`.
fn record_to_string(
    record: &fastq::Record,
    start: usize,
    end: usize,
    total_segments: usize,
    segment_idx: usize,
    update_mods: bool,
) -> Result<String, ModTagError> {
    // The fastq parser splits the header on the first space, which may fall
    // inside a tag value, so SAM tags written as tab-separated fields
    // (`samtools fastq -T MM,ML,...`) can end up spread across the id and the
    // description. Reassemble the original header before picking it apart on
    // tabs, so that every tag is seen wherever the space happened to land.
    let full = match record.desc() {
        Some(d) => Cow::Owned(format!("{} {}", record.id(), d)),
        None => Cow::Borrowed(record.id()),
    };
    // Tags follow the read name as tab-separated fields. Without a tab the
    // whole header is a plain name, optionally followed by a description.
    let (name_field, tags) = match full.split_once('\t') {
        Some((name_field, tags)) => (name_field, Some(tags)),
        None => (full.as_ref(), None),
    };
    // The segment suffix belongs on the read name, not on a description or tag.
    let (name, name_rest) = match name_field.split_once(' ') {
        Some((name, rest)) => (name, Some(rest)),
        None => (name_field, None),
    };

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
