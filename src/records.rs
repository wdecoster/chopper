use std::io::{BufWriter, Write};

use bio::io::fastq;

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
    pub fn new(record: &fastq::Record, start: usize, end: usize, total_segments: usize, segment_idx: usize) -> Self {
        let record = record_to_string(&record, start, end, total_segments, segment_idx);
        
        WritableRecord {
            record,
        }
    }

    /// Writes the record to the provided buffer for stdout output.
    pub fn write_on_buffer<W: Write>(&self, buf: &mut BufWriter<W>) -> Result<usize, std::io::Error> {
        buf.write(self.record.as_bytes())
    }

}

/// Converts a `fastq::record` into a valid FASTQ string within the range `[start..end]`.
fn record_to_string(record: &fastq::Record, start: usize, end: usize, total_segments: usize, segment_idx: usize) -> String {
    // Use a single formatted string with one allocation for the header
    let header = if total_segments > 1 {
        // Add suffix for multiple segments
        match record.desc() {
            Some(d) => format!("@{}_segment_{} {}", record.id(), segment_idx + 1, d),
            None => format!("@{}_segment_{}", record.id(), segment_idx + 1),
        }
    } else {
        // Single segment, use original header
        match record.desc() {
            Some(d) => format!("@{} {}", record.id(), d),
            None => format!("@{}", record.id()),
        }
    };
    
    // Apply the trimming to both sequence and quality data
    let seq_slice = &record.seq()[start..end];
    let qual_slice = &record.qual()[start..end];
    
    format!(
        "{}\n{}\n+\n{}\n",
        header,
        unsafe { std::str::from_utf8_unchecked(seq_slice) },
        unsafe { std::str::from_utf8_unchecked(qual_slice) }
    )
}

#[cfg(test)]
mod tests {
    use bio::io::fastq;

    use crate::records::record_to_string;
    
    #[test]
    fn test_completed_record_to_string() {
        let record = fastq::Record::with_attrs(
            "10-bases",
            None, 
            b"AAAAAAAAAA",
            b"IIIIIIIIII");

        let start = 0;
        let end = 10;
        let total_segments = 1;
        let segment_idx = 0;

        let expected = String::from(
            "@10-bases\nAAAAAAAAAA\n+\nIIIIIIIIII\n");

        let actual = record_to_string(&record, start, end, total_segments, segment_idx);

        assert_eq!(expected, actual);
    }

    #[test]
    fn test_record_to_string_one_segment() {
        let record = fastq::Record::with_attrs(
            "10-bases",
            None, 
            b"TTAAAAAATT",
            b"KKIIIIIIKK");

        let start = 2;
        let end = 8;
        let total_segments = 1;
        let segment_idx = 0;

        let expected = String::from(
            "@10-bases\nAAAAAA\n+\nIIIIII\n");

        let actual = record_to_string(&record, start, end, total_segments, segment_idx);

        assert_eq!(expected, actual);
    }

    #[test]
    #[should_panic]
    fn test_record_to_string_with_no_valid_segment() {
        let record = fastq::Record::with_attrs(
            "10-bases",
            None, 
            b"TTAAAAAATT",
            b"KKIIIIIIKK");

        let start = 8;
        let end = 2;
        let total_segments = 1;
        let segment_idx = 0;


        let _ = record_to_string(&record, start, end, total_segments, segment_idx);
    }

    #[test]
    fn test_record_to_string_multiple_segments() {
        let record = fastq::Record::with_attrs(
            "10-bases",
            None, 
            b"TTAAAAAATT",
            b"KKIIIIIIKK");

        let start = 2;
        let end = 8;
        let total_segments = 2;
        let segment_idx = 1;

        let expected = String::from(
            "@10-bases_segment_2\nAAAAAA\n+\nIIIIII\n");

        let actual = record_to_string(&record, start, end, total_segments, segment_idx);

        assert_eq!(expected, actual);
    }

    #[test]
    fn test_record_to_string_multiple_segments_with_desc() {
        let record = fastq::Record::with_attrs(
            "10-bases",
            Some("description"), 
            b"TTAAAAAATT",
            b"KKIIIIIIKK");

        let start = 2;
        let end = 8;
        let total_segments = 2;
        let segment_idx = 1;

        let expected = String::from(
            "@10-bases_segment_2 description\nAAAAAA\n+\nIIIIII\n");

        let actual = record_to_string(&record, start, end, total_segments, segment_idx);

        assert_eq!(expected, actual);
    }
}
