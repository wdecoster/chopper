use std::io::{BufWriter, StdoutLock, Write};

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
    pub fn new(record: fastq::Record, start: usize, end: usize) -> Self {
        let record = WritableRecord::record_to_string(&record, start, end);
        
        WritableRecord {
            record,
        }
    }

    /// Writes the record to the provided buffer for stdout output.
    pub fn write_on_buffer(&self, buf: &mut BufWriter<StdoutLock>) -> Result<usize, std::io::Error> {
        buf.write(self.record.as_bytes())
    }

    /// Converts a `fastq::record` into a valid FASTQ string within the range `[start..end]`.
    fn record_to_string(record: &fastq::Record, start: usize, end: usize) -> String {
        // Use a single formatted string with one allocation for the header
        let header = match record.desc() {
            Some(d) => format!("@{} {}", record.id(), d),
            None => format!("@{}", record.id()),
        };
        
        // Apply the trimming to both sequence and quality data
        let seq_slice = &record.seq()[start..end];
        let qual_slice = &record.qual()[start..end];
        
        format!(
            "{}\n{}\n+\n{}",
            header,
            unsafe { std::str::from_utf8_unchecked(seq_slice) },
            unsafe { std::str::from_utf8_unchecked(qual_slice) }
        )
    }
}
