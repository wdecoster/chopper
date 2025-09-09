use bio::io::fastq;

pub struct WritableRecord {
    record: fastq::Record,
    start: usize,
    end: usize,
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
        WritableRecord {
            record,
            start,
            end,
        }
    }
    
    /// Write a record to stdout
    pub fn write_record(&self) {
        // Use a single formatted string with one allocation for the header
        let header = match self.record.desc() {
            Some(d) => format!("@{} {}", self.record.id(), d),
            None => format!("@{}", self.record.id()),
        };
        
        // Apply the trimming to both sequence and quality data
        let seq_slice = &self.record.seq()[self.start..self.end];
        let qual_slice = &self.record.qual()[self.start..self.end];
        
        // Use a single print to minimize syscalls
        println!(
            "{}\n{}\n+\n{}",
            header,
            unsafe { std::str::from_utf8_unchecked(seq_slice) },
            unsafe { std::str::from_utf8_unchecked(qual_slice) }
        );
    }
}
