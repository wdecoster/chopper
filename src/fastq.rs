//! A FASTQ reader and record, in place of the `bio` crate's.
//!
//! On well-formed input it yields exactly what `bio::io::fastq` did: the header
//! with trailing whitespace removed, and sequence and qualities with each
//! line's trailing whitespace (so also `\r`) removed. It differs on input that
//! `bio` mishandled:
//!
//! - A sequence and quality string of different lengths is an error. `bio`
//!   accepted it, and slicing the qualities with sequence offsets then panicked.
//! - Wrapped (multi-line) records pair qualities to the sequence by length,
//!   not by counting lines. A sequence on one line (the usual case) takes
//!   exactly one quality line, so a short one is reported for that read
//!   instead of swallowing the next record.
//! - A read without bases (`@r\n\n+\n\n`) is a valid, empty record, which the
//!   length filter drops, rather than an "incomplete record" error.
//! - Blank lines between records and at the end of the input are skipped.
//! - Errors name the line they occurred on.
//!
//! Unlike `bio`, the header is kept whole rather than split into an id and a
//! description at the first space, so it never has to be reassembled.

use std::fmt;
use std::io::BufRead;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Record {
    header: String,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

impl Record {
    /// Builds a record from an id and an optional description, which are
    /// joined by a space into the header, as they are on a FASTQ header line.
    pub fn with_attrs(id: &str, desc: Option<&str>, seq: &[u8], qual: &[u8]) -> Self {
        let header = match desc {
            Some(desc) => format!("{id} {desc}"),
            None => id.to_string(),
        };
        Record {
            header,
            seq: seq.to_vec(),
            qual: qual.to_vec(),
        }
    }

    /// The header line without the leading `@`: the read name, then anything
    /// after it (a description, or SAM tags as tab-separated fields).
    pub fn header(&self) -> &str {
        &self.header
    }

    pub fn seq(&self) -> &[u8] {
        &self.seq
    }

    /// Phred+33 encoded qualities, as many as there are bases.
    pub fn qual(&self) -> &[u8] {
        &self.qual
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct Error {
    /// 1-based line number in the input.
    pub line: usize,
    pub message: String,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "line {}: {}", self.line, self.message)
    }
}

impl std::error::Error for Error {}

pub struct Reader<R> {
    inner: R,
    line: Vec<u8>,
    line_number: usize,
}

impl<R: BufRead> Reader<R> {
    pub fn new(inner: R) -> Self {
        Reader {
            inner,
            line: Vec::new(),
            line_number: 0,
        }
    }

    pub fn records(self) -> Records<R> {
        Records {
            reader: self,
            failed: false,
        }
    }

    /// Reads the next line into `self.line`, returning false at end of input.
    fn next_line(&mut self) -> Result<bool, Error> {
        self.line.clear();
        let n = self
            .inner
            .read_until(b'\n', &mut self.line)
            .map_err(|e| self.error(format!("cannot read input: {e}")))?;
        if n > 0 {
            self.line_number += 1;
        }
        Ok(n > 0)
    }

    fn error(&self, message: impl Into<String>) -> Error {
        Error {
            line: self.line_number,
            message: message.into(),
        }
    }

    /// Reads the next record, or returns `None` at the end of the input.
    pub fn read(&mut self) -> Option<Result<Record, Error>> {
        // Skip blank lines before a record, including at the end of the input.
        loop {
            match self.next_line() {
                Ok(false) => return None,
                Ok(true) if self.line.trim_ascii().is_empty() => continue,
                Ok(true) => break,
                Err(e) => return Some(Err(e)),
            }
        }
        Some(self.read_record())
    }

    fn read_record(&mut self) -> Result<Record, Error> {
        if self.line[0] != b'@' {
            return Err(self.error("expected '@' at the start of a record"));
        }
        let header = std::str::from_utf8(&self.line[1..])
            .map_err(|_| self.error("the header is not valid UTF-8"))?
            .trim_end()
            .to_string();
        let start = self.line_number;
        let truncated = |reader: &Self| {
            reader.error(format!("the record starting on line {start} is truncated"))
        };

        let mut seq = Vec::new();
        let mut seq_lines = 0;
        loop {
            if !self.next_line()? {
                return Err(truncated(self));
            }
            if self.line.first() == Some(&b'+') {
                break;
            }
            seq.extend_from_slice(self.line.trim_ascii_end());
            seq_lines += 1;
        }
        if seq_lines == 0 {
            return Err(self.error(format!(
                "the record starting on line {start} has no sequence line"
            )));
        }

        let mut qual = Vec::with_capacity(seq.len());
        loop {
            if !self.next_line()? {
                return Err(truncated(self));
            }
            qual.extend_from_slice(self.line.trim_ascii_end());
            if seq_lines == 1 || qual.len() >= seq.len() {
                break;
            }
        }
        if qual.len() != seq.len() {
            return Err(self.error(format!(
                "read {} has {} bases but {} quality values",
                header.split(['\t', ' ']).next().unwrap_or(&header),
                seq.len(),
                qual.len()
            )));
        }

        Ok(Record { header, seq, qual })
    }
}

/// Iterates over the records, ending after the first error: once a record
/// is malformed, where the next one starts is anyone's guess.
pub struct Records<R> {
    reader: Reader<R>,
    failed: bool,
}

impl<R: BufRead> Iterator for Records<R> {
    type Item = Result<Record, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.failed {
            return None;
        }
        let next = self.reader.read();
        self.failed = matches!(next, Some(Err(_)));
        next
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(input: &str) -> Vec<Result<Record, Error>> {
        Reader::new(input.as_bytes()).records().collect()
    }

    fn record(header: &str, seq: &str, qual: &str) -> Result<Record, Error> {
        Ok(Record::with_attrs(
            header,
            None,
            seq.as_bytes(),
            qual.as_bytes(),
        ))
    }

    fn error_on(line: usize, input: &str) -> String {
        match parse(input).pop() {
            Some(Err(e)) => {
                assert_eq!(e.line, line, "{e}");
                e.message
            }
            other => panic!("expected an error, got {other:?}"),
        }
    }

    #[test]
    fn reads_four_line_records() {
        assert_eq!(
            parse("@r1 some description\nACGT\n+\nIIII\n@r2\nGG\n+r2\n!!\n"),
            [
                record("r1 some description", "ACGT", "IIII"),
                record("r2", "GG", "!!")
            ]
        );
    }

    #[test]
    fn keeps_the_header_whole() {
        // SAM tags after a tab, a tag value with a space, repeated spaces.
        assert_eq!(
            parse("@r1\tMM:Z:C+m?,1;\tCO:Z:a b  c\nACGT\n+\nIIII\n"),
            [record("r1\tMM:Z:C+m?,1;\tCO:Z:a b  c", "ACGT", "IIII")]
        );
    }

    #[test]
    fn strips_line_endings_and_trailing_whitespace_like_bio() {
        assert_eq!(
            parse("@r1 d \t\r\nACGT \r\n+\r\nIIII\t\r\n"),
            [record("r1 d", "ACGT", "IIII")]
        );
        assert_eq!(parse("@r1\nACGT\n+\nIIII"), [record("r1", "ACGT", "IIII")]);
    }

    #[test]
    fn quality_lines_may_start_with_at_or_plus() {
        assert_eq!(
            parse("@r1\nACGT\n+\n@+II\n@r2\nA\n+\n+\n"),
            [record("r1", "ACGT", "@+II"), record("r2", "A", "+")]
        );
    }

    #[test]
    fn wrapped_records_are_paired_by_length() {
        // Qualities wrapped differently from the sequence.
        assert_eq!(
            parse("@r1\nACGT\nAC\n+\nIII\nIII\n@r2\nA\n+\nI\n"),
            [record("r1", "ACGTAC", "IIIIII"), record("r2", "A", "I")]
        );
    }

    #[test]
    fn a_read_without_bases_is_an_empty_record() {
        let records = parse("@empty\n\n+\n\n@r2\nA\n+\nI\n");
        assert_eq!(records, [record("empty", "", ""), record("r2", "A", "I")]);
        assert!(records[0].as_ref().unwrap().is_empty());
    }

    #[test]
    fn blank_lines_between_and_after_records_are_skipped() {
        assert_eq!(
            parse("\n@r1\nA\n+\nI\n\n\n@r2\nC\n+\nI\n\n"),
            [record("r1", "A", "I"), record("r2", "C", "I")]
        );
        assert_eq!(parse(""), []);
        assert_eq!(parse("\n\n"), []);
    }

    #[test]
    fn mismatched_lengths_are_an_error() {
        let message = error_on(4, "@r1 desc\nACGT\n+\nII\n");
        assert_eq!(message, "read r1 has 4 bases but 2 quality values");
        // Too many qualities on the last quality line.
        error_on(4, "@r1\nACGT\n+\nIIIIII\n");
    }

    #[test]
    fn malformed_records_are_errors_with_their_line() {
        assert!(error_on(1, "r1\nACGT\n+\nIIII\n").contains("expected '@'"));
        assert!(error_on(5, "@r1\nA\n+\nI\nr2\n").contains("expected '@'"));
        assert!(error_on(3, "@r1\nACGT\n+\n").contains("truncated"));
        assert!(error_on(2, "@r1\nACGT\n").contains("truncated"));
        assert!(error_on(5, "@r1\nAC\nGT\n+\nII\n").contains("truncated"));
        assert!(error_on(2, "@r1\n+\nII\n").contains("no sequence line"));
        // A short quality line does not swallow the next record.
        assert!(error_on(4, "@r1\nACGT\n+\nII\n@r2\nA\n+\nI\n").contains("2 quality values"));
    }

    #[test]
    fn a_header_that_is_not_utf8_is_an_error() {
        let input = b"@r1\xff\nACGT\n+\nIIII\n";
        let records: Vec<_> = Reader::new(&input[..]).records().collect();
        assert_eq!(
            records,
            [Err(Error {
                line: 1,
                message: "the header is not valid UTF-8".into()
            })]
        );
    }
}
