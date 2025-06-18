use bio::io::fastq;
use super::phred_score_to_probability;

/// A trait for implementing custom read trimming strategies.
///
/// Implementing this trait allows you to define how a read should be trimmed  
/// according to a specific algorithm or criteria. This is useful when working  
/// with different trimming approaches.
///
/// The `trim` method should return `Some((start, end))` if a valid sub-read  
/// is identified, or `None` if the read should be discarded.
pub trait TrimStrategy {
    /// Performs the trimming process based on the defined strategy.
    ///
    /// Returns `Some((start, end))` with the indices of the trimmed region  
    /// if a valid sub-read is identified, or `None` if the read should be discarded.
    fn trim(&self, record: &fastq::Record) -> Option<(usize, usize)>;
}

/// A trimming strategy that applies a modified Mott algorithm to extract the highest-quality sub-read.
///
/// This strategy identifies the sub-read with the lowest cumulative error probability,
/// following an approach similar to Phred quality trimming. It is useful for removing
/// low-quality bases typically found at the ends of reads.
///
/// When used, it returns `Some((start, end))` with the indices of the highest-quality
/// sub-read if one is found; otherwise, returns None.
pub struct HighestQualityTrimStrategy {
    cutoff: f64,
}

impl HighestQualityTrimStrategy {
    pub fn new(cutoff: f64) -> HighestQualityTrimStrategy  {
        HighestQualityTrimStrategy { cutoff }
    }
}

impl TrimStrategy for HighestQualityTrimStrategy {
    fn trim(&self, record: &fastq::Record) -> Option<(usize, usize)> {
        let mut best_start = usize::MAX;
        let mut best_end = usize::MAX;
        let mut best_cumulative_error = 0.0;
        let mut best_length = 0;

        let mut current_start = 0;
        let mut current_end;
        let mut current_cumulative_error = -1.0;
        for (i, phred_qual) in record.qual().iter().enumerate() {
            let prob_error = self.cutoff - phred_score_to_probability(*phred_qual);
            if current_cumulative_error < 0.0 {
                current_cumulative_error = 0.0;
                current_start = i;
            } 

            current_cumulative_error += prob_error;
            current_end = i;
            
            if best_cumulative_error < current_cumulative_error ||
                (best_cumulative_error == current_cumulative_error && best_length < current_end - current_start + 1) {
                best_start = current_start;
                best_end = current_end;
                best_cumulative_error = current_cumulative_error;
                best_length = current_end - current_start + 1;
            }
        }

        if best_start == usize::MAX {
            None
        } else {
            Some((best_start, best_end+1))
        }
    }
}

#[cfg(test)]
mod tests {
    use bio::io::fastq;
    use super::*;
    
    #[test]
    fn highest_quality_sub_read_test() {
        let cases: [(fastq::Record, f64, Option<(usize, usize)>); 6] = [
            (
                fastq::Record::with_attrs(
                    &"id-01",
                    None,
                    b"AAAAAAAAAAAAAAATTTAA",
                    b"&#3-G27C:(@G7B55+C4I"
                ), 0.01, Some((10, 16)), // cutoff: Q20
            ), ( // TODO: calcular
                fastq::Record::with_attrs(
                    &"id-01",
                    None,
                    b"TTTTTTTTTTTTTTTTTTTT",
                    b"77%'24)FAF9@=94'%054"
                ), 0.19952623149688797, Some((0, 20)), // cutoff: Q7
            ), (
                fastq::Record::with_attrs(
                    &"id-01",
                    None,
                    b"AAAAAAAAAAAAAAATTTTA",
                    b"'8$-BF2!C;+59->H@91#"
                ), 0.03162277660168379, Some((11, 19)), // cutoff: Q15
            ), (
                fastq::Record::with_attrs(
                    &"id-01",
                    None,
                    b"AAAAAAAAAAAAAAAAAAAA",
                    b"%,42$CH*#0+0C6=0,*6/"
                ), 0.0001, None, // cutoff: Q40
            ), (
                fastq::Record::with_attrs(
                    &"id-01",
                    None,
                    b"AAAAAAAAAAAAAAAAAAAT",
                    b"-------------------J"
                ), 0.0001, Some((19,20)), // cutoff: Q40
            ), (
                fastq::Record::with_attrs(
                    &"id-01",
                    None,
                    b"TAAAAAAAAAAAAAAAAAAA",
                    b"I-------------------"
                ), 0.0001, Some((0,1)), // cutoff: Q40
            )
        ];

        for (read, cutoff, expected) in cases {
            let trim_strategy = HighestQualityTrimStrategy::new(cutoff);
            assert_eq!(trim_strategy.trim(&read) , expected);
        }
    }
}
