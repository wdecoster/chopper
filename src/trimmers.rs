use bio::io::fastq;
use super::phred_score_to_probability;

/// A trait for implementing custom read trimming strategies.
///
/// Implementing this trait allows you to define how a read should be trimmed  
/// according to a specific algorithm or criteria. This is useful when working  
/// with different trimming approaches.
///
/// The `trim` method should return `Some((start, end))` if a valid read segment
/// is identified, or `None` if the read should be discarded.
pub trait TrimStrategy: Send + Sync {
    /// Performs the trimming process based on the defined strategy.
    ///
    /// Returns `Some((start, end))` with the indices of the trimmed region  
    /// if a valid read segment is identified, or `None` if the read should be discarded.
    fn trim(&self, record: &fastq::Record) -> Option<(usize, usize)>;
}

/// A trimming strategy that applies a modified Mott algorithm to extract the highest-quality read
/// segment.
///
/// This strategy identifies the read segment with the lowest cumulative error probability,
/// following an approach similar to Phred quality trimming. It is useful for removing
/// low-quality bases typically found at the ends of reads.
///
/// When used, it returns `Some((start, end))` with the indices of the highest-quality
/// read segment if one is found; otherwise, returns None.
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
        let mut current_cumulative_error = -1.0;
        for (i, phred_qual) in record.qual().iter().enumerate() {
            // The FASTQ Phred quality score must be converted from its ASCII representation
            // before calling phred_score_to_probability
            let prob_error = self.cutoff - phred_score_to_probability(*phred_qual - 33);
            if current_cumulative_error < 0.0 {
                current_cumulative_error = 0.0;
                current_start = i;
            } 

            current_cumulative_error += prob_error;
            
            if best_cumulative_error < current_cumulative_error ||
                (best_cumulative_error == current_cumulative_error && best_length < i - current_start + 1) {
                best_start = current_start;
                best_end = i;
                best_cumulative_error = current_cumulative_error;
                best_length = i - current_start + 1;
            }
        }

        if best_start == usize::MAX {
            None
        } else {
            Some((best_start, best_end+1))
        }
    }
}

/// This strategy trims low-quality bases from both ends of the read 
/// until it reaches a base with a quality score ≥ `cutoff` (Q-score).
///
/// Returns `Some((start, end))` with the indices of the valid read segment
/// if one is found; otherwise, returns `None`.
pub struct TrimByQualityStrategy {
    cutoff: u8,
}

impl TrimByQualityStrategy {
    pub fn new(cutoff: u8) -> TrimByQualityStrategy {
        TrimByQualityStrategy { cutoff }
    }
}

impl TrimStrategy for TrimByQualityStrategy {
    fn trim(&self, record: &fastq::Record) -> Option<(usize, usize)> {
        let read_len = record.seq().len();
        let trim_threshold = self.cutoff;

        // Quality-based trimming
        let quals = record.qual();
        
        // Find first position from start with quality >= threshold
        let mut trim_start = 0;
        while trim_start < read_len && (quals[trim_start] - 33) < trim_threshold {
            trim_start += 1;
        }
        
        // Find first position from end with quality >= threshold
        let mut trim_end = read_len;
        while trim_end > trim_start && (quals[trim_end - 1] - 33) < trim_threshold {
            trim_end -= 1;
        }

        if trim_end <= trim_start {
            return None;
        }
        
        Some((trim_start, trim_end))
    }
}

/// This strategy removes a fixed number of bases from both ends of the read.
pub struct FixedCropStrategy {
    head_crop: usize,
    tail_crop: usize,
}

impl FixedCropStrategy {
    pub fn new(head_crop: usize, tail_crop: usize) -> FixedCropStrategy {
        FixedCropStrategy { head_crop, tail_crop }
    }
}

impl TrimStrategy for FixedCropStrategy {
    fn trim(&self, record: &fastq::Record) -> Option<(usize, usize)> {
        if record.seq().len() - self.tail_crop <= self.head_crop {
            return None;
        }
        Some((self.head_crop, record.seq().len() - self.tail_crop))
    }
}

#[cfg(test)]
mod tests {
    use bio::io::fastq;
    use super::*;

    fn get_reads() -> [fastq::Record; 6] {
        [
            fastq::Record::with_attrs(
                &"id-01",
                None,
                b"AAAAAAAAAAAAAAATTTAA",
                b"&#3-G27C:(@G7B55+C4I"
            ),
            fastq::Record::with_attrs(
                &"id-01",
                None,
                b"TTTTTTTTTTTTTTTTTTTT",
                b"77%'24)FAF9@=94'%054"
            ),
            fastq::Record::with_attrs(
                &"id-01",
                None,
                b"AAAAAAAAAAAAAAATTTTA",
                b"'8$-BF2!C;+59->H@91#"
            ),
            fastq::Record::with_attrs(
                &"id-01",
                None,
                b"AAAAAAAAAAAAAAAAAAAA",
                b"%,42$CH*#0+0C6=0,*6/"
            ),
            fastq::Record::with_attrs(
                &"id-01",
                None,
                b"AAAAAAAAAAAAAAAAAAAT",
                b"-------------------J"
            ),
            fastq::Record::with_attrs(
                &"id-01",
                None,
                b"TAAAAAAAAAAAAAAAAAAA",
                b"I-------------------"
            ),
        ]
    }
    
    #[test]
    fn highest_quality_read_test() {
        let cases: [(f64, Option<(usize, usize)>); 6] = [
            (0.01, Some((10, 16))), // cutoff: Q20)
            (0.19952623149688797, Some((0, 20))), // cutoff: Q7
            (0.03162277660168379, Some((11, 19))), // cutoff: Q15
            (0.0001, None), // cutoff: Q40
            (0.0001, Some((19,20))), // cutoff: Q40
            (0.0001, Some((0,1)))  // cutoff: Q40
        ];

        let reads = get_reads();

        for ((cutoff, expected), read) in cases.iter().zip(reads) {
            let trim_strategy = HighestQualityTrimStrategy::new(*cutoff);
            assert_eq!(trim_strategy.trim(&read) , *expected);
        }
    }

    #[test]
    fn trim_by_quality_strategy_test() {
        let cases: [(u8, Option<(usize, usize)>); 6] = [
            ( 20, Some((4, 20))),
            ( 7, Some((0, 20))),
            (15, Some((1, 19))),
            (40, None),
            (40, Some((19,20))),
            (40, Some((0,1)))
        ];

        let reads = get_reads();

        for ((cutoff, expected), read) in cases.iter().zip(reads) {
            let trim_strategy = TrimByQualityStrategy::new(*cutoff);
            assert_eq!(trim_strategy.trim(&read) , *expected);
        }
    }

    #[test]
    fn fixed_crop_strategy_test() {
        let cases: [(usize, usize, Option<(usize, usize)>); 6] = [
            (5, 3, Some((5, 17))),
            (4, 4, Some((4, 16))),
            (1, 1, Some((1, 19))),
            (15, 15, None),
            (19, 0, Some((19,20))),
            (0, 19, Some((0,1)))
        ];

        let reads = get_reads();

        for ((head_crop, tail_crop, expected), read) in cases.iter().zip(reads) {
            let trim_strategy = FixedCropStrategy::new(*head_crop, *tail_crop);
            assert_eq!(trim_strategy.trim(&read) , *expected);
        }
    }
}
