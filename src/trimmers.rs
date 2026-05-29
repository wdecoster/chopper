use super::phred_score_to_probability;
use bio::io::fastq;

/// A trait for implementing custom read trimming strategies.
///
/// Implementing this trait allows you to define how a read should be trimmed  
/// according to a specific algorithm or criteria. This is useful when working  
/// with different trimming approaches.
///
/// The `trim` method should return a vector of `(start, end)` tuples representing
/// valid read segments, or an empty vector if the read should be discarded.
pub trait TrimStrategy: Send + Sync {
    /// Performs the trimming process based on the defined strategy.
    ///
    /// Returns a vector of `(start, end)` tuples with the indices of valid read segments.
    /// Returns an empty vector if the read should be discarded.
    fn trim(&self, record: &fastq::Record) -> Vec<(usize, usize)>;
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
    pub fn new(cutoff: f64) -> HighestQualityTrimStrategy {
        HighestQualityTrimStrategy { cutoff }
    }
}

impl TrimStrategy for HighestQualityTrimStrategy {
    fn trim(&self, record: &fastq::Record) -> Vec<(usize, usize)> {
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

            if best_cumulative_error < current_cumulative_error
                || (best_cumulative_error == current_cumulative_error
                    && best_length < i - current_start + 1)
            {
                best_start = current_start;
                best_end = i;
                best_cumulative_error = current_cumulative_error;
                best_length = i - current_start + 1;
            }
        }

        if best_start == usize::MAX {
            vec![]
        } else {
            vec![(best_start, best_end + 1)]
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
    fn trim(&self, record: &fastq::Record) -> Vec<(usize, usize)> {
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
            return vec![];
        }

        vec![(trim_start, trim_end)]
    }
}

/// This strategy removes a fixed number of bases from both ends of the read.
pub struct FixedCropStrategy {
    head_crop: usize,
    tail_crop: usize,
}

impl FixedCropStrategy {
    pub fn new(head_crop: usize, tail_crop: usize) -> FixedCropStrategy {
        FixedCropStrategy {
            head_crop,
            tail_crop,
        }
    }
}

impl TrimStrategy for FixedCropStrategy {
    fn trim(&self, record: &fastq::Record) -> Vec<(usize, usize)> {
        if record.seq().len().saturating_sub(self.tail_crop) <= self.head_crop {
            return vec![];
        }
        vec![(self.head_crop, record.seq().len() - self.tail_crop)]
    }
}

/// This strategy splits reads by low-quality segments and outputs the high-quality
/// parts on the left and right of low-quality regions, provided they pass the length filter.
///
/// The strategy treats only *runs* of consecutive bases below the quality threshold
/// that are at least `window` long as separators. Shorter dips below the threshold are
/// tolerated and kept inside the surrounding segment. Leading and trailing low-quality
/// bases are trimmed off, and each high-quality segment must reach `min_length` to be
/// returned as a separate trimmed read.
pub struct SplitByLowQualityStrategy {
    cutoff: u8,
    min_length: usize,
    window: usize,
}

impl SplitByLowQualityStrategy {
    pub fn new(cutoff: u8, min_length: usize, window: usize) -> SplitByLowQualityStrategy {
        // A window of 0 would never split; treat it as 1 (split on any bad base).
        let window = window.max(1);
        SplitByLowQualityStrategy {
            cutoff,
            min_length,
            window,
        }
    }
}

impl TrimStrategy for SplitByLowQualityStrategy {
    fn trim(&self, record: &fastq::Record) -> Vec<(usize, usize)> {
        let quals = record.qual();
        let mut segments = Vec::new();

        // `segment_start` is the first high-quality base of the current segment,
        // `last_good` is the most recent high-quality base seen, and `bad_run`
        // counts consecutive low-quality bases since the last good base.
        let mut segment_start: Option<usize> = None;
        let mut last_good: Option<usize> = None;
        let mut bad_run: usize = 0;

        let push_segment = |start: usize, end: usize, segments: &mut Vec<(usize, usize)>| {
            if end - start >= self.min_length {
                segments.push((start, end));
            }
        };

        for (i, &q) in quals.iter().enumerate() {
            let quality = q - 33; // Convert from ASCII to actual Phred score

            if quality >= self.cutoff {
                // High quality base: (re)open or extend the current segment. Any
                // pending short run of bad bases is tolerated and stays inside it.
                if segment_start.is_none() {
                    segment_start = Some(i);
                }
                last_good = Some(i);
                bad_run = 0;
            } else {
                // Low quality base: only split once the run reaches the window size.
                bad_run += 1;
                if bad_run >= self.window {
                    if let (Some(start), Some(lg)) = (segment_start, last_good) {
                        push_segment(start, lg + 1, &mut segments);
                    }
                    segment_start = None;
                    last_good = None;
                }
            }
        }

        // Handle the final segment, trimming any trailing tolerated bad bases.
        if let (Some(start), Some(lg)) = (segment_start, last_good) {
            push_segment(start, lg + 1, &mut segments);
        }

        segments
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::fastq;

    fn get_reads() -> [fastq::Record; 6] {
        [
            fastq::Record::with_attrs(
                "id-01",
                None,
                b"AAAAAAAAAAAAAAATTTAA",
                b"&#3-G27C:(@G7B55+C4I",
            ),
            fastq::Record::with_attrs(
                "id-01",
                None,
                b"TTTTTTTTTTTTTTTTTTTT",
                b"77%'24)FAF9@=94'%054",
            ),
            fastq::Record::with_attrs(
                "id-01",
                None,
                b"AAAAAAAAAAAAAAATTTTA",
                b"'8$-BF2!C;+59->H@91#",
            ),
            fastq::Record::with_attrs(
                "id-01",
                None,
                b"AAAAAAAAAAAAAAAAAAAA",
                b"%,42$CH*#0+0C6=0,*6/",
            ),
            fastq::Record::with_attrs(
                "id-01",
                None,
                b"AAAAAAAAAAAAAAAAAAAT",
                b"-------------------J",
            ),
            fastq::Record::with_attrs(
                "id-01",
                None,
                b"TAAAAAAAAAAAAAAAAAAA",
                b"I-------------------",
            ),
        ]
    }

    #[test]
    fn highest_quality_read_test() {
        let cases: [(f64, Vec<(usize, usize)>); 6] = [
            (0.01, vec![(10, 16)]),                // cutoff: Q20)
            (0.19952623149688797, vec![(0, 20)]),  // cutoff: Q7
            (0.03162277660168379, vec![(11, 19)]), // cutoff: Q15
            (0.0001, vec![]),                      // cutoff: Q40
            (0.0001, vec![(19, 20)]),              // cutoff: Q40
            (0.0001, vec![(0, 1)]),                // cutoff: Q40
        ];

        let reads = get_reads();

        for ((cutoff, expected), read) in cases.iter().zip(reads) {
            let trim_strategy = HighestQualityTrimStrategy::new(*cutoff);
            assert_eq!(trim_strategy.trim(&read), *expected);
        }
    }

    #[test]
    fn trim_by_quality_strategy_test() {
        let cases: [(u8, Vec<(usize, usize)>); 6] = [
            (20, vec![(4, 20)]),
            (7, vec![(0, 20)]),
            (15, vec![(1, 19)]),
            (40, vec![]),
            (40, vec![(19, 20)]),
            (40, vec![(0, 1)]),
        ];

        let reads = get_reads();

        for ((cutoff, expected), read) in cases.iter().zip(reads) {
            let trim_strategy = TrimByQualityStrategy::new(*cutoff);
            assert_eq!(trim_strategy.trim(&read), *expected);
        }
    }

    #[test]
    fn fixed_crop_strategy_test() {
        let cases = [
            (5usize, 3usize, vec![(5usize, 17usize)]),
            (30, 4, vec![]),
            (1, 1, vec![(1, 19)]),
            (15, 30, vec![]),
            (19, 0, vec![(19, 20)]),
            (0, 19, vec![(0, 1)]),
        ];

        let reads = get_reads();

        for ((head_crop, tail_crop, expected), read) in cases.iter().zip(reads) {
            let trim_strategy = FixedCropStrategy::new(*head_crop, *tail_crop);
            assert_eq!(trim_strategy.trim(&read), *expected);
        }
    }

    #[test]
    fn split_by_low_quality_strategy_test() {
        // Test cases: (cutoff, min_length, expected_segments)
        let cases = [
            // Read 1: Cutoff Q20, min_length 3
            (20u8, 3usize, vec![(6usize, 9usize), (10, 16)]),
            // Read 2: Cutoff Q7, min_length 3
            (7, 3, vec![(4, 15), (17, 20)]),
            // Read 3: Cutoff Q15, min_length 3
            (15, 3, vec![(4, 7), (14, 19)]),
            // Read 4: Cutoff Q40, min_length 3
            (40, 3, vec![]),
            // Read 4: Same read with Q40 but min_length 1
            (40, 1, vec![(19, 20)]),
            // Read 5: Same with Q40 but min_length 1
            (40, 1, vec![(0, 1)]),
            // Read 1: Test with Q10 threshold, min_length 5
            (10, 5, vec![(2, 9), (10, 20)]),
            // Read 2: Test with Q10 threshold, min_length 3
            (10, 3, vec![(7, 15), (17, 20)]),
        ];

        let reads = get_reads();

        for ((cutoff, min_length, expected), read) in cases.iter().zip(reads) {
            let trim_strategy = SplitByLowQualityStrategy::new(*cutoff, *min_length, 1);
            let result = trim_strategy.trim(&read);
            assert_eq!(
                result,
                *expected,
                "Failed for cutoff={}, min_length={}, read_id={}",
                cutoff,
                min_length,
                read.id()
            );
        }
    }

    #[test]
    fn split_by_low_quality_multiple_segments_test() {
        // Create a read with distinct high and low quality regions
        // Pattern: HHH LLL HHH LLL HHH (H=high quality, L=low quality)
        let test_read = fastq::Record::with_attrs(
            "test_split",
            None,
            b"AAATTTAAATTTAAATTTAAA", // 21 bases
            b"III###III###III###III", // Q40 and Q2 alternating
        );

        let trim_strategy = SplitByLowQualityStrategy::new(10, 2, 1); // Q10 cutoff, min 2 bases
        let result = trim_strategy.trim(&test_read);

        // Should get multiple segments: (0,3), (6,9), (12,15), (18,21)
        let expected = vec![(0, 3), (6, 9), (12, 15), (18, 21)];
        assert_eq!(result, expected);
    }

    #[test]
    fn split_by_low_quality_edge_cases_test() {
        // Test empty read
        let empty_read = fastq::Record::with_attrs("empty", None, b"", b"");
        let trim_strategy = SplitByLowQualityStrategy::new(10, 1, 1);
        assert_eq!(trim_strategy.trim(&empty_read), vec![]);

        // Test read with all low quality
        let low_qual_read = fastq::Record::with_attrs(
            "low_qual", None, b"AAAA", b"####", // All Q2
        );
        assert_eq!(trim_strategy.trim(&low_qual_read), vec![]);

        // Test read with all high quality
        let high_qual_read = fastq::Record::with_attrs(
            "high_qual",
            None,
            b"AAAA",
            b"IIII", // All Q40
        );
        assert_eq!(trim_strategy.trim(&high_qual_read), vec![(0, 4)]);
    }

    #[test]
    fn split_by_low_quality_window_test() {
        // Two high-quality blocks separated by a single low-quality base,
        // then by a run of three low-quality bases.
        // Pattern (Q40=I high, Q2=# low): III # IIII ### III  -> 14 bases
        let read =
            fastq::Record::with_attrs("window_test", None, b"AAAAAAAAAAAAAA", b"III#IIII###III");

        // Window 1: any single bad base splits -> three segments.
        let w1 = SplitByLowQualityStrategy::new(10, 1, 1);
        assert_eq!(w1.trim(&read), vec![(0, 3), (4, 8), (11, 14)]);

        // Window 3: the single bad base at index 3 is tolerated, so the first
        // two blocks merge; only the run of three bad bases splits the read.
        let w3 = SplitByLowQualityStrategy::new(10, 1, 3);
        assert_eq!(w3.trim(&read), vec![(0, 8), (11, 14)]);

        // Window 4: no bad run reaches length 4, so the whole read stays together
        // (trailing/leading bad bases trimmed; here it spans the good bases).
        let w4 = SplitByLowQualityStrategy::new(10, 1, 4);
        assert_eq!(w4.trim(&read), vec![(0, 14)]);
    }
}
