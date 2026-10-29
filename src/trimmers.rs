use super::phred_score_to_probability;
use crate::fastq;

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

/// Removes a poly(T) head and a poly(A) tail from the outer ends of a read.
///
/// The tails are found as cutadapt's `--poly-a` does: scanning inward from the
/// end of the read, each matching base scores +1 and any other base -2, and the
/// tail extends to wherever the score peaks, provided at most 20% of it are
/// other bases. The scoring tolerates the odd sequencing error inside a tail,
/// and the -2 keeps it from running on into the transcript. The tail must reach
/// the end of the read: adapters and primers have to be removed first.
pub struct PolyTailTrimmer {
    /// Shorter runs are left alone: a read ends in a few A's by chance.
    min_length: usize,
}

/// At most one in five bases of a tail may be something else.
const POLY_TAIL_MAX_ERROR_DIVISOR: usize = 5;

impl PolyTailTrimmer {
    pub fn new(min_length: usize) -> Self {
        PolyTailTrimmer { min_length }
    }

    /// Narrows the outer ends of a read's kept segments: a poly(T) head off
    /// the start of the first segment and a poly(A) tail off the end of the
    /// last. Inner segment ends (from splitting) are not touched, as they are
    /// not read ends. Segments with nothing left are dropped.
    pub fn trim_outer_ends(&self, seq: &[u8], segments: &mut Vec<(usize, usize)>) {
        if let Some(first) = segments.first_mut() {
            first.0 += self.run_length(seq[first.0..first.1].iter(), b"TtUu");
        }
        if let Some(last) = segments.last_mut() {
            last.1 -= self.run_length(seq[last.0..last.1].iter().rev(), b"Aa");
        }
        segments.retain(|&(start, end)| start < end);
    }

    /// Length of the run of `bases` at the start of `seq`, allowing for errors.
    fn run_length<'a>(&self, seq: impl Iterator<Item = &'a u8>, bases: &[u8]) -> usize {
        let (mut score, mut best_score, mut best_length, mut errors) = (0i64, 0i64, 0, 0);
        for (i, base) in seq.enumerate() {
            if bases.contains(base) {
                score += 1;
            } else {
                score -= 2;
                errors += 1;
            }
            let length = i + 1;
            if score > best_score && errors * POLY_TAIL_MAX_ERROR_DIVISOR <= length {
                best_score = score;
                best_length = length;
            }
        }
        if best_length >= self.min_length {
            best_length
        } else {
            0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fastq;

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
                read.header()
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

    fn poly_trim(seq: &[u8], min_length: usize) -> Vec<(usize, usize)> {
        let mut segments = vec![(0, seq.len())];
        PolyTailTrimmer::new(min_length).trim_outer_ends(seq, &mut segments);
        segments
    }

    #[test]
    fn exact_poly_tails_are_removed_from_both_ends() {
        let seq = b"TTTTTTTTTTTTGCATGCATCGAAAAAAAAAAAA";
        assert_eq!(poly_trim(seq, 10), vec![(12, 22)]);
        // U counts as T, lowercase as uppercase.
        assert_eq!(
            poly_trim(b"uuUUUUUUUUUUGCAUGCAUCGaaAAAAAAAAAA", 10),
            vec![(12, 22)]
        );
    }

    #[test]
    fn errors_inside_a_tail_are_tolerated() {
        // A substitution, an insertion and a deleted stretch inside each tail.
        let seq = b"TTTTGTTTTTTTCTTTTGCATGCATCGAAAAAAAAGAAAAACAAAAAA";
        assert_eq!(poly_trim(seq, 10), vec![(17, 27)]);
    }

    #[test]
    fn short_runs_are_left_alone() {
        let seq = b"TTTTGCATGCATCGAAAAAAAAA"; // 4 T's and 9 A's
        assert_eq!(poly_trim(seq, 10), vec![(0, seq.len())]);
        assert_eq!(poly_trim(seq, 4), vec![(4, 14)]);
    }

    #[test]
    fn a_tail_does_not_run_on_into_a_rich_sequence() {
        // The -2 for other bases stops the tail at the transcript, even when
        // that is A-rich; at most the A's directly before the tail go with it.
        let seq = b"GCATGCACAGACATAGCAAAAAAAAAAAAAAA";
        assert_eq!(poly_trim(seq, 10), vec![(0, 17)]);
    }

    #[test]
    fn a_read_that_is_all_tail_is_dropped() {
        assert_eq!(poly_trim(b"AAAAAAAAAAAAAAAAAAAA", 10), vec![]);
        assert_eq!(poly_trim(b"TTTTTTTTTTTTAAAAAAAAAAAA", 10), vec![]);
    }

    #[test]
    fn only_the_outer_ends_of_split_segments_are_trimmed() {
        let seq = b"TTTTTTTTTTTTGCATGCAAAAAAAAAAAATTTTTTTTTTTTGCATGCAAAAAAAAAAAA";
        let mut segments = vec![(0, 30), (30, seq.len())];
        PolyTailTrimmer::new(10).trim_outer_ends(seq, &mut segments);
        // The A's ending the first segment and the T's starting the second stay.
        assert_eq!(segments, vec![(12, 30), (30, 48)]);
    }

    /// Deterministic xorshift, so the simulations below need no extra crate.
    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn base(&mut self) -> u8 {
            b"ACGT"[(self.next() % 4) as usize]
        }
        fn chance(&mut self) -> f64 {
            (self.next() % 1_000_000) as f64 / 1_000_000.0
        }
        /// ONT-like errors: half deletions, a quarter each substitutions and insertions.
        fn noisy(&mut self, seq: &[u8], error_rate: f64) -> Vec<u8> {
            let mut out = Vec::with_capacity(seq.len());
            for &b in seq {
                let r = self.chance();
                if r < error_rate * 0.5 {
                    continue;
                } else if r < error_rate * 0.75 {
                    let mut other = self.base();
                    while other == b {
                        other = self.base();
                    }
                    out.push(other);
                } else if r < error_rate {
                    out.push(b);
                    out.push(self.base());
                } else {
                    out.push(b);
                }
            }
            out
        }
    }

    #[test]
    fn simulated_noisy_tails_are_removed_almost_entirely() {
        let mut rng = Rng(42);
        for error_rate in [0.02, 0.05, 0.10, 0.15] {
            let (mut tail_left, mut body_lost, n) = (0, 0, 2000);
            for _ in 0..n {
                let body: Vec<u8> = (0..300).map(|_| rng.base()).collect();
                let head = rng.noisy(&[b'T'; 30], error_rate);
                let body = rng.noisy(&body, error_rate);
                let tail = rng.noisy(&[b'A'; 60], error_rate);
                let seq = [&head[..], &body[..], &tail[..]].concat();
                let (start, end) = (head.len(), head.len() + body.len());

                let trimmed = poly_trim(&seq, 10);
                assert_eq!(trimmed.len(), 1, "a read with a body is kept");
                let (s, e) = trimmed[0];
                tail_left += start.saturating_sub(s) + e.saturating_sub(end);
                body_lost += s.saturating_sub(start) + end.saturating_sub(e);
            }
            // Tail bases left, and body bases removed along with a tail, per read.
            // Body bases are lost when the body happens to end in A or start
            // with T: those cannot be told apart from the tail.
            let left = tail_left as f64 / n as f64;
            let lost = body_lost as f64 / n as f64;
            assert!(
                left < 1.5,
                "error rate {error_rate}: {left} tail bases left per read"
            );
            assert!(
                lost < 1.5,
                "error rate {error_rate}: {lost} body bases lost per read"
            );
        }
    }

    #[test]
    fn reads_without_tails_are_left_alone() {
        let mut rng = Rng(7);
        let mut trimmed = 0;
        for _ in 0..5000 {
            let seq: Vec<u8> = (0..300).map(|_| rng.base()).collect();
            if poly_trim(&seq, 10) != vec![(0, seq.len())] {
                trimmed += 1;
            }
        }
        // A run of ten-ish A's or T's at a read end is rare by chance.
        assert!(
            trimmed <= 5,
            "{trimmed} of 5000 tailless reads were trimmed"
        );
    }
}
