use bio::io::fastq;
use clap::{ Parser, ValueEnum};
use minimap2::*;
use rayon::prelude::*;
use std::error::Error;
use std::io::Read;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use trimmers::*;
use utils::file_reader;

mod trimmers;
mod utils;

fn parse_usize_or_inf(s: &str) -> Result<usize, String> {
    if s.to_uppercase() == "INF" {
        Ok(usize::MAX)
    } else {
        s.parse().map_err(|e| format!("Invalid number: {e}"))
    }
}

fn parse_gc_value(s: &str) -> Result<f64, String> {
    let f: f64 = s.parse()
        .map_err(|e| format!("Invalid number: {e}"))?;
    if 0.0 <= f && f <= 1.0 {
        Ok(f)
    } else {
        Err(format!("Value out of range. Expected a value from 0 to 1."))
    }
}

// The arguments end up in the Cli struct
#[derive(Parser, Debug)]
#[clap(author, version, about="Filtering and trimming of fastq files. Reads on stdin and writes to stdout.", long_about = None)]
struct Cli {
    /// Sets a minimum Phred average quality score
    #[arg(short = 'q', long = "quality", value_parser, default_value_t = 0.0, help_heading = "Filtering Options")]
    minqual: f64,

    /// Sets a maximum Phred average quality score
    #[arg(long, value_parser, default_value_t = 1000.0, help_heading = "Filtering Options")]
    maxqual: f64,

    /// Sets a minimum read length
    #[arg(short = 'l', long, value_parser, default_value_t = 1, help_heading = "Filtering Options")]
    minlength: usize,

    /// Sets a maximum read length
    #[arg(long, value_parser = parse_usize_or_inf, default_value = "INF", help_heading = "Filtering Options")]
    maxlength: usize,

    /// Filter min GC content
    #[arg(long, value_parser = parse_gc_value, help_heading = "Filtering Options")]
    mingc: Option<f64>,

    /// Filter max GC content
    #[arg(long, value_parser = parse_gc_value, help_heading = "Filtering Options")]
    maxgc: Option<f64>,

    /// Filter contaminants against a fasta
    #[arg(short, long, value_parser, help_heading = "Filtering Options")]
    contam: Option<String>,

    /// Select the trimming strategy to apply to the reads.
    #[arg(long="trim-approach", value_parser, help_heading = "Trimming Options")]
    trim_approach: Option<TrimApproach>,

    /// Set the minimum quality score (Q-score) threshold for trimming low-quality bases from read ends.
    /// Required when using the `trim-by-quality` or `best-subread` trimming approaches.
    #[arg(long, value_parser, help_heading = "Trimming Options")]
    cutoff: Option<u8>,

    /// Trim N bases from the start of each read.
    /// Required only when using the `fixed-crop` trimming approach.
    #[arg(long, value_parser, default_value_t = 0, help_heading = "Trimming Options")]
    headcrop: usize,

    /// Trim N bases from the end of each read.
    /// Required only when using the `fixed-crop` trimming approach.
    #[arg(long, value_parser, default_value_t = 0, help_heading = "Trimming Options")]
    tailcrop: usize,

    /// Use N parallel threads
    #[arg(short, long, value_parser, default_value_t = 4, help_heading = "Setup Options")]
    threads: usize,

    /// Input filename [default: read from stdin]
    #[arg(short = 'i', long = "input", value_parser, help_heading = "Setup Options")]
    input: Option<String>,

    /// Output the opposite of the normal results
    #[arg(long, help_heading = "Setup Options")]
    inverse: bool,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum TrimApproach {
    /// Remove a fixed number of bases from both ends of the read.
    /// Requires setting both --headcrop and --tailcrop.
    FixedCrop,
    /// Trim low-quality bases from the ends of the read until reaching
    /// a base with quality ≥ --cutoff.
    TrimByQuality,
    /// Extract the highest-quality subread based on --cutoff, trimming
    /// low-quality bases from both ends.
    BestSubread
}

fn is_file(pathname: &str) -> Result<(), String> {
    let path = PathBuf::from(pathname);
    if path.is_file() {
        Ok(())
    } else {
        Err(format!("Input file {} is invalid", path.display()))
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Cli::parse();
    
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .expect("Error: Unable to build threadpool");

    let mut reader = file_reader(args.input.as_ref())?;
    filter(&mut reader, args);

    Ok(())
}

/// Validates that all required arguments are set for the selected trimming approach
/// and return a `Option<Arc<TrimStrategy>>`.
/// Exits the program with an error message if required arguments are missing.
fn build_trimming_approach(args: &Cli) -> Option<Arc<dyn TrimStrategy>> {
    if let Some(trim_approach) = args.trim_approach {
        match trim_approach {
            TrimApproach::FixedCrop => {
                if args.headcrop == 0 && args.tailcrop == 0 {
                    eprintln!(
                        "Error: When using the 'fixed-crop' trimming approach, at least one of --headcrop or --tailcrop must be greater than 0."
                    );
                    std::process::exit(1);
                }
                Some(Arc::new(FixedCropStrategy::new(args.headcrop, args.tailcrop)))
            },
            TrimApproach::TrimByQuality => {
                if let Some(cutoff)= args.cutoff {
                    Some(Arc::new(TrimByQualitytrategy::new(cutoff)))
                } else {
                    eprintln!(
                        "Error: When using the 'trim-by-quality' trimming approache, the --cutoff parameter must be set."
                    );
                    std::process::exit(1);
                }

            },
            TrimApproach::BestSubread => {
                if let Some(cutoff) = args.cutoff {
                    Some(Arc::new(HighestQualityTrimStrategy::new(phred_score_to_probability(cutoff + 33))))
                } else {
                    eprintln!(
                        "Error: When using the 'best-subread' trimming approache, the --cutoff parameter must be set."
                    );
                    std::process::exit(1);
                }
            }
        }
    } else {
        None
    }
}

/// This function filters fastq on stdin based on quality, maxlength and minlength
/// and applies trimming before writting to stdout
fn filter<T>(input: &mut T, args: Cli)
where
    T: Read + std::marker::Send,
{
    let aligner_option = if let Some(ref fas) = args.contam {
        is_file(fas)
            .unwrap_or_else(|_| panic!("Fasta file for filtering contaminants is invalid",));

        Some(setup_contamination_filter(fas, &args.threads))
    } else { None };

    let trimmer_strategy = build_trimming_approach(&args);
    
    let total_reads_ = Arc::new(AtomicUsize::new(0));
    let output_reads_ = Arc::new(AtomicUsize::new(0));

    fastq::Reader::new(input)
        .records()
        .par_bridge()
        .for_each(|record| {
            let record = record.expect("ERROR: problem parsing fastq record");
            total_reads_.fetch_add(1, Ordering::SeqCst);

            if record.is_empty() {
                return;
            }

            let valid_qual = is_valid_quality(&record, args.minqual, args.maxqual);
            let valid_len = is_valid_length(&record, args.minlength, args.maxlength);

            // If a GC content filter is set, validate the GC content; otherwise, assume it is valid.
            let valid_gc_p = if args.mingc.is_some() || args.maxgc.is_some() {
                is_valid_gc_percent(&record, args.mingc.unwrap_or(0.0), args.maxgc.unwrap_or(1.0))
            } else { true };

            // If a contaminants filter is set, validate the reads; otherwise, assume it is valid
            let is_not_contam = if let Some(ref aligner) = aligner_option {
                !is_contamination(&record.seq(), aligner)
            } else { true };

            if (valid_gc_p && valid_len && valid_qual && is_not_contam) ^ args.inverse {
                let pos_option = if let Some(ref trimmer) = trimmer_strategy {
                    trimmer.trim(&record)
                } else {
                    Some((0, record.seq().len()))
                };

                if let Some((start, end)) = pos_option {
                    write_record(&record, start, end);
                    output_reads_.fetch_add(1, Ordering::SeqCst);
                }
            }
        });

    let output_reads = output_reads_.load(Ordering::SeqCst);
    let total_reads = total_reads_.load(Ordering::SeqCst);
    eprintln!("Kept {output_reads} reads out of {total_reads} reads");
}

/// Write a record to stdout
fn write_record(record: &fastq::Record, start_pos: usize, end_pos: usize) {
    // Use a single formatted string with one allocation for the header
    let header = match record.desc() {
        Some(d) => format!("@{} {}", record.id(), d),
        None => format!("@{}", record.id()),
    };
    
    // Apply the trimming to both sequence and quality data
    let seq_slice = &record.seq()[start_pos..end_pos];
    let qual_slice = &record.qual()[start_pos..end_pos];
    
    // Use a single print to minimize syscalls
    println!(
        "{}\n{}\n+\n{}",
        header,
        unsafe { std::str::from_utf8_unchecked(seq_slice) },
        unsafe { std::str::from_utf8_unchecked(qual_slice) }
    );
}

/// This function calculates the average quality of a read, and does this correctly
/// First the Phred scores are converted to probabilities using `phred_quality_to_probability` and summed
/// and then divided by the number of bases/scores and converted to Phred again -10*log10(average)
fn ave_qual(quals: &[u8]) -> f64 {
    let probability_sum = quals.iter()
        .fold(0.0, |sum, &q| sum + phred_score_to_probability(q));
    
    (probability_sum / quals.len() as f64).log10() * -10.0
}

/// This function convert a Phred score to a probability
/// using the formula: 10^(-q/10)
fn phred_score_to_probability(phred: u8) -> f64 {
    10_f64.powf(((phred -33) as f64) / -10.0)
}

fn setup_contamination_filter(contam_fasta: &str, threads: &usize) -> Arc<Aligner<Built>> {
    Arc::new(Aligner::builder()
        .preset(Preset::LrHq)
        .with_index_threads(*threads)
        .with_cigar()
        .with_index(contam_fasta, None)
        .expect("Unable to build index"))
}

fn is_valid_quality(record: &fastq::Record, min_qual: f64, max_qual: f64) -> bool {
    let avg_qual = ave_qual(record.qual());
    min_qual <= avg_qual && avg_qual <= max_qual
}

fn is_valid_length(record: &fastq::Record, min_len: usize, max_len: usize) -> bool {
    let len = record.seq().len();
    min_len <= len && len <= max_len
}

fn is_valid_gc_percent(record: &fastq::Record, min_gc_p: f64, max_gc_p: f64) -> bool {
    let gc_percent = cal_gc(&record.seq());
    min_gc_p <= gc_percent && gc_percent <= max_gc_p
}


// Checks if a sequence is a contaminant, and returns true if so
// A sequence is considered a contaminant if there is an alignment of at least 90% between query and  target
fn is_contamination(readseq: &&[u8], contam: &Arc<Aligner<Built>>) -> bool {
    let alignment = contam
        .map(readseq, false, false, None, None, None)
        .expect("Unable to align");
    let query_len = readseq.len() as f64;
    if alignment.is_empty() {
        false
    } else {
        alignment
            .iter()
            .any(|a| (a.match_len as f64) >= (query_len * 0.9))
    }
}

fn cal_gc(readseq: &[u8]) -> f64 {
    let gc_count = readseq.iter().fold(0, |count, &base| {
        count + match base {
            b'G' | b'g' | b'C' | b'c' => 1,
            _ => 0,
        }
    });
    
    (gc_count as f64) / (readseq.len() as f64)
}

#[test]
fn test_ave_qual() {
    // Original test values need to be adjusted by adding 33 to each value
    assert_eq!(ave_qual(&[10+33]), 10.0);
    assert!((ave_qual(&[10+33, 11+33, 12+33]) - 10.923583702678473) < 0.00001);
    assert!((ave_qual(&[10+33, 11+33, 12+33, 20+33, 30+33, 40+33, 50+33]) - 14.408827647036087) < 0.00001);
    assert!(
        (ave_qual(&[
            17+33, 19+33, 11+33, 5+33, 3+33, 19+33, 22+33, 24+33, 20+33, 22+33, 30+33, 31+33, 32+33, 20+33, 21+33, 30+33, 28+33, 10+33, 13+33, 12+33, 18+33, 18+33,
            18+33, 19+33, 24+33, 25+33, 35+33, 33+33, 34+33, 35+33, 34+33, 27+33, 29+33, 25+33, 21+33, 18+33, 19+33, 12+33, 14+33, 15+33, 24+33, 26+33, 24+33, 7+33,
            12+33, 17+33, 17+33, 19+33, 17+33, 8+33, 14+33, 15+33, 13+33, 15+33, 9+33, 3+33, 4+33, 23+33, 23+33, 29+33, 23+33, 10+33, 29+33, 30+33, 31+33, 27+33, 25+33,
            14+33, 2+33, 13+33, 19+33, 14+33, 13+33, 13+33, 3+33, 2+33, 10+33, 17+33, 19+33, 25+33, 27+33, 20+33, 19+33, 11+33, 5+33, 7+33, 8+33, 8+33, 5+33, 2+33, 10+33,
            12+33, 16+33, 18+33, 16+33, 14+33, 12+33, 15+33, 2+33, 3+33, 11+33, 10+33, 15+33, 17+33, 17+33, 16+33, 13+33, 18+33, 26+33, 26+33, 23+33, 25+33, 23+33,
            18+33, 16+33, 33+33, 30+33, 26+33, 26+33, 21+33, 23+33, 8+33, 8+33, 11+33, 11+33, 6+33, 14+33, 19+33, 22+33, 20+33, 20+33, 18+33, 17+33, 20+33, 23+33,
            24+33, 28+33, 28+33, 28+33, 21+33, 20+33, 25+33, 27+33, 37+33, 28+33, 36+33, 29+33, 24+33, 27+33, 16+33, 18+33, 12+33, 8+33, 5+33, 3+33, 4+33, 6+33, 5+33,
            4+33, 4+33, 2+33, 10+33, 12+33, 6+33, 9+33, 9+33, 15+33, 16+33, 11+33, 10+33, 8+33, 8+33, 4+33, 3+33, 5+33, 4+33, 6+33, 15+33, 10+33, 9+33, 8+33, 7+33, 12+33, 4+33,
            5+33, 11+33, 12+33, 17+33, 13+33, 11+33, 17+33, 16+33, 4+33, 4+33, 5+33, 5+33, 12+33, 18+33, 17+33, 21+33
        ]) - 10.017407548271677)
            < 0.00001
    )
}

#[test]
fn test_filter() {
    filter(
        &mut std::fs::File::open("test-data/test.fastq").unwrap(),
        Cli {
            minlength: 100,
            maxlength: 100000,
            minqual: 5.0,
            maxqual: 200.0,
            trim_approach: Some(TrimApproach::FixedCrop),
            cutoff: None,
            headcrop: 10,
            tailcrop: 10,
            threads: 1,
            contam: None,
            inverse: false,
            input: None,
            mingc: Some(0.0),
            maxgc: Some(1.0),
        },
    );
}

#[test]
fn test_filter_with_trim_by_quality_approach() {
    filter(
        &mut std::fs::File::open("test-data/test.fastq").unwrap(),
        Cli {
            minlength: 100,
            maxlength: 100000,
            minqual: 5.0,
            maxqual: 200.0,
            trim_approach: Some(TrimApproach::TrimByQuality),
            cutoff: Some(10),
            headcrop: 0,
            tailcrop: 0,
            threads: 1,
            contam: None,
            inverse: false,
            input: None,
            mingc: Some(0.0),
            maxgc: Some(1.0),
        },
    );
}

#[test]
fn test_filter_with_best_subread_approach() {
    filter(
        &mut std::fs::File::open("test-data/test.fastq").unwrap(),
        Cli {
            minlength: 100,
            maxlength: 100000,
            minqual: 5.0,
            maxqual: 200.0,
            trim_approach: Some(TrimApproach::BestSubread),
            cutoff: Some(10),
            headcrop: 0,
            tailcrop: 0,
            threads: 1,
            contam: None,
            inverse: false,
            input: None,
            mingc: Some(0.0),
            maxgc: Some(1.0),
        },
    );
}

#[test]
fn test_contam() {
    let t: usize = 8;
    let aligner = setup_contamination_filter("test-data/random_contam.fa", &t);
    let rec = fastq::Reader::new(std::fs::File::open("test-data/test.fastq").unwrap())
        .records()
        .next()
        .unwrap()
        .unwrap();
    assert!(is_contamination(&rec.seq(), &aligner));
}

#[test]
fn test_no_contam() {
    let t: usize = 8;
    let aligner = setup_contamination_filter("test-data/random_contam.fa", &t);
    let rec = fastq::Reader::new(std::fs::File::open("test-data/other-test.fastq").unwrap())
        .records()
        .next()
        .unwrap()
        .unwrap();
    assert!(!is_contamination(&rec.seq(), &aligner));
}

#[test]
fn test_filter_with_contam() {
    filter(
        &mut std::fs::File::open("test-data/test.fastq").unwrap(),
        Cli {
            minlength: 100,
            maxlength: 100000,
            minqual: 5.0,
            maxqual: 100.0,
            trim_approach: Some(TrimApproach::FixedCrop),
            cutoff: None,
            headcrop: 10,
            tailcrop: 10,
            threads: 1,
            contam: Some("test-data/random_contam.fa".to_owned()),
            inverse: false,
            input: None,
            mingc: Some(0.0),
            maxgc: Some(1.0),
        },
    );
}

#[test]
fn test_record_qual_len() {
    fastq::Reader::new(std::fs::File::open("test-data/test.fastq").unwrap())
        .records()
        .for_each(|record| {
            let record = record.unwrap();
            if !record.is_empty() {
                let read_len = record.seq().len();
                let quals = record.qual();
                assert_eq!(
                    read_len,
                    quals.len(),
                    "length read doesn't equal length qual"
                );
            }
        })
}

#[test]
fn test_quals() {
    let rec = fastq::Reader::new(std::fs::File::open("test-data/test.fastq").unwrap())
        .records()
        .next()
        .unwrap()
        .unwrap();
    let quals = &rec.qual()[0..100]
        .iter()
        .map(|i| i - 33)
        .collect::<Vec<u8>>();
    assert_eq!(
        quals,
        &vec![
            17, 19, 11, 5, 3, 19, 22, 24, 20, 22, 30, 31, 32, 20, 21, 30, 28, 10, 13, 12, 18, 18,
            18, 19, 24, 25, 35, 33, 34, 35, 34, 27, 29, 25, 21, 18, 19, 12, 14, 15, 24, 26, 24, 7,
            12, 17, 17, 19, 17, 8, 14, 15, 13, 15, 9, 3, 4, 23, 23, 29, 23, 10, 29, 30, 31, 27, 25,
            14, 2, 13, 19, 14, 13, 13, 3, 2, 10, 17, 19, 25, 27, 20, 19, 11, 5, 7, 8, 8, 5, 2, 10,
            12, 16, 18, 16, 14, 12, 15, 2, 3
        ],
        "quals not as expected!"
    )
}

#[test]
fn phred_score_to_probability_test() {
    let cases: [(u8, f64); 4] = [
        (20 + 33, 0.01), // Q20
        (30 + 33, 0.001), // Q30
        (15 + 33, 0.03162277660168379), // Q15
        (25 + 33, 0.0031622776601683794), // Q25
    ];

    for (phred, prob) in cases {
        assert_eq!(phred_score_to_probability(phred), prob);
    }
}
