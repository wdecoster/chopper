// based on https://docs.rs/bio/0.32.0/bio/io/fastq/index.html#read-and-write
use bio::io::fastq;
use clap::Parser;
use minimap2::*;
use rayon::prelude::*;
use std::io::{self, Read};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

// The arguments end up in the Cli struct
#[derive(Parser, Debug)]
#[clap(author, version, about="Filtering and trimming of fastq files. Reads on stdin and writes to stdout.", long_about = None)]
struct Cli {
    /// Sets a minimum Phred average quality score
    #[arg(short = 'q', long = "quality", value_parser, default_value_t = 0.0)]
    minqual: f64,

    /// Sets a maximum Phred average quality score
    #[arg(long, value_parser, default_value_t = 1000.0)]
    maxqual: f64,

    /// Sets a minimum read length
    #[arg(short = 'l', long, value_parser, default_value_t = 1)]
    minlength: usize,

    /// Sets a maximum read length
    // Default is largest i32. Better would be to explicitly use Inf, but couldn't figure it out.
    #[arg(long, value_parser, default_value_t = 2147483647)]
    maxlength: usize,

    /// Trim N nucleotides from the start of a read
    #[arg(long, value_parser, default_value_t = 0)]
    headcrop: usize,

    /// Trim N nucleotides from the end of a read
    #[arg(long, value_parser, default_value_t = 0)]
    tailcrop: usize,

    /// Use N parallel threads
    #[arg(short, long, value_parser, default_value_t = 4)]
    threads: usize,

    /// Filter contaminants against a fasta
    #[arg(short, long, value_parser)]
    contam: Option<String>,
}

fn is_file(pathname: &str) -> Result<(), String> {
    let path = PathBuf::from(pathname);
    if path.is_file() {
        Ok(())
    } else {
        Err(format!("Input file {} is invalid", path.display()))
    }
}

fn main() {
    let args = Cli::parse();
    filter(&mut io::stdin(), args);
}

/// This function filters fastq on stdin based on quality, maxlength and minlength
/// and applies trimming before writting to stdout
fn filter<T>(input: &mut T, args: Cli)
where
    T: Read + std::marker::Send,
{
    match args.contam {
        Some(ref fas) => {
            is_file(fas)
                .unwrap_or_else(|_| panic!("Fasta file for filtering contaminants is invalid",));
            let mut total_reads = 0;
            let mut output_reads = 0;
            let aligner = setup_contamination_filter(fas);
            fastq::Reader::new(input).records().for_each(|record| {
                let record = record.unwrap();
                total_reads += 1;
                if !record.is_empty() {
                    let read_len = record.seq().len();
                    // If a read is shorter than what is to be cropped the read is dropped entirely (filtered out)
                    if args.headcrop + args.tailcrop < read_len {
                        let average_quality =
                            ave_qual(&record.qual().iter().map(|i| i - 33).collect::<Vec<u8>>());
                        if average_quality >= args.minqual
                            && average_quality <= args.maxqual
                            && read_len >= args.minlength
                            && read_len <= args.maxlength
                            && !is_contamination(&record.seq(), &aligner)
                        {
                            write_record(record, &args, read_len);
                            output_reads += 1;
                        }
                    }
                }
            });
            eprintln!("Kept {output_reads} reads out of {total_reads} reads");
        }

        None => {
            let total_reads_ = Arc::new(AtomicUsize::new(0));
            let output_reads_ = Arc::new(AtomicUsize::new(0));
            rayon::ThreadPoolBuilder::new()
                .num_threads(args.threads)
                .build_global()
                .unwrap();
            fastq::Reader::new(input)
                .records()
                .par_bridge()
                .for_each(|record| {
                    let record = record.unwrap();
                    total_reads_.fetch_add(1, Ordering::SeqCst);
                    if !record.is_empty() {
                        let read_len = record.seq().len();
                        // If a read is shorter than what is to be cropped the read is dropped entirely (filtered out)
                        if args.headcrop + args.tailcrop < read_len {
                            let average_quality = ave_qual(
                                &record.qual().iter().map(|i| i - 33).collect::<Vec<u8>>(),
                            );
                            if average_quality >= args.minqual
                                && average_quality <= args.maxqual
                                && read_len >= args.minlength
                                && read_len <= args.maxlength
                            {
                                write_record(record, &args, read_len);
                                output_reads_.fetch_add(1, Ordering::SeqCst);
                            }
                        }
                    }
                });
            let output_reads = output_reads_.load(Ordering::SeqCst);
            let total_reads = total_reads_.load(Ordering::SeqCst);
            eprintln!("Kept {output_reads} reads out of {total_reads} reads");
        }
    }
}

fn write_record(record: fastq::Record, args: &Cli, read_len: usize) {
    // Check if a description attribute is present, taken from the bio-rust code to format fastq
    let header = match record.desc() {
        Some(d) => format!("{} {}", record.id(), d),
        None => record.id().to_owned(),
    };
    // Print out the records passing the filters, applying trimming on seq and qual
    // Could consider to use unsafe `from_utf8_unchecked`
    println!(
        "@{}\n{}\n+\n{}",
        header,
        std::str::from_utf8(&record.seq()[args.headcrop..read_len - args.tailcrop]).unwrap(),
        std::str::from_utf8(&record.qual()[args.headcrop..read_len - args.tailcrop]).unwrap()
    );
}

/// This function calculates the average quality of a read, and does this correctly
/// First the Phred scores are converted to probabilities (10^(q)/-10) and summed
/// and then divided by the number of bases/scores and converted to Phred again -10*log10(average)
fn ave_qual(quals: &[u8]) -> f64 {
    let probability_sum = quals
        .iter()
        .map(|q| 10_f64.powf((*q as f64) / -10.0))
        .sum::<f64>();
    (probability_sum / quals.len() as f64).log10() * -10.0
}

fn setup_contamination_filter(contam_fasta: &str) -> Aligner {
    Aligner::builder()
        .with_threads(8)
        .map_ont()
        .with_index(contam_fasta, None)
        .expect("Unable to build index")
}

// Checks if a sequence is a contaminant, and returns false if so
// A sequence is considered a contaminant if there is an alignment of at least 90% between query and  target
fn is_contamination(readseq: &&[u8], contam: &Aligner) -> bool {
    let alignment = contam
        .map(readseq, false, false, None, None)
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

#[test]
fn test_ave_qual() {
    assert_eq!(ave_qual(&[10]), 10.0);
    assert!((ave_qual(&[10, 11, 12]) - 10.923583702678473) < 0.00001);
    assert!((ave_qual(&[10, 11, 12, 20, 30, 40, 50]) - 14.408827647036087) < 0.00001);
    assert!(
        (ave_qual(&[
            17, 19, 11, 5, 3, 19, 22, 24, 20, 22, 30, 31, 32, 20, 21, 30, 28, 10, 13, 12, 18, 18,
            18, 19, 24, 25, 35, 33, 34, 35, 34, 27, 29, 25, 21, 18, 19, 12, 14, 15, 24, 26, 24, 7,
            12, 17, 17, 19, 17, 8, 14, 15, 13, 15, 9, 3, 4, 23, 23, 29, 23, 10, 29, 30, 31, 27, 25,
            14, 2, 13, 19, 14, 13, 13, 3, 2, 10, 17, 19, 25, 27, 20, 19, 11, 5, 7, 8, 8, 5, 2, 10,
            12, 16, 18, 16, 14, 12, 15, 2, 3, 11, 10, 15, 17, 17, 16, 13, 18, 26, 26, 23, 25, 23,
            18, 16, 33, 30, 26, 26, 21, 23, 8, 8, 11, 11, 6, 14, 19, 22, 20, 20, 18, 17, 20, 23,
            24, 28, 28, 28, 21, 20, 25, 27, 37, 28, 36, 29, 24, 27, 16, 18, 12, 8, 5, 3, 4, 6, 5,
            4, 4, 2, 10, 12, 6, 9, 9, 15, 16, 11, 10, 8, 8, 4, 3, 5, 4, 6, 15, 10, 9, 8, 7, 12, 4,
            5, 11, 12, 17, 13, 11, 17, 16, 4, 4, 5, 5, 12, 18, 17, 21
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
            headcrop: 10,
            tailcrop: 10,
            threads: 1,
            contam: None,
        },
    );
}

#[test]
fn test_contam() {
    let aligner = setup_contamination_filter("test-data/random_contam.fa");
    let rec = fastq::Reader::new(std::fs::File::open("test-data/test.fastq").unwrap())
        .records()
        .next()
        .unwrap()
        .unwrap();
    assert!(is_contamination(&rec.seq(), &aligner));
}

#[test]
fn test_no_contam() {
    let aligner = setup_contamination_filter("test-data/random_contam.fa");
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
            headcrop: 10,
            tailcrop: 10,
            threads: 1,
            contam: Some("test-data/random_contam.fa".to_owned()),
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
