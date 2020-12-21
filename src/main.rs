// based on https://docs.rs/bio/0.32.0/bio/io/fastq/index.html#read-and-write

extern crate bio;
extern crate clap;

use clap::{Arg, App};
use std::io;
use bio::io::fastq;
use bio::io::fastq::FastqRead;


fn main() {
    let config = get_args();
    filter(config.minqual, config.minlen, config.maxlen, config.headcrop, config.tailcrop)
}

struct Config {
    minqual: f64,
    minlen: usize,
    maxlen: usize,
    headcrop: usize,
    tailcrop: usize,
}

impl Config {
    // Create a config object from matches returned by get_args: clap arguments
    // These values are checked by is_int to be valid integers
    // and therefore calling unwrap() should be okay.
    fn new(matches: clap::ArgMatches) -> Config {
        let minqual: f64 = matches.value_of("quality").unwrap().parse().unwrap();
        let minlen: usize = matches.value_of("minlength").unwrap().parse().unwrap();
        let maxlen: usize = matches.value_of("maxlength").unwrap().parse().unwrap();
        let headcrop: usize = matches.value_of("headcrop").unwrap().parse().unwrap();
        let tailcrop: usize = matches.value_of("headcrop").unwrap().parse().unwrap();
        Config { minqual, minlen, maxlen, headcrop, tailcrop }
    }
}

fn get_args() -> Config {
    let matches = App::new("nanofilt")
                      .version("0.1")
                      .author("Wouter De Coster")
                      .about("Filtering and trimming of fastq files. Reads on stdin and writes to stdout.")
                      .after_help("EXAMPLE:\n\tgunzip -c reads.fastq.gz | nanofilt -q 10 -l 500 | gzip > filtered_reads.fastq.gz")
                      .arg(Arg::with_name("quality")
                           .short("q")
                           .long("quality")
                           .help("Sets a minimum Phred average quality score")
                           .takes_value(true)
                           .default_value("0")
                           .validator(is_int))
                       .arg(Arg::with_name("minlength")
                            .short("l")
                            .long("minlength")
                            .help("Sets a minimum read length")
                            .takes_value(true)
                            .default_value("1")
                            .validator(is_int))
                       .arg(Arg::with_name("maxlength")
                            .long("maxlength")
                            .help("Sets a maximum read length")
                            .takes_value(true)
                            // Default is largest i32. Better would be to explicitly use Inf, but couldn't figure it out.
                            .default_value("2147483647")
                            .validator(is_int))
                       .arg(Arg::with_name("headcrop")
                            .long("headcrop")
                            .help("Trim N nucleotides from the start of a read")
                            .takes_value(true)
                            .default_value("0")
                            .validator(is_int))
                       .arg(Arg::with_name("tailcrop")
                            .long("tailcrop")
                            .help("Trim N nucleotides from the end of a read")
                            .takes_value(true)
                            .default_value("0")
                            .validator(is_int))
                      .get_matches();
    Config::new(matches)
}

// Function to check if the supplied (string) argument on the command line are actually integers.
fn is_int(v: String) -> Result<(), String> {
    match v.parse::<i32>() {
        Ok(_i) => Ok(()),
        Err(_e) => Err(String::from("The value should be a positive integer!"))
    }
}
/// This function filters fastq on stdin based on quality, maxlength and minlength
/// and applies trimming before writting to stdout
fn filter(minqual: f64, minlen: usize, maxlen: usize, headcrop: usize, tailcrop:usize) {
    let mut reader = fastq::Reader::new(io::stdin());
    let mut record = fastq::Record::new();

    while let Ok(()) = reader.read(&mut record) {
        if record.is_empty() {
            break;
        }
        let read_len = record.seq().len();
        // If a read is shorter than what is to be cropped the read is dropped entirely (filtered out)
        if headcrop + tailcrop > read_len { continue }
        let average_quality = ave_qual(record.qual());
        if average_quality >= minqual && read_len >= minlen && read_len <= maxlen {
            // Check if a description attribute is present, taken from the bio-rust code to format fastq
            let header = match record.desc() {
                Some(d) => format!("{} {}", record.id(), d),
                None => record.id().to_owned(),
            };
            // Print out the records passing the filters, applying trimming on seq and qual
            // Could consider to use unsafe `from_utf8_unchecked`
            println!("@{}\n{}\n+\n{}",
                header,
                std::str::from_utf8(&record.seq()[headcrop..read_len-tailcrop]).unwrap(),
                std::str::from_utf8(&record.qual()[headcrop..read_len-tailcrop]).unwrap());
        }
    }

}

/// This function calculates the average quality of a read, and does this correctly
/// First the Phred scores are converted to probabilities (10^(q-33)/-10) and summed
/// and then divided by the number of bases/scores and converted to Phred again -10*log10(average)
fn ave_qual(quals: &[u8]) -> f64 {
    let probability_sum = quals.iter().map(|q| 10_f64.powf((*q as f64 - 33.0) / -10.0)).sum::<f64>();
    (probability_sum / quals.len() as f64).log10() * -10.0
}

// FEATURES TO ADD
// Write test for ave_qual
// write integration tests
// package
