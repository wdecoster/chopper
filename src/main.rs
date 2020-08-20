// based on https://docs.rs/bio/0.32.0/bio/io/fastq/index.html#read-and-write

extern crate bio;
use std::io;
use bio::io::fastq;
use bio::io::fastq::FastqRead;


fn main() {
    let mut reader = fastq::Reader::new(io::stdin());
    let mut writer = fastq::Writer::new(io::stdout());
    let mut record = fastq::Record::new();

    while let Ok(()) = reader.read(&mut record) {
        if record.is_empty() {
            break;
        }
        let average_quality = ave_qual(record.qual());
        let read_len = record.seq().len();

        if average_quality >= minqual && read_len >= minlen && read_len <= maxlen {
            writer.write_record(&record).unwrap_or_else(|error| {
                panic!("Problem writing to stdout: {:?}", error);
            })
        }
    }

}

fn ave_qual(quals: &[u8]) -> f64 {
    let probability_sum = quals.iter().map(|q| 10_f64.powf((*q as f64 - 33.0) / -10.0)).sum::<f64>();
    (probability_sum / quals.len() as f64).log10() * -10.0
}

