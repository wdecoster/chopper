// based on https://docs.rs/bio/0.32.0/bio/io/fastq/index.html#read-and-write

use std::io;
use bio::io::fastq;
use bio::io::fastq::FastqRead;


fn main() {
    let mut reader = fastq::Reader::new(io::stdin());
    let mut writer = fastq::Writer::new(io::stdout());
    let mut record = fastq::Record::new();

    while let Ok(()) = reader.read(&mut record) {
        if record.is_empty() {
            let check = record.check();
            break;
        }

        let mut sum_qual = record.qual().iter().sum::<u8>() as f64;

        if (sum_qual / record.seq().len() as f64 - 33.0) > 30.0 {
            writer.write_record(&record);
        }
    }

}
