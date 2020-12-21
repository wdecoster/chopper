# nanofilt
Rust implementation of [NanoFilt](https://github.com/wdecoster/nanofilt), originally written in Python. This tool filters and trims a fastq file as was originally intended for long read sequencing such as PacBio or ONT.  
Filtering is done on average read quality and minimal or maximal read length, and applying a headcrop (start of read) and tailcrop (end of read) while printing the reads passing the filter.

Compared to the Python implementation the scope is to deliver the same results, almost the same functionality, at much faster execution times. At the moment this tool does not support filtering using a sequencing_summary file or filtering on GC content. If those features are of interest then please reach out.  
As this is my first rust project I welcome all feedback!

## Installation

## Usage
Reads on stdin and writes to stdout.

```
FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
        --headcrop <headcrop>      Trim N nucleotides from the start of a read [default: 0]
        --maxlength <maxlength>    Sets a maximum read length [default: 2147483647]
    -l, --minlength <minlength>    Sets a minimum read length [default: 1]
    -q, --quality <quality>        Sets a minimum Phred average quality score [default: 0]
        --tailcrop <tailcrop>      Trim N nucleotides from the end of a read [default: 0]
```

EXAMPLE:  
	`gunzip -c reads.fastq.gz | nanofilt -q 10 -l 500 | gzip > filtered_reads.fastq.gz`

## Citation
If you use this tool, please consider citing our [publication](https://academic.oup.com/bioinformatics/article/34/15/2666/4934939).