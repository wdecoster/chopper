# chopper

Rust implementation of [NanoFilt](https://github.com/wdecoster/nanofilt)+[NanoLyse](https://github.com/wdecoster/nanolyse), both originally written in Python. This tool,  intended for long read sequencing such as PacBio or ONT, filters and trims a fastq file.  
Filtering is done on average read quality and minimal or maximal read length, and applying a headcrop (start of read) and tailcrop (end of read) while printing the reads passing the filter.

Compared to the Python implementation the scope is to deliver the same results, almost the same functionality, at much faster execution times. At the moment this tool does not support filtering using a sequencing_summary file or filtering on GC content. If those features are of interest then please reach out.  
As this is my first rust project I welcome all feedback!

## Installation

## Usage

Reads on stdin and writes to stdout.

```text
FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
        --headcrop      Trim N nucleotides from the start of a read [default: 0]
        --maxlength     Sets a maximum read length [default: 2147483647]
    -l, --minlength     Sets a minimum read length [default: 1]
    -q, --quality       Sets a minimum Phred average quality score [default: 0]
        --tailcrop      Trim N nucleotides from the end of a read [default: 0]
        --threads       Number of parallel threads to use [default: 4]
```

EXAMPLE:  
 `gunzip -c reads.fastq.gz | chopper -q 10 -l 500 | gzip > filtered_reads.fastq.gz`

## Citation

If you use this tool, please consider citing our [publication](https://academic.oup.com/bioinformatics/article/34/15/2666/4934939).
