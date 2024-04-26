# chopper

Rust implementation of [NanoFilt](https://github.com/wdecoster/nanofilt)+[NanoLyse](https://github.com/wdecoster/nanolyse), both originally written in Python. This tool, intended for long read sequencing such as PacBio or ONT, filters and trims a fastq file.  
Filtering is done on average read quality and minimal or maximal read length, and applying a headcrop (start of read) and tailcrop (end of read) while printing the reads passing the filter.

Compared to the Python implementation the scope is to deliver the same results, almost the same functionality, at much faster execution times. At the moment this tool does not support filtering using a sequencing_summary file. If those features are of interest then please reach out.  

## Installation

Preferably, for most users, download a ready-to-use binary for your system to add directory on your $PATH from the [releases](https://github.com/wdecoster/chopper/releases).  
You may have to change the file permissions to execute it with `chmod +x chopper`

Alternatively, use conda to install  
`conda install -c bioconda chopper`

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
        --contam        Fasta file with reference to check potential contaminants against [default None]
    -i, --input         Input filename [default: read from stdin]
        --maxgc         Sets a maximum GC content [default: 1.0]
        --mingc         Sets a minimum GC content [default: 0.0]
```

EXAMPLES:

```bash
gunzip -c reads.fastq.gz | chopper -q 10 -l 500 | gzip > filtered_reads.fastq.gz
chopper -q 10 -l 500 -i reads.fastq > filtered_reads.fastq
chopper -q 10 -l 500 -i reads.fastq.gz | gzip > filtered_reads.fastq.gz
```

Note that the tool may be substantially slower in the third example above, and piping while decompressing is recommended (as in the first example). 

## CITATION

If you use this tool, please consider citing our [publication](https://academic.oup.com/bioinformatics/article/39/5/btad311/7160911).
