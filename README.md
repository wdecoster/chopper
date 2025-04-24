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
Usage: chopper [OPTIONS]

Options:
  -q, --quality <MINQUAL>      Sets a minimum Phred average quality score [default: 0]
      --maxqual <MAXQUAL>      Sets a maximum Phred average quality score [default: 1000]
  -l, --minlength <MINLENGTH>  Sets a minimum read length [default: 1]
      --maxlength <MAXLENGTH>  Sets a maximum read length
      --headcrop <HEADCROP>    Trim N nucleotides from the start of a read [default: 0]
      --tailcrop <TAILCROP>    Trim N nucleotides from the end of a read [default: 0]
  -t, --threads <THREADS>      Use N parallel threads [default: 4]
  -c, --contam <CONTAM>        Filter contaminants against a fasta
      --inverse                Output the opposite of the normal results
  -i, --input <INPUT>          Input filename [default: read from stdin]
      --maxgc <MAXGC>          Filter max GC content [default: 1]
      --mingc <MINGC>          Filter min GC content [default: 0]
      --trim <TRIM>            Set a q-score cutoff to trim read ends
  -h, --help                   Print help
  -V, --version                Print version
```

EXAMPLES:

```bash
gunzip -c reads.fastq.gz | chopper -q 10 -l 500 | gzip > filtered_reads.fastq.gz
chopper -q 10 -l 500 -i reads.fastq > filtered_reads.fastq
chopper -q 10 -l 500 -i reads.fastq.gz | gzip > filtered_reads.fastq.gz
```

## CITATION

If you use this tool, please consider citing our [publication](https://academic.oup.com/bioinformatics/article/39/5/btad311/7160911).
