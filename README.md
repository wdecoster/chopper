# chopper

Rust implementation of [NanoFilt](https://github.com/wdecoster/nanofilt)+[NanoLyse](https://github.com/wdecoster/nanolyse), both originally written in Python. This tool, intended for long read sequencing such as PacBio or ONT, filters and trims a fastq file.  

Filtering is based on average read quality, minimum or maximum read length, and GC content percentage.

On the other hand, trimming is performed using one of three approaches:

* Fixed cropping (head crop at the start of the read and tail crop at the end),

* Quality-based trimming using a threshold,

* Extracting the highest-quality sub-read based on a quality score.

Reads that pass the filters are printed to standard output (STDOUT).

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
  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version

Filtering Options:
  -q, --quality <MINQUAL>
          Sets a minimum Phred average quality score
          
          [default: 0]

      --maxqual <MAXQUAL>
          Sets a maximum Phred average quality score
          
          [default: 1000]

  -l, --minlength <MINLENGTH>
          Sets a minimum read length
          
          [default: 1]

      --maxlength <MAXLENGTH>
          Sets a maximum read length
          
          [default: INF]

      --mingc <MINGC>
          Filter min GC content

      --maxgc <MAXGC>
          Filter max GC content

  -c, --contam <CONTAM>
          Filter contaminants against a fasta

Trimming Options:
      --trim-approach <TRIM_APPROACH>
          Select the trimming strategy to apply to the reads

          Possible values:
          - fixed-crop:      Remove a fixed number of bases from both ends of the read. Requires setting both --headcrop and --tailcrop
          - trim-by-quality: Trim low-quality bases from the ends of the read until reaching a base with quality ≥ --cutoff
          - best-read-segment:    Extract the highest-quality read segment based on --cutoff, trimming low-quality bases from both ends

      --cutoff <CUTOFF>
          Set the minimum quality score (Q-score) threshold for trimming low-quality bases from read ends. Required when using the `trim-by-quality` or `best-read-segment` trimming approaches

      --headcrop <HEADCROP>
          Trim N bases from the start of each read. Required only when using the `fixed-crop` trimming approach
          
          [default: 0]

      --tailcrop <TAILCROP>
          Trim N bases from the end of each read. Required only when using the `fixed-crop` trimming approach
          
          [default: 0]

Setup Options:
  -t, --threads <THREADS>
          Use N parallel threads
          
          [default: 4]

  -i, --input <INPUT>
          Input filename [default: read from stdin]

      --inverse
          Output the opposite of the normal results
```

## Examples:

```bash
gunzip -c reads.fastq.gz | chopper -q 10 -l 500 | gzip > filtered_reads.fastq.gz
chopper -q 10 -l 500 -i reads.fastq > filtered_reads.fastq
chopper -q 10 -l 500 -i reads.fastq.gz | gzip > filtered_reads.fastq.gz
```

## Citation

If you use this tool, please consider citing our [publication](https://academic.oup.com/bioinformatics/article/39/5/btad311/7160911).
