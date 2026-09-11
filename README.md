# chopper

Rust implementation of [NanoFilt](https://github.com/wdecoster/nanofilt)+[NanoLyse](https://github.com/wdecoster/nanolyse), both originally written in Python. This tool, intended for long read sequencing such as PacBio or ONT, filters and trims a fastq file.  

Filtering is based on average read quality, minimum or maximum read length, and GC content percentage.

On the other hand, trimming is performed using one of four approaches:

* Fixed cropping (head crop at the start of the read and tail crop at the end),

* Quality-based trimming using a threshold,

* Extracting the highest-quality sub-read based on a quality score,

* Splitting reads by low-quality segments and outputting high-quality parts.

Reads that pass the filters are printed to standard output (STDOUT).

Base modification tags (`MM`/`ML`) carried over from a modified-basecalled BAM can
be kept in step with the trimming with `--update-mods`, see
[Base modification tags](#base-modification-tags).

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
          - split-by-low-quality: Split reads by low-quality segments and output high-quality parts on the left and right, provided they pass the length filter

      --cutoff <CUTOFF>
          Set the minimum quality score (Q-score) threshold for trimming low-quality bases from read ends. Required when using the `trim-by-quality`, `best-read-segment`, or `split-by-low-quality` trimming approaches

      --split-window <SPLIT_WINDOW>
          Minimum number of consecutive bases below --cutoff required to split a read. Shorter low-quality runs are tolerated and kept inside the segment. Only used with the `split-by-low-quality` trimming approach
          
          [default: 1]

      --headcrop <HEADCROP>
          Trim N bases from the start of each read. Required only when using the `fixed-crop` trimming approach
          
          [default: 0]

      --tailcrop <TAILCROP>
          Trim N bases from the end of each read. Required only when using the `fixed-crop` trimming approach
          
          [default: 0]

      --update-mods
          Recompute the base modification tags (MM, ML and MN) written by `samtools fastq -T MM,ML,MN` so that they match the trimmed read. Only those tags are corrected; other position- or quality-dependent tags (qs, ns, ts, du, ...) are passed through unchanged and go stale. Reads whose tags are malformed or do not describe their sequence are reported as an error rather than written out with wrong coordinates

Setup Options:
  -t, --threads <THREADS>
          Use N parallel threads
          
          [default: 4]

  -i, --input <INPUT>
          Input filename [default: read from stdin]

      --inverse
          Output the opposite of the normal results
```

## Examples

```bash
gunzip -c reads.fastq.gz | chopper -q 10 -l 500 | gzip > filtered_reads.fastq.gz
chopper -q 10 -l 500 -i reads.fastq > filtered_reads.fastq
chopper -q 10 -l 500 -i reads.fastq.gz | gzip > filtered_reads.fastq.gz

# Trim low-quality bases from read ends
chopper --trim-approach trim-by-quality --cutoff 15 -i reads.fastq > trimmed_reads.fastq

# Extract the highest-quality segment from each read
chopper --trim-approach best-read-segment --cutoff 15 -l 100 -i reads.fastq > best_segments.fastq

# Split reads by low-quality segments and output high-quality parts
chopper --trim-approach split-by-low-quality --cutoff 15 -l 50 -i reads.fastq > split_reads.fastq

# Only split when at least 5 consecutive bases fall below the cutoff (tolerate shorter dips)
chopper --trim-approach split-by-low-quality --cutoff 15 --split-window 5 -l 50 -i reads.fastq > split_reads.fastq

# Keep base modification tags in step with the trimming
samtools fastq -T MM,ML,MN reads.bam \
  | chopper --trim-approach fixed-crop --headcrop 20 --tailcrop 20 --update-mods \
  > trimmed_reads.fastq
```

## Base modification tags

Modified basecalls from dorado are carried into fastq as SAM tags:

```bash
samtools fastq -T MM,ML,MN reads.bam > reads.fastq
```

`MM` stores the position of each modified base as a count of canonical bases to
skip, so trimming or splitting a read silently invalidates it: the coordinates
now point at the wrong bases, or past the end of the read entirely. Passing
`--update-mods` makes chopper recompute `MM`, `ML` and `MN` for whatever part of
the read it keeps, for every trimming approach. Calls that fall outside the kept
range are dropped along with their `ML` probabilities. It is off by default, and
without it the tags are written through untouched, as before.

**Only `MM`, `ML` and `MN` are corrected.** Other tags that trimming also
invalidates are passed through unchanged and will be stale, including the mean
quality `qs`, the signal offsets `ns` and `ts`, and the duration `du`. If those
matter downstream, recompute or drop them yourself.

Reads whose modification tags cannot be rewritten are a fatal error rather than
a warning, because writing them through would attach modification coordinates to
a sequence they no longer describe. This happens when the tags are malformed, when
`MM` and `ML` disagree on how many calls there are, or when the tags describe a
different sequence than the read carries (for example because the read was
already trimmed by another tool without updating them). Rerun without
`--update-mods` to write such reads through unchanged.

Note that `--update-mods` only fixes the tags. It does not make chopper aware of
modification calls when deciding what to trim.

## Performance

chopper's filtering is fast, and in practice the runtime is often dominated by
(de)compressing gzipped FASTQ rather than by the filtering itself. Since chopper
reads from stdin and writes to stdout, the recommended way to run it is to pipe
data through [`pigz`](https://zlib.net/pigz/):

```bash
pigz -dc reads.fastq.gz | chopper -q 10 -l 500 | pigz > filtered_reads.fastq.gz
```

This helps in two ways:

- **Output compression is genuinely parallelised.** `pigz` compresses the
  filtered output across multiple cores, which is usually much faster than a
  single-threaded `gzip`.
- **Decompression runs in its own process.** Note that gzip decompression is
  inherently single-threaded (a gzip stream must be inflated sequentially, so
  `pigz` cannot parallelise it either). The benefit of piping is that the
  decompression happens in a separate process, overlapping with chopper's
  filtering on other cores. chopper itself already uses the fast `zlib-ng`
  backend for decompression, so reading a `.gz` file directly with `-i` is also
  efficient.

## Citation

If you use this tool, please consider citing our [publication](https://academic.oup.com/bioinformatics/article/39/5/btad311/7160911).
