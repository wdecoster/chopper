# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
make ci                 # fmt-check + clippy + test, what CI runs
make all                # fmt, clippy, test, build
cargo test              # 72 unit tests + 19 end-to-end tests in tests/formats.rs
cargo test modtags      # one module
cargo test untrimmed_range_is_unchanged   # one test by name substring
cargo test -- --ignored # the 4 #[ignore]d end-to-end tests (see below)
make musl               # static linux binary, as shipped in releases
```

Clippy is run as `cargo clippy --all-targets --all-features -- -D warnings`; warnings fail the build.

The four `#[ignore]`d tests in `src/main.rs` call `filter()` directly, which writes to real stdout and asserts nothing. They are smoke tests only — run them explicitly, and do not add assertions expecting captured output.

## Architecture

### Everything is a `(start, end)` range over the original record

The single most important invariant: a read is never copied or mutated while being filtered or trimmed. Filters inspect `&fastq::Record`, trimmers return `Vec<(usize, usize)>` of half-open ranges into that record, and only `record_to_string` / `record_to_bam` (`src/records.rs`) ever slice. Multiple ranges mean the read was split into segments.

This is why a new trimming approach only needs a `TrimStrategy` impl, and why `--update-mods` covers all four approaches for free: everything downstream already speaks in offsets into the untrimmed sequence.

### Pipeline

`main` → `filter` → dispatches on `--threads` to `sequential_filter` or `parallel_filter` → `get_valid_segments` per record → `WritableRecord::new` per segment → writer.

`get_valid_segments` (`src/main.rs`) is the heart of it, and its ordering matters in two ways:

- Filters run **cheapest-first with early bail-out** (length → quality → GC → contamination, the last being a minimap2 alignment). The early returns are only sound in normal mode; with `--inverse` a read that fails a filter may still be kept, so every filter must be evaluated and the decision deferred to the final XOR. Preserve that distinction when touching this function.
- Filters apply to the **whole untrimmed read**. Trimming happens afterwards, and only `--minlength` is re-checked per segment. Mean quality is never recomputed per segment, so a read passing `-q 10` overall yields segments that are kept regardless of their own quality. This surprises users (see issue #69) — it is current intended behaviour, not a bug to fix incidentally.

### Parallelism and output order

`parallel_filter` uses rayon `par_bridge` over the record iterator plus a pool of crossbeam channels: `for_each_init` pops one `Sender` per worker, and a dedicated writer thread multiplexes the receivers with `Select`.

**Output order is nondeterministic whenever `--threads > 1.`** Any test or before/after comparison of output bytes must either use `-t 1` or compare order-insensitively.

### FASTQ parsing (`src/fastq.rs`)

chopper has its own FASTQ reader; it replaced the `bio` crate, which brought in about 60 crates for this one use. On well-formed input it must yield exactly what `bio` did (trailing whitespace trimmed from every line, CRLF accepted), and the v0.14.1 digest test guards that. On malformed input it is deliberately stricter or more lenient: mismatched sequence/quality lengths are an error (bio accepted them, then trimming panicked), wrapped records pair qualities by length, a one-line sequence takes exactly one quality line, empty reads and blank lines are accepted, and errors carry line numbers and exit 1 rather than panicking.

The header is kept whole (`Record::header()`). SAM tags from `samtools fastq -T MM,ML,MN` follow the read name as tab-separated fields, and a tag value may contain spaces (`CO:Z:some comment`), so `split_header` (`src/records.rs`) splits on the first tab for the tags and then on the first space for the name. The `_segment_N` suffix goes on the read name, never at the end of the line: writing it there silently corrupts the last tag.

### Base modification tags (`src/modtags.rs`)

Rewrites `MM`/`ML`/`MN` for a subsequence, behind `--update-mods`. Invariants worth knowing before editing, each checked against htslib's `sam_mods.c`:

- FASTQ is always in original read orientation, so the canonical base is counted left-to-right and the `+`/`-` strand character is metadata carried through unchanged. (In BAM this is driven by `BAM_FREVERSE`, which has no FASTQ equivalent.)
- The canonical base `N` matches **every** base, not just literal `N` characters (`freq[15] = l_qseq` in htslib).
- `ML` is optional; `MM` alone is valid.
- A list with multiple codes (`C+mh`) has an `ML` stride > 1, and `ML` values across lists are concatenated in list order.

Design contract: tags that cannot be rewritten are **fatal**, not a warning — passing them through would attach modification coordinates to a sequence they no longer describe, which is the bug the flag exists to prevent. Do not soften this to a warning without discussion.

Changes here should be validated against htslib (via pysam) rather than against chopper's own output — decode the original and the trimmed read independently and check the surviving calls land at the expected positions with the expected probabilities.

### Output writing

Use `write_all`, never `write`. Records at or above the `BufWriter`'s 8 KiB capacity bypass the buffer and reach stdout in a single `write` call that may consume only part of it — that is any read over roughly 4 kb, routine here. The fatal paths (`abort_on_record_error`, `abort_on_input_error`, `abort_on_write_error`) use `eprintln!` + `process::exit`, not `Result` propagation; note that exiting from a rayon worker discards whatever is still buffered, so the output of a failed run is deliberately incomplete.

### Input

`src/utils.rs` sniffs gzip/bzip2/xz by magic bytes (falling back to file extension) and returns a boxed `BufRead`; stdin is peeked for gzip only. No `-i` means stdin, and a terminal on stdin is a fatal error.

### Unaligned BAM (`src/formats.rs`)

BAM is BGZF, so it arrives through the gzip decoder and is recognised by the decompressed `BAM\1` magic. Each BAM record becomes an `InputRecord`: a `fastq::Record` view (phred+33 qualities, tags as SAM text in the id after a tab) that filters and trimmers use unchanged, plus the original `RecordBuf`. Which tags go into the view depends on the output (`TagText`): all of them for FASTQ, only `MM`/`ML`/`MN`/`Mm`/`Ml` for BAM with `--update-mods`, none otherwise. BAM output clones the original record so every other tag keeps its type and order, and the input header is passed through with a chained `@PG`. FASTQ output must stay byte-identical to earlier releases; `tests/formats.rs` checks it against `test-data/regression/v0.14.1.sha256`. Validate BAM changes against htslib (pysam); `test-data/ubam/make_fixtures.py` regenerates the fixtures.

## Releases

1. Bump `version` in `Cargo.toml` **and** refresh `Cargo.lock` — `publish.yml` builds with `--locked` and fails otherwise.
2. Push to master, wait for CI.
3. `git tag -a vX.Y.Z` with the **release notes as the tag body**. The `release` job in `.github/workflows/publish.yml` reads `%(contents:body)` from the annotated tag and creates the GitHub release with it; the build matrix then uploads linux, linux-musl and macos binaries into that release. A lightweight tag falls back to GitHub's generated notes.
4. Push the tag to trigger the build.

Releases must not be bare — put the effort into the tag annotation rather than editing the release afterwards.
