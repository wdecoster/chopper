//! Input and output formats: FASTQ and unaligned BAM.
//!
//! Filters and trimmers only ever see a `fastq::Record`. A BAM record is turned
//! into one on input and kept alongside it, so that BAM output starts from the
//! original record and its typed tags, flags and header survive the round trip.
//! Only the read name, sequence and qualities are replaced, plus the base
//! modification tags under `--update-mods`.

use std::fmt;
use std::io::{self, BufRead, Stdout, StdoutLock, Write};
use std::num::NonZeroUsize;

use bio::io::fastq;
use clap::ValueEnum;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::{
    self as sam,
    alignment::{
        io::Write as _,
        record::data::field::Tag,
        record_buf::{
            data::field::{value::Array, Value},
            Data,
        },
        RecordBuf,
    },
    header::record::value::{
        map::{self, program, Program},
        Map,
    },
};

use crate::records::WritableRecord;
use crate::utils::peek;

/// The decompressed contents of every BAM file start with these bytes.
const BAM_MAGIC: &[u8] = b"BAM\x01";

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum Format {
    Fastq,
    /// Unaligned BAM
    Bam,
}

/// A read as it comes off the input, before it is converted for the filters.
/// The conversion happens on the worker threads rather than in the reader.
pub enum RawRecord {
    Fastq(fastq::Record),
    Bam(RecordBuf),
}

/// One input read. `fastq` is what the filters and trimmers see; for BAM input
/// it is built from `bam`, which is kept for writing BAM output.
pub struct InputRecord {
    pub fastq: fastq::Record,
    pub bam: Option<RecordBuf>,
}

/// Which BAM tags are copied as SAM text into the header of the FASTQ view.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TagText {
    /// FASTQ output: every tag, as `samtools fastq -T '*'` would write them.
    All,
    /// BAM output with `--update-mods`: just the tags that get rewritten. The
    /// rest are carried over from the original BAM record untouched.
    ModsOnly,
    /// BAM output: no tags, the original record keeps them.
    None,
}

impl TagText {
    pub fn new(output: Format, update_mods: bool) -> Self {
        match (output, update_mods) {
            (Format::Fastq, _) => TagText::All,
            (Format::Bam, true) => TagText::ModsOnly,
            (Format::Bam, false) => TagText::None,
        }
    }
}

/// An opened input, with its format detected from its contents.
pub enum Input {
    Fastq(Box<dyn BufRead + Send>),
    Bam(bam::io::Reader<Box<dyn BufRead + Send>>, Box<sam::Header>),
}

impl Input {
    /// Detects the format of an already decompressed stream. BAM is recognised
    /// by its magic number; everything else is parsed as FASTQ, as before.
    pub fn detect(reader: Box<dyn BufRead + Send>) -> io::Result<Self> {
        let (magic, reader) = peek(reader, BAM_MAGIC.len())?;
        let reader: Box<dyn BufRead + Send> = Box::new(io::BufReader::new(reader));
        if magic == BAM_MAGIC {
            let mut reader = bam::io::Reader::from(reader);
            let header = reader.read_header()?;
            Ok(Input::Bam(reader, Box::new(header)))
        } else {
            Ok(Input::Fastq(reader))
        }
    }

    pub fn format(&self) -> Format {
        match self {
            Input::Fastq(_) => Format::Fastq,
            Input::Bam(..) => Format::Bam,
        }
    }

    pub fn header(&self) -> Option<&sam::Header> {
        match self {
            Input::Fastq(_) => None,
            Input::Bam(_, header) => Some(header),
        }
    }

    pub fn into_records(self) -> Box<dyn Iterator<Item = RawRecord> + Send> {
        match self {
            Input::Fastq(reader) => Box::new(fastq::Reader::new(reader).records().map(|record| {
                RawRecord::Fastq(record.expect("ERROR: problem parsing fastq record"))
            })),
            Input::Bam(mut reader, header) => Box::new(std::iter::from_fn(move || {
                let mut record = RecordBuf::default();
                match reader.read_record_buf(&header, &mut record) {
                    Ok(0) => None,
                    Ok(_) => Some(RawRecord::Bam(record)),
                    Err(e) => {
                        eprintln!("Error: problem parsing BAM record: {e}");
                        std::process::exit(1);
                    }
                }
            })),
        }
    }
}

impl InputRecord {
    pub fn new(raw: RawRecord, tag_text: TagText) -> Result<Self, String> {
        match raw {
            RawRecord::Fastq(fastq) => Ok(InputRecord { fastq, bam: None }),
            RawRecord::Bam(bam) => Ok(InputRecord {
                fastq: bam_to_fastq(&bam, tag_text)?,
                bam: Some(bam),
            }),
        }
    }

    /// The read name, for error messages.
    pub fn name(&self) -> &str {
        let id = self.fastq.id();
        id.split('\t').next().unwrap_or(id)
    }
}

/// Builds the FASTQ view of a BAM record that the filters and trimmers work on.
fn bam_to_fastq(record: &RecordBuf, tag_text: TagText) -> Result<fastq::Record, String> {
    let name = record
        .name()
        .map(|n| n.to_string())
        .unwrap_or_else(|| "*".to_string());
    // An aligned record may be reverse complemented and carries a CIGAR and
    // position that trimming would invalidate. Neither can be handled here.
    let flags = record.flags();
    if !flags.is_unmapped() {
        return Err(format!(
            "read {name} is aligned; only unaligned BAM is supported as input"
        ));
    }
    // A secondary or supplementary record repeats (part of) a read that is
    // also in the file as a primary record, and would be written out twice.
    if flags.is_secondary() || flags.is_supplementary() {
        return Err(format!(
            "read {name} is a secondary or supplementary record (flag {}); \
             only unaligned BAM is supported as input",
            u16::from(flags)
        ));
    }
    // A reverse-strand record is stored reverse complemented, and its MM tag
    // counts bases on the original strand, which trimming here cannot follow.
    if flags.is_reverse_complemented() {
        return Err(format!(
            "read {name} is flagged as reverse complemented (flag {}); \
             only unaligned BAM is supported as input",
            u16::from(flags)
        ));
    }
    // Mates are filtered and trimmed one by one, which would leave records
    // flagged as paired whose mate was dropped or renamed into segments.
    if flags.is_segmented() {
        return Err(format!(
            "read {name} is one of a read pair (flag {}); \
             paired reads are not supported, as filtering would separate mates",
            u16::from(flags)
        ));
    }
    let seq = record.sequence().as_ref();
    let qual = record.quality_scores().as_ref();
    if qual.len() != seq.len() {
        return Err(format!(
            "read {name} has no base qualities, which chopper needs for filtering and trimming"
        ));
    }
    let qual = qual
        .iter()
        .map(|&q| {
            if q <= 93 {
                Ok(q + 33)
            } else {
                Err(format!(
                    "read {name} has a base quality of {q}, above the maximum of 93"
                ))
            }
        })
        .collect::<Result<Vec<u8>, String>>()?;

    let mut id = name;
    let fields = match tag_text {
        TagText::All => tags_to_text(record.data(), |_| true)?,
        TagText::ModsOnly => tags_to_text(record.data(), is_mod_tag)?,
        TagText::None => vec![],
    };
    for field in fields {
        id.push('\t');
        id.push_str(&field);
    }
    Ok(fastq::Record::with_attrs(&id, None, seq, &qual))
}

/// The tags `--update-mods` rewrites. `Mm`/`Ml` are the pre-standard spelling
/// that `modtags` also accepts.
pub fn is_mod_tag(tag: &Tag) -> bool {
    matches!(tag.as_ref(), b"MM" | b"ML" | b"MN" | b"Mm" | b"Ml")
}

/// The fixed fields of an unmapped SAM record without a name, sequence or
/// qualities, so that a line of SAM tags can be parsed or written on its own.
const EMPTY_SAM_FIELDS: usize = 11;
const EMPTY_SAM_PREFIX: &str = "*\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*";

/// Formats the selected tags as SAM text fields (`TG:T:value`), in their
/// original order. This goes through noodles' SAM writer so that the text is
/// exactly what the SAM spec prescribes; floats are written in their shortest
/// form that reads back to the same value.
fn tags_to_text(data: &Data, keep: impl Fn(&Tag) -> bool) -> Result<Vec<String>, String> {
    let data: Data = data
        .iter()
        .filter(|(tag, _)| keep(tag))
        .map(|(tag, value)| (tag, value.clone()))
        .collect();
    if data.is_empty() {
        return Ok(vec![]);
    }
    let record = RecordBuf::builder().set_data(data).build();
    let mut line = Vec::new();
    sam::io::Writer::new(&mut line)
        .write_alignment_record(&sam::Header::default(), &record)
        .map_err(|e| format!("cannot write tags as SAM text: {e}"))?;
    let line = String::from_utf8(line).map_err(|e| e.to_string())?;
    Ok(line
        .trim_end_matches('\n')
        .split('\t')
        .skip(EMPTY_SAM_FIELDS)
        .map(str::to_string)
        .collect())
}

/// Parses SAM text tag fields (`TG:T:value`) into typed BAM tags, in order.
/// Integers get the smallest type that holds them, as htslib does when it
/// converts SAM.
pub fn text_to_tags<S: AsRef<str>>(fields: &[S]) -> Result<Data, String> {
    if fields.is_empty() {
        return Ok(Data::default());
    }
    let error = |e: &dyn fmt::Display| {
        let fields: Vec<&str> = fields.iter().map(AsRef::as_ref).collect();
        format!("'{}' are not valid SAM tags: {e}", fields.join("\t"))
    };
    // The SAM spec allows a B array without values (`ML:B:C`), which is what
    // `--update-mods` writes for a segment that kept no calls, but noodles
    // does not parse one. Those are built here; noodles parses the rest.
    let empty: Vec<Option<(Tag, Value)>> = fields.iter().map(|f| empty_array(f.as_ref())).collect();
    let mut line = String::from(EMPTY_SAM_PREFIX);
    for (field, empty) in fields.iter().zip(&empty) {
        if empty.is_none() {
            line.push('\t');
            line.push_str(field.as_ref());
        }
    }
    let parse = || -> io::Result<RecordBuf> {
        let record = sam::Record::try_from(line.as_bytes())?;
        RecordBuf::try_from_alignment_record(&sam::Header::default(), &record)
    };
    let parsed = parse().map_err(|e| error(&e))?;
    let mut parsed = parsed.data().iter();
    let data: Data = empty
        .into_iter()
        .map(|empty| match empty {
            Some(field) => Some(field),
            // noodles reads every SAM integer as int32; narrow it like htslib.
            None => parsed.next().map(|(tag, value)| match value {
                Value::Int32(n) => (tag, Value::try_from(i64::from(*n)).expect("fits in i32")),
                _ => (tag, value.clone()),
            }),
        })
        .collect::<Option<_>>()
        .ok_or_else(|| error(&"a field is missing"))?;
    // A repeated tag collapses into one entry, which would lose a value.
    if data.len() != fields.len() {
        return Err(error(&"a tag appears more than once"));
    }
    Ok(data)
}

/// Parses a B array field without values, e.g. `ML:B:C`.
fn empty_array(field: &str) -> Option<(Tag, Value)> {
    let b = field.as_bytes();
    if b.len() != 6 || &b[2..5] != b":B:" {
        return None;
    }
    let array = match b[5] {
        b'c' => Array::Int8(vec![]),
        b'C' => Array::UInt8(vec![]),
        b's' => Array::Int16(vec![]),
        b'S' => Array::UInt16(vec![]),
        b'i' => Array::Int32(vec![]),
        b'I' => Array::UInt32(vec![]),
        b'f' => Array::Float(vec![]),
        _ => return None,
    };
    Some((Tag::from([b[0], b[1]]), Value::Array(array)))
}

/// The header for BAM output: the input's header when it was BAM, a minimal
/// one otherwise, with a `@PG` line for this run appended to the chain.
///
/// `order_kept` says whether the records come out in input order with their
/// names unchanged. Otherwise a sort or grouping order the input declared no
/// longer holds, and the header says `SO:unsorted` instead.
pub fn output_header(input: Option<&sam::Header>, order_kept: bool) -> io::Result<sam::Header> {
    let mut header = match input {
        Some(header) => header.clone(),
        None => {
            let hd = Map::<map::Header>::builder()
                .insert(map::header::tag::SORT_ORDER, "unknown")
                .build()
                .map_err(io::Error::other)?;
            sam::Header::builder().set_header(hd).build()
        }
    };
    let command_line: Vec<String> = std::env::args().collect();
    let pg = Map::<Program>::builder()
        .insert(program::tag::NAME, "chopper")
        .insert(program::tag::VERSION, env!("CARGO_PKG_VERSION"))
        .insert(program::tag::COMMAND_LINE, command_line.join(" "))
        .build()
        .map_err(io::Error::other)?;
    add_program(header.programs_mut(), pg);

    if !order_kept {
        if let Some(hd) = header.header_mut() {
            let fields = hd.other_fields_mut();
            if fields.contains_key(&map::header::tag::SORT_ORDER) {
                fields.insert(map::header::tag::SORT_ORDER, "unsorted".into());
            }
            fields.shift_remove(&map::header::tag::GROUP_ORDER);
            fields.shift_remove(&map::header::tag::SUBSORT_ORDER);
        }
    }
    Ok(header)
}

/// Adds chopper's `@PG` line, chained after the last program(s) when the
/// existing chain allows it.
///
/// A `PP` that names a missing program (common after merging or reheadering)
/// makes noodles panic while it looks for the end of the chain, and a cyclic
/// chain has no end at all. In either case the line is added without a `PP`,
/// so a damaged header does not stop a run that does not depend on it.
fn add_program(programs: &mut sam::header::Programs, pg: Map<Program>) {
    let ids = programs.as_ref();
    let chain_intact = ids.values().all(|map| {
        map.other_fields()
            .get(&program::tag::PREVIOUS_PROGRAM_ID)
            .is_none_or(|pp| ids.contains_key(pp))
    });
    if chain_intact && programs.add("chopper", pg.clone()).is_ok() {
        return;
    }
    let ids = programs.as_mut();
    let mut id = String::from("chopper");
    let mut n = 0;
    while ids.contains_key(id.as_bytes()) {
        n += 1;
        id = format!("chopper.{n}");
    }
    ids.insert(id.into(), pg);
}

/// Where kept reads go: FASTQ text or BGZF-compressed BAM, both on stdout.
pub enum OutputWriter {
    Fastq(io::BufWriter<StdoutLock<'static>>),
    /// BGZF compression runs on its own worker threads, which need to own the
    /// stream, so this writes to `Stdout` rather than to a lock on it.
    Bam(
        bam::io::Writer<bgzf::io::MultithreadedWriter<Stdout>>,
        Box<sam::Header>,
    ),
}

impl OutputWriter {
    /// For BAM, writes the header straight away, so that a run that keeps no
    /// reads still produces a valid, empty BAM file.
    ///
    /// BGZF compression gets a quarter of `threads` (at least one) on top of
    /// the filtering threads, rather than as many again: compressing is
    /// cheaper than filtering, and the total stays close to what `--threads`
    /// asks for.
    pub fn new(format: Format, header: Option<&sam::Header>, threads: usize) -> io::Result<Self> {
        match format {
            Format::Fastq => Ok(OutputWriter::Fastq(io::BufWriter::new(io::stdout().lock()))),
            Format::Bam => {
                let header = header.expect("BAM output needs a header").clone();
                let workers = NonZeroUsize::new(threads / 4).unwrap_or(NonZeroUsize::MIN);
                let mut writer = bam::io::Writer::from(
                    bgzf::io::MultithreadedWriter::with_worker_count(workers, io::stdout()),
                );
                writer.write_header(&header)?;
                Ok(OutputWriter::Bam(writer, Box::new(header)))
            }
        }
    }

    pub fn write(&mut self, record: &WritableRecord) -> io::Result<()> {
        match (self, record) {
            (OutputWriter::Fastq(buf), WritableRecord::Fastq(_)) => record.write_on_buffer(buf),
            (OutputWriter::Bam(writer, header), WritableRecord::Bam(record)) => {
                writer.write_alignment_record(header, record)
            }
            _ => unreachable!("records are built in the output format"),
        }
    }

    /// Flushes everything and, for BAM, writes the BGZF end-of-file marker.
    pub fn finish(self) -> io::Result<()> {
        match self {
            OutputWriter::Fastq(mut buf) => buf.flush(),
            // Stdout is line buffered and the BGZF end-of-file block holds no
            // newline, so without this flush a write error on the last bytes
            // would surface only at exit, where it is ignored.
            OutputWriter::Bam(writer, _) => writer.into_inner().finish()?.flush(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::alignment::{
        record::Flags,
        record_buf::{QualityScores, Sequence},
    };

    fn ubam_record(data: Data) -> RecordBuf {
        RecordBuf::builder()
            .set_name("read1")
            .set_flags(Flags::UNMAPPED)
            .set_sequence(Sequence::from(b"ACGTC".to_vec()))
            .set_quality_scores(QualityScores::from(vec![10, 20, 30, 40, 50]))
            .set_data(data)
            .build()
    }

    #[test]
    fn bam_record_becomes_fastq_with_all_tags_as_text() {
        let data: Data = [
            (Tag::from(*b"qs"), Value::Float(12.5)),
            (Tag::from(*b"MM"), Value::from("C+m?,0;")),
            (Tag::from(*b"ML"), Value::Array(Array::UInt8(vec![200]))),
            (Tag::from(*b"RG"), Value::from("run1_model")),
        ]
        .into_iter()
        .collect();
        let fastq = bam_to_fastq(&ubam_record(data), TagText::All).unwrap();
        assert_eq!(
            fastq.id(),
            "read1\tqs:f:12.5\tMM:Z:C+m?,0;\tML:B:C,200\tRG:Z:run1_model"
        );
        assert_eq!(fastq.seq(), b"ACGTC");
        assert_eq!(fastq.qual(), b"+5?IS");
    }

    #[test]
    fn mods_only_keeps_just_the_modification_tags() {
        let data: Data = [
            (Tag::from(*b"qs"), Value::Float(12.5)),
            (Tag::from(*b"MM"), Value::from("C+m?,0;")),
            (Tag::from(*b"MN"), Value::UInt16(5)),
        ]
        .into_iter()
        .collect();
        let fastq = bam_to_fastq(&ubam_record(data), TagText::ModsOnly).unwrap();
        assert_eq!(fastq.id(), "read1\tMM:Z:C+m?,0;\tMN:i:5");
        let fastq = bam_to_fastq(&ubam_record(Data::default()), TagText::None).unwrap();
        assert_eq!(fastq.id(), "read1");
    }

    #[test]
    fn aligned_bam_records_are_rejected() {
        let mut record = ubam_record(Data::default());
        *record.flags_mut() = Flags::empty();
        assert!(bam_to_fastq(&record, TagText::All)
            .unwrap_err()
            .contains("aligned"));
    }

    #[test]
    fn bam_records_without_qualities_are_rejected() {
        let mut record = ubam_record(Data::default());
        *record.quality_scores_mut() = QualityScores::default();
        assert!(bam_to_fastq(&record, TagText::All)
            .unwrap_err()
            .contains("no base qualities"));
    }

    #[test]
    fn tags_round_trip_through_text() {
        let data: Data = [
            (Tag::from(*b"qs"), Value::Float(0.1)),
            (Tag::from(*b"ch"), Value::UInt16(1234)),
            (
                Tag::from(*b"mv"),
                Value::Array(Array::Int8(vec![5, 1, 0, -1])),
            ),
            (Tag::from(*b"CO"), Value::from("a comment with spaces")),
            (Tag::from(*b"tp"), Value::Character(b'A')),
        ]
        .into_iter()
        .collect();
        let text = tags_to_text(&data, |_| true).unwrap();
        let back = text_to_tags(&text).unwrap();
        assert_eq!(back.get(&Tag::from(*b"qs")), Some(&Value::Float(0.1)));
        assert_eq!(
            back.get(&Tag::from(*b"mv")),
            Some(&Value::Array(Array::Int8(vec![5, 1, 0, -1])))
        );
        assert_eq!(
            back.get(&Tag::from(*b"CO")),
            Some(&Value::from("a comment with spaces"))
        );
        // SAM text has a single integer type, so the width is not kept.
        assert_eq!(back.get(&Tag::from(*b"ch")), Some(&Value::UInt16(1234)));
    }

    #[test]
    fn invalid_tag_text_is_an_error() {
        assert!(text_to_tags(&["not-a-tag"]).is_err());
        assert!(text_to_tags(&["XY:i:1", "XY:i:2"]).is_err());
    }

    #[test]
    fn empty_arrays_keep_their_place() {
        // What --update-mods writes for a segment that kept no calls.
        let data = text_to_tags(&["MM:Z:C+m?;", "ML:B:C", "MN:i:11"]).unwrap();
        let tags: Vec<Tag> = data.keys().collect();
        assert_eq!(
            tags,
            [Tag::from(*b"MM"), Tag::from(*b"ML"), Tag::from(*b"MN")]
        );
        assert_eq!(
            data.get(&Tag::from(*b"ML")),
            Some(&Value::Array(Array::UInt8(vec![])))
        );
        assert_eq!(data.get(&Tag::from(*b"MN")), Some(&Value::UInt8(11)));
    }

    #[test]
    fn output_header_keeps_input_and_chains_a_program_line() {
        let input = sam::Header::builder()
            .add_program("basecaller", Map::<Program>::default())
            .add_comment("kept")
            .build();
        let header = output_header(Some(&input), true).unwrap();
        assert_eq!(header.comments(), input.comments());
        let chopper = header.programs().as_ref().get(&b"chopper"[..]).unwrap();
        assert_eq!(
            chopper
                .other_fields()
                .get(&program::tag::PREVIOUS_PROGRAM_ID)
                .map(|v| v.as_ref()),
            Some(&b"basecaller"[..])
        );
    }
}
