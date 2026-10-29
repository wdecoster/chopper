//! End-to-end tests of the input and output formats, running the built binary.
//!
//! FASTQ output is checked against digests of the v0.14.1 release; uBAM
//! output is read back with noodles and compared field by field against the
//! fixtures in `test-data/ubam`, which were written by htslib (see
//! `test-data/ubam/make_fixtures.py`).

use std::io::Write;
use std::process::{Command, Output, Stdio};

use noodles_bam as bam;
use noodles_sam::{
    self as sam,
    alignment::{
        io::Write as _,
        record::{data::field::Tag, Flags},
        record_buf::{data::field::Value, RecordBuf},
    },
};
use sha2::{Digest, Sha256};

const UBAM: &str = "test-data/ubam/reads.bam";
const TAGGED_FASTQ: &str = "test-data/ubam/reads.fastq";

fn chopper(args: &[&str], stdin: Option<&[u8]>) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_chopper"))
        .args(args)
        .stdin(if stdin.is_some() {
            Stdio::piped()
        } else {
            Stdio::null()
        })
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("chopper runs");
    // Feed stdin from its own thread: writing all of it before reading stdout
    // deadlocks once both pipe buffers are full.
    let feeder = stdin.map(|input| {
        let mut pipe = child.stdin.take().unwrap();
        let input = input.to_vec();
        std::thread::spawn(move || pipe.write_all(&input))
    });
    let output = child.wait_with_output().unwrap();
    if let Some(feeder) = feeder {
        feeder.join().unwrap().unwrap();
    }
    output
}

/// Runs chopper and returns its stdout, failing the test if it did not succeed.
fn run(args: &[&str]) -> Vec<u8> {
    let out = chopper(args, None);
    assert!(
        out.status.success(),
        "chopper {args:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

fn read_bam(bytes: &[u8]) -> (sam::Header, Vec<RecordBuf>) {
    let mut reader = bam::io::Reader::new(bytes);
    let header = reader.read_header().expect("valid BAM header");
    let records = reader
        .record_bufs(&header)
        .collect::<Result<_, _>>()
        .expect("valid BAM records");
    (header, records)
}

fn fixture_bam() -> (sam::Header, Vec<RecordBuf>) {
    read_bam(&std::fs::read(UBAM).unwrap())
}

fn tag(name: &[u8; 2]) -> Tag {
    Tag::from(*name)
}

/// The value of an integer tag regardless of its width, since SAM text (and so
/// FASTQ) has a single integer type. Every integer in the fixtures fits in i32.
fn widen(value: &Value) -> Value {
    let int = |n: i64| Value::Int32(i32::try_from(n).expect("fixture integers fit in i32"));
    match value {
        Value::Int8(n) => int(i64::from(*n)),
        Value::UInt8(n) => int(i64::from(*n)),
        Value::Int16(n) => int(i64::from(*n)),
        Value::UInt16(n) => int(i64::from(*n)),
        Value::Int32(n) => int(i64::from(*n)),
        Value::UInt32(n) => int(i64::from(*n)),
        v => v.clone(),
    }
}

fn sha256_hex(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

/// Sorts FASTQ records the way `paste - - - - | LC_ALL=C sort | tr '\t' '\n'`
/// does, to compare multi-threaded output, whose order is not deterministic.
fn sort_fastq(bytes: &[u8]) -> Vec<u8> {
    let text = String::from_utf8(bytes.to_vec()).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    let mut records: Vec<String> = lines.chunks(4).map(|r| r.join("\t")).collect();
    records.sort();
    records
        .iter()
        .map(|r| r.replace('\t', "\n") + "\n")
        .collect::<String>()
        .into_bytes()
}

#[test]
fn fastq_output_is_identical_to_v0_14_1() {
    let digests = std::fs::read_to_string("test-data/regression/v0.14.1.sha256").unwrap();
    let test_fastq = std::fs::read("test-data/test.fastq").unwrap();
    let mut checked = 0;
    for line in digests.lines().filter(|l| !l.starts_with('#')) {
        let mut fields = line.splitn(3, "  ");
        let (expected, mode, args) = (
            fields.next().unwrap(),
            fields.next().unwrap(),
            fields.next().unwrap(),
        );
        let args: Vec<&str> = args.split(' ').collect();
        let output = match mode {
            "exact" => run(&args),
            "sorted" => sort_fastq(&run(&args)),
            "stdin" => {
                let out = chopper(&args, Some(&test_fastq));
                assert!(out.status.success());
                out.stdout
            }
            _ => panic!("unknown mode {mode}"),
        };
        let actual = sha256_hex(&output);
        assert_eq!(actual, expected, "output of chopper {args:?} changed");
        checked += 1;
    }
    assert_eq!(checked, 14);
}

#[test]
fn ubam_to_ubam_without_trimming_keeps_every_record_unchanged() {
    let (input_header, input) = fixture_bam();
    let (header, output) = read_bam(&run(&["-t", "1", "-i", UBAM]));

    // Every field of every record, tags with their exact types and order.
    assert_eq!(output, input);

    // The header is passed through, with chopper chained after the last @PG.
    assert_eq!(header.read_groups(), input_header.read_groups());
    assert_eq!(header.comments(), input_header.comments());
    assert_eq!(header.header(), input_header.header());
    let programs = header.programs().as_ref();
    assert_eq!(
        programs.keys().map(|k| k.to_string()).collect::<Vec<_>>(),
        ["basecaller", "samtools", "chopper"]
    );
    let pg = &programs[&b"chopper"[..]];
    let field = |t: &[u8; 2]| {
        pg.other_fields()
            .iter()
            .find(|(k, _)| k.as_ref() == t)
            .map(|(_, v)| v.to_string())
    };
    assert_eq!(field(b"PP").as_deref(), Some("samtools"));
    assert_eq!(field(b"PN").as_deref(), Some("chopper"));
    assert!(field(b"CL")
        .unwrap()
        .ends_with("-t 1 -i test-data/ubam/reads.bam"));
}

#[test]
fn ubam_to_ubam_trimming_rewrites_only_sequence_and_mod_tags() {
    let crop = [
        "-t",
        "1",
        "--trim-approach",
        "fixed-crop",
        "--headcrop",
        "7",
        "--tailcrop",
        "5",
        "--update-mods",
    ];
    let (_, input) = fixture_bam();
    let (_, output) = read_bam(&run(&[&["-i", UBAM][..], &crop].concat()));
    // The same trimming on the FASTQ fixture, as the reference for MM/ML/MN.
    let fastq = String::from_utf8(run(&[&["-i", TAGGED_FASTQ][..], &crop].concat())).unwrap();
    let fastq_headers: Vec<&str> = fastq.lines().step_by(4).collect();

    assert_eq!(output.len(), input.len());
    for ((out, inp), fastq_header) in output.iter().zip(&input).zip(fastq_headers) {
        let len = inp.sequence().len();
        assert_eq!(out.name(), inp.name());
        assert_eq!(out.flags(), inp.flags());
        assert_eq!(
            out.sequence().as_ref(),
            &inp.sequence().as_ref()[7..len - 5]
        );
        assert_eq!(
            out.quality_scores().as_ref(),
            &inp.quality_scores().as_ref()[7..len - 5]
        );
        // Same tags, same order; all but the modification tags untouched.
        let keys = |r: &RecordBuf| r.data().keys().collect::<Vec<_>>();
        assert_eq!(keys(out), keys(inp));
        for (t, value) in inp.data().iter() {
            if ![tag(b"MM"), tag(b"ML"), tag(b"MN")].contains(&t) {
                assert_eq!(out.data().get(&t), Some(value));
            }
        }
        if let Some(Value::String(mm)) = out.data().get(&tag(b"MM")) {
            assert!(fastq_header.contains(&format!("\tMM:Z:{mm}")));
        }
        if let Some(mn) = out.data().get(&tag(b"MN")) {
            assert_eq!(widen(mn), Value::Int32((len - 12) as i32));
        }
    }
}

#[test]
fn ubam_to_fastq_writes_all_tags_like_samtools_fastq() {
    // reads.fastq is reads.bam as htslib formats every tag (`samtools fastq -T '*'`).
    let expected = std::fs::read(TAGGED_FASTQ).unwrap();
    let out = run(&["-t", "1", "-i", UBAM, "--output-format", "fastq"]);
    assert_eq!(
        String::from_utf8(out).unwrap(),
        String::from_utf8(expected).unwrap()
    );
}

#[test]
fn ubam_to_fastq_matches_fastq_to_fastq_when_trimming() {
    let trim = [
        "-t",
        "1",
        "--trim-approach",
        "split-by-low-quality",
        "--cutoff",
        "8",
        "-l",
        "10",
        "--update-mods",
    ];
    let from_bam = run(&[&["-i", UBAM, "--output-format", "fastq"][..], &trim].concat());
    let from_fastq = run(&[&["-i", TAGGED_FASTQ][..], &trim].concat());
    assert!(!from_bam.is_empty());
    assert_eq!(from_bam, from_fastq);
}

#[test]
fn fastq_to_ubam_types_the_header_tags() {
    let (_, expected) = fixture_bam();
    let (header, output) = read_bam(&run(&[
        "-t",
        "1",
        "-i",
        TAGGED_FASTQ,
        "--output-format",
        "bam",
    ]));

    // A minimal header: no @RG/@CO to pass through from FASTQ.
    assert!(header.read_groups().is_empty());
    assert_eq!(
        header
            .programs()
            .as_ref()
            .keys()
            .map(|k| k.to_string())
            .collect::<Vec<_>>(),
        ["chopper"]
    );

    assert_eq!(output.len(), expected.len());
    for (out, exp) in output.iter().zip(&expected) {
        assert_eq!(out.name(), exp.name());
        assert_eq!(out.flags(), exp.flags());
        assert_eq!(out.sequence(), exp.sequence());
        assert_eq!(out.quality_scores(), exp.quality_scores());
        // Same tags and values in the same order; only integer widths may
        // differ, since FASTQ text does not carry them.
        let tags = |r: &RecordBuf| {
            r.data()
                .iter()
                .map(|(t, v)| (t, widen(v)))
                .collect::<Vec<_>>()
        };
        assert_eq!(tags(out), tags(exp));
    }
}

#[test]
fn fastq_description_is_kept_as_a_comment_tag() {
    let dir = env!("CARGO_TARGET_TMPDIR");
    let path = format!("{dir}/description.fastq");
    std::fs::write(&path, "@r1 runid=abc ch=5\nACGTACGT\n+\nIIIIIIII\n").unwrap();
    let (_, records) = read_bam(&run(&["-i", &path, "--output-format", "bam"]));
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].name().unwrap(), "r1");
    assert_eq!(
        records[0].data().get(&tag(b"CO")),
        Some(&Value::from("runid=abc ch=5"))
    );
}

#[test]
fn split_segments_are_named_like_fastq_output() {
    let args = [
        "-t",
        "1",
        "--trim-approach",
        "split-by-low-quality",
        "--cutoff",
        "8",
        "-l",
        "10",
    ];
    let (_, records) = read_bam(&run(&[&["-i", UBAM][..], &args].concat()));
    let fastq = String::from_utf8(run(&[&["-i", TAGGED_FASTQ][..], &args].concat())).unwrap();
    let names: Vec<String> = fastq
        .lines()
        .step_by(4)
        .map(|h| h[1..].split('\t').next().unwrap().to_string())
        .collect();
    let bam_names: Vec<String> = records
        .iter()
        .map(|r| r.name().unwrap().to_string())
        .collect();
    assert!(bam_names.iter().any(|n| n.contains("_segment_")));
    assert_eq!(bam_names, names);
}

#[test]
fn ubam_on_stdin_is_detected() {
    let bam = std::fs::read(UBAM).unwrap();
    let out = chopper(&["-t", "1"], Some(&bam));
    assert!(out.status.success());
    let (_, records) = read_bam(&out.stdout);
    assert_eq!(records, fixture_bam().1);
}

#[test]
fn multithreaded_ubam_output_holds_the_same_records() {
    let sorted = |bytes: Vec<u8>| {
        let (_, mut records) = read_bam(&bytes);
        records.sort_by(|a, b| a.name().cmp(&b.name()));
        records
    };
    let args = ["--trim-approach", "trim-by-quality", "--cutoff", "15"];
    assert_eq!(
        sorted(run(&[&["-t", "4", "-i", UBAM][..], &args].concat())),
        sorted(run(&[&["-t", "1", "-i", UBAM][..], &args].concat()))
    );
}

#[test]
fn ubam_output_without_reads_is_a_valid_empty_bam() {
    let (header, records) = read_bam(&run(&["-i", UBAM, "-l", "100000"]));
    assert!(records.is_empty());
    assert!(header.programs().as_ref().contains_key(&b"chopper"[..]));
}

#[test]
fn aligned_bam_is_rejected() {
    let out = chopper(&["-i", "test-data/ubam/aligned.bam"], None);
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("read mapped1 is aligned"), "{stderr}");
}

/// Writes a BAM file into the test's temporary directory and returns its path.
fn write_bam(name: &str, header: &sam::Header, records: &[RecordBuf]) -> String {
    let path = format!("{}/{name}", env!("CARGO_TARGET_TMPDIR"));
    let mut writer = bam::io::Writer::new(std::fs::File::create(&path).unwrap());
    writer.write_header(header).unwrap();
    for record in records {
        writer.write_alignment_record(header, record).unwrap();
    }
    writer.try_finish().unwrap();
    path
}

/// The fixture's header with its text edited, for headers noodles would not build.
fn fixture_header_with(edit: impl Fn(String) -> String) -> sam::Header {
    let (header, _) = fixture_bam();
    let mut text = Vec::new();
    sam::io::Writer::new(&mut text)
        .write_header(&header)
        .unwrap();
    edit(String::from_utf8(text).unwrap())
        .parse()
        .expect("edited header parses")
}

fn program_ids(header: &sam::Header) -> Vec<String> {
    header
        .programs()
        .as_ref()
        .keys()
        .map(|k| k.to_string())
        .collect()
}

#[test]
fn a_dangling_previous_program_does_not_stop_bam_output() {
    // `PP` naming a program that is not in the header makes noodles panic
    // while chaining a new @PG line.
    let header = fixture_header_with(|h| h.replace("PP:basecaller", "PP:missing"));
    let path = write_bam("dangling_pp.bam", &header, &fixture_bam().1);
    let (out_header, records) = read_bam(&run(&["-t", "1", "-i", &path]));
    assert_eq!(records, fixture_bam().1);
    assert_eq!(
        program_ids(&out_header),
        ["basecaller", "samtools", "chopper"]
    );
    let chopper = &out_header.programs().as_ref()[&b"chopper"[..]];
    assert!(!chopper
        .other_fields()
        .iter()
        .any(|(k, _)| k.as_ref() == b"PP"));
}

#[test]
fn sort_order_is_only_kept_when_output_order_is() {
    let header = fixture_header_with(|h| h.replace("SO:unknown", "SO:queryname\tGO:query"));
    let path = write_bam("queryname.bam", &header, &fixture_bam().1);
    let sort_fields = |args: &[&str]| {
        let (header, _) = read_bam(&run(&[&["-i", &path][..], args].concat()));
        let hd = header.header().unwrap().other_fields();
        let get = |t: &[u8; 2]| {
            hd.iter()
                .find(|(k, _)| k.as_ref() == t)
                .map(|(_, v)| v.to_string())
        };
        (get(b"SO"), get(b"GO"))
    };
    let kept = (Some("queryname".to_string()), Some("query".to_string()));
    let dropped = (Some("unsorted".to_string()), None);
    assert_eq!(sort_fields(&["-t", "1"]), kept);
    assert_eq!(sort_fields(&["-t", "4"]), dropped);
    assert_eq!(
        sort_fields(&[
            "-t",
            "1",
            "--trim-approach",
            "split-by-low-quality",
            "--cutoff",
            "8"
        ]),
        dropped
    );
}

#[test]
fn secondary_supplementary_and_reverse_records_are_rejected() {
    let (header, records) = fixture_bam();
    for (flag, expected) in [
        (0x4 | 0x100, "secondary or supplementary"),
        (0x4 | 0x800, "secondary or supplementary"),
        (0x4 | 0x10, "reverse complemented"),
    ] {
        let mut record = records[0].clone();
        *record.flags_mut() = Flags::from(flag);
        let path = write_bam(&format!("flag_{flag}.bam"), &header, &[record]);
        let out = chopper(&["-i", &path], None);
        assert_eq!(out.status.code(), Some(1), "flag {flag}");
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(stderr.contains(expected), "flag {flag}: {stderr}");
    }
}

#[test]
fn a_read_name_bam_cannot_hold_is_reported_with_the_read() {
    let dir = env!("CARGO_TARGET_TMPDIR");
    let long = "r".repeat(300);
    for (file, name) in [("long_name.fastq", long.as_str()), ("empty_name.fastq", "")] {
        let path = format!("{dir}/{file}");
        std::fs::write(&path, format!("@{name}\nACGT\n+\nIIII\n")).unwrap();
        let out = chopper(&["-i", &path, "--output-format", "bam"], None);
        assert_eq!(out.status.code(), Some(1));
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(stderr.contains("is not a valid BAM read name"), "{stderr}");
        assert!(!stderr.contains("failed to write output"), "{stderr}");
    }
}

#[test]
fn gzipped_fastq_is_read_from_a_file_and_from_stdin() {
    let plain = std::fs::read(TAGGED_FASTQ).unwrap();
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gz.write_all(&plain).unwrap();
    let gz = gz.finish().unwrap();
    let path = format!("{}/reads.fastq.gz", env!("CARGO_TARGET_TMPDIR"));
    std::fs::write(&path, &gz).unwrap();

    let expected = run(&["-t", "1", "-i", TAGGED_FASTQ]);
    assert_eq!(run(&["-t", "1", "-i", &path]), expected);
    let out = chopper(&["-t", "1"], Some(&gz));
    assert!(out.status.success());
    assert_eq!(out.stdout, expected);
}

#[test]
fn legacy_mm_ml_spelling_is_rewritten_in_place_in_bam() {
    // Pre-standard dorado wrote Mm/Ml; rename the fixture's tags to those.
    let (header, records) = fixture_bam();
    let legacy: Vec<RecordBuf> = records
        .iter()
        .map(|r| {
            let mut r = r.clone();
            let data = r
                .data()
                .iter()
                .map(|(t, v)| match t.as_ref() {
                    b"MM" => (tag(b"Mm"), v.clone()),
                    b"ML" => (tag(b"Ml"), v.clone()),
                    _ => (t, v.clone()),
                })
                .collect();
            *r.data_mut() = data;
            r
        })
        .collect();
    let path = write_bam("legacy_mods.bam", &header, &legacy);
    let crop = [
        "-t",
        "1",
        "--trim-approach",
        "fixed-crop",
        "--headcrop",
        "7",
        "--tailcrop",
        "5",
        "--update-mods",
    ];
    let (_, from_legacy) = read_bam(&run(&[&["-i", &path][..], &crop].concat()));
    let (_, from_standard) = read_bam(&run(&[&["-i", UBAM][..], &crop].concat()));
    for (legacy, standard) in from_legacy.iter().zip(&from_standard) {
        assert!(legacy.data().get(&tag(b"MM")).is_none());
        assert_eq!(legacy.data().len(), standard.data().len());
        assert_eq!(
            legacy.data().get(&tag(b"Mm")),
            standard.data().get(&tag(b"MM"))
        );
        assert_eq!(
            legacy.data().get(&tag(b"Ml")),
            standard.data().get(&tag(b"ML"))
        );
    }
}

#[test]
fn multithreaded_ubam_split_with_update_mods_matches_one_thread() {
    let sorted = |bytes: Vec<u8>| {
        let (_, mut records) = read_bam(&bytes);
        records.sort_by(|a, b| a.name().cmp(&b.name()));
        records
    };
    let args = [
        "-i",
        UBAM,
        "--trim-approach",
        "split-by-low-quality",
        "--cutoff",
        "8",
        "-l",
        "10",
        "--update-mods",
    ];
    let one = sorted(run(&[&["-t", "1"][..], &args].concat()));
    assert!(one.len() > 6);
    assert_eq!(sorted(run(&[&["-t", "4"][..], &args].concat())), one);
}
