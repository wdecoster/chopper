"""Generates the uBAM test fixtures in this directory.

    python3 test-data/ubam/make_fixtures.py

Writes, from the same synthetic reads:

* reads.bam    an unaligned BAM shaped like dorado's output: a header with
               @RG/@PG/@CO lines and typed tags on every record (MM/ML/MN,
               qs, ch, du, mv, RG, ...), one read without tags.
* reads.fastq  those reads as `samtools fastq -T '*'` writes them: every tag
               as a tab-separated SAM field after the read name.
* aligned.bam  a single mapped read, which chopper must refuse.

The output is deterministic (fixed seed), and is written by htslib through
pysam rather than by chopper, so the tests compare chopper against an
independent implementation.
"""

import array
import os
import random

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
rng = random.Random(42)


def random_read(length):
    seq = "".join(rng.choice("ACGT") for _ in range(length))
    qual = [rng.randint(2, 40) for _ in range(length)]
    return seq, qual


def mm_for(seq, base, codes, every):
    """An MM/ML pair calling every `every`-th occurrence of `base`."""
    n = seq.count(base) if base != "N" else len(seq)
    called = list(range(0, n, every))
    skips, last = [], -1
    for occ in called:
        skips.append(occ - last - 1)
        last = occ
    mm = f"{base}+{codes}?" + "".join(f",{s}" for s in skips) + ";"
    ml = [rng.randint(0, 255) for _ in range(len(called) * len(codes))]
    return mm, ml


def make_reads():
    reads = []
    for i, (length, mods) in enumerate(
        [
            (180, ("C", "m", 3)),
            (240, ("C", "mh", 2)),
            (90, ("A", "a", 4)),
            (300, ("C", "m", 1)),
            (60, None),
            (150, ("C", "m", 5)),
        ]
    ):
        seq, qual = random_read(length)
        tags = [("qs", 12.5, "f"), ("du", 0.5, "f"), ("ch", 1234 + i, "i")]
        if mods is not None:
            mm, ml = mm_for(seq, *mods)
            tags += [("MM", mm, "Z"), ("ML", ml, "C"), ("MN", length, "i")]
        tags += [
            ("mv", [5, 1, 0, 1, -1], "c"),
            ("RG", "run1_model", "Z"),
        ]
        if i == 3:
            tags.append(("CO", "a comment with spaces", "Z"))
        if mods is None:
            tags = []
        reads.append((f"read{i + 1}", seq, qual, tags))
    return reads


def header():
    return {
        "HD": {"VN": "1.6", "SO": "unknown"},
        "RG": [{"ID": "run1_model", "PL": "ONT", "DS": "fixture"}],
        "PG": [
            {"ID": "basecaller", "PN": "dorado", "VN": "1.0.0", "CL": "dorado basecaller"},
            {"ID": "samtools", "PN": "samtools", "PP": "basecaller", "VN": "1.21"},
        ],
        "CO": ["chopper uBAM test fixture"],
    }


def segment(name, seq, qual, tags, hdr):
    seg = pysam.AlignedSegment(hdr)
    seg.query_name = name
    seg.flag = 4
    seg.reference_id = -1
    seg.reference_start = -1
    seg.mapping_quality = 0
    seg.query_sequence = seq
    seg.query_qualities = pysam.qualitystring_to_array("".join(chr(q + 33) for q in qual))
    for tag, value, ty in tags:
        if ty in ("C", "c"):
            # B arrays: the array's typecode sets the subtype (uint8 / int8).
            seg.set_tag(tag, array.array("B" if ty == "C" else "b", value))
        else:
            seg.set_tag(tag, value, value_type=ty)
    return seg


def main():
    hdr = pysam.AlignmentHeader.from_dict(header())
    reads = make_reads()

    with pysam.AlignmentFile(os.path.join(HERE, "reads.bam"), "wb", header=hdr) as out:
        for read in reads:
            out.write(segment(*read, hdr))

    # samtools fastq -T '*': the SAM text of every tag, after the name.
    with open(os.path.join(HERE, "reads.fastq"), "w") as out, pysam.AlignmentFile(
        os.path.join(HERE, "reads.bam"), "rb", check_sq=False
    ) as bam:
        for seg in bam:
            tags = seg.to_string().split("\t")[11:]
            head = "\t".join([seg.query_name] + tags)
            qual = "".join(chr(q + 33) for q in seg.query_qualities)
            out.write(f"@{head}\n{seg.query_sequence}\n+\n{qual}\n")

    aligned = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    )
    with pysam.AlignmentFile(os.path.join(HERE, "aligned.bam"), "wb", header=aligned) as out:
        seq, qual = random_read(50)
        seg = pysam.AlignedSegment(aligned)
        seg.query_name = "mapped1"
        seg.flag = 0
        seg.reference_id = 0
        seg.reference_start = 100
        seg.mapping_quality = 60
        seg.cigarstring = "50M"
        seg.query_sequence = seq
        seg.query_qualities = pysam.qualitystring_to_array("".join(chr(q + 33) for q in qual))
        out.write(seg)


if __name__ == "__main__":
    main()
