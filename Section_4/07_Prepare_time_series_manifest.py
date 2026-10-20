#!/usr/bin/env python3
# ============================================================
# 07_Prepare_time_series_manifest.py
# Tool    : standard library
# Input   : Clean-read FASTQ directory (searched recursively)
#             /path/to/WWTP_filtering/clean_data/zfz-reads/DC/b/fastq
# Output  : 07_samples/00_manifest/DC_read_pairs.tsv
# ============================================================

import argparse
import csv
import re
from pathlib import Path

SAMPLE_PATTERN = re.compile(r"(d([3-8])((?:0[1-9])|(?:1[0-2]))b)", re.IGNORECASE)

MATE_PATTERN = re.compile(
    r"(?:^|[._-])R?([12])(?:[._-]?001)"
    r"(?=\.(?:fastq|fq)(?:\.gz)?$)",
    re.IGNORECASE
)

FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--read-root", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    read_root = Path(args.read_root).resolve()
    output = Path(args.output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    if not read_root.is_dir():
        raise FileNotFoundError(read_root)

    pairs = {}
    unrecognized = []

    for path in sorted(read_root.rglob("*")):
        if not path.is_file():
            continue

        if not path.name.lower().endswith(FASTQ_SUFFIXES):
            continue

        sample_match = SAMPLE_PATTERN.search(path.name)
        if sample_match is None:
            unrecognized.append(str(path))
            continue

        sample = sample_match.group(1).lower()
        suffix = path.name[sample_match.end():]

        mate_match = MATE_PATTERN.search(suffix)
        if mate_match is None:
            unrecognized.append(str(path))
            continue

        mate = int(mate_match.group(1))
        pairs.setdefault(sample, {}).setdefault(mate, []).append(path.resolve())

    rows = []

    for sample, mates in pairs.items():
        if 1 not in mates or 2 not in mates:
            print(f"WARNING: incomplete pair for {sample}: {mates}")
            continue

        if len(mates[1]) != 1 or len(mates[2]) != 1:
            raise RuntimeError(
                f"Multiple FASTQs detected for {sample}: "
                f"R1={mates[1]}, R2={mates[2]}"
            )

        year = 2010 + int(sample[1])
        month = int(sample[2:4])

        rows.append({
            "sample": sample,
            "year_month": f"{year:04d}-{month:02d}",
            "timepoint_month": (year - 2013) * 12 + (month - 1),
            "R1": str(mates[1][0]),
            "R2": str(mates[2][0])
        })

    rows.sort(key=lambda x: (x["timepoint_month"], x["sample"]))

    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample", "year_month", "timepoint_month", "R1", "R2"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Complete paired samples: {len(rows)}")
    print(f"Unrecognized FASTQ files: {len(unrecognized)}")
    print(f"Output: {output}")

    if rows:
        print(f"First sample: {rows[0]['sample']} {rows[0]['year_month']}")
        print(f"Last sample:  {rows[-1]['sample']} {rows[-1]['year_month']}")


if __name__ == "__main__":
    main()
