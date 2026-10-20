#!/usr/bin/env python3
# ============================================================
# 11_Fix_pair_ids.py
# Tool    : standard library
# Input   : Filtered ChronoStrain read table
#             07_samples/03_filter/filtered_DC_chronostrain_reads.csv
# Output  : 07_samples/04_reads_pairfixed/<original fastq names>
#           07_samples/04_reads_pairfixed/filtered_DC_chronostrain_reads_pairfixed.csv
#           07_samples/04_reads_pairfixed/pairing_summary.tsv
# ============================================================

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


def open_text(path: Path, mode: str):
    if path.name.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def rewrite_fastq(src: Path, dst: Path, required_suffix: str) -> set:
    """Add /1 or /2 to the first token of each FASTQ header.

    Returns the set of mate-prefix identifiers, i.e. read identifiers
    with the /1 or /2 suffix removed.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(dst.stem + ".tmp" + dst.suffix)

    mate_prefixes = set()
    n_records = 0

    header_pattern = re.compile(r"^@(\S+)(.*?)(\r?\n)$")

    with open_text(src, "rt") as fin, open_text(tmp, "wt") as fout:
        while True:
            header = fin.readline()
            if header == "":
                break

            sequence = fin.readline()
            plus = fin.readline()
            quality = fin.readline()

            if not sequence or not plus or not quality:
                raise RuntimeError(
                    f"Truncated FASTQ record in {src}, near record {n_records + 1}"
                )

            match = header_pattern.match(header)
            if match is None:
                raise RuntimeError(f"Invalid FASTQ header in {src}: {header.rstrip()}")

            old_id = match.group(1)
            remainder = match.group(2)
            newline = match.group(3)

            if old_id.endswith("/1") or old_id.endswith("/2"):
                if not old_id.endswith(required_suffix):
                    raise RuntimeError(
                        f"Unexpected mate suffix in {src}: {old_id}, "
                        f"expected {required_suffix}"
                    )
                new_id = old_id
            else:
                new_id = old_id + required_suffix

            mate_prefix = new_id[:-2]

            if mate_prefix in mate_prefixes:
                raise RuntimeError(f"Duplicated read identifier in {src}: {new_id}")

            mate_prefixes.add(mate_prefix)

            fout.write(f"@{new_id}{remainder}{newline}")
            fout.write(sequence)
            fout.write(plus)
            fout.write(quality)

            n_records += 1

    tmp.replace(dst)

    print(f"{src.name} -> {dst.name}: {n_records} reads, suffix {required_suffix}")

    return mate_prefixes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-manifest", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--output-manifest-name", required=True)
    args = parser.parse_args()

    input_manifest = Path(args.input_manifest).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    output_manifest = output_dir / args.output_manifest_name
    pairing_summary = output_dir / "pairing_summary.tsv"

    rows_out = []
    sample_ids = defaultdict(dict)

    with input_manifest.open(newline="", encoding="utf-8") as handle:
        for line_number, row in enumerate(csv.reader(handle), start=1):
            if not row:
                continue

            if len(row) != 6:
                raise RuntimeError(
                    f"Manifest line {line_number}: expected 6 columns, "
                    f"found {len(row)}"
                )

            timepoint, sample, depth, path_str, read_type, quality = row

            src = Path(path_str)
            if not src.is_absolute():
                src = input_manifest.parent / src
            src = src.resolve()

            if not src.is_file():
                raise FileNotFoundError(src)

            if read_type == "paired_1":
                suffix = "/1"
            elif read_type == "paired_2":
                suffix = "/2"
            else:
                raise RuntimeError(
                    f"Unexpected read type in paired analysis: {read_type}"
                )

            dst = output_dir / src.name

            sample_ids[sample][read_type] = rewrite_fastq(src, dst, suffix)

            rows_out.append([timepoint, sample, depth, dst.name, read_type, quality])

    with output_manifest.open("w", newline="", encoding="utf-8") as handle:
        csv.writer(handle).writerows(rows_out)

    with pairing_summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "sample", "forward_reads", "reverse_reads",
            "matched_pairs", "forward_only", "reverse_only"
        ])

        for sample in sorted(sample_ids):
            forward = sample_ids[sample].get("paired_1", set())
            reverse = sample_ids[sample].get("paired_2", set())
            matched = forward & reverse

            writer.writerow([
                sample, len(forward), len(reverse), len(matched),
                len(forward - reverse), len(reverse - forward)
            ])

    print()
    print(f"Output manifest: {output_manifest}")
    print(f"Pairing summary: {pairing_summary}")
    print(f"Manifest rows: {len(rows_out)}")


if __name__ == "__main__":
    main()
