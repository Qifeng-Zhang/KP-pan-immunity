#!/usr/bin/env python3
# ============================================================
# 09_Build_chronostrain_reads.py
# Tool    : standard library
# Input   : Time-series manifest
#             07_samples/00_manifest/DC_read_pairs.tsv
#           Read counts table, two columns: sample <tab> read_pairs
#             07_samples/01_read_counts/DC_read_counts.tsv
# Output  : 07_samples/02_reads_table/DC_chronostrain_reads.tsv
#             headerless: timepoint, sample, depth, path, read type,
#                         quality format
#           07_samples/02_reads_table/DC_sample_metadata.tsv
# ============================================================

import argparse
import csv
import re
from pathlib import Path


def normalize_sample(sample: str) -> str:
    return sample.strip().lower()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--metadata-output", required=True)
    args = parser.parse_args()

    pairs_path = Path(args.pairs).resolve()
    counts_path = Path(args.counts).resolve()
    output_path = Path(args.output).resolve()
    metadata_path = Path(args.metadata_output).resolve()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)

    read_counts = {}

    with counts_path.open(encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")

        for line_number, row in enumerate(reader, start=1):
            if not row or not any(field.strip() for field in row):
                continue

            if len(row) < 2:
                raise ValueError(
                    f"Counts table line {line_number} has fewer than "
                    f"two columns: {row}"
                )

            sample = normalize_sample(row[0])

            if line_number == 1 and sample in {"sample", "sample_id", "sampleid"}:
                continue

            try:
                read_count = int(row[1].strip().replace(",", ""))
            except ValueError as exc:
                raise ValueError(
                    f"Invalid read count at line {line_number}: {row[1]}"
                ) from exc

            if read_count <= 0:
                raise ValueError(
                    f"Non-positive read count for {sample}: {read_count}"
                )

            if sample in read_counts:
                raise ValueError(f"Duplicated sample in counts table: {sample}")

            read_counts[sample] = read_count

    # --- Build the ChronoStrain read table ---------------------------------
    manifest_rows = []
    metadata_rows = []
    pair_samples = set()

    sample_regex = re.compile(r"^d([3-8])((?:0[1-9])|(?:1[0-2]))b$", re.IGNORECASE)

    with pairs_path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required_columns = {"sample", "year_month", "timepoint_month", "R1", "R2"}

        if reader.fieldnames is None:
            raise RuntimeError("Pair table has no header.")

        missing_columns = required_columns - set(reader.fieldnames)
        if missing_columns:
            raise RuntimeError(
                f"Pair table is missing columns: {sorted(missing_columns)}"
            )

        for row in reader:
            sample = normalize_sample(row["sample"])
            pair_samples.add(sample)

            if sample not in read_counts:
                raise KeyError(
                    f"Sample is present in pair table but absent from "
                    f"counts table: {sample}"
                )

            if sample_regex.fullmatch(sample) is None:
                raise ValueError(f"Unexpected DC sample name: {sample}")

            r1 = Path(row["R1"]).expanduser().resolve()
            r2 = Path(row["R2"]).expanduser().resolve()

            if not r1.is_file():
                raise FileNotFoundError(f"R1 does not exist: {r1}")
            if not r2.is_file():
                raise FileNotFoundError(f"R2 does not exist: {r2}")

            timepoint = float(row["timepoint_month"])
            read_count = read_counts[sample]

            manifest_rows.extend([
                [f"{timepoint:.1f}", sample, str(read_count),
                 str(r1), "paired_1", "fastq-sanger"],
                [f"{timepoint:.1f}", sample, str(read_count),
                 str(r2), "paired_2", "fastq-sanger"]
            ])

            metadata_rows.append({
                "sample": sample,
                "year_month": row["year_month"],
                "timepoint_month": f"{timepoint:.1f}",
                "read_pairs": read_count,
                "R1": str(r1),
                "R2": str(r2)
            })

    extra_count_samples = sorted(set(read_counts) - pair_samples)
    if extra_count_samples:
        raise RuntimeError(
            "Samples are present in counts table but absent from pair table: "
            + ", ".join(extra_count_samples)
        )

    manifest_rows.sort(
        key=lambda row: (float(row[0]), row[1], 0 if row[4] == "paired_1" else 1)
    )
    metadata_rows.sort(
        key=lambda row: (float(row["timepoint_month"]), row["sample"])
    )

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        csv.writer(handle, delimiter="\t").writerows(manifest_rows)

    with metadata_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample", "year_month", "timepoint_month",
                        "read_pairs", "R1", "R2"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(metadata_rows)

    print(f"Samples in count table: {len(read_counts)}")
    print(f"Samples in pair table:  {len(pair_samples)}")
    print(f"ChronoStrain rows:      {len(manifest_rows)}")
    print(f"First sample:           {metadata_rows[0]['sample']}")
    print(f"Last sample:            {metadata_rows[-1]['sample']}")
    print(f"Output:                 {output_path}")


if __name__ == "__main__":
    main()
