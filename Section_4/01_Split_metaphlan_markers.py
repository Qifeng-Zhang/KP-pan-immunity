#!/usr/bin/env python3
# ============================================================
# 01_Split_metaphlan_markers.py
# Tool    : Biopython
# Input   : MetaPhlAn marker FASTA of the target species
#             /path/to/Strain/mpa4_marker/t__SGB12599.fna
# Output  : /path/to/Strain/chronostrain_thauera/01_markers/split/<marker>.fna
#           /path/to/Strain/chronostrain_thauera/01_markers/marker_seeds_mpa4.tsv
# ============================================================

import argparse
import re
from pathlib import Path
from Bio import SeqIO


def sanitize(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--seed-tsv", required=True)
    args = parser.parse_args()

    input_fasta = Path(args.input).resolve()
    output_dir = Path(args.output_dir).resolve()
    seed_tsv = Path(args.seed_tsv).resolve()

    output_dir.mkdir(parents=True, exist_ok=True)
    seed_tsv.parent.mkdir(parents=True, exist_ok=True)

    records = list(SeqIO.parse(str(input_fasta), "fasta"))
    if not records:
        raise RuntimeError(f"No marker sequences found in {input_fasta}")

    seen = set()

    with seed_tsv.open("w", encoding="utf-8") as table:
        for idx, record in enumerate(records, start=1):
            marker_name = f"MPA_{idx:04d}_{sanitize(record.id)}"

            if marker_name in seen:
                raise RuntimeError(f"Duplicated marker name: {marker_name}")
            seen.add(marker_name)

            output_fasta = output_dir / f"{marker_name}.fna"

            record.id = marker_name
            record.name = marker_name
            record.description = ""

            SeqIO.write([record], str(output_fasta), "fasta")
            table.write(f"{marker_name}\t{output_fasta}\n")

    print(f"Markers written: {len(records)}")
    print(f"Seed table: {seed_tsv}")


if __name__ == "__main__":
    main()
