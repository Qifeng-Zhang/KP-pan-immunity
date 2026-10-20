#!/usr/bin/env python3
# ============================================================
# 04_Build_reference_catalog.py
# Tool    : Biopython
# Input   : MAG directory
#             /path/to/Strain/thauera_bin_90
#           MAG-to-lineage table (99% ANI lineages), two columns:
#             MAG <tab> lineage99
#             data/thauera_MAG_lineages.tsv
# Output  : 03_reference/pseudochromosomes/<accession>.fna
#           03_reference/thauera_reference_catalog.tsv
#           03_reference/thauera_accession_map.tsv
# ============================================================

import argparse
import csv
import re
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mag-dir", required=True)
    parser.add_argument("--lineage-table", required=True)
    parser.add_argument("--pseudo-dir", required=True)
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--separator-length", type=int, default=100)
    parser.add_argument("--expected-mags", type=int, default=15)
    args = parser.parse_args()

    mag_dir = Path(args.mag_dir).resolve()
    lineage_table = Path(args.lineage_table).resolve()
    pseudo_dir = Path(args.pseudo_dir).resolve()
    catalog_path = Path(args.catalog).resolve()
    mapping_path = Path(args.mapping).resolve()

    pseudo_dir.mkdir(parents=True, exist_ok=True)
    catalog_path.parent.mkdir(parents=True, exist_ok=True)
    mapping_path.parent.mkdir(parents=True, exist_ok=True)

    lineages = {}
    with lineage_table.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            lineages[row["MAG"]] = row["lineage99"]

    if len(lineages) != args.expected_mags:
        raise RuntimeError(
            f"Expected {args.expected_mags} MAG entries, found {len(lineages)}"
        )

    catalog_rows = []
    mapping_rows = []

    for mag_name, lineage in lineages.items():
        mag_path = mag_dir / mag_name

        if not mag_path.is_file():
            raise FileNotFoundError(f"Missing MAG: {mag_path}")

        records = [
            record for record in SeqIO.parse(str(mag_path), "fasta")
            if len(record.seq) > 0
        ]

        if not records:
            raise RuntimeError(f"No sequence found in {mag_path}")

        biological_length = sum(len(record.seq) for record in records)
        separator = "N" * args.separator_length
        pseudo_sequence = separator.join(
            str(record.seq).upper() for record in records
        )

        mag_stem = mag_path.stem
        accession = safe_name(f"TAA_{lineage}_{mag_stem}")
        pseudo_path = pseudo_dir / f"{accession}.fna"

        pseudo_record = SeqRecord(
            Seq(pseudo_sequence),
            id=accession,
            name=accession,
            description=f"source_MAG={mag_name};lineage99={lineage}"
        )
        SeqIO.write([pseudo_record], str(pseudo_path), "fasta")

        catalog_rows.append({
            "Accession": accession,
            "Genus": "Thauera",
            "Species": "aminoaromatica",
            "Strain": mag_stem,
            "ChromosomeLen": biological_length,
            "SeqPath": str(pseudo_path),
            "GFF": "None"
        })

        mapping_rows.append({
            "Accession": accession,
            "MAG": mag_name,
            "lineage99": lineage,
            "n_contigs": len(records),
            "biological_length": biological_length,
            "pseudo_length": len(pseudo_sequence),
            "SeqPath": str(pseudo_path)
        })

    with catalog_path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=["Accession", "Genus", "Species", "Strain",
                        "ChromosomeLen", "SeqPath", "GFF"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(catalog_rows)

    with mapping_path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=["Accession", "MAG", "lineage99", "n_contigs",
                        "biological_length", "pseudo_length", "SeqPath"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(mapping_rows)

    print(f"Reference entries: {len(catalog_rows)}")
    print(f"Catalog: {catalog_path}")
    print(f"Mapping: {mapping_path}")


if __name__ == "__main__":
    main()
