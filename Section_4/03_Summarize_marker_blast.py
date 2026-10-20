#!/usr/bin/env python3
# ============================================================
# 03_Summarize_marker_blast.py
# Tool    : Biopython
# Input   : Marker FASTA
#             /path/to/Strain/mpa4_marker/t__SGB12599.fna
#           BLAST output directory
#             /path/to/Strain/chronostrain_thauera/02_marker_qc/blast
# Output  : marker_coverage_by_MAG.tsv
#           marker_prevalence_across_MAGs.tsv
#           marker_presence_matrix.tsv
# ============================================================

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--marker-fasta", required=True)
    parser.add_argument("--blast-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    marker_fasta = Path(args.marker_fasta).resolve()
    blast_dir = Path(args.blast_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    markers = [record.id for record in SeqIO.parse(str(marker_fasta), "fasta")]
    marker_set = set(markers)

    blast_files = sorted(blast_dir.glob("*.tsv"))
    if not blast_files:
        raise RuntimeError(f"No BLAST TSV files found in {blast_dir}")

    mags = [path.stem for path in blast_files]

    hits = {mag: defaultdict(set) for mag in mags}

    for path in blast_files:
        mag = path.stem

        with path.open(encoding="utf-8") as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if len(row) < 2:
                    continue

                marker, subject = row[0], row[1]

                if marker in marker_set:
                    hits[mag][marker].add(subject)

    # 1. Per-MAG coverage
    with (output_dir / "marker_coverage_by_MAG.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "MAG", "total_markers", "markers_present",
            "marker_fraction", "markers_with_multiple_subject_hits"
        ])

        for mag in mags:
            present = {m for m in markers if len(hits[mag][m]) > 0}
            multi = {m for m in markers if len(hits[mag][m]) > 1}

            writer.writerow([
                mag, len(markers), len(present),
                len(present) / len(markers) if markers else 0,
                len(multi)
            ])

    # 2. Per-marker prevalence
    with (output_dir / "marker_prevalence_across_MAGs.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow([
            "marker", "MAGs_present", "MAG_fraction",
            "MAGs_with_multiple_subject_hits"
        ])

        for marker in markers:
            n_present = sum(len(hits[mag][marker]) > 0 for mag in mags)
            n_multi = sum(len(hits[mag][marker]) > 1 for mag in mags)

            writer.writerow([marker, n_present, n_present / len(mags), n_multi])

    # 3. Presence/absence matrix
    with (output_dir / "marker_presence_matrix.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["marker"] + mags)

        for marker in markers:
            writer.writerow(
                [marker] + [int(len(hits[mag][marker]) > 0) for mag in mags]
            )

    print(f"Markers: {len(markers)}")
    print(f"MAGs: {len(mags)}")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
