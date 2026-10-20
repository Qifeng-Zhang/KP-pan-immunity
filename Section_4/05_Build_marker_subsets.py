#!/usr/bin/env python3
# ============================================================
# 05_Build_marker_subsets.py
# Tool    : standard library
# Input   : Marker prevalence table
#             02_marker_qc/marker_prevalence_across_MAGs.tsv
#           Full seed table
#             01_markers/marker_seeds_mpa4.tsv
# Output  : 04_marker_subsets/markers_<label>.ids
#           04_marker_subsets/marker_seeds_<label>.tsv
#           04_marker_subsets/marker_id_map_<label>.tsv
# ============================================================

import argparse
import csv
import re
from pathlib import Path


def remove_mpa_prefix(marker_id: str) -> str:
    """MPA_0001_SGB123_markerA -> SGB123_markerA"""
    return re.sub(r"^MPA_[0-9]+_", "", marker_id, flags=re.IGNORECASE)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prevalence", required=True)
    parser.add_argument("--seed-table", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--label", default="U80")
    parser.add_argument("--min-presence", type=int, default=14)
    args = parser.parse_args()

    prevalence_path = Path(args.prevalence).resolve()
    seed_path = Path(args.seed_table).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Read marker prevalence -------------------------------------------
    prevalence = {}

    with prevalence_path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"marker", "MAGs_present"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Missing columns in prevalence table: {sorted(missing)}"
            )

        for row in reader:
            marker = row["marker"].strip()
            if marker in prevalence:
                raise ValueError(f"Duplicated marker in prevalence table: {marker}")
            prevalence[marker] = int(row["MAGs_present"])

    # --- Read seed table ---------------------------------------------------
    seed_rows = []

    with seed_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 2:
                raise ValueError(
                    f"Seed table line {line_number} has fewer than "
                    f"two columns: {line}"
                )

            full_id = fields[0].strip()
            seed_rows.append({
                "full_id": full_id,
                "raw_id": remove_mpa_prefix(full_id),
                "line": line
            })

    # --- Match prevalence IDs against seed IDs -----------------------------
    def find_seed(marker_id: str):
        matches = [
            row for row in seed_rows
            if marker_id == row["full_id"] or marker_id == row["raw_id"]
        ]

        if len(matches) == 0:
            raise KeyError(
                f"No seed-table entry matched prevalence marker: {marker_id}"
            )

        if len(matches) > 1:
            raise ValueError(
                f"Multiple seed-table entries matched {marker_id}: "
                f"{[row['full_id'] for row in matches]}"
            )

        return matches[0]

    wanted_ids = [
        marker for marker, n_present in prevalence.items()
        if n_present >= args.min_presence
    ]

    selected = []
    seen_full_ids = set()

    for marker_id in wanted_ids:
        seed = find_seed(marker_id)

        if seed["full_id"] in seen_full_ids:
            raise ValueError(f"Duplicated selected seed: {seed['full_id']}")

        seen_full_ids.add(seed["full_id"])
        selected.append((marker_id, seed))

    ids_path = output_dir / f"markers_{args.label}.ids"
    seed_out = output_dir / f"marker_seeds_{args.label}.tsv"
    map_out = output_dir / f"marker_id_map_{args.label}.tsv"

    with ids_path.open("w", encoding="utf-8") as out:
        for original_id, _ in selected:
            out.write(f"{original_id}\n")

    with seed_out.open("w", encoding="utf-8") as out:
        for _, seed in selected:
            out.write(seed["line"] + "\n")

    with map_out.open("w", encoding="utf-8") as out:
        out.write("prevalence_id\tseed_id\n")
        for original_id, seed in selected:
            out.write(f"{original_id}\t{seed['full_id']}\n")

    print(f"{args.label}: selected {len(selected)} markers")
    print(f"  IDs:   {ids_path}")
    print(f"  Seeds: {seed_out}")
    print(f"  Map:   {map_out}")


if __name__ == "__main__":
    main()
