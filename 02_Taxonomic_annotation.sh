#!/usr/bin/env bash
# ============================================================
# 02_Taxonomic_annotation.sh
# Purpose : Species-level taxonomic assignment of the assembled
#           isolate genomes.
# Tool    : GTDB-Tk v2.3.2 (classify_wf)
# Input   : /path/to/WGS/assemblies/*.fasta
# Output  : /path/to/WGS/gtdbtk
# ============================================================
set -euo pipefail

GENOME_DIR=/path/to/WGS/assemblies
OUT_DIR=/path/to/WGS/gtdbtk
THREADS=80

mkdir -p "${OUT_DIR}"

gtdbtk classify_wf \
    --genome_dir "${GENOME_DIR}" \
    --out_dir "${OUT_DIR}" \
    --extension fasta \
    --skip_ani_screen \
    --cpus "${THREADS}"
