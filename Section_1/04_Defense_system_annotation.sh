#!/usr/bin/env bash
# ============================================================
# 04_Defense_system_annotation.sh
# Purpose : Annotate antiviral defense systems in the predicted
#           proteins of each isolate genome.
# Tool    : DefenseFinder v2.0.0
# Input   : /path/to/WGS/prodigal/protein/*.protein.fa
# Output  : /path/to/WGS/defensefinder/gembase/<sample>/
# ============================================================
set -euo pipefail

PROTEIN_DIR=/path/to/WGS/prodigal/protein
OUT_DIR=/path/to/WGS/defensefinder/gembase
WORKERS=100

mkdir -p "${OUT_DIR}"

for protein in "${PROTEIN_DIR}"/*.protein.fa; do
    sample=$(basename "${protein}" .protein.fa)

    defense-finder run \
        -o "${OUT_DIR}/${sample}" \
        -w "${WORKERS}" \
        --db-type gembase \
        "${protein}"
done
