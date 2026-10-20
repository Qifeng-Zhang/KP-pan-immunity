#!/usr/bin/env bash
# ============================================================
# 03_Gene_prediction.sh
# Purpose : Predict open reading frames in the assembled genomes.
# Tool    : Prodigal v2.6.3 (single-genome mode)
# Input   : /path/to/WGS/assemblies/*.fasta
# Output  : /path/to/WGS/prodigal/orf/*.gff
#           /path/to/WGS/prodigal/protein/*.protein.fa
#           /path/to/WGS/prodigal/DNA/*.DNA.fa
# ============================================================
set -euo pipefail

GENOME_DIR=/path/to/WGS/assemblies
OUT_DIR=/path/to/WGS/prodigal

mkdir -p "${OUT_DIR}"/{orf,protein,DNA}

for genome in "${GENOME_DIR}"/*.fasta; do
    sample=$(basename "${genome}" .fasta)

    prodigal \
        -i "${genome}" \
        -o "${OUT_DIR}/orf/${sample}.gff" \
        -a "${OUT_DIR}/protein/${sample}.protein.fa" \
        -d "${OUT_DIR}/DNA/${sample}.DNA.fa" \
        -f gff \
        -p single
done
