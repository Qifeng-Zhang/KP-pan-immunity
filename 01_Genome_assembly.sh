#!/usr/bin/env bash
# ============================================================
# 01_Genome_assembly.sh
# Purpose : Assemble whole-genome sequencing reads of the bacterial
#           isolates (one assembly per isolate).
# Tool    : SPAdes 4.1.0 (--isolate mode)
# Input   : Paired-end reads per isolate
#             /path/to/WGS/raw_reads/<sample>/<sample>_1.fq.gz
#             /path/to/WGS/raw_reads/<sample>/<sample>_2.fq.gz
# Output  : /path/to/WGS/assemblies/<sample>/
# Manuscript: Supplementary Note 1 (isolation and identification of strains)
# ============================================================
set -euo pipefail

RAW_DIR=/path/to/WGS/raw_reads
OUT_DIR=/path/to/WGS/assemblies
THREADS=30

mkdir -p "${OUT_DIR}"

for read1 in "${RAW_DIR}"/*/*_1.fq.gz; do
    sample=$(basename "${read1}" _1.fq.gz)
    read2="$(dirname "${read1}")/${sample}_2.fq.gz"

    spades.py \
        --isolate \
        -t "${THREADS}" \
        -1 "${read1}" \
        -2 "${read2}" \
        -o "${OUT_DIR}/${sample}"
done
