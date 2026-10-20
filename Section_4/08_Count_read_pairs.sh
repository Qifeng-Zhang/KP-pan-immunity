#!/usr/bin/env bash
# ============================================================
# 08_Count_read_pairs.sh
# Tool    : standard shell utilities
# Input   : Time-series manifest
#             07_samples/00_manifest/DC_read_pairs.tsv
# Output  : 07_samples/01_read_counts/DC_read_counts.tsv
#             headerless, two columns: sample <tab> read_pairs
# ============================================================
set -euo pipefail

ROOT=/path/to/Strain
WORK=${ROOT}/chronostrain_thauera

MANIFEST=${WORK}/07_samples/00_manifest/DC_read_pairs.tsv
COUNTS_DIR=${WORK}/07_samples/01_read_counts
COUNTS=${COUNTS_DIR}/DC_read_counts.tsv

mkdir -p "${COUNTS_DIR}"
: > "${COUNTS}"

tail -n +2 "${MANIFEST}" | awk -F'\t' '{print $1 "\t" $4}' | while IFS=$'\t' read -r sample r1; do
    if [[ "${r1}" == *.gz ]]; then
        n_lines=$(zcat "${r1}" | wc -l)
    else
        n_lines=$(wc -l < "${r1}")
    fi

    n_pairs=$((n_lines / 4))
    printf '%s\t%s\n' "${sample}" "${n_pairs}" >> "${COUNTS}"
done

column -t -s $'\t' "${COUNTS}"
