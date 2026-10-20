#!/usr/bin/env bash
# ============================================================
# 10_Filter_reads.sh
# Tool    : ChronoStrain
# Input   : ChronoStrain read table
#             07_samples/02_reads_table/DC_chronostrain_reads.tsv
#           ChronoStrain configuration
#             config/chronostrain.ini
# Output  : 07_samples/03_filter/filtered_DC_chronostrain_reads.tsv
# ============================================================
set -euo pipefail

ROOT=/path/to/Strain
WORK=${ROOT}/chronostrain_thauera

CONFIG=${WORK}/config/chronostrain.ini
READS=${WORK}/07_samples/02_reads_table/DC_chronostrain_reads.tsv
OUTDIR=${WORK}/07_samples/03_filter

mkdir -p "${OUTDIR}"

chronostrain -c "${CONFIG}" filter \
  -r "${READS}" \
  -o "${OUTDIR}"
