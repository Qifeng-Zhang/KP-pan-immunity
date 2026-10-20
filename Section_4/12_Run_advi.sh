#!/usr/bin/env bash
# ============================================================
# 12_Run_advi.sh
# Tool    : ChronoStrain
# Input   : ChronoStrain configuration
#             config/chronostrain.ini
#           Clustered database entries
#             06_clusters/clusters_U80_t0998.txt
#           Mate-corrected ChronoStrain read table
#             07_samples/04_reads_pairfixed/filtered_DC_chronostrain_reads_pairfixed.csv
# Output  : 07_samples/05_inference/U80_t0998/
#             posterior samples and ELBO trace
# ============================================================
set -euo pipefail


ROOT=/path/to/Strain
WORK=${ROOT}/chronostrain_thauera

CONFIG=${WORK}/config/chronostrain.ini
CLUSTERS=${WORK}/06_clusters/clusters_U80_t0998.txt
READS=${WORK}/07_samples/04_reads_pairfixed/filtered_DC_chronostrain_reads_pairfixed.csv
OUTDIR=${WORK}/07_samples/05_inference/U80_t0998

mkdir -p "${OUTDIR}"

chronostrain -c "${CONFIG}" advi \
  -r "${READS}" \
  -o "${OUTDIR}" \
  -s "${CLUSTERS}" \
  --plot-elbo
