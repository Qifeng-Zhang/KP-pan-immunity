#!/usr/bin/env bash
# ============================================================
# 02_Marker_prevalence_QC.sh
# Tool    : BLASTN
# Input   : Marker FASTA
#             /path/to/Strain/mpa4_marker/t__SGB12599.fna
#           MAGs of the species-level cluster
#             /path/to/Strain/thauera_bin_90/*.fa
# Output  : /path/to/Strain/chronostrain_thauera/02_marker_qc/blast/<MAG>.tsv
# ============================================================
set -euo pipefail

ROOT=/path/to/Strain
WORK=${ROOT}/chronostrain_thauera

MAG_DIR=${ROOT}/thauera_bin_90
MARKERS=${ROOT}/mpa4_marker/t__SGB12599.fna
BLAST_DIR=${WORK}/02_marker_qc/blast

mkdir -p "${BLAST_DIR}"

for mag in "${MAG_DIR}"/*.fa; do
    mag_id=$(basename "${mag}" .fa)

    blastn \
      -task blastn \
      -query "${MARKERS}" \
      -subject "${mag}" \
      -perc_identity 90 \
      -qcov_hsp_perc 80 \
      -dust no \
      -max_hsps 1 \
      -outfmt '6 qseqid sseqid pident length qlen qcovhsp evalue bitscore' \
      > "${BLAST_DIR}/${mag_id}.tsv"
done
