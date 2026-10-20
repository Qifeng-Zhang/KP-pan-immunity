#!/usr/bin/env bash
# ============================================================
# 01_Full_length_16S.sh
# Purpose : Full-length 16S rRNA gene amplicon profiling of the
#           synthetic communities (PacBio CCS reads).
# Tool    : QIIME 2 2024.10
# Input   : CCS reads and a QIIME 2 single-end manifest
#             /path/to/16S/manifest_single.txt
# Output  : /path/to/16S/qza/         (QIIME 2 artifacts)
#           /path/to/16S/exported/tables  (TSV feature tables)
# ============================================================
set -euo pipefail

# mamba activate qiime2-amplicon-2024.10

WORK_DIR=/path/to/16S
MANIFEST=${WORK_DIR}/manifest_single.txt

QZA_DIR=${WORK_DIR}/qza
EXPORT_DIR=${WORK_DIR}/exported

THREADS=80
RAREFY_DEPTH=10000

mkdir -p "${QZA_DIR}" "${EXPORT_DIR}"

# --- Import single-end sequences -----------------------------------------
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "${MANIFEST}" \
  --output-path "${QZA_DIR}/demux.qza" \
  --input-format SingleEndFastqManifestPhred33V2

# --- Sequencing summary -----------------
qiime demux summarize \
  --i-data "${QZA_DIR}/demux.qza" \
  --o-visualization "${QZA_DIR}/demux.qzv"

# --- DADA2 denoising ------------------------------------------
qiime dada2 denoise-ccs \
  --i-demultiplexed-seqs "${QZA_DIR}/demux.qza" \
  --p-trunc-len 0 \
  --p-trim-left 0 \
  --p-max-ee 2 \
  --p-trunc-q 2 \
  --p-min-len 1000 \
  --p-max-len 1600 \
  --p-min-fold-parent-over-abundance 3.5 \
  --p-n-threads "${THREADS}" \
  --o-representative-sequences "${QZA_DIR}/rep-seqs_dada2.qza" \
  --o-table "${QZA_DIR}/table_dada2.qza" \
  --o-denoising-stats "${QZA_DIR}/stats_dada2.qza"

qiime metadata tabulate \
  --m-input-file "${QZA_DIR}/stats_dada2.qza" \
  --o-visualization "${QZA_DIR}/stats_dada2.qzv"

# --- Filter low-abundance features ---------------------------------------
qiime feature-table filter-features \
  --i-table "${QZA_DIR}/table_dada2.qza" \
  --p-min-samples 1 \
  --p-min-frequency 500 \
  --o-filtered-table "${QZA_DIR}/table_filtered.qza"

qiime feature-table filter-seqs \
  --i-data "${QZA_DIR}/rep-seqs_dada2.qza" \
  --i-table "${QZA_DIR}/table_filtered.qza" \
  --o-filtered-data "${QZA_DIR}/rep-seqs_filtered.qza"

# --- Export filtered table -----------------------------------------------
qiime tools export \
  --input-path "${QZA_DIR}/table_filtered.qza" \
  --output-path "${EXPORT_DIR}/table_filtered"

biom convert \
  -i "${EXPORT_DIR}/table_filtered/feature-table.biom" \
  -o "${EXPORT_DIR}/table_filtered/feature-table_filtered.txt" \
  --to-tsv

# --- Rarefy to an even sequencing depth ----------------------------------
qiime feature-table rarefy \
  --i-table "${QZA_DIR}/table_filtered.qza" \
  --p-sampling-depth "${RAREFY_DEPTH}" \
  --o-rarefied-table "${QZA_DIR}/table_rarefied.qza"

qiime feature-table filter-seqs \
  --i-data "${QZA_DIR}/rep-seqs_filtered.qza" \
  --i-table "${QZA_DIR}/table_rarefied.qza" \
  --o-filtered-data "${QZA_DIR}/rep-seqs_rarefied.qza"

# --- Export rarefied table -----------------------------------------------
qiime tools export \
  --input-path "${QZA_DIR}/table_rarefied.qza" \
  --output-path "${EXPORT_DIR}/table_rarefied"

biom convert \
  -i "${EXPORT_DIR}/table_rarefied/feature-table.biom" \
  -o "${EXPORT_DIR}/table_rarefied/feature-table_rarefied.txt" \
  --to-tsv
