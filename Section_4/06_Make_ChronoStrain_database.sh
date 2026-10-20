#!/usr/bin/env bash
# ============================================================
# 06_Make_ChronoStrain_database.sh
# Tool    : ChronoStrain
# Input   : Conserved-marker seed table
#             04_marker_subsets/marker_seeds_U80.tsv
#           Reference catalog (pseudo-chromosomes)
#             03_reference/thauera_reference_catalog.tsv
# Output  : config/chronostrain.ini
#           05_database/thauera_mpa4.json
#           06_clusters/clusters_U80_t0998.txt
# ============================================================
set -euo pipefail

ROOT=/path/to/Strain
WORK=${ROOT}/chronostrain_thauera

SEEDS=${WORK}/04_marker_subsets/marker_seeds_U80.tsv
CATALOG=${WORK}/03_reference/thauera_reference_catalog.tsv

DB_JSON=${WORK}/05_database/thauera_mpa4.json
DB_DATA=${WORK}/05_database/data
BLAST_TMP=${WORK}/05_database/blast_tmp
CACHE_DIR=${WORK}/cache
CONFIG=${WORK}/config/chronostrain.ini

CLUSTER_THRESHOLD=0.998

mkdir -p "${DB_DATA}" "${BLAST_TMP}" "${CACHE_DIR}" \
         "${WORK}/06_clusters" "${WORK}/config" "${WORK}/logs"

# --- Write configuration ---------------------------------------------------
cat > "${CONFIG}" <<EOF
[Model]
USE_QUALITY_SCORES=True
FRAG_LEN_NB_N=11186530.399297457
FRAG_LEN_NB_P=0.9999865398828846
EXTRA_STRAIN=False
SICS_DOF_1=0.1
SICS_SCALE_1=100.0
SICS_DOF=0.1
SICS_SCALE=1.0
INV_TEMPERATURE=1
NUM_CORES=80
SPARSE_MATRICES=True
INSERTION_LL_1=-12.78589
DELETION_LL_1=-12.18627
INSERTION_LL_2=-12.56274
DELETION_LL_2=-12.22627
MIN_OVERLAP_RATIO=0.5
CACHE_DIR=${CACHE_DIR}

[Entrez]
ENABLED=False
EMAIL=

[Engine]
PRNG_Key=314159
DTYPE=float32

[Database]
DB_PARSER_CLASS=chronostrain.database.JSONParser
DB_DATA_DIR=${DB_DATA}

[Database.ParserArgs]
ENTRIES_FILE=${DB_JSON}
marker_max_len=1000000

[ExternalTools]
PAIRWISE_ALN_BACKEND=bwa
PAIRWISE_ALN_BAM=True
EOF

# --- Build database --------------------------------------------------------
chronostrain -c "${CONFIG}" make-db \
  -m "${SEEDS}" \
  -r "${CATALOG}" \
  -b thauera_mpa4 \
  -bd "${BLAST_TMP}" \
  -o "${DB_JSON}" \
  --threads 64

# --- Cluster database entries into strain lineages -------------------------
chronostrain -c "${CONFIG}" cluster-db \
  -i "${DB_JSON}" \
  -o "${WORK}/06_clusters/clusters_U80_t0998.txt" \
  -t "${CLUSTER_THRESHOLD}"
