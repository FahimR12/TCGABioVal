#!/bin/bash
# Run the TCGABioVal GBM preprocessing and BN learning pipeline.
#
# Stages run in sequence:
#   01  R/01_ingest.R         — read STAR count files from GDCdata (fixed input, no network)
#   02  R/02_qc.R             — MAD-based sample QC flagging + diagnostic plots
#   03  R/03_normalise.R      — gene filtering, TMM+logCPM, IDH covariate parsing
#   04  R/04_feature_select.R — top-variance genes ∪ curated EGFR reference set
#   05  R/05_bn_bnlearn.R     — Hill Climbing BN + bootstrap edge stability
#
# Submit (all stages):
#   sbatch /users/frazai/fastscratch/TCGABioVal/run_pipeline_gbm.sbatch
#
# Skip ingestion if 01_raw_counts.rds already exists:
#   sbatch --export=ALL,SKIP_INGEST=1 run_pipeline_gbm.sbatch
#
# Run a specific stage range (e.g. re-run from normalisation):
#   sbatch --export=ALL,FROM_STAGE=03 run_pipeline_gbm.sbatch
#
# Run a single stage:
#   sbatch --export=ALL,ONLY_STAGE=05 run_pipeline_gbm.sbatch
#
# Optional: tag outputs with a run identifier (written to results/):
#   sbatch --export=ALL,RUN_ID=pilot_v1 run_pipeline_gbm.sbatch

#SBATCH -J tcga-bioval-gbm
#SBATCH -p nodes
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH -D /users/frazai/fastscratch/TCGABioVal
#SBATCH -o /users/frazai/fastscratch/TCGABioVal/logs/%x-%j.out
#SBATCH -e /users/frazai/fastscratch/TCGABioVal/logs/%x-%j.err

set -euo pipefail

mkdir -p /users/frazai/fastscratch/TCGABioVal/logs

# =============================================================================
# Environment
# =============================================================================
module purge
module load r/4.5.2-gcc14.2.0

# If the R module is unavailable, fall back to the conda smc_r_env:
# source /users/frazai/.bashrc
# conda activate smc_r_env

# Prevent nested parallelism from over-subscribing cores
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Non-interactive plot backend (no display on HPC)
export MPLBACKEND=Agg

# =============================================================================
# Config — point to Barkla paths
# =============================================================================
CONFIG="config/config.yml"

# Override the GDCdata root for Barkla (different from Windows dev path).
# The pipeline reads this env var if set, overriding config.yml.
# Ensure this points to the same GDCdata used in BNPipeline on Barkla.
export TCGABIOVAL_GDCDATA_ROOT="/users/frazai/fastscratch/BNPipeline/GDCdata"

# Optional clinical TSV override (set if you have a pre-downloaded clinical file)
# export TCGABIOVAL_CLINICAL_TSV="/users/frazai/fastscratch/BNPipeline/data/clinical_gbm.csv"

# =============================================================================
# Run ID
# =============================================================================
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
echo "[pipeline] job=${SLURM_JOB_ID:-local} run_id=${RUN_ID} config=${CONFIG}"
echo "[pipeline] gdcdata=${TCGABIOVAL_GDCDATA_ROOT}"
echo "[pipeline] host=$(hostname) cpus=${SLURM_CPUS_PER_TASK:-8}"

# =============================================================================
# Stage selection
# =============================================================================
# Determine which stages to run. Modes (mutually exclusive, checked in order):
#   ONLY_STAGE  — run exactly one stage (e.g. ONLY_STAGE=05)
#   FROM_STAGE  — run from a stage onward (e.g. FROM_STAGE=03 runs 03 04 05)
#   SKIP_INGEST — run 02 03 04 05 (stage 01 already done)
#   default     — run all stages
if [[ -n "${ONLY_STAGE:-}" ]]; then
  STAGES="${ONLY_STAGE}"
  echo "[pipeline] mode: single stage ${ONLY_STAGE}"
elif [[ -n "${FROM_STAGE:-}" ]]; then
  # Build stage list from FROM_STAGE to 05
  all_stages=(01 02 03 04 05)
  STAGES=""
  in_range=0
  for s in "${all_stages[@]}"; do
    [[ "$s" == "${FROM_STAGE}" ]] && in_range=1
    [[ $in_range -eq 1 ]] && STAGES="${STAGES} ${s}"
  done
  STAGES="${STAGES## }"   # trim leading space
  echo "[pipeline] mode: from stage ${FROM_STAGE} — running: ${STAGES}"
elif [[ -n "${SKIP_INGEST:-}" ]]; then
  STAGES="02 03 04 05"
  echo "[pipeline] mode: skip ingest — running: ${STAGES}"
else
  STAGES="01 02 03 04 05"
  echo "[pipeline] mode: full pipeline — running: ${STAGES}"
fi

# =============================================================================
# Helper: run a pipeline stage via the master orchestrator
# =============================================================================
run_stage() {
  local stage_id="$1"
  local stage_label="$2"
  local t_start=$(date +%s)

  echo ""
  echo "[pipeline] ── Stage ${stage_id}: ${stage_label} $(printf '─%.0s' {1..40})"
  echo "[pipeline] start: $(date '+%Y-%m-%d %H:%M:%S')"

  Rscript scripts/run_pipeline.R "${stage_id}"
  local exit_code=$?

  local t_end=$(date +%s)
  local elapsed=$(( t_end - t_start ))
  local mins=$(( elapsed / 60 ))
  local secs=$(( elapsed % 60 ))

  if [[ $exit_code -ne 0 ]]; then
    echo "[pipeline] ERROR: Stage ${stage_id} failed after ${mins}m${secs}s (exit ${exit_code})"
    exit $exit_code
  fi

  echo "[pipeline] done:  Stage ${stage_id} in ${mins}m${secs}s"
}

# =============================================================================
# Stage labels (for readable log output)
# =============================================================================
declare -A STAGE_LABELS
STAGE_LABELS[01]="Ingest STAR count files"
STAGE_LABELS[02]="Sample QC (MAD-based flagging)"
STAGE_LABELS[03]="Normalise + IDH covariate"
STAGE_LABELS[04]="Feature selection (variance + EGFR)"
STAGE_LABELS[05]="BN learning — bnlearn Hill Climbing + bootstrap"

# =============================================================================
# Run selected stages
# =============================================================================
PIPELINE_START=$(date +%s)

for stage in $STAGES; do
  run_stage "${stage}" "${STAGE_LABELS[$stage]:-unknown}"
done

# =============================================================================
# Summary
# =============================================================================
PIPELINE_END=$(date +%s)
TOTAL_ELAPSED=$(( PIPELINE_END - PIPELINE_START ))
TOTAL_MINS=$(( TOTAL_ELAPSED / 60 ))
TOTAL_SECS=$(( TOTAL_ELAPSED % 60 ))

echo ""
echo "[pipeline] ════════════════════════════════════════════════════════"
echo "[pipeline] Pipeline complete — ${TOTAL_MINS}m${TOTAL_SECS}s"
echo "[pipeline] Run ID: ${RUN_ID}"
echo "[pipeline] Outputs: data/processed/ and results/"
echo "[pipeline] ════════════════════════════════════════════════════════"
