# =============================================================================
# run_pipeline.R — Master orchestrator
# =============================================================================
# Run the full pipeline (or individual stages) from the project root.
#
# Usage:
#   Rscript scripts/run_pipeline.R              # run all stages
#   Rscript scripts/run_pipeline.R 01           # run stage 01 only
#   Rscript scripts/run_pipeline.R 02 03        # run stages 02 and 03
#   Rscript scripts/run_pipeline.R --from 03    # run stages 03 onward
# =============================================================================

# Ensure working directory is project root
# (handles being called from scripts/ or project root)
if (basename(getwd()) == "scripts") setwd("..")

source(file.path("R", "00_utils.R"))

# =============================================================================
# Stage registry
# =============================================================================

stages <- list(
  "01" = list(
    name   = "Ingest TCGA data",
    script = file.path("R", "01_ingest.R"),
    func   = "ingest_tcga"
  ),
  "02" = list(
    name   = "Sample QC",
    script = file.path("R", "02_qc.R"),
    func   = "run_qc"
  ),
  "03" = list(
    name   = "Normalise + IDH covariate",
    script = file.path("R", "03_normalise.R"),
    func   = "run_normalise"
  ),
  "04" = list(
    name   = "Feature selection",
    script = file.path("R", "04_feature_select.R"),
    func   = "run_feature_select"
  ),
  "05" = list(
    name   = "BN learning (bnlearn HC)",
    script = file.path("R", "05_bn_bnlearn.R"),
    func   = "run_bn_bnlearn"
  )
)

# =============================================================================
# CLI argument parsing
# =============================================================================

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) == 0) {
    return(names(stages))  # run all
  }

  if (args[1] == "--from") {
    start <- args[2]
    all_ids <- names(stages)
    idx <- match(start, all_ids)
    if (is.na(idx)) stop("Unknown stage: ", start)
    return(all_ids[idx:length(all_ids)])
  }

  # Specific stages
  for (s in args) {
    if (!(s %in% names(stages))) stop("Unknown stage: ", s)
  }
  args
}

# =============================================================================
# Run
# =============================================================================

main <- function() {
  cfg <- load_config()
  run_stages <- parse_args()

  log_info(strrep("=", 60))
  log_info("TCGABioVal Pipeline")
  log_info("Cohort: ", cfg$cohort$project)
  log_info("Stages: ", paste(run_stages, collapse = ", "))
  log_info(strrep("=", 60))

  pipeline_start <- Sys.time()

  for (stage_id in run_stages) {
    stage <- stages[[stage_id]]
    log_info("")
    log_info(">>> Stage ", stage_id, ": ", stage$name)

    # Source the stage script (defines the function)
    source(stage$script, local = TRUE)

    # Call the stage function with the config
    result <- do.call(stage$func, list(cfg = cfg))
  }

  elapsed <- difftime(Sys.time(), pipeline_start, units = "mins")
  log_info("")
  log_info(strrep("=", 60))
  log_info("Pipeline complete — ", round(as.numeric(elapsed), 1), " minutes")
  log_info(strrep("=", 60))
}

if (sys.nframe() == 0) {
  main()
}
