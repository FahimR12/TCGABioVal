# =============================================================================
# 00_utils.R — Shared utility functions
# =============================================================================
# Loaded by every pipeline stage. Provides:
#   - Config loading and validation
#   - Path resolution (relative to project root)
#   - Logging helpers
#   - Common checks
# =============================================================================

suppressPackageStartupMessages({
  library(yaml)
})

# Null-coalescing operator (not in base R)
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

#' Load and validate the pipeline config
#'
#' @param config_path Path to config.yml (default: config/config.yml from
#'   project root)
#' @return Named list of config values
load_config <- function(config_path = NULL) {
  if (is.null(config_path)) {
    config_path <- file.path(find_project_root(), "config", "config.yml")
  }
  stopifnot(
    "Config file not found" = file.exists(config_path)
  )
  cfg <- yaml::read_yaml(config_path)
  validate_config(cfg)
  cfg
}

#' Minimal validation — check required top-level keys exist
validate_config <- function(cfg) {
  required_keys <- c("paths", "cohort", "qc", "normalisation",
                     "clinical", "feature_selection", "bn_learn")
  missing <- setdiff(required_keys, names(cfg))
  if (length(missing) > 0) {
    stop("Config is missing required sections: ",
         paste(missing, collapse = ", "))
  }
  # Validate GDCdata root exists
  gdc_root <- cfg$paths$gdcdata_root
  if (!dir.exists(gdc_root)) {
    warning("GDCdata root does not exist: ", gdc_root,
            "\n  Update paths.gdcdata_root in config.yml")
  }
  invisible(cfg)
}

# ---------------------------------------------------------------------------
# Project root detection
# ---------------------------------------------------------------------------

#' Walk up from the working directory until we find config/config.yml
#' (or a .git directory). Returns the project root as an absolute path.
find_project_root <- function(start = getwd()) {
  dir <- normalizePath(start, mustWork = FALSE)
  while (TRUE) {
    if (file.exists(file.path(dir, "config", "config.yml")) ||
        dir.exists(file.path(dir, ".git"))) {
      return(dir)
    }
    parent <- dirname(dir)
    if (parent == dir) {
      stop("Could not find project root. ",
           "Run scripts from within the TCGABioVal directory.")
    }
    dir <- parent
  }
}

# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

#' Resolve a path relative to the project root
#' @param ... Path components (passed to file.path)
project_path <- function(...) {
  file.path(find_project_root(), ...)
}

#' Build the cohort-specific GDCdata path
#' e.g. C:/Users/fahim/Desktop/scripts/BNPipeline/GDCdata/TCGA-GBM
gdcdata_cohort_path <- function(cfg) {
  file.path(cfg$paths$gdcdata_root, cfg$cohort$project)
}

#' Build a processed-data output path, creating the directory if needed
processed_path <- function(cfg, ...) {
  p <- project_path(cfg$paths$processed_dir, ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

#' Build a results output path, creating the directory if needed
results_path <- function(cfg, ...) {
  p <- project_path(cfg$paths$results_dir, ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

#' Simple timestamped logging to console
#' @param ... Message parts (concatenated with paste0)
#' @param level One of "INFO", "WARN", "ERROR"
log_msg <- function(..., level = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0(...)
  message(sprintf("[%s] %s: %s", ts, level, msg))
}

log_info  <- function(...) log_msg(..., level = "INFO")
log_warn  <- function(...) log_msg(..., level = "WARN")
log_error <- function(...) log_msg(..., level = "ERROR")

#' Log the start of a pipeline stage
log_stage_start <- function(stage_name) {
  log_info(strrep("=", 60))
  log_info("Starting: ", stage_name)
  log_info(strrep("=", 60))
}

#' Log stage completion with elapsed time
log_stage_end <- function(stage_name, start_time) {
  elapsed <- difftime(Sys.time(), start_time, units = "secs")
  log_info("Completed: ", stage_name,
           " (", round(as.numeric(elapsed), 1), "s)")
  log_info(strrep("-", 60))
}

# ---------------------------------------------------------------------------
# Common checks
# ---------------------------------------------------------------------------

#' Assert that required packages are installed
check_packages <- function(pkgs) {
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required packages: ", paste(missing, collapse = ", "),
         "\nInstall with: install.packages(c('",
         paste(missing, collapse = "', '"), "'))")
  }
}

#' Assert a data.frame has expected columns
assert_columns <- function(df, expected, df_name = "data") {
  missing <- setdiff(expected, colnames(df))
  if (length(missing) > 0) {
    stop(df_name, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

#' Save an R object with a log message
save_output <- function(obj, path, obj_name = deparse(substitute(obj))) {
  saveRDS(obj, path)
  log_info("Saved ", obj_name, " → ", basename(path),
           " (", format(file.size(path), big.mark = ","), " bytes)")
}
