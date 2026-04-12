# =============================================================================
# 03_normalise.R — Gene filtering, normalisation, and IDH covariate handling
# =============================================================================
# Takes raw counts + QC flags, filters low-expression genes, normalises with
# edgeR (TMM + logCPM), and parses the IDH mutation status covariate.
#
# Inputs:
#   - data/processed/01_raw_counts.rds
#   - data/processed/01_clinical.rds
#   - data/processed/02_qc_metrics.rds
#
# Outputs:
#   - data/processed/03_expr_normalised.rds  : normalised expression matrix
#   - data/processed/03_clinical_final.rds   : clinical with IDH + QC merged
#   - data/processed/03_dge_object.rds       : edgeR DGEList (for diagnostics)
#   - data/processed/03_normalise_summary.txt
# =============================================================================

source(file.path("R", "00_utils.R"))

check_packages(c("edgeR", "limma"))

suppressPackageStartupMessages({
  library(edgeR)
})

# =============================================================================
# Gene filtering
# =============================================================================

#' Filter low-expression genes
#'
#' Keeps genes with CPM >= min_cpm in at least (min_sample_fraction × n_samples)
#' samples. This is the standard edgeR filterByExpr logic, but explicit for
#' transparency.
#'
#' @param counts Raw count matrix (genes × samples)
#' @param cfg    Pipeline config
#' @return Logical vector (TRUE = keep) with names = gene IDs
filter_genes <- function(counts, cfg) {
  filt_cfg <- cfg$normalisation$gene_filter
  min_cpm <- filt_cfg$min_cpm
  min_frac <- filt_cfg$min_sample_fraction
  min_samples <- ceiling(min_frac * ncol(counts))

  log_info("Gene filter: CPM >= ", min_cpm,
           " in >= ", min_frac * 100, "% of samples",
           " (n >= ", min_samples, ")")

  # Compute CPM from raw counts (no normalisation yet)
  lib_sizes <- colSums(counts)
  cpm_mat <- t(t(counts) / lib_sizes) * 1e6

  keep <- rowSums(cpm_mat >= min_cpm) >= min_samples

  log_info("Genes: ", nrow(counts), " → ", sum(keep),
           " (removed ", sum(!keep), " low-expression genes)")

  keep
}

# =============================================================================
# Normalisation
# =============================================================================

#' Normalise expression using edgeR
#'
#' @param counts Filtered count matrix
#' @param cfg    Pipeline config
#' @return list with:
#'   - expr: normalised expression matrix (genes × samples)
#'   - dge:  the DGEList object (for downstream diagnostics)
normalise_expression <- function(counts, cfg) {
  norm_cfg <- cfg$normalisation

  log_info("Normalisation method: ", norm_cfg$method)
  log_info("Transform: ", norm_cfg$transform)

  # Build DGEList
  dge <- DGEList(counts = counts)

  # Calculate normalisation factors
  dge <- calcNormFactors(dge, method = norm_cfg$method)

  log_info("Norm factors — range: [",
           signif(min(dge$samples$norm.factors), 3), ", ",
           signif(max(dge$samples$norm.factors), 3), "]",
           ", median: ", signif(median(dge$samples$norm.factors), 3))

  # Apply transform
  if (norm_cfg$transform == "logCPM") {
    expr <- cpm(dge, log = TRUE, prior.count = 1)
    log_info("logCPM: prior.count = 1")
  } else if (norm_cfg$transform == "CPM") {
    expr <- cpm(dge, log = FALSE)
  } else {
    stop("Unknown transform: ", norm_cfg$transform)
  }

  log_info("Expression matrix: ", nrow(expr), " genes × ", ncol(expr),
           " samples")

  list(expr = expr, dge = dge)
}

# =============================================================================
# IDH covariate parsing
# =============================================================================

#' Parse IDH mutation status from clinical metadata
#'
#' Implements the column detection priority cascade and regex normalisation
#' from our previous discussion. Key points:
#'   - Tries multiple column names in priority order
#'   - Guards against "unmutated" matching the "mut" pattern
#'   - Distinguishes genuinely untested from confirmed wildtype
#'
#' @param clinical data.frame of clinical metadata
#' @param idh_cfg  The clinical.covariates.idh_status section of config
#' @return character vector: "mutant", "wildtype", or NA (length = nrow)
parse_idh_status <- function(clinical, idh_cfg) {
  log_info("Parsing IDH mutation status")

  # Step 1: Find the source column
  source_col <- NULL
  for (col_name in idh_cfg$source_columns) {
    if (col_name %in% colnames(clinical)) {
      source_col <- col_name
      break
    }
  }

  if (is.null(source_col)) {
    log_warn("No IDH status column found in clinical data. ",
             "Tried: ", paste(idh_cfg$source_columns, collapse = ", "))
    log_warn("Available columns: ",
             paste(head(colnames(clinical), 30), collapse = ", "))
    return(rep(NA_character_, nrow(clinical)))
  }

  log_info("Using IDH column: '", source_col, "'")
  raw_values <- tolower(trimws(as.character(clinical[[source_col]])))
  log_info("Raw value distribution: ",
           paste(names(table(raw_values, useNA = "ifany")), "=",
                 table(raw_values, useNA = "ifany"), collapse = ", "))

  # Step 2: Classify each value
  result <- rep(NA_character_, length(raw_values))

  # First, check exclusion patterns (guard against "unmutated" etc.)
  is_excluded <- rep(FALSE, length(raw_values))
  for (pat in idh_cfg$exclude_patterns) {
    is_excluded <- is_excluded | grepl(pat, raw_values, ignore.case = TRUE)
  }
  if (any(is_excluded)) {
    log_info("Exclusion guard matched ", sum(is_excluded),
             " values (e.g. 'unmutated')")
  }

  # Classify mutant (but only if not excluded)
  for (pat in idh_cfg$mutant_patterns) {
    matches <- grepl(pat, raw_values, ignore.case = TRUE) & !is_excluded
    result[matches & is.na(result)] <- "mutant"
  }

  # Classify wildtype (excluded values that match mutant patterns are
  # reclassified as wildtype — "unmutated" = wildtype)
  for (pat in idh_cfg$wildtype_patterns) {
    matches <- grepl(pat, raw_values, ignore.case = TRUE)
    result[matches & is.na(result)] <- "wildtype"
  }
  # Excluded values (e.g. "unmutated") → wildtype
  result[is_excluded & is.na(result)] <- "wildtype"

  # Step 3: Handle NAs
  # In TCGA-GBM, NA often means "not tested" (older samples predate routine
  # IDH testing), NOT "confirmed wildtype"
  if (isTRUE(idh_cfg$treat_na_as_unknown)) {
    log_info("NAs treated as genuinely unknown (not assumed wildtype)")
  }

  n_mut <- sum(result == "mutant", na.rm = TRUE)
  n_wt  <- sum(result == "wildtype", na.rm = TRUE)
  n_na  <- sum(is.na(result))

  log_info("IDH status: mutant = ", n_mut, ", wildtype = ", n_wt,
           ", unknown/NA = ", n_na)

  if (n_mut < 15) {
    log_warn("Only ", n_mut, " IDH-mutant samples — ",
             "class imbalance may affect downstream BN modelling. ",
             "Consider stratification or excluding IDH as a BN node.")
  }

  result
}

# =============================================================================
# Main normalisation function
# =============================================================================

run_normalise <- function(cfg) {
  log_stage_start("03 — Normalise")
  start_time <- Sys.time()

  # Load data
  counts   <- readRDS(processed_path(cfg, "01_raw_counts.rds"))
  clinical <- readRDS(processed_path(cfg, "01_clinical.rds"))
  qc       <- readRDS(processed_path(cfg, "02_qc_metrics.rds"))

  # -------------------------------------------------------------------------
  # Optionally exclude flagged samples (manual decision point)
  # -------------------------------------------------------------------------
  # By default, we keep all samples. Uncomment the next block to drop flagged
  # samples automatically. The conservative approach is to review QC plots
  # first and edit this section based on what you see.
  #
  # drop_samples <- qc$sample_id[qc$flag_any]
  # if (length(drop_samples) > 0) {
  #   log_info("Dropping ", length(drop_samples), " QC-flagged samples")
  #   keep <- !(colnames(counts) %in% drop_samples)
  #   counts   <- counts[, keep]
  #   clinical <- clinical[keep, ]
  # }

  log_info("Proceeding with ", ncol(counts), " samples (flagged samples retained)")

  # -------------------------------------------------------------------------
  # Filter genes
  # -------------------------------------------------------------------------
  keep_genes <- filter_genes(counts, cfg)
  counts_filt <- counts[keep_genes, ]

  # -------------------------------------------------------------------------
  # Normalise
  # -------------------------------------------------------------------------
  norm_result <- normalise_expression(counts_filt, cfg)

  # -------------------------------------------------------------------------
  # Parse IDH covariate
  # -------------------------------------------------------------------------
  idh_status <- parse_idh_status(clinical, cfg$clinical$covariates$idh_status)
  clinical$idh_status <- idh_status

  # Merge QC flags into clinical
  clinical <- merge(clinical, qc[, c("sample_id", "flag_any")],
                    by.x = 1, by.y = "sample_id", all.x = TRUE)

  # -------------------------------------------------------------------------
  # Save outputs
  # -------------------------------------------------------------------------
  save_output(norm_result$expr,
              processed_path(cfg, "03_expr_normalised.rds"),
              "normalised expression")

  save_output(clinical,
              processed_path(cfg, "03_clinical_final.rds"),
              "clinical metadata (final)")

  save_output(norm_result$dge,
              processed_path(cfg, "03_dge_object.rds"),
              "DGEList object")

  write_normalise_summary(norm_result, clinical, cfg)

  log_stage_end("03 — Normalise", start_time)
  invisible(list(expr = norm_result$expr, clinical = clinical))
}

write_normalise_summary <- function(norm_result, clinical, cfg) {
  expr <- norm_result$expr
  out_path <- processed_path(cfg, "03_normalise_summary.txt")

  lines <- c(
    "========================================",
    "03_normalise.R — Summary",
    paste("Date:", Sys.time()),
    "========================================",
    "",
    paste("Normalisation:", cfg$normalisation$method, "+", cfg$normalisation$transform),
    paste("Genes (after filter):", nrow(expr)),
    paste("Samples:", ncol(expr)),
    "",
    "Expression range:",
    paste("  min:", signif(min(expr), 4)),
    paste("  max:", signif(max(expr), 4)),
    paste("  median:", signif(median(expr), 4)),
    "",
    "IDH status distribution:",
    capture.output(print(table(clinical$idh_status, useNA = "ifany")))
  )

  writeLines(lines, out_path)
  log_info("Wrote normalisation summary → ", basename(out_path))
}

# =============================================================================
# Run if called directly
# =============================================================================
if (sys.nframe() == 0) {
  cfg <- load_config()
  run_normalise(cfg)
}
