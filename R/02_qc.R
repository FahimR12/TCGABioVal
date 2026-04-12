# =============================================================================
# 02_qc.R — MAD-based sample QC with flagging
# =============================================================================
# Computes per-sample QC metrics and flags outliers using median ± k × MAD.
# Samples are FLAGGED, not dropped. Dropping decisions are deferred to manual
# review of PCA/MDS/correlation alongside these flags.
#
# Inputs:
#   - data/processed/01_raw_counts.rds
#   - data/processed/01_clinical.rds
#
# Outputs:
#   - data/processed/02_qc_metrics.rds     : per-sample QC metrics + flags
#   - data/processed/02_qc_thresholds.rds  : computed thresholds (for record)
#   - results/02_qc_diagnostic_plots.pdf   : visual diagnostics
#   - data/processed/02_qc_summary.txt     : human-readable summary
# =============================================================================

source(file.path("R", "00_utils.R"))

check_packages(c("matrixStats"))

# =============================================================================
# QC metric computation
# =============================================================================

#' Compute per-sample QC metrics from a raw count matrix
#'
#' @param counts Integer matrix (genes × samples)
#' @param gene_info Optional data.frame with gene_id + gene_name columns.
#'   If provided and contains gene_type or gene_name, mitochondrial genes
#'   are identified automatically.
#' @return data.frame with one row per sample
compute_qc_metrics <- function(counts, gene_info = NULL) {
  log_info("Computing per-sample QC metrics")

  n_samples <- ncol(counts)
  n_genes   <- nrow(counts)
  log_info("Input: ", n_genes, " genes × ", n_samples, " samples")

  # Library size (total counts per sample)
  lib_size <- colSums(counts)

  # Number of detected genes (count > 0)
  n_detected <- colSums(counts > 0)

  # Mitochondrial percentage
  mito_genes <- identify_mito_genes(rownames(counts), gene_info)
  if (length(mito_genes) > 0) {
    mito_counts <- colSums(counts[mito_genes, , drop = FALSE])
    pct_mito <- 100 * mito_counts / lib_size
    log_info("Identified ", length(mito_genes), " mitochondrial genes")
  } else {
    pct_mito <- rep(NA_real_, n_samples)
    log_warn("No mitochondrial genes identified — pct_mito will be NA. ",
             "Check that gene names use 'MT-' prefix or gene_info is provided.")
  }

  data.frame(
    sample_id  = colnames(counts),
    lib_size   = lib_size,
    n_detected = n_detected,
    pct_mito   = pct_mito,
    stringsAsFactors = FALSE
  )
}

#' Identify mitochondrial genes by name pattern
identify_mito_genes <- function(gene_ids, gene_info = NULL) {
  # Try matching on gene names if gene_info is available
  if (!is.null(gene_info) && "gene_name" %in% colnames(gene_info)) {
    mito_idx <- grepl("^MT-", gene_info$gene_name, ignore.case = TRUE)
    if (sum(mito_idx) > 0) {
      # Return the corresponding gene_ids (rownames of the count matrix)
      return(gene_info$gene_id[mito_idx])
    }
  }
  # Fall back to matching on gene_ids directly
  mito <- gene_ids[grepl("^MT-", gene_ids, ignore.case = TRUE)]
  mito
}

# =============================================================================
# MAD-based outlier flagging
# =============================================================================

#' Flag outlier samples using MAD-based thresholds
#'
#' @param qc_metrics data.frame from compute_qc_metrics()
#' @param qc_config  The qc section of the pipeline config
#' @return list with:
#'   - metrics: qc_metrics with added flag columns
#'   - thresholds: data.frame of computed thresholds per metric
flag_outliers <- function(qc_metrics, qc_config) {
  log_info("Flagging outliers (MAD threshold = ", qc_config$mad_threshold, ")")

  k <- qc_config$mad_threshold
  thresholds <- list()

  for (metric_name in names(qc_config$metrics)) {
    metric_cfg <- qc_config$metrics[[metric_name]]
    if (!isTRUE(metric_cfg$enabled)) next

    values <- qc_metrics[[metric_name]]
    if (all(is.na(values))) {
      log_warn("Skipping ", metric_name, " — all values are NA")
      qc_metrics[[paste0("flag_", metric_name)]] <- FALSE
      next
    }

    med <- median(values, na.rm = TRUE)
    mad_val <- mad(values, na.rm = TRUE)

    # Compute threshold based on direction
    if (metric_cfg$direction == "lower") {
      threshold <- med - k * mad_val
      flagged <- !is.na(values) & values < threshold
      threshold_label <- paste0("< ", signif(threshold, 4))
    } else if (metric_cfg$direction == "upper") {
      threshold <- med + k * mad_val
      flagged <- !is.na(values) & values > threshold
      threshold_label <- paste0("> ", signif(threshold, 4))
    } else {
      # Both tails
      lower <- med - k * mad_val
      upper <- med + k * mad_val
      flagged <- !is.na(values) & (values < lower | values > upper)
      threshold_label <- paste0("< ", signif(lower, 4),
                                " or > ", signif(upper, 4))
    }

    qc_metrics[[paste0("flag_", metric_name)]] <- flagged

    thresholds[[metric_name]] <- data.frame(
      metric    = metric_name,
      median    = med,
      mad       = mad_val,
      k         = k,
      direction = metric_cfg$direction,
      threshold = threshold_label,
      n_flagged = sum(flagged),
      stringsAsFactors = FALSE
    )

    log_info("  ", metric_name, ": median = ", signif(med, 4),
             ", MAD = ", signif(mad_val, 4),
             ", threshold ", threshold_label,
             ", flagged = ", sum(flagged))
  }

  # Composite flag: flagged on ANY metric
  flag_cols <- grep("^flag_", colnames(qc_metrics), value = TRUE)
  qc_metrics$flag_any <- rowSums(qc_metrics[flag_cols], na.rm = TRUE) > 0

  n_any <- sum(qc_metrics$flag_any)
  log_info("Total samples flagged (any metric): ", n_any, " / ",
           nrow(qc_metrics))

  thresholds_df <- do.call(rbind, thresholds)

  list(metrics = qc_metrics, thresholds = thresholds_df)
}

# =============================================================================
# Diagnostic plots
# =============================================================================

#' Generate QC diagnostic plots
#'
#' @param qc_result Output of flag_outliers()
#' @param out_path  Path for the output PDF
plot_qc_diagnostics <- function(qc_result, out_path) {
  log_info("Generating QC diagnostic plots")

  metrics <- qc_result$metrics
  thresholds <- qc_result$thresholds

  pdf(out_path, width = 10, height = 8)
  on.exit(dev.off())

  # Set up multi-panel layout
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  # 1. Library size distribution
  plot_metric_hist(metrics, "lib_size", thresholds,
                   main = "Library Size",
                   xlab = "Total counts")

  # 2. Genes detected
  plot_metric_hist(metrics, "n_detected", thresholds,
                   main = "Genes Detected",
                   xlab = "Number of genes (count > 0)")

  # 3. Mitochondrial percentage
  if (!all(is.na(metrics$pct_mito))) {
    plot_metric_hist(metrics, "pct_mito", thresholds,
                     main = "Mitochondrial %",
                     xlab = "% mitochondrial reads")
  }

  # 4. Lib size vs genes detected (scatter, coloured by flag)
  col <- ifelse(metrics$flag_any, "red", "grey40")
  plot(metrics$lib_size, metrics$n_detected,
       col = col, pch = 16, cex = 0.7,
       xlab = "Library size", ylab = "Genes detected",
       main = "Library Size vs Genes Detected")
  legend("bottomright",
         legend = c("Pass", "Flagged"),
         col = c("grey40", "red"), pch = 16, cex = 0.8)

  log_info("Saved QC plots → ", basename(out_path))
}

#' Helper: histogram for a single QC metric with threshold line
plot_metric_hist <- function(metrics, metric_name, thresholds, ...) {
  values <- metrics[[metric_name]]
  flagged <- metrics[[paste0("flag_", metric_name)]]
  if (all(is.na(values))) return(invisible(NULL))

  hist(values, breaks = 30, col = "steelblue", border = "white", ...)

  # Add threshold line
  if (metric_name %in% thresholds$metric) {
    thr_row <- thresholds[thresholds$metric == metric_name, ]
    # Parse threshold value from the label
    thr_str <- thr_row$threshold
    thr_val <- as.numeric(gsub("[<> ]", "", strsplit(thr_str, " or ")[[1]]))
    abline(v = thr_val, col = "red", lwd = 2, lty = 2)
  }

  # Mark flagged count
  if (any(flagged, na.rm = TRUE)) {
    mtext(paste0(sum(flagged, na.rm = TRUE), " flagged"),
          side = 3, adj = 1, col = "red", cex = 0.8)
  }
}

# =============================================================================
# Main QC function
# =============================================================================

run_qc <- function(cfg) {
  log_stage_start("02 — Sample QC")
  start_time <- Sys.time()

  # Load ingested data
  counts   <- readRDS(processed_path(cfg, "01_raw_counts.rds"))
  clinical <- readRDS(processed_path(cfg, "01_clinical.rds"))

  # Try to load gene_info if available
  gene_info_path <- processed_path(cfg, "01_gene_info.rds")
  gene_info <- if (file.exists(gene_info_path)) readRDS(gene_info_path) else NULL

  # Compute metrics
  qc_metrics <- compute_qc_metrics(counts, gene_info)

  # Flag outliers
  qc_result <- flag_outliers(qc_metrics, cfg$qc)

  # Merge flags into clinical metadata
  clinical_with_qc <- merge(clinical, qc_result$metrics,
                            by.x = 1, by.y = "sample_id",
                            all.x = TRUE)

  # Save outputs
  save_output(qc_result$metrics,
              processed_path(cfg, "02_qc_metrics.rds"),
              "QC metrics + flags")

  save_output(qc_result$thresholds,
              processed_path(cfg, "02_qc_thresholds.rds"),
              "QC thresholds")

  # Generate diagnostic plots
  plot_qc_diagnostics(
    qc_result,
    results_path(cfg, "02_qc_diagnostic_plots.pdf")
  )

  # Write human-readable summary
  write_qc_summary(qc_result, cfg)

  log_stage_end("02 — Sample QC", start_time)
  invisible(qc_result)
}

write_qc_summary <- function(qc_result, cfg) {
  metrics <- qc_result$metrics
  thresholds <- qc_result$thresholds
  out_path <- processed_path(cfg, "02_qc_summary.txt")

  lines <- c(
    "========================================",
    "02_qc.R — Sample QC Summary",
    paste("Date:", Sys.time()),
    paste("MAD threshold:", cfg$qc$mad_threshold),
    "========================================",
    "",
    "Thresholds:",
    capture.output(print(thresholds, row.names = FALSE)),
    "",
    paste("Total samples:", nrow(metrics)),
    paste("Flagged (any metric):", sum(metrics$flag_any)),
    paste("Passing:", sum(!metrics$flag_any)),
    "",
    "Flagged samples:",
    if (any(metrics$flag_any)) {
      flagged <- metrics[metrics$flag_any, ]
      capture.output(print(flagged[, c("sample_id", grep("flag_",
                     colnames(flagged), value = TRUE))], row.names = FALSE))
    } else {
      "  (none)"
    }
  )

  writeLines(lines, out_path)
  log_info("Wrote QC summary → ", basename(out_path))
}

# =============================================================================
# Run if called directly
# =============================================================================
if (sys.nframe() == 0) {
  cfg <- load_config()
  run_qc(cfg)
}
