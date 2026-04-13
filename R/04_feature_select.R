# =============================================================================
# 04_feature_select.R — Variance filtering + EGFR pathway gene set
# =============================================================================
# Selects genes for BN learning using a two-pronged strategy:
#   1. Top-N genes by variance (data-driven)
#   2. Union with curated EGFR reference gene set (biology-driven)
#
# This ensures the BN input is small enough to learn from (~150–200 genes)
# while guaranteeing the EGFR pathway members are present for validation.
#
# Inputs:
#   - data/processed/03_expr_normalised.rds
#
# Outputs:
#   - data/processed/04_expr_bn_input.rds       : expression matrix for BN
#   - data/processed/04_selected_genes.rds       : gene selection metadata
#   - data/processed/04_feature_select_summary.txt
# =============================================================================

source(file.path("R", "00_utils.R"))

# =============================================================================
# Gene selection
# =============================================================================

#' Select genes for BN learning
#'
#' @param expr Normalised expression matrix (genes × samples)
#' @param cfg  Pipeline config
#' @return list with:
#'   - expr_selected: subsetted expression matrix
#'   - gene_meta: data.frame tracking why each gene was included
select_features <- function(expr, cfg) {
  log_info("Feature selection strategy: ", cfg$feature_selection$strategy)

  fs_cfg <- cfg$feature_selection

  # -------------------------------------------------------------------------
  # Step 1: Compute per-gene variance
  # -------------------------------------------------------------------------
  gene_var <- apply(expr, 1, var)
  gene_var_sorted <- sort(gene_var, decreasing = TRUE)

  log_info("Gene variance — range: [",
           signif(min(gene_var), 3), ", ",
           signif(max(gene_var), 3), "]",
           ", median: ", signif(median(gene_var), 3))

  # -------------------------------------------------------------------------
  # Step 2: Top-N by variance
  # -------------------------------------------------------------------------
  top_n <- fs_cfg$top_variance_n
  top_var_genes <- names(gene_var_sorted)[1:min(top_n, length(gene_var_sorted))]
  log_info("Top ", top_n, " variance genes selected")

  # -------------------------------------------------------------------------
  # Step 3: EGFR reference gene set — find matches in the expression data
  # -------------------------------------------------------------------------
  egfr_ref <- fs_cfg$egfr_reference_genes

  # Match by gene name — expression matrix rownames could be Ensembl IDs
  # or gene symbols. Try both.
  if (any(grepl("^ENSG", rownames(expr)))) {
    log_info("Rownames appear to be Ensembl IDs — ",
             "will need gene_info for symbol mapping")
    # Try to load gene_info for mapping
    gene_info_path <- processed_path(cfg, "01_gene_info.rds")
    if (file.exists(gene_info_path)) {
      gene_info <- readRDS(gene_info_path)
      # Build a symbol → ensembl map
      symbol_map <- setNames(gene_info$gene_id, gene_info$gene_name)
      egfr_ensembl <- symbol_map[egfr_ref]
      egfr_ensembl <- egfr_ensembl[!is.na(egfr_ensembl)]
      egfr_in_data <- intersect(egfr_ensembl, rownames(expr))
      log_info("EGFR gene set: ", length(egfr_ref), " symbols → ",
               length(egfr_in_data), " mapped to Ensembl IDs in data")
    } else {
      log_warn("No gene_info available for symbol → Ensembl mapping. ",
               "EGFR gene set matching may be incomplete.")
      egfr_in_data <- intersect(egfr_ref, rownames(expr))
    }
  } else {
    # Rownames are gene symbols (or at least not Ensembl)
    egfr_in_data <- intersect(egfr_ref, rownames(expr))
    log_info("EGFR gene set: ", length(egfr_ref), " genes, ",
             length(egfr_in_data), " found in expression data")
  }

  egfr_missing <- setdiff(egfr_ref,
                           if (exists("egfr_ensembl")) names(egfr_ensembl[egfr_ensembl %in% egfr_in_data])
                           else egfr_in_data)
  if (length(egfr_missing) > 0) {
    log_warn("EGFR genes NOT in data: ",
             paste(egfr_missing, collapse = ", "))
  }

  # -------------------------------------------------------------------------
  # Step 4: Union
  # -------------------------------------------------------------------------
  all_selected <- union(top_var_genes, egfr_in_data)

  # Cap at max_genes if needed
  if (length(all_selected) > fs_cfg$max_genes) {
    log_warn("Selected ", length(all_selected), " genes > max (",
             fs_cfg$max_genes, "). Trimming low-variance non-EGFR genes.")
    # Keep all EGFR genes, trim from the variance set
    non_egfr <- setdiff(all_selected, egfr_in_data)
    non_egfr_ranked <- non_egfr[order(gene_var[non_egfr], decreasing = TRUE)]
    n_keep <- fs_cfg$max_genes - length(egfr_in_data)
    all_selected <- union(egfr_in_data, non_egfr_ranked[1:n_keep])
  }

  log_info("Final gene set: ", length(all_selected), " genes",
           " (", length(intersect(all_selected, egfr_in_data)),
           " from EGFR reference)")

  # -------------------------------------------------------------------------
  # Build gene metadata
  # -------------------------------------------------------------------------
  gene_meta <- data.frame(
    gene = all_selected,
    variance = gene_var[all_selected],
    variance_rank = match(all_selected, names(gene_var_sorted)),
    in_egfr_ref = all_selected %in% egfr_in_data,
    selection_reason = ifelse(
      all_selected %in% egfr_in_data & all_selected %in% top_var_genes,
      "both",
      ifelse(all_selected %in% egfr_in_data, "egfr_reference", "high_variance")
    ),
    stringsAsFactors = FALSE
  )
  gene_meta <- gene_meta[order(gene_meta$variance_rank), ]

  # Subset expression matrix
  expr_selected <- expr[all_selected, ]

  list(expr_selected = expr_selected, gene_meta = gene_meta)
}

# =============================================================================
# Main
# =============================================================================

run_feature_select <- function(cfg) {
  log_stage_start("04 — Feature Selection")
  start_time <- Sys.time()

  expr <- readRDS(processed_path(cfg, "03_expr_normalised.rds"))

  result <- select_features(expr, cfg)

  # Save
  save_output(result$expr_selected,
              processed_path(cfg, "04_expr_bn_input.rds"),
              "BN input expression matrix")

  save_output(result$gene_meta,
              processed_path(cfg, "04_selected_genes.rds"),
              "gene selection metadata")

  # Summary
  out_path <- processed_path(cfg, "04_feature_select_summary.txt")
  meta <- result$gene_meta
  lines <- c(
    "========================================",
    "04_feature_select.R — Summary",
    paste("Date:", Sys.time()),
    "========================================",
    "",
    paste("Strategy:", cfg$feature_selection$strategy),
    paste("Total genes selected:", nrow(meta)),
    paste("  - High variance only:", sum(meta$selection_reason == "high_variance")),
    paste("  - EGFR reference only:", sum(meta$selection_reason == "egfr_reference")),
    paste("  - Both:", sum(meta$selection_reason == "both")),
    "",
    "Top 20 genes by variance:",
    capture.output(print(head(meta, 20), row.names = FALSE)),
    "",
    "EGFR reference genes in final set:",
    paste(" ", meta$gene[meta$in_egfr_ref], collapse = "\n")
  )
  writeLines(lines, out_path)
  log_info("Wrote feature selection summary → ", basename(out_path))

  # Export CSV for SMC and Python-based BN methods
  export_bn_csv(result$expr_selected, result$gene_meta, cfg)

  log_stage_end("04 — Feature Selection", start_time)
  invisible(result)
}

# =============================================================================
# CSV export for SMC + external methods
# =============================================================================

#' Export BN-input expression matrix as CSV (samples × genes, gene symbols)
#'
#' This is the primary handoff to:
#'   - Your SMC method (Python)
#'   - NOTEARS / DAGMA (Python)
#'   - Any tool that reads tabular input
#'
#' Output: data/processed/04_bn_input.csv
#'   Rows    = samples (TCGA barcodes or UUIDs)
#'   Columns = gene symbols (or Ensembl IDs where symbol unavailable)
#'   Values  = logCPM normalised expression (continuous, not counts)
#'
#' @param expr_selected  genes × samples expression matrix (Ensembl ID rownames)
#' @param gene_meta      data.frame from select_features() with gene column
#' @param cfg            pipeline config
export_bn_csv <- function(expr_selected, gene_meta, cfg) {
  log_info("Exporting BN input CSV (samples × genes, gene symbols)")

  # -------------------------------------------------------------------------
  # Step 1: Map Ensembl IDs → gene symbols using gene_info from stage 01
  # -------------------------------------------------------------------------
  gene_info_path <- processed_path(cfg, "01_gene_info.rds")

  if (file.exists(gene_info_path)) {
    gene_info <- readRDS(gene_info_path)

    # Build Ensembl → symbol lookup
    sym_lookup <- setNames(gene_info$gene_name, gene_info$gene_id)

    # Map rownames; fall back to Ensembl ID if no symbol
    mapped_symbols <- sym_lookup[rownames(expr_selected)]
    fallback_mask  <- is.na(mapped_symbols) | nchar(mapped_symbols) == 0
    mapped_symbols[fallback_mask] <- rownames(expr_selected)[fallback_mask]

    n_mapped   <- sum(!fallback_mask)
    n_fallback <- sum(fallback_mask)
    log_info("  Symbol mapping: ", n_mapped, " mapped, ",
             n_fallback, " using Ensembl ID fallback")
  } else {
    log_warn("  gene_info not found — using Ensembl IDs as column names")
    mapped_symbols <- rownames(expr_selected)
  }

  # Make column names safe (no duplicates, no spaces)
  col_names <- make.unique(mapped_symbols, sep = "_")

  # -------------------------------------------------------------------------
  # Step 2: Transpose to samples × genes and write CSV
  # -------------------------------------------------------------------------
  df_out <- as.data.frame(t(expr_selected))
  colnames(df_out) <- col_names

  # Build patient_id column from sample names.
  # After GDC API resolution in stage 01, colnames are full TCGA aliquot
  # barcodes (e.g. TCGA-02-0003-01A-01R-0177-01).  Truncate to 12 chars
  # to get the patient-level ID (TCGA-02-0003) matching bn_ready_matrix.csv.
  # If names are still UUIDs (API unavailable), keep them unchanged.
  sample_names <- rownames(df_out)
  patient_ids <- ifelse(
    grepl("^TCGA-", sample_names),
    substr(sample_names, 1, 12),
    sample_names   # UUID fallback — stage 01 API resolution did not run
  )
  n_truncated <- sum(grepl("^TCGA-", sample_names))
  if (n_truncated > 0) {
    log_info("  Converted ", n_truncated, " aliquot barcodes → 12-char patient IDs")
  } else {
    log_warn("  Sample names are UUIDs — re-run stage 01 to resolve to TCGA barcodes")
  }
  df_out <- cbind(patient_id = patient_ids, df_out)

  csv_path <- processed_path(cfg, "04_bn_input.csv")
  write.csv(df_out, csv_path, row.names = FALSE)

  file_mb <- round(file.size(csv_path) / 1e6, 2)
  log_info("  Wrote ", nrow(df_out), " samples × ", ncol(df_out) - 1,
           " genes → ", basename(csv_path), " (", file_mb, " MB)")

  # Also write the gene metadata with symbols as a reference
  gene_meta$symbol <- mapped_symbols[match(gene_meta$gene, rownames(expr_selected))]
  meta_path <- processed_path(cfg, "04_selected_genes_symbols.csv")
  write.csv(gene_meta, meta_path, row.names = FALSE)
  log_info("  Gene reference table → ", basename(meta_path))

  invisible(csv_path)
}

# =============================================================================
# Run if called directly
# =============================================================================
if (sys.nframe() == 0) {
  cfg <- load_config()
  run_feature_select(cfg)
}
