# =============================================================================
# 07_bn_genie3.R — Feature importance network: GENIE3
# =============================================================================
# GENIE3 (GEne Network Inference with Ensemble of trees) learns a weighted
# directed network using Random Forests / Extra-Trees. Each gene's expression
# is regressed on all other genes; edge weights reflect predictor importance.
#
# NOTE on interpretation:
#   GENIE3 is NOT a causal method — it produces a relevance network.
#   It is included as a standard benchmark comparison (it performs well on
#   DREAM challenges) but edges are not directionally interpretable the same
#   way as BN or PC edges. Include in comparisons with this caveat noted.
#
# Inputs:
#   - data/processed/04_expr_bn_input.rds
#
# Outputs (in results/genie3/):
#   - weight_matrix.rds      : p × p importance weight matrix
#   - edge_list.csv          : ranked edges above threshold
#   - egfr_neighbourhood.rds
#   - summary.txt
# =============================================================================

source(file.path("R", "00_utils.R"))

check_packages(c("GENIE3"))

suppressPackageStartupMessages({
  library(GENIE3)
})

# =============================================================================
# GENIE3 runner
# =============================================================================

#' Run GENIE3 on an expression matrix
#'
#' @param expr       genes × samples matrix (GENIE3 expects this orientation)
#' @param n_trees    Number of trees per gene (default 1000 for final; 100 for fast)
#' @param n_cores    Number of parallel cores (1 = no parallelism)
#' @return weight matrix (p × p), rows = regulators, cols = targets
run_genie3 <- function(expr, n_trees = 500, n_cores = 1) {
  log_info("[GENIE3] ", nrow(expr), " genes × ", ncol(expr), " samples")
  log_info("[GENIE3] n_trees = ", n_trees, ", n_cores = ", n_cores)

  # GENIE3 expects genes × samples (matches our internal format)
  # Input should be a numeric matrix
  expr_mat <- as.matrix(expr)

  weight_mat <- GENIE3(
    exprMatrix  = expr_mat,
    nTrees      = n_trees,
    nCores      = n_cores,
    verbose     = TRUE
  )

  log_info("[GENIE3] Weight matrix: ", nrow(weight_mat), " × ", ncol(weight_mat))
  weight_mat
}

# =============================================================================
# Weight matrix → ranked edge list
# =============================================================================

#' Convert GENIE3 weight matrix to a ranked edge list
#'
#' @param weight_mat   p × p weight matrix from GENIE3()
#' @param name_map     Named vector: safe_name → original_name
#' @param threshold    Minimum weight to include an edge (default: top 5% of weights)
#' @return data.frame sorted by weight descending
weight_matrix_to_edges <- function(weight_mat, name_map = NULL,
                                    threshold = NULL) {
  p <- nrow(weight_mat)

  # Flatten the matrix (exclude diagonal)
  edges <- data.frame(
    from   = rep(rownames(weight_mat), times = p),
    to     = rep(colnames(weight_mat), each  = p),
    weight = as.vector(weight_mat),
    stringsAsFactors = FALSE
  )

  # Remove self-edges (GENIE3 sets diagonal to 0, but be explicit)
  edges <- edges[edges$from != edges$to, ]

  # Apply threshold — default: keep edges above 5th percentile
  if (is.null(threshold)) {
    threshold <- quantile(edges$weight, 0.95)
    log_info("[GENIE3] Auto-threshold (top 5%): ", signif(threshold, 3))
  }
  edges <- edges[edges$weight >= threshold, ]

  # Sort descending
  edges <- edges[order(-edges$weight), ]

  # Map to original names if provided
  if (!is.null(name_map)) {
    edges$from_original <- name_map[edges$from] %||% edges$from
    edges$to_original   <- name_map[edges$to]   %||% edges$to
  }

  log_info("[GENIE3] Edges above threshold: ", nrow(edges))
  edges
}

# =============================================================================
# EGFR neighbourhood
# =============================================================================

#' Extract EGFR neighbourhood from GENIE3 weight matrix
extract_egfr_genie3 <- function(weight_mat, name_map = NULL,
                                  hops = 1, top_n_per_gene = 10) {
  gene_names <- rownames(weight_mat)

  # Find EGFR
  target <- grep("^EGFR$", gene_names, value = TRUE, ignore.case = TRUE)
  if (length(target) == 0 && !is.null(name_map)) {
    target <- names(name_map)[name_map == "EGFR"]
    target <- intersect(target, gene_names)
  }
  if (length(target) == 0) {
    log_warn("  EGFR not found in GENIE3 weight matrix")
    return(NULL)
  }
  target <- target[1]
  original_target <- if (!is.null(name_map)) name_map[target] %||% target else target

  # Top regulators of EGFR (high weight in weight_mat[regulator, EGFR])
  egfr_col   <- weight_mat[, target]
  regulators <- names(sort(egfr_col, decreasing = TRUE))[1:min(top_n_per_gene,
                                                                 length(egfr_col))]
  regulators <- setdiff(regulators, target)

  # Top targets of EGFR (high weight in weight_mat[EGFR, target])
  egfr_row <- weight_mat[target, ]
  targets  <- names(sort(egfr_row, decreasing = TRUE))[1:min(top_n_per_gene,
                                                               length(egfr_row))]
  targets <- setdiff(targets, target)

  remap <- function(x) {
    if (!is.null(name_map)) name_map[x] %||% x else x
  }

  log_info("  EGFR top regulators (in): ", paste(head(remap(regulators), 5),
                                                   collapse = ", "))
  log_info("  EGFR top targets (out):   ", paste(head(remap(targets), 5),
                                                   collapse = ", "))

  list(
    target       = original_target,
    regulators   = remap(regulators),   # genes that predict EGFR
    targets      = remap(targets),      # genes that EGFR predicts
    neighbourhood = remap(unique(c(target, regulators, targets)))
  )
}

# =============================================================================
# Main
# =============================================================================

run_bn_genie3 <- function(cfg) {
  log_stage_start("07 — GENIE3 (importance-based network)")
  start_time <- Sys.time()

  expr <- readRDS(processed_path(cfg, "04_expr_bn_input.rds"))
  log_info("Input: ", nrow(expr), " genes × ", ncol(expr), " samples")

  out_dir <- results_path(cfg, "genie3")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Name mapping
  original_names <- rownames(expr)
  safe_names     <- make.names(rownames(expr), unique = TRUE)
  name_map       <- setNames(original_names, safe_names)
  rownames(expr) <- safe_names

  # Determine n_trees from config (default 500; reduce for faster testing)
  n_trees <- cfg$bn_learn$genie3_n_trees %||% 500
  n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
  log_info("Using ", n_cores, " core(s) (SLURM_CPUS_PER_TASK)")

  # Run GENIE3
  weight_mat <- run_genie3(expr, n_trees = n_trees, n_cores = n_cores)
  saveRDS(weight_mat, file.path(out_dir, "weight_matrix.rds"))

  # Edge list
  edge_df <- weight_matrix_to_edges(weight_mat, name_map)
  write.csv(edge_df, file.path(out_dir, "edge_list.csv"), row.names = FALSE)

  # EGFR neighbourhood
  egfr_nhd <- extract_egfr_genie3(weight_mat, name_map,
                                    hops = cfg$validation$egfr_neighbourhood$hops)
  if (!is.null(egfr_nhd)) {
    saveRDS(egfr_nhd, file.path(out_dir, "egfr_neighbourhood.rds"))
  }

  # Summary
  lines <- c(
    "========================================",
    "07_bn_genie3.R — Summary",
    paste("Date:", Sys.time()),
    "========================================",
    "",
    "NOTE: GENIE3 is NOT a causal method.",
    "Edges reflect feature importance, not causal direction.",
    "",
    paste("Genes:", nrow(expr)),
    paste("Samples:", ncol(expr)),
    paste("n_trees:", n_trees),
    paste("Edges above threshold:", nrow(edge_df)),
    "",
    "Top 10 edges by weight:",
    capture.output(print(head(edge_df[,
      intersect(c("from_original", "to_original", "weight"), colnames(edge_df))
    ], 10), row.names = FALSE))
  )
  if (!is.null(egfr_nhd)) {
    lines <- c(lines,
      "",
      "EGFR top regulators:",
      paste(" ", head(egfr_nhd$regulators, 10), collapse = "\n"),
      "EGFR top targets:",
      paste(" ", head(egfr_nhd$targets, 10), collapse = "\n")
    )
  }
  writeLines(lines, file.path(out_dir, "summary.txt"))
  log_info("Wrote summary → ", file.path(out_dir, "summary.txt"))

  log_stage_end("07 — GENIE3", start_time)
  invisible(list(weight_mat = weight_mat, edge_df = edge_df, egfr = egfr_nhd))
}

# =============================================================================
# Run if called directly
# =============================================================================
if (sys.nframe() == 0) {
  cfg <- load_config()
  run_bn_genie3(cfg)
}
