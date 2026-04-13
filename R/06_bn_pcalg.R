# =============================================================================
# 06_bn_pcalg.R — Constraint-based BN learning: PC algorithm + GES
# =============================================================================
# Runs two causal discovery methods from the pcalg package:
#
#   PC  (Peter-Clark) — constraint-based; uses conditional independence tests
#                       to orient edges. Returns a CPDAG (equivalence class).
#   GES (Greedy Equivalence Search) — score-based search over CPDAGs.
#                       Typically faster than PC on large graphs.
#
# Both return a completed partially directed acyclic graph (CPDAG), not a
# fully oriented DAG. For comparison with bnlearn/SMC, we:
#   - Report the CPDAG skeleton (undirected edges)
#   - Orient as many edges as the algorithm can (directed edges)
#   - Extract the EGFR neighbourhood from whatever orientation is available
#
# Inputs:
#   - data/processed/04_expr_bn_input.rds
#
# Outputs (in results/pcalg/{pc,ges}/):
#   - cpdag.rds            : raw pcalg output object
#   - edge_list.csv        : directed + undirected edges
#   - egfr_neighbourhood.rds
#   - summary.txt
#
# Also:
#   - results/pcalg/comparison_edges.csv
# =============================================================================

source(file.path("R", "00_utils.R"))

check_packages(c("pcalg", "graph", "RBGL"))

suppressPackageStartupMessages({
  library(pcalg)
  library(graph)
})

# =============================================================================
# PC algorithm
# =============================================================================

#' Run the PC algorithm
#'
#' @param data  samples × genes data.frame (scaled expression)
#' @param alpha Significance level for conditional independence tests
#' @return pcAlgo object
run_pc <- function(data, alpha = 0.01) {
  log_info("[PC] alpha = ", alpha, ", n = ", nrow(data),
           ", p = ", ncol(data))

  n <- nrow(data)
  p <- ncol(data)

  # Sufficient statistics for Gaussian CI tests (correlation matrix)
  suff_stat <- list(C = cor(data), n = n)

  # PC algorithm with Fisher-Z conditional independence test
  pc_fit <- pc(
    suffStat  = suff_stat,
    indepTest = gaussCItest,
    alpha     = alpha,
    p         = p,
    labels    = colnames(data),
    verbose   = FALSE
  )

  log_info("[PC] Done — edges in CPDAG: ",
           sum(as(pc_fit@graph, "matrix") != 0))
  pc_fit
}

# =============================================================================
# GES
# =============================================================================

#' Run GES (Greedy Equivalence Search)
#'
#' @param data  samples × genes data.frame
#' @return GES result list
run_ges <- function(data) {
  log_info("[GES] n = ", nrow(data), ", p = ", ncol(data))

  # GES score — BIC for Gaussian data
  score <- new("GaussL0penObsScore", data = data)

  ges_fit <- ges(score, verbose = FALSE)

  n_edges <- sum(as(ges_fit$essgraph, "matrix") != 0)
  log_info("[GES] Done — edges in CPDAG: ", n_edges)
  ges_fit
}

# =============================================================================
# Shared helpers for pcalg output
# =============================================================================

#' Extract edge list from a pcalg CPDAG/PDAG as a data.frame
#'
#' pcalg encodes:
#'   amat[i, j] = 1, amat[j, i] = 0  →  i → j  (directed)
#'   amat[i, j] = 1, amat[j, i] = 1  →  i — j  (undirected)
#'
#' @param cpdag_obj  A pcAlgo object (PC) or the $essgraph slot (GES)
#' @param node_names Character vector of variable names
#' @return data.frame with columns: from, to, type ("directed" or "undirected")
cpdag_to_edge_list <- function(cpdag_obj, node_names) {
  # Extract adjacency matrix — handle both pcAlgo and EssGraph
  if (inherits(cpdag_obj, "pcAlgo")) {
    amat <- as(cpdag_obj@graph, "matrix")
  } else if (inherits(cpdag_obj, "EssGraph")) {
    amat <- as(cpdag_obj, "matrix")
  } else {
    amat <- as.matrix(cpdag_obj)
  }

  rownames(amat) <- node_names
  colnames(amat) <- node_names

  edges <- list()
  n <- nrow(amat)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (amat[i, j] == 1) {
        type <- if (amat[j, i] == 1) "undirected" else "directed"
        # Avoid duplicate undirected edges (i-j and j-i)
        if (type == "undirected" && j < i) next
        edges[[length(edges) + 1]] <- data.frame(
          from = node_names[i],
          to   = node_names[j],
          type = type,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(edges) == 0) {
    return(data.frame(from = character(), to = character(),
                       type = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, edges)
}

#' Extract EGFR neighbourhood from a pcalg adjacency matrix
extract_egfr_pcalg <- function(cpdag_obj, node_names, hops = 1) {
  # Find EGFR node
  target_idx <- which(grepl("^EGFR$", node_names, ignore.case = TRUE))
  if (length(target_idx) == 0) {
    log_warn("  EGFR not found in node list")
    return(NULL)
  }
  target_idx  <- target_idx[1]
  target_name <- node_names[target_idx]

  # Adjacency matrix
  if (inherits(cpdag_obj, "pcAlgo")) {
    amat <- as(cpdag_obj@graph, "matrix")
  } else {
    amat <- as(cpdag_obj, "matrix")
  }

  # Neighbours (any edge to/from EGFR, directed or undirected)
  neighbourhood_idx <- which(amat[target_idx, ] != 0 | amat[, target_idx] != 0)
  parents_idx   <- which(amat[, target_idx] == 1 & amat[target_idx, ] == 0)
  children_idx  <- which(amat[target_idx, ] == 1 & amat[, target_idx] == 0)
  adjacent_idx  <- which(amat[target_idx, ] == 1 & amat[, target_idx] == 1)

  log_info("  EGFR: ", length(parents_idx), " parents, ",
           length(children_idx), " children, ",
           length(adjacent_idx), " undirected neighbours")

  list(
    target    = target_name,
    parents   = node_names[parents_idx],
    children  = node_names[children_idx],
    adjacent  = node_names[adjacent_idx],   # undirected connections
    neighbourhood = node_names[c(target_idx, neighbourhood_idx)]
  )
}

# =============================================================================
# Main runner per method
# =============================================================================

run_single_pcalg <- function(expr, method, cfg) {
  out_dir <- results_path(cfg, "pcalg", method)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Transpose + scale: pcalg expects samples × genes, standardised
  data <- as.data.frame(t(expr))
  node_names_orig <- colnames(data)
  # Safe names for pcalg
  safe_names <- make.names(colnames(data), unique = TRUE)
  name_map   <- setNames(node_names_orig, safe_names)
  colnames(data) <- safe_names

  data_scaled <- as.data.frame(scale(data))

  # Run the algorithm
  fit <- switch(method,
    pc  = run_pc(data_scaled, alpha = cfg$bn_learn$pc_alpha %||% 0.01),
    ges = run_ges(data_scaled),
    stop("Unknown pcalg method: ", method)
  )

  saveRDS(fit, file.path(out_dir, "cpdag.rds"))

  # Extract edge list (map back to original names)
  cpdag_obj <- if (method == "pc") fit else fit$essgraph
  edge_df <- cpdag_to_edge_list(cpdag_obj, safe_names)

  # Map to original names
  if (nrow(edge_df) > 0) {
    edge_df$from_original <- name_map[edge_df$from]
    edge_df$to_original   <- name_map[edge_df$to]
  }
  write.csv(edge_df, file.path(out_dir, "edge_list.csv"), row.names = FALSE)
  log_info("[", toupper(method), "] ",
           sum(edge_df$type == "directed"), " directed + ",
           sum(edge_df$type == "undirected"), " undirected edges")

  # EGFR neighbourhood
  egfr_nhd <- extract_egfr_pcalg(cpdag_obj, safe_names,
                                   hops = cfg$validation$egfr_neighbourhood$hops)
  if (!is.null(egfr_nhd)) {
    # Remap to original names
    remap <- function(x) name_map[x] %||% x
    egfr_nhd$parents      <- remap(egfr_nhd$parents)
    egfr_nhd$children     <- remap(egfr_nhd$children)
    egfr_nhd$adjacent     <- remap(egfr_nhd$adjacent)
    egfr_nhd$neighbourhood <- remap(egfr_nhd$neighbourhood)
    saveRDS(egfr_nhd, file.path(out_dir, "egfr_neighbourhood.rds"))
  }

  # Summary
  lines <- c(
    "========================================",
    paste0("pcalg — ", toupper(method)),
    paste("Date:", Sys.time()),
    "========================================",
    "",
    paste("Nodes:", ncol(data)),
    paste("Samples:", nrow(data)),
    paste("Directed edges:", sum(edge_df$type == "directed")),
    paste("Undirected edges:", sum(edge_df$type == "undirected")),
    ""
  )
  if (!is.null(egfr_nhd)) {
    lines <- c(lines,
      "--- EGFR neighbourhood ---",
      paste("Parents:   ", paste(egfr_nhd$parents,   collapse = ", ")),
      paste("Children:  ", paste(egfr_nhd$children,  collapse = ", ")),
      paste("Undirected:", paste(egfr_nhd$adjacent,  collapse = ", "))
    )
  }
  writeLines(lines, file.path(out_dir, "summary.txt"))

  list(fit = fit, edge_df = edge_df, egfr = egfr_nhd)
}

# =============================================================================
# Main
# =============================================================================

run_bn_pcalg <- function(cfg) {
  log_stage_start("06 — BN Learning (pcalg: PC + GES)")
  start_time <- Sys.time()

  expr <- readRDS(processed_path(cfg, "04_expr_bn_input.rds"))
  log_info("Input: ", nrow(expr), " genes × ", ncol(expr), " samples")

  methods <- c("pc", "ges")
  results <- list()

  for (m in methods) {
    log_info("")
    log_info(strrep("-", 50))
    log_info("Method: ", toupper(m))
    log_info(strrep("-", 50))

    results[[m]] <- tryCatch(
      run_single_pcalg(expr, m, cfg),
      error = function(e) {
        log_error("pcalg method ", toupper(m), " failed: ", conditionMessage(e))
        NULL
      }
    )
  }

  log_stage_end("06 — BN Learning (pcalg)", start_time)
  invisible(results)
}

# =============================================================================
# Run if called directly
# =============================================================================
if (sys.nframe() == 0) {
  cfg <- load_config()
  run_bn_pcalg(cfg)
}
