# =============================================================================
# 05_bn_bnlearn.R — Bayesian Network learning with bnlearn (Hill Climbing)
# =============================================================================
# First BN method in the comparison pipeline. Learns a DAG from the
# normalised expression matrix using bnlearn's Hill Climbing algorithm,
# then assesses edge stability via bootstrapping.
#
# Inputs:
#   - data/processed/04_expr_bn_input.rds
#
# Outputs:
#   - results/05_bn_hc_dag.rds              : learned bn object
#   - results/05_bn_hc_bootstrap.rds        : bootstrap strength object
#   - results/05_bn_hc_averaged.rds         : averaged network (stable edges)
#   - results/05_bn_hc_edge_list.csv        : edge list with bootstrap freq
#   - results/05_egfr_neighbourhood.rds     : EGFR local graph extraction
#   - results/05_bn_hc_summary.txt
# =============================================================================

source(file.path("R", "00_utils.R"))

check_packages(c("bnlearn"))

suppressPackageStartupMessages({
  library(bnlearn)
})

# =============================================================================
# BN structure learning
# =============================================================================

#' Learn a DAG using Hill Climbing
#'
#' @param expr Expression matrix (genes × samples) — will be transposed
#'   to samples × genes for bnlearn
#' @param cfg  Pipeline config
#' @return bn object (learned DAG)
learn_bn_hc <- function(expr, cfg) {
  bn_cfg <- cfg$bn_learn
  log_info("Algorithm: ", bn_cfg$algorithm)
  log_info("Score: ", bn_cfg$score)
  log_info("Restarts: ", bn_cfg$restart, ", Perturb: ", bn_cfg$perturb)

  # bnlearn expects samples × variables
  data <- as.data.frame(t(expr))

  # Sanitise column names — bnlearn doesn't like special characters
  original_names <- colnames(data)
  safe_names <- make.names(colnames(data), unique = TRUE)
  colnames(data) <- safe_names
  name_map <- setNames(original_names, safe_names)

  log_info("Learning DAG from ", nrow(data), " samples × ",
           ncol(data), " genes")

  dag <- hc(
    data,
    score    = bn_cfg$score,
    restart  = bn_cfg$restart,
    perturb  = bn_cfg$perturb,
    max.iter = bn_cfg$max_iter
  )

  n_edges <- nrow(arcs(dag))
  log_info("Learned DAG: ", length(nodes(dag)), " nodes, ",
           n_edges, " edges")

  # Store the name mapping as an attribute for later use
  attr(dag, "name_map") <- name_map

  dag
}

# =============================================================================
# Bootstrap stability assessment
# =============================================================================

#' Assess edge stability via nonparametric bootstrap
#'
#' @param expr Expression matrix (genes × samples)
#' @param cfg  Pipeline config
#' @return list with:
#'   - strength: data.frame of edge bootstrap frequencies
#'   - averaged: bn object with only stable edges
bootstrap_stability <- function(expr, cfg) {
  bn_cfg <- cfg$bn_learn
  stab_cfg <- bn_cfg$stability

  data <- as.data.frame(t(expr))
  colnames(data) <- make.names(colnames(data), unique = TRUE)

  log_info("Bootstrap stability: ", stab_cfg$n_bootstrap, " resamples")
  log_info("Edge threshold: ", stab_cfg$edge_threshold)

  boot_strength <- boot.strength(
    data,
    R         = stab_cfg$n_bootstrap,
    algorithm = bn_cfg$algorithm,
    algorithm.args = list(
      score   = bn_cfg$score,
      restart = bn_cfg$restart,
      perturb = bn_cfg$perturb
    )
  )

  # Averaged network — keep edges above threshold
  avg_net <- averaged.network(boot_strength,
                               threshold = stab_cfg$edge_threshold)

  n_stable <- nrow(arcs(avg_net))
  log_info("Stable edges (freq >= ", stab_cfg$edge_threshold, "): ", n_stable)

  # Summary of bootstrap strengths
  log_info("Bootstrap strength distribution:")
  log_info("  min: ", signif(min(boot_strength$strength), 3))
  log_info("  median: ", signif(median(boot_strength$strength), 3))
  log_info("  mean: ", signif(mean(boot_strength$strength), 3))
  log_info("  max: ", signif(max(boot_strength$strength), 3))

  list(strength = boot_strength, averaged = avg_net)
}

# =============================================================================
# EGFR neighbourhood extraction
# =============================================================================

#' Extract the local neighbourhood around EGFR
#'
#' @param dag    bn object (or averaged network)
#' @param hops   Number of hops (1 = parents + children)
#' @param target Node name for EGFR (auto-detected)
#' @return list with:
#'   - target: the EGFR node name
#'   - parents: parent nodes
#'   - children: child nodes
#'   - neighbourhood: all nodes within `hops` distance
#'   - subgraph_edges: edge list for the subgraph
extract_egfr_neighbourhood <- function(dag, hops = 1, target = NULL) {
  all_nodes <- nodes(dag)

  # Auto-detect EGFR node (might be sanitised to EGFR or egfr etc.)
  if (is.null(target)) {
    target <- grep("^EGFR$", all_nodes, value = TRUE, ignore.case = TRUE)
    if (length(target) == 0) {
      target <- grep("EGFR", all_nodes, value = TRUE, ignore.case = TRUE)
    }
    if (length(target) == 0) {
      log_warn("EGFR node not found in DAG. Available nodes: ",
               paste(head(all_nodes, 20), collapse = ", "))
      return(NULL)
    }
    target <- target[1]
  }

  log_info("Extracting neighbourhood for: ", target, " (hops = ", hops, ")")

  # Get parents and children
  parents  <- bnlearn::parents(dag, target)
  children <- bnlearn::children(dag, target)

  log_info("  Parents:  ", if (length(parents) > 0) paste(parents, collapse = ", ") else "(none)")
  log_info("  Children: ", if (length(children) > 0) paste(children, collapse = ", ") else "(none)")

  # Expand to k-hop neighbourhood if requested
  neighbourhood <- unique(c(target, parents, children))
  if (hops >= 2) {
    for (h in 2:hops) {
      new_nodes <- character()
      for (node in neighbourhood) {
        new_nodes <- c(new_nodes,
                       bnlearn::parents(dag, node),
                       bnlearn::children(dag, node))
      }
      neighbourhood <- unique(c(neighbourhood, new_nodes))
    }
    log_info("  ", hops, "-hop neighbourhood: ", length(neighbourhood), " nodes")
  }

  # Extract edges within the neighbourhood
  all_arcs <- arcs(dag)
  sub_arcs <- all_arcs[all_arcs[, 1] %in% neighbourhood &
                        all_arcs[, 2] %in% neighbourhood, , drop = FALSE]

  list(
    target       = target,
    parents      = parents,
    children     = children,
    neighbourhood = neighbourhood,
    subgraph_edges = sub_arcs,
    n_edges      = nrow(sub_arcs)
  )
}

# =============================================================================
# Edge list export
# =============================================================================

#' Export edge list with bootstrap frequencies
export_edge_list <- function(dag, boot_result, name_map, out_path) {
  edges <- as.data.frame(arcs(dag), stringsAsFactors = FALSE)
  colnames(edges) <- c("from", "to")

  # Add bootstrap strength if available
  if (!is.null(boot_result)) {
    strength_df <- boot_result$strength
    edges <- merge(edges, strength_df,
                   by.x = c("from", "to"),
                   by.y = c("from", "to"),
                   all.x = TRUE)
  }

  # Map back to original gene names
  if (!is.null(name_map)) {
    edges$from_original <- name_map[edges$from]
    edges$to_original   <- name_map[edges$to]
  }

  # Sort by strength (descending)
  if ("strength" %in% colnames(edges)) {
    edges <- edges[order(-edges$strength), ]
  }

  write.csv(edges, out_path, row.names = FALSE)
  log_info("Exported ", nrow(edges), " edges → ", basename(out_path))

  edges
}

# =============================================================================
# Main
# =============================================================================

run_bn_bnlearn <- function(cfg) {
  log_stage_start("05 — BN Learning (bnlearn HC)")
  start_time <- Sys.time()

  expr <- readRDS(processed_path(cfg, "04_expr_bn_input.rds"))

  # -------------------------------------------------------------------------
  # Learn structure
  # -------------------------------------------------------------------------
  log_info("--- Single-run structure learning ---")
  dag <- learn_bn_hc(expr, cfg)
  save_output(dag, results_path(cfg, "05_bn_hc_dag.rds"), "HC DAG")

  # -------------------------------------------------------------------------
  # Bootstrap stability
  # -------------------------------------------------------------------------
  log_info("--- Bootstrap stability assessment ---")
  boot_result <- bootstrap_stability(expr, cfg)
  save_output(boot_result$strength,
              results_path(cfg, "05_bn_hc_bootstrap.rds"),
              "bootstrap strengths")
  save_output(boot_result$averaged,
              results_path(cfg, "05_bn_hc_averaged.rds"),
              "averaged network")

  # -------------------------------------------------------------------------
  # Export edge list
  # -------------------------------------------------------------------------
  name_map <- attr(dag, "name_map")
  edge_df <- export_edge_list(
    boot_result$averaged, boot_result,
    name_map,
    results_path(cfg, "05_bn_hc_edge_list.csv")
  )

  # -------------------------------------------------------------------------
  # EGFR neighbourhood extraction
  # -------------------------------------------------------------------------
  egfr_nhd <- extract_egfr_neighbourhood(
    boot_result$averaged,
    hops = cfg$validation$egfr_neighbourhood$hops
  )
  if (!is.null(egfr_nhd)) {
    save_output(egfr_nhd,
                results_path(cfg, "05_egfr_neighbourhood.rds"),
                "EGFR neighbourhood")
  }

  # -------------------------------------------------------------------------
  # Summary
  # -------------------------------------------------------------------------
  write_bn_summary(dag, boot_result, egfr_nhd, cfg)

  log_stage_end("05 — BN Learning (bnlearn HC)", start_time)
  invisible(list(dag = dag, boot = boot_result, egfr = egfr_nhd))
}

write_bn_summary <- function(dag, boot_result, egfr_nhd, cfg) {
  out_path <- results_path(cfg, "05_bn_hc_summary.txt")
  avg <- boot_result$averaged

  lines <- c(
    "========================================",
    "05_bn_bnlearn.R — Summary",
    paste("Date:", Sys.time()),
    "========================================",
    "",
    "Algorithm: Hill Climbing (bnlearn::hc)",
    paste("Score:", cfg$bn_learn$score),
    paste("Restarts:", cfg$bn_learn$restart),
    "",
    "--- Single-run DAG ---",
    paste("Nodes:", length(nodes(dag))),
    paste("Edges:", nrow(arcs(dag))),
    "",
    "--- Averaged Network (bootstrap) ---",
    paste("Bootstrap resamples:", cfg$bn_learn$stability$n_bootstrap),
    paste("Edge threshold:", cfg$bn_learn$stability$edge_threshold),
    paste("Stable edges:", nrow(arcs(avg))),
    ""
  )

  if (!is.null(egfr_nhd)) {
    lines <- c(lines,
      "--- EGFR Neighbourhood ---",
      paste("Target node:", egfr_nhd$target),
      paste("Parents:", paste(egfr_nhd$parents, collapse = ", ")),
      paste("Children:", paste(egfr_nhd$children, collapse = ", ")),
      paste("1-hop neighbourhood size:", length(egfr_nhd$neighbourhood)),
      paste("Subgraph edges:", egfr_nhd$n_edges)
    )
  } else {
    lines <- c(lines, "EGFR neighbourhood: NOT FOUND IN DAG")
  }

  writeLines(lines, out_path)
  log_info("Wrote BN summary → ", basename(out_path))
}

# =============================================================================
# Run if called directly
# =============================================================================
if (sys.nframe() == 0) {
  cfg <- load_config()
  run_bn_bnlearn(cfg)
}
