# =============================================================================
# 05_bn_bnlearn.R — BN structure learning: HC, Tabu, MMHC (bnlearn)
# =============================================================================
# Runs all three bnlearn score-based algorithms on the same expression matrix
# so that results are directly comparable. Each algorithm:
#   - Learns a single-run DAG
#   - Assesses edge stability via bootstrap resampling
#   - Extracts the EGFR 1-hop neighbourhood
#   - Writes outputs to results/bnlearn/<algorithm>/
#
# Inputs:
#   - data/processed/04_expr_bn_input.rds
#
# Outputs per algorithm (in results/bnlearn/{hc,tabu,mmhc}/):
#   - dag.rds              : single-run learned bn object
#   - bootstrap.rds        : boot.strength data.frame
#   - averaged.rds         : averaged network (stable edges only)
#   - edge_list.csv        : edges with bootstrap frequency
#   - egfr_neighbourhood.rds
#   - summary.txt
#
# Also writes:
#   - results/bnlearn/comparison_edges.csv : edge overlap across all 3 methods
# =============================================================================

source(file.path("R", "00_utils.R"))

check_packages(c("bnlearn"))

suppressPackageStartupMessages({
  library(bnlearn)
})

# =============================================================================
# Single algorithm runner
# =============================================================================

#' Run one bnlearn algorithm end-to-end
#'
#' @param expr      genes × samples matrix (transposed internally)
#' @param algo      algorithm id: "hc", "tabu", or "mmhc"
#' @param cfg       pipeline config
#' @return list(dag, boot, egfr)
run_single_bnlearn <- function(expr, algo, cfg) {
  bn_cfg  <- cfg$bn_learn
  stab_cfg <- bn_cfg$stability

  out_dir <- results_path(cfg, "bnlearn", algo)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  data <- as.data.frame(t(expr))

  # Sanitise names — bnlearn rejects special characters
  original_names <- colnames(data)
  safe_names     <- make.names(colnames(data), unique = TRUE)
  colnames(data) <- safe_names
  name_map <- setNames(original_names, safe_names)

  log_info("[", toupper(algo), "] Learning DAG from ",
           nrow(data), " samples × ", ncol(data), " genes")

  # -------------------------------------------------------------------------
  # Single-run structure learning
  # -------------------------------------------------------------------------
  dag <- switch(algo,
    hc = hc(
      data,
      score    = bn_cfg$score,
      restart  = bn_cfg$restart,
      perturb  = bn_cfg$perturb,
      max.iter = bn_cfg$max_iter
    ),
    tabu = tabu(
      data,
      score    = bn_cfg$score,
      tabu     = bn_cfg$tabu_length %||% 10L,
      max.iter = bn_cfg$max_iter
    ),
    mmhc = mmhc(data),
    stop("Unknown bnlearn algorithm: ", algo)
  )

  n_edges <- nrow(arcs(dag))
  log_info("[", toupper(algo), "] Learned DAG: ",
           length(nodes(dag)), " nodes, ", n_edges, " edges")

  attr(dag, "name_map") <- name_map
  saveRDS(dag, file.path(out_dir, "dag.rds"))

  # -------------------------------------------------------------------------
  # Bootstrap stability
  # -------------------------------------------------------------------------
  log_info("[", toupper(algo), "] Bootstrap: ",
           stab_cfg$n_bootstrap, " resamples")

  algo_args <- switch(algo,
    hc   = list(score = bn_cfg$score, restart = bn_cfg$restart,
                perturb = bn_cfg$perturb),
    tabu = list(score = bn_cfg$score,
                tabu = bn_cfg$tabu_length %||% 10L),
    mmhc = list()
  )

  boot_strength <- boot.strength(
    data,
    R              = stab_cfg$n_bootstrap,
    algorithm      = algo,
    algorithm.args = algo_args
  )

  avg_net <- averaged.network(boot_strength, threshold = stab_cfg$edge_threshold)
  log_info("[", toupper(algo), "] Stable edges (freq >= ",
           stab_cfg$edge_threshold, "): ", nrow(arcs(avg_net)))

  saveRDS(boot_strength, file.path(out_dir, "bootstrap.rds"))
  saveRDS(avg_net,       file.path(out_dir, "averaged.rds"))

  # -------------------------------------------------------------------------
  # Edge list CSV
  # -------------------------------------------------------------------------
  edge_df <- export_edge_list_bnlearn(avg_net, boot_strength, name_map,
                                       file.path(out_dir, "edge_list.csv"))

  # -------------------------------------------------------------------------
  # EGFR neighbourhood
  # -------------------------------------------------------------------------
  egfr_nhd <- extract_egfr_neighbourhood(
    avg_net,
    hops   = cfg$validation$egfr_neighbourhood$hops,
    name_map = name_map
  )
  if (!is.null(egfr_nhd)) {
    saveRDS(egfr_nhd, file.path(out_dir, "egfr_neighbourhood.rds"))
  }

  # -------------------------------------------------------------------------
  # Per-algorithm summary
  # -------------------------------------------------------------------------
  write_bnlearn_algo_summary(dag, avg_net, boot_strength, egfr_nhd,
                              algo, stab_cfg, out_dir)

  list(dag = dag, boot = boot_strength, avg = avg_net,
       egfr = egfr_nhd, name_map = name_map, edge_df = edge_df)
}

# =============================================================================
# Cross-method edge comparison
# =============================================================================

#' Build a long-format comparison table: which edges appear in which methods
build_edge_comparison <- function(results_by_algo) {
  # Collect all edges seen across any method
  all_edges <- do.call(rbind, lapply(names(results_by_algo), function(algo) {
    df <- results_by_algo[[algo]]$edge_df
    if (is.null(df) || nrow(df) == 0) return(NULL)
    # Use original (unmapped) names if available
    from_col <- if ("from_original" %in% colnames(df)) "from_original" else "from"
    to_col   <- if ("to_original" %in% colnames(df)) "to_original" else "to"
    data.frame(
      from      = df[[from_col]],
      to        = df[[to_col]],
      algorithm = algo,
      strength  = df$strength %||% NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(all_edges)) return(NULL)

  # Pivot: one row per edge, one column per algorithm's bootstrap strength
  edge_keys <- unique(paste(all_edges$from, "->", all_edges$to))
  algos <- names(results_by_algo)

  comparison <- do.call(rbind, lapply(edge_keys, function(key) {
    parts <- strsplit(key, " -> ")[[1]]
    row <- data.frame(from = parts[1], to = parts[2], stringsAsFactors = FALSE)
    for (a in algos) {
      m <- all_edges[all_edges$algorithm == a &
                       paste(all_edges$from, "->", all_edges$to) == key, ]
      row[[paste0("strength_", a)]] <- if (nrow(m) > 0) m$strength[1] else 0
    }
    row$n_methods <- sum(sapply(algos, function(a) {
      (row[[paste0("strength_", a)]] %||% 0) >= 0.01
    }))
    row
  }))

  # Sort: edges present in most methods first, then by HC strength
  comparison <- comparison[order(-comparison$n_methods,
                                  -comparison$strength_hc %||% 0), ]
  comparison
}

# =============================================================================
# Helpers
# =============================================================================

#' Export edge list with bootstrap frequencies
export_edge_list_bnlearn <- function(dag, boot_result, name_map, out_path) {
  if (nrow(arcs(dag)) == 0) {
    log_warn("No edges in averaged network — edge list will be empty")
    write.csv(data.frame(from = character(), to = character(),
                          strength = numeric(), direction = numeric()),
              out_path, row.names = FALSE)
    return(invisible(NULL))
  }

  edges <- as.data.frame(arcs(dag), stringsAsFactors = FALSE)
  colnames(edges) <- c("from", "to")

  # Merge bootstrap strength
  strength_df <- as.data.frame(boot_result, stringsAsFactors = FALSE)
  edges <- merge(edges, strength_df,
                 by = c("from", "to"), all.x = TRUE)

  # Map back to original gene names
  if (!is.null(name_map)) {
    edges$from_original <- name_map[edges$from]
    edges$to_original   <- name_map[edges$to]
  }

  edges <- edges[order(-edges$strength), ]
  write.csv(edges, out_path, row.names = FALSE)
  log_info("  Exported ", nrow(edges), " edges → ", basename(out_path))
  invisible(edges)
}

#' Extract the local neighbourhood around EGFR from a bnlearn DAG
#'
#' Works with sanitised node names by using name_map to find EGFR.
extract_egfr_neighbourhood <- function(dag, hops = 1, target = NULL,
                                        name_map = NULL) {
  all_nodes <- nodes(dag)

  # Find EGFR — first try safe names, then map back from original names
  if (is.null(target)) {
    target <- grep("^EGFR$", all_nodes, value = TRUE, ignore.case = TRUE)
    if (length(target) == 0 && !is.null(name_map)) {
      # Find safe name whose original is EGFR
      egfr_safe <- names(name_map)[name_map == "EGFR"]
      target <- intersect(egfr_safe, all_nodes)
    }
    if (length(target) == 0) {
      log_warn("EGFR node not found in DAG")
      return(NULL)
    }
    target <- target[1]
  }

  original_target <- if (!is.null(name_map)) name_map[target] else target
  log_info("  EGFR node: '", target, "' (original: '", original_target, "')")

  parents_nodes  <- bnlearn::parents(dag, target)
  children_nodes <- bnlearn::children(dag, target)

  neighbourhood <- unique(c(target, parents_nodes, children_nodes))
  if (hops >= 2) {
    for (h in seq_len(hops - 1)) {
      new_nodes <- unlist(lapply(neighbourhood, function(n) {
        c(bnlearn::parents(dag, n), bnlearn::children(dag, n))
      }))
      neighbourhood <- unique(c(neighbourhood, new_nodes))
    }
  }

  sub_arcs <- arcs(dag)
  sub_arcs <- sub_arcs[sub_arcs[, 1] %in% neighbourhood &
                          sub_arcs[, 2] %in% neighbourhood, , drop = FALSE]

  # Map to original names for readability
  to_original <- function(x) {
    if (!is.null(name_map)) name_map[x] %||% x else x
  }

  list(
    target           = original_target,
    parents          = to_original(parents_nodes),
    children         = to_original(children_nodes),
    neighbourhood    = to_original(neighbourhood),
    subgraph_edges   = sub_arcs,
    n_edges          = nrow(sub_arcs)
  )
}

write_bnlearn_algo_summary <- function(dag, avg, boot, egfr, algo,
                                        stab_cfg, out_dir) {
  lines <- c(
    "========================================",
    paste0("bnlearn — ", toupper(algo)),
    paste("Date:", Sys.time()),
    "========================================",
    "",
    paste("Single-run edges:", nrow(arcs(dag))),
    paste("Bootstrap resamples:", stab_cfg$n_bootstrap),
    paste("Edge threshold:", stab_cfg$edge_threshold),
    paste("Stable edges:", nrow(arcs(avg))),
    ""
  )
  if (!is.null(egfr)) {
    lines <- c(lines,
      "--- EGFR 1-hop neighbourhood ---",
      paste("Parents: ", paste(egfr$parents, collapse = ", ")),
      paste("Children:", paste(egfr$children, collapse = ", ")),
      paste("Neighbourhood size:", length(egfr$neighbourhood)),
      paste("Subgraph edges:", egfr$n_edges)
    )
  } else {
    lines <- c(lines, "EGFR: not found in DAG")
  }
  writeLines(lines, file.path(out_dir, "summary.txt"))
}

# =============================================================================
# Main
# =============================================================================

run_bn_bnlearn <- function(cfg) {
  log_stage_start("05 — BN Learning (bnlearn: HC, Tabu, MMHC)")
  start_time <- Sys.time()

  expr <- readRDS(processed_path(cfg, "04_expr_bn_input.rds"))
  log_info("Input: ", nrow(expr), " genes × ", ncol(expr), " samples")

  algorithms <- c("hc", "tabu", "mmhc")
  results <- list()

  for (algo in algorithms) {
    log_info("")
    log_info(strrep("-", 50))
    log_info("Algorithm: ", toupper(algo))
    log_info(strrep("-", 50))

    results[[algo]] <- tryCatch(
      run_single_bnlearn(expr, algo, cfg),
      error = function(e) {
        log_error("Algorithm ", toupper(algo), " failed: ", conditionMessage(e))
        NULL
      }
    )
  }

  # Cross-method edge comparison
  valid_results <- Filter(Negate(is.null), results)
  if (length(valid_results) > 1) {
    comparison <- build_edge_comparison(valid_results)
    if (!is.null(comparison)) {
      comp_path <- results_path(cfg, "bnlearn", "comparison_edges.csv")
      dir.create(dirname(comp_path), recursive = TRUE, showWarnings = FALSE)
      write.csv(comparison, comp_path, row.names = FALSE)
      log_info("Edge comparison table → ", basename(comp_path),
               " (", nrow(comparison), " unique edges)")

      # Log edges agreed on by all methods
      n_agree_all <- sum(comparison$n_methods == length(valid_results))
      log_info("Edges in ALL ", length(valid_results), " methods: ", n_agree_all)
    }
  }

  log_stage_end("05 — BN Learning (bnlearn)", start_time)
  invisible(results)
}

# =============================================================================
# Run if called directly
# =============================================================================
if (sys.nframe() == 0) {
  cfg <- load_config()
  run_bn_bnlearn(cfg)
}
