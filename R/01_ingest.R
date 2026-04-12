# =============================================================================
# 01_ingest.R — Read TCGA expression + clinical data from existing GDCdata
# =============================================================================
# This script treats the BNPipeline/GDCdata directory as a FIXED INPUT.
# It does NOT fetch from GDC. It reads whatever TCGAbiolinks downloaded
# and produces a clean SummarizedExperiment (or list) for downstream stages.
#
# Outputs (saved to data/processed/):
#   - 01_raw_counts.rds    : genes × samples integer count matrix
#   - 01_clinical.rds      : sample-level clinical metadata data.frame
#   - 01_ingest_summary.txt: human-readable summary of what was loaded
# =============================================================================

source(file.path("R", "00_utils.R"))

check_packages(c("SummarizedExperiment", "TCGAbiolinks", "GenomicRanges"))

suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

# =============================================================================
# Main ingestion function
# =============================================================================

#' Ingest TCGA data from a pre-downloaded GDCdata directory
#'
#' @param cfg Pipeline config (from load_config())
#' @return A list with components: counts, clinical, gene_info
ingest_tcga <- function(cfg) {
  log_stage_start("01 — Ingest TCGA data")
  start_time <- Sys.time()

  project   <- cfg$cohort$project
  gdc_root  <- cfg$paths$gdcdata_root
  cohort_dir <- gdcdata_cohort_path(cfg)

  log_info("Project:    ", project)
  log_info("GDCdata:    ", gdc_root)
  log_info("Cohort dir: ", cohort_dir)

  # -------------------------------------------------------------------------
  # Strategy: try multiple ingestion approaches in order of preference
  # -------------------------------------------------------------------------
  result <- NULL

  # Approach 1: Load a saved SummarizedExperiment / RangedSummarizedExperiment
  #   TCGAbiolinks::GDCprepare() often saves these as .rda/.rds files
  result <- try_load_se(cohort_dir, project)

  # Approach 2: Use TCGAbiolinks query + prepare on existing local files
  if (is.null(result)) {
    result <- try_tcgabiolinks_prepare(gdc_root, project)
  }

  # Approach 3: Read raw count files directly from the GDC directory structure
  if (is.null(result)) {
    result <- try_read_raw_files(cohort_dir, project)
  }

  if (is.null(result)) {
    stop("Could not ingest data from ", cohort_dir,
         "\n  Check that GDCdata was downloaded correctly in BNPipeline.")
  }

  # -------------------------------------------------------------------------
  # Filter to primary tumour samples
  # -------------------------------------------------------------------------
  result <- filter_sample_types(result, cfg$cohort$sample_types)

  # -------------------------------------------------------------------------
  # Save outputs
  # -------------------------------------------------------------------------
  save_output(result$counts,
              processed_path(cfg, "01_raw_counts.rds"),
              "raw count matrix")

  save_output(result$clinical,
              processed_path(cfg, "01_clinical.rds"),
              "clinical metadata")

  if (!is.null(result$gene_info)) {
    save_output(result$gene_info,
                processed_path(cfg, "01_gene_info.rds"),
                "gene info")
  }

  write_ingest_summary(result, cfg)

  log_stage_end("01 — Ingest", start_time)
  invisible(result)
}

# =============================================================================
# Ingestion approaches
# =============================================================================

#' Approach 1: Look for a pre-saved SE object
try_load_se <- function(cohort_dir, project) {
  # Common patterns for saved SE objects
  candidates <- c(
    list.files(cohort_dir, pattern = "\\.rds$", full.names = TRUE,
               recursive = TRUE),
    list.files(cohort_dir, pattern = "\\.rda$", full.names = TRUE,
               recursive = TRUE)
  )
  # Also check one level up (some setups save SE alongside GDCdata/)
  parent_candidates <- list.files(
    dirname(cohort_dir),
    pattern = sprintf("(%s|se|summarized).*\\.(rds|rda)$", project),
    full.names = TRUE, ignore.case = TRUE
  )
  candidates <- c(candidates, parent_candidates)

  if (length(candidates) == 0) {
    log_info("Approach 1 (saved SE): no .rds/.rda files found")
    return(NULL)
  }

  for (f in candidates) {
    log_info("Trying to load: ", basename(f))
    obj <- tryCatch({
      if (grepl("\\.rds$", f)) readRDS(f) else {
        env <- new.env(); load(f, envir = env); as.list(env)[[1]]
      }
    }, error = function(e) NULL)

    if (inherits(obj, "SummarizedExperiment") ||
        inherits(obj, "RangedSummarizedExperiment")) {
      log_info("Loaded SummarizedExperiment from: ", basename(f))
      return(se_to_list(obj))
    }
  }
  log_info("Approach 1: no valid SE objects found in candidates")
  NULL
}

#' Approach 2: Use TCGAbiolinks to prepare from local GDCdata
try_tcgabiolinks_prepare <- function(gdc_root, project) {
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    log_info("Approach 2 (TCGAbiolinks): package not installed, skipping")
    return(NULL)
  }

  log_info("Approach 2: attempting TCGAbiolinks::GDCprepare() on local files")

  tryCatch({
    # Set GDCdata directory so TCGAbiolinks looks in the right place
    query <- TCGAbiolinks::GDCquery(
      project = project,
      data.category = "Transcriptome Profiling",
      data.type     = "Gene Expression Quantification",
      workflow.type  = "STAR - Counts"
    )

    se <- TCGAbiolinks::GDCprepare(
      query     = query,
      directory = gdc_root   # points to existing GDCdata root
    )

    log_info("TCGAbiolinks::GDCprepare() succeeded")
    se_to_list(se)
  }, error = function(e) {
    log_warn("Approach 2 failed: ", conditionMessage(e))
    NULL
  })
}

#' Approach 3: Read raw STAR count files from GDC directory structure
try_read_raw_files <- function(cohort_dir, project) {
  log_info("Approach 3: scanning for raw count files in GDC directory structure")

  # GDC stores files as: TCGA-GBM/<uuid>/<filename>.counts or .tsv
  count_files <- list.files(
    cohort_dir,
    pattern = "\\.(counts|tsv|txt)(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE
  )
  # Filter to likely STAR/HTSeq count files (exclude logs, junctions, etc.)
  count_files <- count_files[!grepl("(junction|log|summary|manifest)",
                                     count_files, ignore.case = TRUE)]

  if (length(count_files) == 0) {
    log_info("Approach 3: no count files found")
    return(NULL)
  }

  log_info("Found ", length(count_files), " candidate count files")

  # Read first file to determine format
  sample_file <- count_files[1]
  header_lines <- readLines(sample_file, n = 10)
  log_info("Sample file header: ", paste(header_lines[1:min(3, length(header_lines))],
                                          collapse = " | "))

  # Determine if STAR (has header with unstranded/stranded columns)
  # or HTSeq (two columns: gene_id \t count)
  is_star <- any(grepl("unstranded|stranded_first|stranded_second|tpm_unstranded",
                        header_lines, ignore.case = TRUE))

  if (is_star) {
    result <- read_star_counts(count_files)
  } else {
    result <- read_htseq_counts(count_files)
  }

  result
}

# =============================================================================
# File-reading helpers
# =============================================================================

#' Read STAR count files (GDC Harmonized format)
#' These have columns: gene_id, gene_name, gene_type,
#'   unstranded, stranded_first, stranded_second, tpm_unstranded, ...
read_star_counts <- function(files) {
  log_info("Reading STAR count files (GDC Harmonized format)")

  counts_list <- lapply(seq_along(files), function(i) {
    if (i %% 50 == 0) log_info("  Reading file ", i, "/", length(files))
    df <- read.delim(files[i], comment.char = "#", stringsAsFactors = FALSE)

    # Skip the first 4 rows if they're summary stats (N_unmapped, etc.)
    if (any(grepl("^N_", df[[1]]))) {
      df <- df[!grepl("^N_", df[[1]]), ]
    }

    # Extract the sample barcode from the directory structure
    # GDC path: .../TCGA-GBM/<uuid>/<filename>
    uuid <- basename(dirname(files[i]))

    # Use unstranded counts by default
    count_col <- if ("unstranded" %in% colnames(df)) "unstranded" else
                 if ("HTSeq...Counts" %in% colnames(df)) "HTSeq...Counts" else
                 colnames(df)[ncol(df)]

    data.frame(
      gene_id = df$gene_id %||% df[[1]],
      count   = as.integer(df[[count_col]]),
      stringsAsFactors = FALSE
    )
  })

  # Get gene info from first file
  first_df <- read.delim(files[1], comment.char = "#", stringsAsFactors = FALSE)
  first_df <- first_df[!grepl("^N_", first_df[[1]]), ]
  gene_info <- data.frame(
    gene_id   = first_df$gene_id %||% first_df[[1]],
    gene_name = first_df$gene_name %||% NA_character_,
    gene_type = first_df$gene_type %||% NA_character_,
    stringsAsFactors = FALSE
  )

  # Build count matrix
  gene_ids <- counts_list[[1]]$gene_id
  count_mat <- do.call(cbind, lapply(counts_list, `[[`, "count"))
  rownames(count_mat) <- gene_ids
  colnames(count_mat) <- basename(dirname(files))  # UUIDs as initial colnames

  log_info("Assembled count matrix: ", nrow(count_mat), " genes × ",
           ncol(count_mat), " samples")

  list(
    counts   = count_mat,
    clinical = data.frame(
      uuid = colnames(count_mat),
      file_path = files,
      stringsAsFactors = FALSE
    ),
    gene_info = gene_info
  )
}

#' Read HTSeq count files (two-column format)
read_htseq_counts <- function(files) {
  log_info("Reading HTSeq count files")

  counts_list <- lapply(files, function(f) {
    df <- read.delim(f, header = FALSE, stringsAsFactors = FALSE,
                     col.names = c("gene_id", "count"))
    # Remove summary lines (__no_feature, __ambiguous, etc.)
    df <- df[!grepl("^__", df$gene_id), ]
    df$count <- as.integer(df$count)
    df
  })

  gene_ids <- counts_list[[1]]$gene_id
  count_mat <- do.call(cbind, lapply(counts_list, `[[`, "count"))
  rownames(count_mat) <- gene_ids
  colnames(count_mat) <- basename(dirname(files))

  log_info("Assembled count matrix: ", nrow(count_mat), " genes × ",
           ncol(count_mat), " samples")

  list(
    counts   = count_mat,
    clinical = data.frame(
      uuid = colnames(count_mat),
      file_path = files,
      stringsAsFactors = FALSE
    ),
    gene_info = data.frame(gene_id = gene_ids, stringsAsFactors = FALSE)
  )
}

# =============================================================================
# Conversion + filtering helpers
# =============================================================================

#' Convert a SummarizedExperiment to our standard list format
se_to_list <- function(se) {
  # Get the count assay — try common names
  assay_names <- SummarizedExperiment::assayNames(se)
  count_assay <- if ("unstranded" %in% assay_names) "unstranded" else
                 if ("HTSeq - Counts" %in% assay_names) "HTSeq - Counts" else
                 if ("raw_count" %in% assay_names) "raw_count" else
                 assay_names[1]

  log_info("Using assay: '", count_assay, "'")

  counts <- SummarizedExperiment::assay(se, count_assay)
  clinical <- as.data.frame(SummarizedExperiment::colData(se))
  gene_info <- as.data.frame(SummarizedExperiment::rowData(se))

  log_info("SE dimensions: ", nrow(counts), " genes × ", ncol(counts),
           " samples")
  log_info("Clinical columns: ", ncol(clinical))

  list(counts = counts, clinical = clinical, gene_info = gene_info)
}

#' Filter to specified sample types (e.g. Primary Tumor only)
filter_sample_types <- function(result, keep_types) {
  clinical <- result$clinical

  # Try common column names for sample type
  type_col <- intersect(
    c("sample_type", "shortLetterCode", "definition",
      "sample_type_id", "Sample.Type"),
    colnames(clinical)
  )

  if (length(type_col) == 0) {
    log_warn("No sample_type column found in clinical data — ",
             "skipping sample type filter. ",
             "Available columns: ",
             paste(head(colnames(clinical), 20), collapse = ", "))
    return(result)
  }

  type_col <- type_col[1]
  log_info("Filtering on column '", type_col, "'")
  log_info("Sample types present: ",
           paste(names(table(clinical[[type_col]])), collapse = ", "))

  keep <- clinical[[type_col]] %in% keep_types
  n_before <- nrow(clinical)
  n_after  <- sum(keep)

  if (n_after == 0) {
    # Try short letter codes (TP = Primary Tumor)
    if (type_col != "shortLetterCode" &&
        "shortLetterCode" %in% colnames(clinical)) {
      keep <- clinical$shortLetterCode == "TP"
      n_after <- sum(keep)
      log_info("Retrying with shortLetterCode == 'TP': ", n_after, " samples")
    }
  }

  if (n_after == 0) {
    log_warn("No samples match filter — keeping all ", n_before, " samples")
    return(result)
  }

  log_info("Sample filter: ", n_before, " → ", n_after,
           " (kept ", paste(keep_types, collapse = "/"), ")")

  result$counts   <- result$counts[, keep, drop = FALSE]
  result$clinical <- clinical[keep, , drop = FALSE]
  result
}

# =============================================================================
# Summary report
# =============================================================================

write_ingest_summary <- function(result, cfg) {
  out_path <- processed_path(cfg, "01_ingest_summary.txt")
  lines <- c(
    "========================================",
    "01_ingest.R — Summary",
    paste("Date:", Sys.time()),
    paste("Cohort:", cfg$cohort$project),
    "========================================",
    "",
    paste("Genes:  ", nrow(result$counts)),
    paste("Samples:", ncol(result$counts)),
    "",
    "Clinical columns:",
    paste(" ", colnames(result$clinical), collapse = "\n"),
    "",
    "Sample barcodes (first 5):",
    paste(" ", head(colnames(result$counts), 5), collapse = "\n")
  )

  if (!is.null(result$gene_info) && "gene_type" %in% colnames(result$gene_info)) {
    gene_type_tab <- sort(table(result$gene_info$gene_type), decreasing = TRUE)
    lines <- c(lines, "", "Gene biotypes (top 10):",
               paste(" ", names(gene_type_tab)[1:min(10, length(gene_type_tab))],
                     ": ", gene_type_tab[1:min(10, length(gene_type_tab))]))
  }

  writeLines(lines, out_path)
  log_info("Wrote ingest summary → ", basename(out_path))
}

# =============================================================================
# Run if called directly
# =============================================================================
if (sys.nframe() == 0) {
  cfg <- load_config()
  result <- ingest_tcga(cfg)
}
