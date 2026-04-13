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
  # Skip if network access is unavailable (saves ~6 minutes of timeouts)
  skip_network <- isTRUE(cfg$ingest$skip_network_approaches)
  if (is.null(result) && !skip_network) {
    result <- try_tcgabiolinks_prepare(gdc_root, project)
  } else if (skip_network) {
    log_info("Approach 2 (TCGAbiolinks): skipped (skip_network_approaches = true)")
  }

  # Approach 3: Read raw count files directly from the GDC directory structure
  if (is.null(result)) {
    result <- try_read_raw_files(cohort_dir, project)
  }

  if (is.null(result)) {
    # Log directory structure to help diagnose what's actually there
    log_error("All ingestion approaches failed for: ", cohort_dir)
    log_error("Directory contents (2 levels):")
    dirs <- list.dirs(cohort_dir, recursive = TRUE)
    dirs <- dirs[nchar(dirs) - nchar(gsub("/|\\\\", "", dirs)) <= 2]
    for (d in head(dirs, 30)) log_error("  ", d)
    stop("Could not ingest data from ", cohort_dir,
         "\n  See the directory listing above.",
         "\n  Update config.yml paths.gdcdata_root if the path is wrong,",
         "\n  or check the BNPipeline download completed successfully.")
  }

  # -------------------------------------------------------------------------
  # Augment with clinical metadata (barcodes, sample type, IDH etc.)
  # -------------------------------------------------------------------------
  result <- augment_clinical(result, cohort_dir, gdc_root, cfg)

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
#' Tries multiple workflow type strings since naming differs across GDC releases.
try_tcgabiolinks_prepare <- function(gdc_root, project) {
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    log_info("Approach 2 (TCGAbiolinks): package not installed, skipping")
    return(NULL)
  }

  log_info("Approach 2: attempting TCGAbiolinks::GDCprepare() on local files")

  # GDC has used different workflow type strings across releases
  workflow_types <- c(
    "STAR - Counts",
    "HTSeq - Counts",
    "HTSeq - FPKM",
    "HTSeq - FPKM-UQ"
  )

  for (wf in workflow_types) {
    log_info("  Trying workflow.type = '", wf, "'")
    result <- tryCatch({
      query <- TCGAbiolinks::GDCquery(
        project       = project,
        data.category = "Transcriptome Profiling",
        data.type     = "Gene Expression Quantification",
        workflow.type = wf
      )
      se <- TCGAbiolinks::GDCprepare(query = query, directory = gdc_root)
      log_info("TCGAbiolinks::GDCprepare() succeeded with workflow '", wf, "'")
      se_to_list(se)
    }, error = function(e) {
      log_info("  '", wf, "' failed: ", conditionMessage(e))
      NULL
    })
    if (!is.null(result)) return(result)
  }

  log_warn("Approach 2: all workflow types failed")
  NULL
}

#' Approach 3: Read raw STAR/HTSeq count files from GDC directory structure
#'
#' GDCdata/ contains multiple data categories (RNA-seq, CNV, SNV, clinical).
#' This function scopes to the Transcriptome_Profiling subdirectory first,
#' then validates each candidate file by header content before reading.
try_read_raw_files <- function(cohort_dir, project) {
  log_info("Approach 3: scanning for RNA-seq count files")

  # -------------------------------------------------------------------------
  # Step 1: Scope the search directory
  # -------------------------------------------------------------------------
  # GDC layout (harmonized): TCGA-GBM/harmonized/Transcriptome_Profiling/
  #                                                 Gene_Expression_Quantification/
  # Older layout: TCGA-GBM/Transcriptome_Profiling/Gene_Expression_Quantification/
  candidate_roots <- c(
    file.path(cohort_dir, "harmonized", "Transcriptome_Profiling",
              "Gene_Expression_Quantification"),
    file.path(cohort_dir, "Transcriptome_Profiling",
              "Gene_Expression_Quantification"),
    file.path(cohort_dir, "harmonized", "Transcriptome_Profiling"),
    file.path(cohort_dir, "Transcriptome_Profiling"),
    cohort_dir   # last resort: search everywhere
  )

  search_root <- NULL
  for (r in candidate_roots) {
    if (dir.exists(r)) {
      search_root <- r
      log_info("  Scoping search to: ", r)
      break
    }
  }
  if (is.null(search_root)) {
    log_info("Approach 3: Transcriptome_Profiling directory not found")
    return(NULL)
  }

  # Log top-level directory contents to help diagnose structure
  top_dirs <- list.dirs(search_root, recursive = FALSE)
  log_info("  Top-level subdirectories (first 5): ",
           paste(basename(head(top_dirs, 5)), collapse = ", "))

  # -------------------------------------------------------------------------
  # Step 2: Find all .tsv/.txt/.counts files in scope
  # -------------------------------------------------------------------------
  all_files <- list.files(
    search_root,
    pattern = "\\.(counts|tsv|txt)(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE
  )

  # Exclude clearly non-count files by filename patterns
  exclude_pat <- "(junction|log|summary|manifest|segment|cnv|maf|vcf|clinical)"
  all_files <- all_files[!grepl(exclude_pat, all_files, ignore.case = TRUE)]

  if (length(all_files) == 0) {
    log_info("Approach 3: no candidate files found under Transcriptome_Profiling")
    return(NULL)
  }

  log_info("  Found ", length(all_files), " candidate files — validating format")

  # -------------------------------------------------------------------------
  # Step 3: Validate file format by inspecting headers
  # -------------------------------------------------------------------------
  # A valid count file must:
  #   (a) Have a first data column that looks like gene IDs (ENSG.../gene names)
  #   (b) Have at least one numeric column
  # CNV segment files have GDC_Aliquot/Chromosome/Start/End in the header.
  count_files <- character(0)
  detected_format <- NULL

  for (f in all_files) {
    header <- tryCatch(readLines(f, n = 6), error = function(e) character(0))
    if (length(header) == 0) next

    header_flat <- paste(header, collapse = "\t")

    # Reject files with CNV/SNV/clinical signatures
    if (grepl("Chromosome|Segment_Mean|GDC_Aliquot|Hugo_Symbol|bcr_patient",
              header_flat, ignore.case = TRUE)) next

    # Accept files that look like STAR gene counts
    if (grepl("unstranded|stranded_first|tpm_unstranded|gene_id.*gene_name",
              header_flat, ignore.case = TRUE)) {
      count_files <- c(count_files, f)
      detected_format <- "star"
      next
    }

    # Accept files with ENSG IDs in first column (HTSeq-style or bare STAR)
    first_col_vals <- gsub("\t.*", "", header)
    if (any(grepl("^(ENSG|__)", first_col_vals))) {
      count_files <- c(count_files, f)
      if (is.null(detected_format)) detected_format <- "htseq"
    }
  }

  if (length(count_files) == 0) {
    log_warn("Approach 3: no valid count files identified after header validation")
    log_warn("  Sample headers from first file (", basename(all_files[1]), "):")
    log_warn("  ", paste(readLines(all_files[1], n = 3), collapse = " | "))
    return(NULL)
  }

  log_info("  Validated ", length(count_files), " count files, format: ",
           detected_format)

  # Use one-file-per-UUID: if multiple files per UUID, take the first
  uuids <- basename(dirname(count_files))
  count_files <- count_files[!duplicated(uuids)]
  log_info("  Unique samples (UUIDs): ", length(count_files))

  # -------------------------------------------------------------------------
  # Step 4: Read based on detected format
  # -------------------------------------------------------------------------
  if (detected_format == "star") {
    read_star_counts(count_files)
  } else {
    read_htseq_counts(count_files)
  }
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
# Clinical data augmentation
# =============================================================================

# =============================================================================
# GDC API UUID resolution
# =============================================================================

#' Resolve file UUIDs to TCGA case barcodes via the GDC REST API
#'
#' Sends a batch POST request to https://api.gdc.cancer.gov/files.
#' Returns a named character vector: file_uuid -> case_submitter_id
#' (e.g. "TCGA-02-0003").
#'
#' Returns NULL silently if httr/jsonlite are not installed or if the
#' network is unavailable. The caller should fall back to UUIDs in that case.
#'
#' @param uuids  Character vector of GDC file UUIDs
#' @return Named character vector or NULL
resolve_uuids_via_gdc_api <- function(uuids) {
  if (!requireNamespace("httr",     quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)) {
    log_warn("  httr/jsonlite not installed — skipping GDC API UUID resolution")
    return(NULL)
  }

  log_info("  Querying GDC API for ", length(uuids), " file UUIDs ...")

  body <- list(
    filters = list(
      op      = "in",
      content = list(field = "file_id", value = as.list(uuids))
    ),
    fields = "file_id,cases.submitter_id",
    size   = length(uuids) + 10
  )

  resp <- tryCatch({
    httr::POST(
      "https://api.gdc.cancer.gov/files",
      httr::content_type_json(),
      body    = jsonlite::toJSON(body, auto_unbox = TRUE),
      httr::timeout(45)
    )
  }, error = function(e) {
    log_warn("  GDC API request error: ", conditionMessage(e))
    NULL
  })

  if (is.null(resp)) return(NULL)

  if (httr::http_error(resp)) {
    log_warn("  GDC API HTTP error: ", httr::status_code(resp))
    return(NULL)
  }

  raw_content <- httr::content(resp, as = "text", encoding = "UTF-8")
  parsed <- tryCatch(
    jsonlite::fromJSON(raw_content, flatten = TRUE),
    error = function(e) {
      log_warn("  GDC API JSON parse error: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(parsed)) return(NULL)

  hits <- parsed$data$hits
  if (is.null(hits) || nrow(hits) == 0) {
    log_warn("  GDC API returned 0 hits")
    return(NULL)
  }

  # cases.submitter_id is a list column when >1 case per file (rare for GBM)
  file_ids <- hits$file_id
  case_ids <- if ("cases.submitter_id" %in% colnames(hits)) {
    vapply(hits$cases.submitter_id, function(x) {
      if (is.null(x) || length(x) == 0) NA_character_
      else as.character(x[[1]][1])
    }, character(1))
  } else {
    rep(NA_character_, nrow(hits))
  }

  lookup <- setNames(case_ids, file_ids)
  lookup
}

#' Augment the minimal clinical stub (uuid + file_path) with real TCGA metadata.
#'
#' Strategy:
#'  1. Extract TCGA barcodes from STAR filenames (often embedded)
#'  2. Look for clinical TSV files in the GDCdata directory
#'  3. Parse sample type from barcode (position 14 = sample type code)
#'
#' @param result  list(counts, clinical, gene_info) from read_star_counts()
#' @param cohort_dir  Path to GDCdata/TCGA-GBM
#' @param gdc_root    Path to GDCdata root
#' @return result with enriched clinical data.frame
augment_clinical <- function(result, cohort_dir, gdc_root, cfg = NULL) {
  clinical <- result$clinical

  # -------------------------------------------------------------------------
  # Step 1: Extract TCGA barcode from the STAR count filename
  # GDC STAR files are named: TCGA-XX-XXXX-YYA-ZZR-NNNN-NN.*.tsv
  # The barcode IS the filename minus the extension suffix
  # -------------------------------------------------------------------------
  file_basenames <- basename(clinical$file_path)
  # Try to extract TCGA barcode pattern from filename
  barcode_match <- regmatches(
    file_basenames,
    regexpr("TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}[A-Z]-[0-9]{2}[A-Z]-[0-9]{4}-[0-9]{2}",
            file_basenames, ignore.case = FALSE)
  )
  n_barcodes <- sum(nchar(barcode_match) > 0)

  if (n_barcodes > 0) {
    log_info("Extracted TCGA barcodes from filenames: ", n_barcodes, "/",
             nrow(clinical), " samples")
    clinical$barcode <- ifelse(nchar(barcode_match) > 0, barcode_match, NA_character_)
  } else {
    log_info("No TCGA barcodes in filenames — trying GDC API UUID resolution")
    clinical$barcode <- NA_character_
    # Attempt to resolve file UUIDs to TCGA barcodes via the GDC REST API
    gdc_map <- resolve_uuids_via_gdc_api(clinical$uuid)
    if (!is.null(gdc_map)) {
      clinical$barcode <- gdc_map[clinical$uuid]
      n_barcodes <- sum(!is.na(clinical$barcode))
      log_info("GDC API resolved: ", n_barcodes, "/", nrow(clinical), " UUIDs to barcodes")
    } else {
      log_warn("GDC API resolution failed — sample IDs will remain as UUIDs")
    }
  }

  # -------------------------------------------------------------------------
  # Step 2: Parse sample type from barcode position 14
  # TCGA barcode: TCGA-XX-XXXX-YY[A-Z]-...
  #   YY = sample type code: 01 = Primary Tumor, 11 = Normal, 06 = Metastatic
  # -------------------------------------------------------------------------
  if (any(!is.na(clinical$barcode))) {
    sample_code <- substr(clinical$barcode, 14, 15)
    clinical$sample_type <- dplyr_recode_sample_type(sample_code)
    clinical$shortLetterCode <- sample_code_to_letter(sample_code)
    type_tab <- table(clinical$sample_type, useNA = "ifany")
    log_info("Sample types from barcode: ",
             paste(names(type_tab), "=", type_tab, collapse = ", "))
  }

  # Rename count matrix columns from UUIDs to barcodes where available
  if (any(!is.na(clinical$barcode))) {
    new_names <- ifelse(!is.na(clinical$barcode),
                        clinical$barcode, clinical$uuid)
    colnames(result$counts) <- new_names
  }

  # -------------------------------------------------------------------------
  # Step 3: Look for a clinical TSV in the GDCdata directory (or config path)
  # -------------------------------------------------------------------------
  cfg_clinical_path <- if (!is.null(cfg)) cfg$paths$clinical_tsv else NULL
  clinical_tsv <- find_clinical_tsv(cohort_dir, gdc_root, cfg_clinical_path)

  if (!is.null(clinical_tsv)) {
    log_info("Loading clinical data from: ", basename(clinical_tsv))
    clin_df <- tryCatch(
      read.delim(clinical_tsv, stringsAsFactors = FALSE, na.strings = c("", "NA", "'--")),
      error = function(e) {
        log_warn("Failed to read clinical TSV: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(clin_df)) {
      log_info("Clinical TSV: ", nrow(clin_df), " rows, ", ncol(clin_df), " columns")
      log_info("Clinical columns (first 20): ",
               paste(head(colnames(clin_df), 20), collapse = ", "))
      clinical <- merge_clinical_tsv(clinical, clin_df)
    }
  } else {
    log_warn("No clinical TSV found in GDCdata — ",
             "IDH status and full metadata will not be available. ",
             "Consider adding clinical data at: ",
             file.path(cohort_dir, "Clinical"))
  }

  result$clinical <- clinical
  result
}

#' Decode TCGA sample type code (positions 14-15 of barcode) to label
dplyr_recode_sample_type <- function(codes) {
  lookup <- c(
    "01" = "Primary Tumor",
    "02" = "Recurrent Tumor",
    "03" = "Primary Blood Derived Cancer",
    "06" = "Metastatic",
    "11" = "Solid Tissue Normal",
    "12" = "Buccal Cell Normal",
    "14" = "Bone Marrow Normal",
    "20" = "Control Analyte"
  )
  result <- lookup[codes]
  result[is.na(result)] <- paste0("Unknown (", codes[is.na(result)], ")")
  unname(result)
}

#' Decode TCGA sample type code to short letter code
sample_code_to_letter <- function(codes) {
  lookup <- c("01" = "TP", "02" = "TR", "06" = "TM", "11" = "NT",
              "03" = "TB", "12" = "NBC")
  result <- lookup[codes]
  result[is.na(result)] <- codes[is.na(result)]
  unname(result)
}

#' Search for clinical TSV files — checks config path first, then GDC locations
find_clinical_tsv <- function(cohort_dir, gdc_root, cfg_path = NULL) {
  # Priority 1: explicit path from config
  if (!is.null(cfg_path) && nchar(cfg_path) > 0 && file.exists(cfg_path)) {
    log_info("Using clinical TSV from config: ", cfg_path)
    return(cfg_path)
  }

  # Priority 2: common locations in GDC download structure
  search_dirs <- c(
    file.path(cohort_dir, "harmonized", "Clinical", "Clinical_Supplement"),
    file.path(cohort_dir, "harmonized", "Clinical"),
    file.path(cohort_dir, "Clinical"),
    file.path(gdc_root),
    dirname(cohort_dir)
  )

  for (d in search_dirs) {
    if (!dir.exists(d)) next
    tsv_files <- list.files(d, pattern = "\\.(tsv|txt|csv)$",
                             full.names = TRUE, recursive = FALSE)
    # Prefer files with "clinical" in the name
    clinical_files <- tsv_files[grepl("clinical|patient|case",
                                       tsv_files, ignore.case = TRUE)]
    if (length(clinical_files) > 0) return(clinical_files[1])
    if (length(tsv_files) > 0) return(tsv_files[1])
  }
  NULL
}

#' Merge clinical TSV into the sample clinical data.frame using barcode matching
merge_clinical_tsv <- function(clinical, clin_df) {
  # Find the barcode/case column in the clinical TSV
  barcode_col <- intersect(
    c("bcr_patient_barcode", "submitter_id", "case_id",
      "bcr_sample_barcode", "cases.submitter_id"),
    colnames(clin_df)
  )

  if (length(barcode_col) == 0) {
    log_warn("Could not find barcode column in clinical TSV. ",
             "Columns: ", paste(head(colnames(clin_df), 15), collapse = ", "))
    return(clinical)
  }

  barcode_col <- barcode_col[1]
  log_info("Joining on clinical column: '", barcode_col, "'")

  # For sample-level barcodes (16 chars), match to case barcode (12 chars)
  if (!is.na(clinical$barcode[1])) {
    case_from_sample <- substr(clinical$barcode, 1, 12)
    case_from_tsv    <- substr(clin_df[[barcode_col]], 1, 12)
    clin_df[[barcode_col]] <- case_from_tsv  # normalise to case level
    clinical$case_barcode <- case_from_sample
    merged <- merge(clinical, clin_df,
                    by.x = "case_barcode", by.y = barcode_col,
                    all.x = TRUE)
  } else {
    merged <- clinical
  }

  n_matched <- sum(!is.na(merged[[setdiff(colnames(clin_df),
                                           barcode_col)[1]]]))
  log_info("Clinical merge: ", n_matched, "/", nrow(clinical),
           " samples matched to clinical TSV")
  merged
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
