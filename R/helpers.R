### Helpers for amRdata live in this script
#########################
# Data curation helpers #
#########################
#' Helps ensure trailing 0s are retained in genome IDs for proper downloading
#' @keywords internal
.id_checker <- function(x) {
  # Taxon IDs are just numbers, genome IDs have decimals, this tells them apart
  grepl("^[0-9]+$", x)
}

#' A helper used in data_curation.R and data_processing.R to ensure exported tables
#' don't lose trailing zeroes. Should be relocated into a common helpers/utilities
#' script later.
#' @keywords internal
.preserve_export_id_text <- function(df) {
  df <- tibble::as_tibble(df)

  id_pattern <- paste0(
    "(^|[._])(",
    "genome(_drug)?_id|taxon_id|",
    "assembly_accession|bioproject_accession|biosample_accession|",
    "refseq_accessions?|genbank_accessions?|sra_accession|pmid|",
    "gene_id|protein_id|domain_id|cluster_id|AccNum|id",
    ")$"
  )

  id_cols <- names(df)[grepl(id_pattern, names(df), ignore.case = TRUE)]
  if (length(id_cols)) {
    df[id_cols] <- lapply(df[id_cols], as.character)
  }

  df
}

#' Helps normalize Docker paths
#' @keywords internal
.docker_path <- function(p) gsub("\\\\", "/", normalizePath(p, mustWork = FALSE))

#' Helps run a shell inside a container, and prefers bash (don't we all?)
#' @keywords internal
.pick_shell <- function(image) {
  chk <- suppressWarnings(system2("docker",
                                  c(
                                    "run", "--rm", image, "sh", "-lc",
                                    "command -v bash >/dev/null || echo NOBASH"
                                  ),
                                  stdout = TRUE, stderr = TRUE
  ))
  if (length(chk) && any(grepl("NOBASH", chk))) "sh" else "bash"
}

# FASTA sanitizer to ensure Panaroo compatibility with BV-BRC CLI downloads
.strip_fasta_preamble <- function(fna_path) {
  if (!file.exists(fna_path)) {
    return(invisible(FALSE))
  }
  txt <- readLines(fna_path, warn = FALSE)
  first <- which(grepl("^\\s*>", txt))[1]
  if (is.na(first)) {
    return(invisible(FALSE))
  }
  if (first > 1L) {
    txt <- txt[first:length(txt)]
    txt[1] <- sub("^\\ufeff", "", txt[1])
    writeLines(txt, fna_path, sep = "\n", useBytes = TRUE)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

# GFF sanitizer to ensure Panaroo compatibility with BV-BRC CLI downloads
.sanitize_gff <- function(gff_path) {
  if (!file.exists(gff_path)) {
    return(invisible(FALSE))
  }
  lines <- readLines(gff_path, warn = FALSE)
  if (length(lines) == 0L) {
    return(invisible(FALSE))
  }
  if (!grepl("^##gff-version\\s*3", lines[1])) {
    lines <- c("##gff-version 3", lines)
  }
  out <- purrr::map_chr(lines, function(line) {
    if (grepl("^#", line)) {
      return(line)
    }
    parts <- strsplit(line, "[\t ]", perl = TRUE)[[1]]
    if (length(parts) >= 9) {
      paste(c(parts[1:8], paste(parts[9:length(parts)], collapse = " ")), collapse = "\t")
    } else {
      line
    }
  })
  writeLines(out, gff_path, sep = "\n", useBytes = TRUE)
  invisible(TRUE)
}

#########################
#   Manifest helpers    #
#########################

#' Returns the basics about a file for manifest logging
#'
#' @param path Character vector of file paths.
#' @param hash Logical. If TRUE, calculate MD5 checksums.
#'
#' @return A list of file records.
#' @keywords internal
.manifest_file_info <- function(path, hash = FALSE) {
  path <- unique(as.character(path))
  path <- path[nzchar(path)]

  if (!length(path)) {
    return(list())
  }

  # See what exists
  purrr::map(path, function(x) {
    exists <- file.exists(x)

    out <- list(
      path = x,
      exists = exists,
      size_bytes = if (exists) file.info(x)$size else NA_real_,
      modified_at = if (exists) as.character(file.info(x)$mtime) else NA_character_
    )

    # Hash what exists, if desired
    if (isTRUE(hash) && exists && !dir.exists(x)) {
      out$md5 <- unname(tools::md5sum(x))
    }

    out
  })
}


#' Capture basic GitHub repo state for manifest provenance
#'
#' @param base_dir Character. Project root.
#'
#' @return A named list.
#' @keywords internal
.manifest_git_info <- function(base_dir = ".") {
  base_dir <- normalizePath(base_dir, mustWork = FALSE)

  # Find Git
  git <- Sys.which("git")

  if (!nzchar(git)) {
    return(list(
      available = FALSE
    ))
  }

  # Run Git through system commands
  run_git <- function(args) {
    tryCatch(
      system2(
        git,
        args = args,
        stdout = TRUE,
        stderr = FALSE
      ),
      error = function(e) character()
    )
  }

  inside <- run_git(c("-C", shQuote(base_dir), "rev-parse", "--is-inside-work-tree"))

  if (!length(inside) || !identical(trimws(inside[[1]]), "true")) {
    return(list(
      available = TRUE,
      repository = FALSE
    ))
  }

  commit <- run_git(c("-C", shQuote(base_dir), "rev-parse", "HEAD"))
  branch <- run_git(c("-C", shQuote(base_dir), "rev-parse", "--abbrev-ref", "HEAD"))
  dirty <- run_git(c("-C", shQuote(base_dir), "status", "--porcelain"))

  list(
    available = TRUE,
    repository = TRUE,
    commit = if (length(commit)) trimws(commit[[1]]) else NA_character_,
    branch = if (length(branch)) trimws(branch[[1]]) else NA_character_,
    dirty = length(dirty) > 0L
  )
}


#' Capture package versions currently loaded in the R session
#'
#' @return Named character vector of package versions.
#' @keywords internal
.manifest_package_versions <- function() {
  pkgs <- sort(loadedNamespaces())

  stats::setNames(
    as.list(
      purrr::map_chr(
        pkgs,
        function(pkg) {
          tryCatch(
            as.character(utils::packageVersion(pkg)),
            error = function(e) NA_character_
          )
        }
      )
    ),
    pkgs
  )
}


#' Generate a unique manifest run identifier
#'
#' @return Character scalar.
#' @keywords internal
.manifest_run_id <- function() {
  paste0(
    "run_",
    format(Sys.time(), "%Y%m%dT%H%M%OS3", tz = "UTC"),
    "_pid",
    Sys.getpid()
  ) |>
    gsub("[^A-Za-z0-9_]", "", x = _)
}


#' Start or load a dataset provenance manifest
#'
#' @param manifest_path Character. Path to the JSON manifest.
#' @param dataset_id Character scalar.
#' @param duckdb_path Character scalar.
#' @param base_dir Character scalar.
#' @param selection Optional named list describing the dataset selection.
#' @param hash_files Logical. Calculate SHA-256 for manifest-recorded files.
#'
#' @return A manifest object with `path` and `run_index`.
#' @keywords internal
.manifest_start <- function(
    manifest_path,
    dataset_id,
    duckdb_path,
    base_dir = ".",
    selection = list(),
    hash_files = FALSE
) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for manifest generation.")
  }

  manifest_path <- normalizePath(
    manifest_path,
    mustWork = FALSE
  )

  dir.create(
    dirname(manifest_path),
    recursive = TRUE,
    showWarnings = FALSE
  )

  manifest <- list(
    schema_version = 1L,
    manifest_created_at = as.character(Sys.time()),
    manifest_updated_at = as.character(Sys.time()),
    dataset_id = dataset_id,
    dataset = list(
      duckdb = duckdb_path,
      selection = selection
    ),
    runs = list()
  )

  run <- list(
    run_id = .manifest_run_id(),
    status = "running",
    started_at = as.character(Sys.time()),
    finished_at = NA_character_,
    command = commandArgs(trailingOnly = FALSE),
    working_directory = getwd(),
    host = as.list(Sys.info()),
    r = list(
      version = R.version.string,
      platform = R.version$platform
    ),
    git = .manifest_git_info(base_dir),
    packages = .manifest_package_versions(),
    stages = list(),
    events = list()
  )

  if (is.null(manifest$runs)) {
    manifest$runs <- list()
  }

  manifest$runs[[length(manifest$runs) + 1L]] <- run
  manifest$manifest_updated_at <- as.character(Sys.time())

  run_index <- length(manifest$runs)

  jsonlite::write_json(
    manifest,
    manifest_path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  structure(
    list(
      manifest = manifest,
      path = manifest_path,
      run_index = run_index,
      hash_files = isTRUE(hash_files)
    ),
    class = "amr_manifest"
  )
}


#' Update a manifest stage
#'
#' @param manifest_state Manifest state returned by [.manifest_start()].
#' @param name Character stage name.
#' @param status Character stage status.
#' @param parameters Optional named list.
#' @param inputs Optional character vector of input paths.
#' @param outputs Optional character vector of output paths.
#' @param tool Optional named list describing the tool.
#' @param metrics Optional named list of metrics.
#' @param message Optional log message.
#'
#' @return Updated manifest state.
#' @keywords internal
.manifest_stage <- function(
    manifest_state,
    name,
    status = "success",
    parameters = list(),
    inputs = character(),
    outputs = character(),
    tool = list(),
    metrics = list(),
    message = NULL
) {
  if (!inherits(manifest_state, "amr_manifest")) {
    stop("Invalid manifest state.")
  }

  stage_index <- which(
    purrr::map_lgl(
      manifest_state$manifest$runs[[manifest_state$run_index]]$stages,
      ~ identical(.x$name, name) && identical(.x$status, "running")
    )
  )

  stage <- list(
    name = name,
    status = status,
    started_at = as.character(Sys.time()),
    parameters = parameters,
    inputs = .manifest_file_info(inputs, hash = manifest_state$hash_files),
    outputs = .manifest_file_info(outputs, hash = manifest_state$hash_files),
    tool = tool,
    metrics = metrics
  )

  if (!is.null(message)) {
    stage$message <- as.character(message)
  }

  if (length(stage_index) == 1L) {
    existing <- manifest_state$manifest$runs[[manifest_state$run_index]]$stages[[stage_index]]

    stage$started_at <- existing$started_at
    stage$finished_at <- if (status != "running") {
      as.character(Sys.time())
    } else {
      NULL
    }

    manifest_state$manifest$runs[[manifest_state$run_index]]$stages[[stage_index]] <- stage
  } else {
    if (status != "running") {
      stage$finished_at <- as.character(Sys.time())
    }

    manifest_state$manifest$runs[[manifest_state$run_index]]$stages <-
      append(
        manifest_state$manifest$runs[[manifest_state$run_index]]$stages,
        list(stage)
      )
  }

  manifest_state$manifest$manifest_updated_at <- as.character(Sys.time())

  jsonlite::write_json(
    manifest_state$manifest,
    manifest_state$path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  manifest_state
}


#' Append a provenance event to the active manifest run
#'
#' @param manifest_state Manifest state returned by [.manifest_start()].
#' @param level Character event level.
#' @param message Character message.
#' @param details Optional named list.
#'
#' @return Updated manifest state.
#' @keywords internal
.manifest_event <- function(
    manifest_state,
    level = "info",
    message,
    details = list()
) {
  manifest_state$manifest$runs[[manifest_state$run_index]]$events <-
    append(
      manifest_state$manifest$runs[[manifest_state$run_index]]$events,
      list(
        list(
          timestamp = as.character(Sys.time()),
          level = level,
          message = message,
          details = details
        )
      )
    )

  manifest_state$manifest$manifest_updated_at <- as.character(Sys.time())

  jsonlite::write_json(
    manifest_state$manifest,
    manifest_state$path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  manifest_state
}


#' Finish an active provenance manifest run
#'
#' @param manifest_state Manifest state returned by [.manifest_start()].
#' @param status Final run status.
#' @param error Optional error message.
#'
#' @return Invisibly returns the final manifest state.
#' @keywords internal
.manifest_finish <- function(
    manifest_state,
    status = "success",
    error = NULL
) {
  manifest_state$manifest$runs[[manifest_state$run_index]]$status <- status
  manifest_state$manifest$runs[[manifest_state$run_index]]$finished_at <-
    as.character(Sys.time())

  # Patching to resolve an indefinite `running` failure state in the manifest
  if (identical(status, "failed")) {
    stages <- manifest_state$manifest$runs[[manifest_state$run_index]]$stages
    running_stage <- which(purrr::map_lgl(stages, ~ identical(.x$status, "running")))

    if (length(running_stage)) {
      stage_error <- if (!is.null(error)) {
        as.character(error)
      } else {
        "Parent run failed before this stage completed."
      }

      for (i in running_stage) {
        stages[[i]]$status <- "failed"
        stages[[i]]$finished_at <- as.character(Sys.time())
        stages[[i]]$error <- stage_error
      }

      manifest_state$manifest$runs[[manifest_state$run_index]]$stages <- stages
    }
  }

  if (!is.null(error)) {
    manifest_state$manifest$runs[[manifest_state$run_index]]$error <- as.character(error)
  }

  manifest_state$manifest$manifest_updated_at <- as.character(Sys.time())

  jsonlite::write_json(
    manifest_state$manifest,
    manifest_state$path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  invisible(manifest_state)
}

# To distinguish multiple manifests in the same bug directory
.manifest_find_latest <- function(duckdb_path) {
  manifest_dir <- dirname(normalizePath(
    duckdb_path,
    mustWork = FALSE
  ))

  manifests <- list.files(
    manifest_dir,
    pattern = "^manifest_.*\\.json$",
    full.names = TRUE
  )

  if (!length(manifests)) {
    return(NULL)
  }

  manifests[which.max(file.info(manifests)$mtime)]
}

#' Resume provenance logging in an existing manifest
#'
#' Loads an existing manifest and appends a new run.
#'
#' @param manifest_path Character. Path to an existing JSON manifest.
#' @param base_dir Character. Project root.
#' @param hash_files Logical. Calculate SHA-256 checksums for manifest-recorded files.
#'
#' @return A manifest object with `path` and `run_index`.
#' @keywords internal
.manifest_resume <- function(
    manifest_path,
    base_dir = ".",
    hash_files = FALSE
) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for manifest generation.")
  }

  manifest_path <- normalizePath(
    manifest_path,
    mustWork = TRUE
  )

  manifest <- jsonlite::read_json(
    manifest_path,
    simplifyVector = FALSE
  )

  if (is.null(manifest$runs)) {
    manifest$runs <- list()
  }

  run <- list(
    run_id = .manifest_run_id(),
    status = "running",
    started_at = as.character(Sys.time()),
    finished_at = NA_character_,
    command = commandArgs(trailingOnly = FALSE),
    working_directory = getwd(),
    host = as.list(Sys.info()),
    r = list(
      version = R.version.string,
      platform = R.version$platform
    ),
    git = .manifest_git_info(base_dir),
    packages = .manifest_package_versions(),
    stages = list(),
    events = list()
  )

  manifest$runs[[length(manifest$runs) + 1L]] <- run
  manifest$manifest_updated_at <- as.character(Sys.time())

  run_index <- length(manifest$runs)

  jsonlite::write_json(
    manifest,
    manifest_path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )

  structure(
    list(
      manifest = manifest,
      path = manifest_path,
      run_index = run_index,
      hash_files = isTRUE(hash_files)
    ),
    class = "amr_manifest"
  )
}

###########################
# Data processing helpers #
###########################

# Map host paths under mounted root to container path
#' .to_container()
#'
#' Used for OS-agnostic mapping of Docker directories and mount paths
#'
#' @keywords internal
#' @examples NULL
.to_container <- function(x, host_root, container_root = "/work") {
  host_root_unix <- .docker_path(host_root)
  x_unix <- .docker_path(x)
  pattern <- paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\\\\\1", host_root_unix))
  sub(pattern, container_root, x_unix)
}

#' Remove pseudogene annotations from Panaroo input GFF files
#'
#' Cleans GFF annotation files of `pseudogene` feature records only. Cleaned
#' GFFs are written to a subdirectory under `output_path` and swapped into the
#' Panaroo input list, leaving the original genome annotations alone.
#'
#' This optional preprocessing step can reduce weird runtime stalls during
#' Panaroo graph construction for some BV-BRC/PATRIC genome annotations that
#' contain troublesome pseudogene features.
#'
#' @param panaroo_input_files Character vector of `"gff fna"` input lines used
#'   by Panaroo.
#' @param output_path Character scalar. Base directory for temporary cleaned
#'   GFF files and audit outputs.
#' @param clean_dir Character scalar. Name of the subdirectory created beneath
#'   `output_path` to store cleaned GFF files. Default `"gff_clean"`.
#'
#' @return A list containing:
#' \itemize{
#'   \item `panaroo_input_files` — rewritten Panaroo input lines pointing to the
#'   cleaned GFF files.
#'   \item `audit` — a tibble summarizing, for each genome, the total number of
#'   annotated features, the number of pseudogenes removed, and the number of
#'   remaining features.
#' }
#'
#' @details
#' This performs lightweight preprocessing only, removing feature records whose
#' third GFF column is exactly `"pseudogene"` and does not otherwise modify
#' annotation coordinates, attributes, or sequence files. FASTA paths are unchanged.
#'
#' @keywords internal
.stripPseudogeneGFFs <- function(panaroo_input_files,
                                 output_path,
                                 clean_dir = "gff_clean") {
  # Normalize our paths
  panaroo_input_files <- as.character(panaroo_input_files)
  output_path <- .docker_path(output_path)

  # Set the directory to place cleaned GFFs into
  clean_root <- file.path(output_path, clean_dir)
  dir.create(clean_root, recursive = TRUE, showWarnings = FALSE)

  # Where the cleaned up Panaroo input and clean audit is stored
  out_lines <- character(length(panaroo_input_files))
  audit <- vector("list", length(panaroo_input_files))

  for (i in seq_along(panaroo_input_files)) {

    # Read the Panaroo gff + fna input lines
    line <- panaroo_input_files[[i]]
    parts <- strsplit(line, "\\s+")[[1]]

    # If you're missing either a gff or an fna file in there somehow
    if (length(parts) < 2L) {
      stop("Broken Panaroo input line: ", line)
    }

    # Read in the files parsed above
    gff_in <- .docker_path(parts[1])
    fna_in <- .docker_path(parts[2])

    # If it didn't read in
    if (!file.exists(gff_in)) {
      stop("Missing GFF file: ", gff_in)
    }
    if (!file.exists(fna_in)) {
      stop("Missing FNA file: ", fna_in)
    }

    # What we're saving out
    gff_out <- file.path(clean_root, basename(gff_in))

    # Read in the GFF lines and fine the comment lined headers
    gff_lines <- readLines(gff_in, warn = FALSE)
    is_header <- startsWith(gff_lines, "#")
    body <- gff_lines[!is_header]

    # If there's nothing in there to parse
    if (length(body) == 0L) {
      writeLines(gff_lines, gff_out, useBytes = TRUE)
      n_total <- 0L
      n_pseudogene <- 0L
      n_kept <- 0L
    } else {
      # Otherwise, find the pseudogene lines and save everything but those
      # Now with 100% more purrr
      fields <- strsplit(body, "\t", fixed = TRUE)
      types <- purrr::map_chr(
        fields,
        \(x) if (length(x) >= 3L) x[[3]] else NA_character_
      )
      keep <- !is.na(types) & types != "pseudogene"

      cleaned <- c(gff_lines[is_header], body[keep])
      writeLines(cleaned, gff_out, useBytes = TRUE)

      # We love stats
      n_total <- length(body)
      n_pseudogene <- sum(!keep, na.rm = TRUE)
      n_kept <- sum(keep, na.rm = TRUE)
    }

    # Record what we did and save it into the audit log
    audit[[i]] <- tibble::tibble(
      gff_in = gff_in,
      gff_out = gff_out,
      n_total_features = n_total,
      n_pseudogene = n_pseudogene,
      n_kept = n_kept
    )

    out_lines[[i]] <- paste(gff_out, fna_in)
  }

  list(
    panaroo_input_files = out_lines,
    audit = dplyr::bind_rows(audit)
  )
}

#' Export a dyad-centric feature table
#'
#' Builds a one-row-per-dyad table linking protein-gene dyads to structural
#' features and HMMER annotations recorded in the dataset provenance manifest.
#' Feature values are deduplicated and combined into semicolon-separated
#' character fields.
#'
#' @param duckdb_path Character. Path to the source dataset DuckDB. Associated
#'   Parquet files and provenance manifest are expected there.
#' @param output_path Character or NULL. Directory where the output Parquet file
#'   will be written. Defaults to the DuckDB directory.
#' @param output_stem Character. Output filename stem. Default
#'   `"dyad_annotations"`.
#' @param feature_scales Character vector of optional feature types to include.
#'   If NULL, includes `struct` plus all HMMER databases recorded in the
#'   manifest.
#' @param verbose Logical. Print progress messages.
#'
#' @return Invisibly returns the path to the generated Parquet file.
#'
#' @keywords internal
.exportDyadAnnotations <- function(
    duckdb_path,
    feature_scales = NULL,
    verbose = TRUE
) {
  duckdb_path <- normalizePath(duckdb_path, mustWork = TRUE)
  parquet_dir <- dirname(duckdb_path)

  manifest_path <- .manifest_find_latest(duckdb_path)

  if (is.null(manifest_path)) {
    stop(
      "No provenance manifest found for: ",
      duckdb_path
    )
  }

  manifest <- jsonlite::read_json(
    manifest_path,
    simplifyVector = FALSE
  )

  # Find latest successful HMMER stage
  hmmer_stage <- NULL

  for (run in rev(manifest$runs %||% list())) {
    stages <- run$stages %||% list()

    matches <- purrr::keep(
      stages,
      ~ identical(.x$name, "hmmer") &&
        identical(.x$status, "success")
    )

    if (length(matches)) {
      hmmer_stage <- matches[[1]]
      break
    }
  }

  if (is.null(hmmer_stage)) {
    stop(
      "No successful HMMER stage found in manifest: ",
      manifest_path
    )
  }

  hmmer_databases <- unique(as.character(
    unlist(
      hmmer_stage$parameters$databases %||% character(),
      use.names = FALSE
    )
  ))

  allowed_features <- c("struct", hmmer_databases)

  if (is.null(feature_scales)) {
    feature_scales <- allowed_features
  } else {
    feature_scales <- unique(as.character(feature_scales))

    unknown_features <- setdiff(
      feature_scales,
      allowed_features
    )

    if (length(unknown_features)) {
      stop(
        "Unsupported feature scale(s): ",
        paste(unknown_features, collapse = ", "),
        ". Available features: ",
        paste(allowed_features, collapse = ", ")
      )
    }
  }

  con <- DBI::dbConnect(
    duckdb::duckdb(),
    dbdir = ":memory:"
  )

  duckdb_temp_dir <- file.path(
    parquet_dir,
    ".duckdb_temp"
  )

  dir.create(
    duckdb_temp_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  DBI::dbExecute(
    con,
    sprintf(
      "SET temp_directory=%s",
      DBI::dbQuoteString(
        con,
        normalizePath(
          duckdb_temp_dir,
          winslash = "/",
          mustWork = TRUE
        )
      )
    )
  )

  on.exit(
    {
      try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
      unlink(
        duckdb_temp_dir,
        recursive = TRUE,
        force = TRUE
      )
    },
    add = TRUE
  )

  parquet_sql <- function(dataset_name) {
    path <- file.path(
      parquet_dir,
      paste0(dataset_name, ".parquet")
    )

    if (!file.exists(path)) {
      return(NULL)
    }

    normalizePath(
      path,
      winslash = "/",
      mustWork = TRUE
    )
  }

  # Initialize using protein-gene dyads
  genome_gene_protein_path <- parquet_sql(
    "genome_gene_protein"
  )

  if (is.null(genome_gene_protein_path)) {
    stop(
      "Required Parquet file not found: ",
      file.path(
        parquet_dir,
        "genome_gene_protein.parquet"
      )
    )
  }

  DBI::dbExecute(
    con,
    sprintf(
      "
      CREATE OR REPLACE VIEW protein_gene AS
      SELECT DISTINCT
        protein_ids AS protein,
        REPLACE(Gene, '~', '.') AS gene
      FROM read_parquet('%s')
      WHERE protein_ids IS NOT NULL
        AND Gene IS NOT NULL
      ",
      genome_gene_protein_path
    )
  )

  DBI::dbExecute(
    con,
    "
    CREATE OR REPLACE VIEW protein_gene_dyad AS
    SELECT DISTINCT
      protein,
      gene,
      CONCAT(protein, '|', gene) AS dyad
    FROM protein_gene
    "
  )

  # Finish initializing with one row per dyad
  feature_select <- character()
  feature_joins <- character()

  # Pangenome graph structural variant ('struct') annotations
  if ("struct" %in% feature_scales) {
    struct_path <- parquet_sql("struct")

    if (is.null(struct_path)) {
      if (isTRUE(verbose)) {
        message(
          "Skipping struct: struct.parquet was not found."
        )
      }
    } else {
      DBI::dbExecute(
        con,
        sprintf(
          "
          CREATE OR REPLACE VIEW struct_genes AS
          SELECT DISTINCT
            s.struct,
            t.gene
          FROM read_parquet('%s') s
          CROSS JOIN UNNEST(
            string_split(s.struct, '.')
          ) AS t(gene)
          WHERE s.value = 1
          ",
          struct_path
        )
      )

      DBI::dbExecute(
        con,
        "
        CREATE OR REPLACE VIEW dyad_struct AS
        SELECT
          pgd.dyad,
          string_agg(
            DISTINCT sg.struct,
            ';'
            ORDER BY sg.struct
          ) AS struct
        FROM protein_gene_dyad pgd
        JOIN struct_genes sg
          ON pgd.gene = sg.gene
        GROUP BY pgd.dyad
        "
      )

      feature_select <- c(
        feature_select,
        "ds.struct"
      )

      feature_joins <- c(
        feature_joins,
        "LEFT JOIN dyad_struct ds ON b.dyad = ds.dyad"
      )
    }
  }

  # HMMER feature annotations
  for (database in intersect(
    hmmer_databases,
    feature_scales
  )) {
    dataset_name <- paste0(
      "protein_",
      database
    )

    hmmer_path <- parquet_sql(dataset_name)

    if (is.null(hmmer_path)) {
      if (isTRUE(verbose)) {
        message(
          "Skipping ",
          database,
          ": ",
          dataset_name,
          ".parquet was not found."
        )
      }
      next
    }

    view_name <- paste0(
      "dyad_",
      make.names(database)
    )

    DBI::dbExecute(
      con,
      sprintf(
        "
        CREATE OR REPLACE VIEW %s AS
        SELECT
          pgd.dyad,
          string_agg(
            DISTINCT h.query_name,
            ';'
            ORDER BY h.query_name
          ) AS feature
        FROM protein_gene_dyad pgd
        JOIN read_parquet('%s') h
          ON pgd.protein = h.protein
        WHERE h.query_name IS NOT NULL
        GROUP BY pgd.dyad
        ",
        view_name,
        hmmer_path
      )
    )

    # Give AMRFinder a better human-readable name (ARG, in this case)
    output_column <- if (identical(database, "AMRFinder")) {
      "ARG"
    } else {
      database
    }

    alias <- paste0(
      "d_",
      make.names(database)
    )

    feature_select <- c(
      feature_select,
      sprintf(
        '%s.feature AS "%s"',
        alias,
        output_column
      )
    )

    feature_joins <- c(
      feature_joins,
      sprintf(
        "LEFT JOIN %s %s ON b.dyad = %s.dyad",
        view_name,
        alias,
        alias
      )
    )
  }

  select_features <- if (length(feature_select)) {
    paste0(
      ",\n      ",
      paste(feature_select, collapse = ",\n      ")
    )
  } else {
    ""
  }

  join_features <- if (length(feature_joins)) {
    paste0(
      "\n    ",
      paste(feature_joins, collapse = "\n    ")
    )
  } else {
    ""
  }

  result_sql <- paste0(
    "
    SELECT
      b.dyad,
      b.protein,
      b.gene",
    select_features,
    "
    FROM protein_gene_dyad b",
    join_features
  )

  if (isTRUE(verbose)) {
    message("Building dyad annotation table.")
  }

  DBI::dbGetQuery(
    con,
    result_sql
  ) |>
    tibble::as_tibble()
}

#########################
#     HMMER helpers     #
#########################
#' Validate if a HMM file has old HMMER3 format
#'
#' @param hmm_file Path to a `.hmm` file.
#'
#' @returns `TRUE` if the file has valid HMMER3 formatting (starts with
#'   `HMMER3/f` and ends with `//`), `FALSE` otherwise.
#'
#' @keywords internal
.isValidHmmFile <- function(hmm_file) {

  lines <- tryCatch(
    readLines(hmm_file, warn = FALSE),
    error = function(e) character(0)
  )

  if (length(lines) == 0) {
    return(FALSE)
  }

  first_line <- trimws(lines[1])
  last_line  <- trimws(tail(lines, 1))

  starts_ok <- grepl("^HMMER3/f", first_line)
  ends_ok   <- identical(last_line, "//")

  starts_ok && ends_ok
}

#' Parsing HMM database to extract profile names, accessions and descriptions
#'
#' @param hmm_file path to the HMM database file (`.hmm`)
#'
#' @returns a tibble
#'
#' @keywords internal
.parse_hmmer_profiles <- function(hmm_file) {

  lines <- readLines(hmm_file, warn = FALSE)

  starts <- c(
    which(grepl("^NAME\\s+", lines)),
    length(lines) + 1L
  )

  blocks <- purrr::map2(
    starts[-length(starts)],
    starts[-1L] - 1L,
    ~ lines[.x:.y]
  )

  extract_field <- function(block, pattern) {

    hit <- stringr::str_subset(block, pattern)

    if (length(hit) == 0) {
      return(NA_character_)
    }

    stringr::str_remove(hit[[1]], pattern)
  }

  purrr::map_dfr(
    blocks,
    ~ tibble::tibble(
      profile_name = extract_field(.x, "^NAME\\s+"),
      profile_accession = extract_field(.x, "^ACC\\s+"),
      profile_description = extract_field(.x, "^DESC\\s+")
    )
  )
}

#' Parse HMMER tabular output into a tibble
#'
#' Reads a HMMER `--domtblout` file and returns a tidy tibble with one row per
#' target-query hit. Comment lines are stripped and the free-text description
#' field is reunited from the remaining whitespace-delimited columns.
#'
#' @param file Path to a HMMER `.tbl` output file produced with `--domtblout`.
#'
#' @return A tibble with 19 columns matching the HMMER per-sequence hit table
#'
#' @references Adapted from the rhmmer package
#'   (<https://github.com/arendsee/rhmmer>).
#'
#' @examples
#' \dontrun{
#' hits <- .parseHMMEROutput("results/Ecoli/protein_chunk_01_COG.tbl")
#' hits |> dplyr::filter(sequence_evalue < 1e-5)
#' }
#'
#' @keywords internal
.parseHMMEROutput <- function(file) {

  # target name         accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
  col_types <- readr::cols(
    protein           = readr::col_character(),  # target name
    protein_accession = readr::col_character(),
    tlen              = readr::col_integer(),

    query_name               = readr::col_character(),  # query name
    query_accession     = readr::col_character(),
    qlen              = readr::col_integer(),

    sequence_evalue   = readr::col_double(),
    sequence_score    = readr::col_double(),
    sequence_bias     = readr::col_double(),

    domain_num        = readr::col_integer(),
    domain_of         = readr::col_integer(),

    c_evalue          = readr::col_double(),
    i_evalue          = readr::col_double(),

    domain_score      = readr::col_double(),
    domain_bias       = readr::col_double(),

    hmm_from          = readr::col_integer(),
    hmm_to            = readr::col_integer(),

    ali_from          = readr::col_integer(),
    ali_to            = readr::col_integer(),

    env_from          = readr::col_integer(),
    env_to            = readr::col_integer(),

    acc               = readr::col_double(),

    target_description = readr::col_character()
  )
  # the line delimiter should always be just "\n", even on Windows
  lines <- readr::read_lines(file, lazy = FALSE, progress = FALSE)

  # drop comment lines
  data_lines <- lines[!grepl("^#", lines)]

  # split: whitespace-separated fields
  split_fields <- strsplit(data_lines, "\\s+", perl = TRUE)

  if (length(split_fields) == 0L) {
    return(readr::read_tsv(I(""), col_names = names(col_types$cols),
                            col_types = col_types, lazy = FALSE, progress = FALSE))
  }

  # count space separated fields
  N <- max(sapply(split_fields, length))

  # Parsing differently to avoid fussy read_tsv() warnings
  txt <- sub(
    pattern = sprintf("(%s).*", paste0(rep("\\S+", N), collapse = " +")),
    replacement = "\\1",
    x = lines,
    perl = TRUE
  ) |>
    gsub(pattern = "  *", replacement = "\t") |>
    paste0(collapse = "\n")

  table <- readr::read_tsv(
    I(txt),
    col_names = names(col_types$cols),
    comment = "#",
    na = "-",
    col_types = col_types,
    lazy = FALSE,
    progress = FALSE
  )

  table
}

#' Write a data frame to a compressed Parquet file
#'
#' @param df A data frame or tibble to write.
#' @param path Output file path (`.parquet` extension).
#'
#' @keywords internal
.write_compressed_parquet <- function(df, path) {
  arrow::write_parquet(
    df,
    path,
    compression = "zstd",
    compression_level = 9,
    use_dictionary = TRUE
  )
}

# Default persistent cache for shared HMMER databases
.defaultHmmerDbDir <- function() {
  file.path(
    tools::R_user_dir("amRdata", "cache"),
    "hmmer"
  )
}

.hmmer_version <- function(docker_image = "staphb/hmmer") {
  output <- system2(
    "docker",
    args = c(
      "run",
      "--rm",
      docker_image,
      "hmmsearch",
      "-h"
    ),
    stdout = TRUE,
    stderr = TRUE
  )

  version <- stringr::str_match(
    paste(output, collapse = "\n"),
    "HMMER ([0-9.]+)"
  )[, 2]

  if (is.na(version)) {
    return(NA_character_)
  }

  version
}
