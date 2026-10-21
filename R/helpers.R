### Helpers for amRdata live in this script


##########################
#    Verbosity helper    #
##########################

#' Print or log a pipeline status message
#'
#' Writes a status message to the console when `verbose = TRUE` and optionally
#' appends the same message to a log file. Logging is independent of console
#' verbosity, so status information can still be retained during quiet runs.
#'
#' @param ... Components passed to `paste0()` to construct the message.
#' @param verbose Logical. If TRUE, print the message to the console.
#'   Default: `FALSE`.
#' @param log_path Character or `NULL`. Optional path to a log file where the
#'   message should be appended.
#'
#' @return Invisibly returns the constructed message.
#' @keywords internal
.amr_status <- function(
    ...,
    verbose = FALSE,
    log_path = NULL
) {
  msg <- paste0(...)

  if (isTRUE(verbose)) {
    message(msg)
  }

  if (!is.null(log_path)) {
    .log_write(log_path, msg)
  }

  invisible(msg)
}

#' Check whether optional progress reporting is available
#'
#' Progress reporting is enabled only when requested by the caller and the
#' optional `progressr` package is installed.
#'
#' @param progress Logical. Whether progress reporting was requested.
#'   Default: `TRUE`.
#'
#' @return A logical scalar indicating whether progress reporting can be used.
#' @keywords internal
.amr_progress_available <- function(progress = TRUE) {
  isTRUE(progress) &&
    requireNamespace(
      "progressr",
      quietly = TRUE
    )
}

#' Select the progress handler used by amRdata
#'
#' Chooses the transient progress display used by amRdata. Status progress is
#' shown as a compact spinner and message, while step progress also reports the
#' current and total number of completed steps.
#'
#' Explicit user configuration through `progressr.handlers` takes precedence.
#' Otherwise, the `progress` backend is preferred when available, followed by
#' `cli`, with `progressr`'s text progress bar used as a fallback.
#'
#' @param type Character. Progress display type. `"status"` is used for
#'   long-running pipeline stages without a meaningful step count; `"steps"`
#'   reports completed steps as `current/total`.
#'
#' @return A `progressr` progression handler.
#' @keywords internal
.amr_progress_handler <- function(type = c("status", "steps")) {
  type <- match.arg(type)

  if (!is.null(getOption("progressr.handlers", NULL))) {
    return(progressr::handlers())
  }

  if (requireNamespace("progress", quietly = TRUE)) {
    if (identical(type, "steps")) {
      return(
        progressr::handler_progress(
          format = ":spin :message :current/:total",
          clear = TRUE
        )
      )
    }

    return(
      progressr::handler_progress(
        format = ":spin :message",
        clear = TRUE
      )
    )
  }

  if (requireNamespace("cli", quietly = TRUE)) {
    if (identical(type, "steps")) {
      return(
        progressr::handler_cli(
          format = "{cli::pb_spin} {cli::pb_status} {cli::pb_current}/{cli::pb_total}",
          format_done = "{cli::pb_status} {cli::pb_current}/{cli::pb_total}",
          clear = TRUE
        )
      )
    }

    return(
      progressr::handler_cli(
        format = "{cli::pb_spin} {cli::pb_status}",
        format_done = "{cli::pb_status}",
        clear = TRUE
      )
    )
  }

  progressr::handler_txtprogressbar(clear = TRUE)
}

#' Evaluate an expression with optional progress reporting
#'
#' @keywords internal
.amr_with_progress <- function(expr,
                               progress = TRUE,
                               type = c("status", "steps")) {
  type <- match.arg(type)
  expr <- substitute(expr)
  env <- parent.frame()

  if (!.amr_progress_available(progress)) {
    return(eval(expr, envir = env))
  }

  progressr::with_progress(eval(expr, envir = env),
                           handlers = .amr_progress_handler(type = type),
                           cleanup = TRUE)
}

#' How much time left on that big download? This will tell you
#'
#' Downloads a file with `httr2`, including retries and an extended timeout for
#' the large reference databases used by amRdata. When requested, `httr2`
#' reports download progress.
#'
#' @param url Character. URL of the file to download.
#' @param destfile Character. Local path where the downloaded file should be
#'   written.
#' @param progress Logical. Show download progress while transferring the file.
#'   Default: `TRUE`.
#'
#' @return Invisibly returns `destfile`.
#' @keywords internal
.amr_download_file <- function(url, destfile, progress = TRUE) {
  req <- httr2::request(url) |>
    httr2::req_retry(
      max_tries = 3L,
      retry_on_failure = TRUE
    ) |>
    httr2::req_timeout(
      max(3600, getOption("timeout"))
    )

  if (isTRUE(progress)) {
    old_options <- options(
      cli.progress_show_after = 0
    )

    on.exit(
      options(old_options),
      add = TRUE
    )

    req <- httr2::req_progress(
      req,
      type = "down"
    )
  }

  httr2::req_perform(
    req,
    path = destfile
  )

  invisible(destfile)
}


#' Create an optional progressor
#'
#' Creates a `progressr` progressor when progress reporting is available.
#' Otherwise, returns a no-op function so callers can issue progress updates
#' without their own package-availability checks.
#'
#' @param steps Integer. Number of progress steps expected.
#' @param progress Logical. Enable progress reporting when available.
#'   Default: `TRUE`.
#' @param label Character. Optional progressor label.
#' @param message Character. Initial progress message.
#' @param auto_finish Logical. If TRUE, automatically finish the progressor
#'   when all steps have been reported. Default: `TRUE`.
#'
#' @return A progressor function, or a no-op function when progress reporting
#'   is unavailable.
#' @keywords internal
.amr_progressor <- function(steps,
                            progress = TRUE,
                            label = NA_character_,
                            message = character(),
                            auto_finish = TRUE) {
  if (!.amr_progress_available(progress)) {
    return(function(...) {
      invisible(NULL)
    })
  }

  progressr::progressor(
    steps = steps,
    label = label,
    message = message,
    auto_finish = auto_finish,
    on_exit = FALSE
  )
}


#' Run one workflow stage with optional progress reporting
#'
#' Evaluates a workflow stage while displaying a tenmporary progress message.
#' Permanent console output is controlled independently by `verbose`, and the
#' stage message can optionally be written to a processing log.
#'
#' @param message Character scalar describing the workflow stage.
#' @param expr Expression to evaluate.
#' @param progress Logical. Show temporary progress while the expression runs.
#'   Default: `TRUE`.
#' @param verbose Logical. Print persistent start and completion messages.
#'   Default: `FALSE`.
#' @param log_path Character or `NULL`. Optional processing log path.
#'
#' @return The value returned by `expr`.
#' @keywords internal
.amr_progress_step <- function(message, expr,
                               progress = TRUE,
                               verbose = FALSE,
                               log_path = NULL) {
  if (missing(message) ||
      length(message) != 1L ||
      !nzchar(as.character(message))) {
    stop("`message` must be a non-empty character scalar.", call. = FALSE)
  }

  expr <- substitute(expr)
  env <- parent.frame()
  message <- as.character(message)

  .amr_status(
    paste0(message, "..."),
    verbose = verbose,
    log_path = log_path
  )

  value <- if (!.amr_progress_available(progress)) {
    eval(expr, envir = env)
  } else {
    .amr_with_progress({
      p <- .amr_progressor(
        steps = 1L,
        progress = progress,
        message = message
      )

      p(amount = 0, message = message)
      value <- eval(expr, envir = env)
      p(message = paste0(message, " complete"))
      value
    }, progress = progress)
  }

  .amr_status(
    paste0(message, " complete."),
    verbose = verbose,
    log_path = log_path
  )

  value
}

##########################
#  Package path helper   #
##########################

.amr_dataset_paths <- function(base_dir, user_bacs) {
  base_dir <- normalizePath(base_dir, mustWork = FALSE)

  dataset_name <- paste(user_bacs, collapse = "__") |>
    stringr::str_replace_all("\\s+", "_") |>
    stringr::str_replace_all("[^A-Za-z0-9._-]", "")

  dataset_id <- .generateDBname(user_bacs)
  root <- file.path(base_dir, "data", dataset_name)

  list(
    root = root,
    genomes = file.path(root, "genomes"),
    panaroo = file.path(root, "panaroo"),
    cdhit = file.path(root, "cd-hit"),
    hmmer = file.path(root, "hmmer"),

    work = file.path(root, "work"),
    working_duckdb = file.path(
      root,
      "work",
      paste0(dataset_id, ".duckdb")
    ),
    panaroo_input = file.path(
      root,
      "work",
      paste0(dataset_id, ".txt")
    ),

    orb = file.path(root, "orb"),
    parquet_duckdb = file.path(
      root,
      "orb",
      paste0(dataset_id, "_parquet.duckdb")
    ),
    processing_log = file.path(
      root,
      "orb",
      paste0(dataset_id, "_processing.log")
    ),

    exports = file.path(root, "exports")
  )
}

.amr_paths_from_duckdb <- function(duckdb_path) {
  duckdb_path <- normalizePath(
    duckdb_path,
    mustWork = FALSE
  )

  work_dir <- dirname(duckdb_path)

  if (!identical(basename(work_dir), "work")) {
    stop(
      "Working DuckDB must be located inside the dataset 'work/' directory."
    )
  }

  root <- dirname(work_dir)

  dataset_id <- tools::file_path_sans_ext(
    basename(duckdb_path)
  )

  list(
    root = root,

    genomes = file.path(root, "genomes"),
    panaroo = file.path(root, "panaroo"),
    cdhit = file.path(root, "cd-hit"),
    hmmer = file.path(root, "hmmer"),

    work = work_dir,
    working_duckdb = duckdb_path,
    panaroo_input = file.path(
      work_dir,
      paste0(dataset_id, ".txt")
    ),

    orb = file.path(root, "orb"),
    parquet_duckdb = file.path(
      root,
      "orb",
      paste0(dataset_id, "_parquet.duckdb")
    ),
    processing_log = file.path(
      root,
      "orb",
      paste0(dataset_id, "_processing.log")
    ),

    exports = file.path(root, "exports")
  )
}


.amr_paths_from_dataset_db <- function(dataset_path) {
  dataset_path <- normalizePath(
    dataset_path,
    mustWork = FALSE
  )

  container_dir <- dirname(dataset_path)
  container_name <- basename(container_dir)

  if (!container_name %in% c("work", "orb")) {
    stop(
      "Dataset DuckDB must be located inside the dataset 'work/' or 'orb/' directory."
    )
  }

  filename <- basename(dataset_path)

  if (
    identical(container_name, "orb") &&
    !grepl("_parquet\\.duckdb$", filename)
  ) {
    stop(
      "ORB DuckDB must use the '<dataset>_parquet.duckdb' filename."
    )
  }

  dataset_id <- tools::file_path_sans_ext(
    filename
  )

  dataset_id <- sub(
    "_parquet$",
    "",
    dataset_id
  )

  root <- dirname(container_dir)

  list(
    root = root,

    genomes = file.path(root, "genomes"),
    panaroo = file.path(root, "panaroo"),
    cdhit = file.path(root, "cd-hit"),
    hmmer = file.path(root, "hmmer"),

    work = file.path(root, "work"),
    working_duckdb = file.path(
      root,
      "work",
      paste0(dataset_id, ".duckdb")
    ),
    panaroo_input = file.path(
      root,
      "work",
      paste0(dataset_id, ".txt")
    ),

    orb = file.path(root, "orb"),
    parquet_duckdb = file.path(
      root,
      "orb",
      paste0(dataset_id, "_parquet.duckdb")
    ),
    processing_log = file.path(
      root,
      "orb",
      paste0(dataset_id, "_processing.log")
    ),

    exports = file.path(root, "exports")
  )
}

.amr_connect_dataset_db <- function(
    dataset_path,
    read_only = FALSE
) {
  dataset_path <- normalizePath(
    dataset_path,
    mustWork = isTRUE(read_only)
  )

  paths <- .amr_paths_from_dataset_db(
    dataset_path
  )

  con <- DBI::dbConnect(
    duckdb::duckdb(),
    dbdir = dataset_path,
    read_only = read_only
  )

  if (identical(
    basename(dirname(dataset_path)),
    "orb"
  )) {
    orb_dir <- normalizePath(
      paths$orb,
      mustWork = TRUE
    )

    tryCatch(
      {
        DBI::dbExecute(
          con,
          paste(
            "SET file_search_path =",
            DBI::dbQuoteString(
              con,
              orb_dir
            )
          )
        )
      },
      error = function(e) {
        try(
          DBI::dbDisconnect(con),
          silent = TRUE
        )

        stop(
          "Could not configure ORB Parquet search path: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }

  con
}


### Ye olde compatibility helper
#' Build the working DuckDB path for a user-bacs selection
#'
#' Compatibility wrapper around `.amr_dataset_paths()`.
#' Places the per-selection working database at:
#'   <base_dir>/data/<dataset>/work/<abbrev>.duckdb
#'
#' New code should prefer `.amr_dataset_paths()` when it needs paths other
#' than the working DuckDB!
#' @keywords internal
.buildDBpath <- function(base_dir = ".", user_bacs) {
  paths <- .amr_dataset_paths(
    base_dir = base_dir,
    user_bacs = user_bacs
  )

  dir.create(
    paths$work,
    recursive = TRUE,
    showWarnings = FALSE
  )

  list(
    db_dir = paths$root,
    db_path = paths$working_duckdb
  )
}

##########################
# CPU allocation helpers #
##########################

#' Resolve requested worker counts againstCPUs available
#' @keywords internal
.resolve_workers <- function(requested = NULL, n_tasks = NULL, warn = TRUE) {
  # No detectCores shenanigans
  available <- as.integer(parallelly::availableCores())

  if (is.null(requested)) {
    workers <- available
  } else {
    if (length(requested) != 1L ||
        is.na(requested) ||
        !is.numeric(requested) ||
        requested < 1 ||
        requested != floor(requested)) {
      stop("`requested` must be a positive number.")
    }

    requested <- as.integer(requested)
    workers <- min(requested, available)

    # Have you requested too many? The code politely figures it out for you
    if (isTRUE(warn) && requested > available) {
      warning(
        sprintf(
          "Requested %d parallel workers, but only %d CPU cores are available ",
          requested, available
        ),
        "to this R process; using ",
        available,
        " workers instead.",
        call. = FALSE
      )
    }
  }

  if (!is.null(n_tasks)) {
    workers <- min(workers, max(1L, as.integer(n_tasks)))
  }

  max(1L, as.integer(workers))
}

#' Set the active future plan for a resolved worker count
#'
#' Sequential when there's only one worker (avoids multisession overhead
#' for a single-threaded task), multisession otherwise. Callers are
#' responsible for restoring the previous plan (e.g. via
#' `old_plan <- future::plan(); on.exit(future::plan(old_plan), add = TRUE)`).
#' @keywords internal
.amr_set_future_plan <- function(n_workers) {
  if (n_workers == 1L) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession, workers = n_workers)
  }

  invisible(NULL)
}

#' Run independent BV-BRC API requests with explicit future plan
#' @keywords internal
.bvbrcFutureMap <- function(.x, .f, num_workers = 8L, ...) {
  if (!length(.x)) {
    return(list())
  }

  n_workers <- .resolve_workers(
    requested = num_workers,
    n_tasks = length(.x)
  )

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  .amr_set_future_plan(n_workers)

  furrr::future_map(
    .x,
    .f,
    ...,
    .options = furrr::furrr_options(seed = TRUE)
  )
}

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

#########################
# BiocFileCache helpers #
#########################

# BFC instance used by amRdata and other amR packages.
# The default BiocFileCache location is intentionally used so the registry is
# package-independent and can be shared across the amR suite
.amr_bfc <- function() {
  if (!requireNamespace("BiocFileCache", quietly = TRUE)) {
    stop("Package 'BiocFileCache' is required for amRdata resource caching.")
  }

  BiocFileCache::BiocFileCache(ask = FALSE)
}

# Find a BFC record by exact resource name.
.amr_bfc_find <- function(bfc, rname) {
  hits <- BiocFileCache::bfcquery(
    bfc,
    query = rname,
    field = "rname",
    exact = TRUE
  )

  if (!nrow(hits)) NULL else hits[1, , drop = FALSE]
}

# Return the BFC-managed path for the shared BV-BRC CLI metadata DuckDB.
#
# If the resource has already been registered but the DuckDB has not yet been
# created, create = TRUE returns the existing reserved path rather than replacing
# the BFC record
.amr_bfc_bvbrc_path <- function(create = FALSE, rname = .amr_bfc_bvbrc_rname()) {
  bfc <- .amr_bfc()
  hit <- .amr_bfc_find(bfc, rname)

  if (!is.null(hit)) {
    path <- as.character(hit$rpath[[1]])

    if (length(path) == 1L &&
        !is.na(path) &&
        nzchar(path)) {
      # BFC stores "relative" resources relative to its cache directory.
      if (identical(as.character(hit$rtype[[1]]), "relative")) {
        path <- file.path(BiocFileCache::bfccache(bfc), path)
      }

      path <- normalizePath(path, mustWork = FALSE)

      if (file.exists(path)) {
        return(normalizePath(path, mustWork = TRUE))
      }

      # The registry entry can exist before the DuckDB itself is written.
      if (isTRUE(create)) {
        return(path)
      }

      return(NULL)
    }

    if (isTRUE(create)) {
      BiocFileCache::bfcremove(bfc, hit$rid[[1]])
    } else {
      return(NULL)
    }
  }

  if (!isTRUE(create)) {
    return(NULL)
  }

  path <- BiocFileCache::bfcnew(
    bfc,
    rname = rname,
    rtype = "relative",
    ext = ".duckdb",
    fname = "exact"
  )

  normalizePath(path, mustWork = FALSE)
}

# Remove the shared BV-BRC CLI metadata DuckDB and its BFC registration
.amr_bfc_remove_bvbrc <- function(rname = .amr_bfc_bvbrc_rname()) {
  bfc <- .amr_bfc()

  hit <- .amr_bfc_find(bfc, rname)

  if (is.null(hit)) {
    return(invisible(list(
      registered = FALSE, path = NULL
    )))
  }

  path <- as.character(hit$rpath[[1]])

  if (length(path) == 1L &&
      !is.na(path) &&
      nzchar(path)) {
    if (identical(as.character(hit$rtype[[1]]), "relative")) {
      path <- file.path(BiocFileCache::bfccache(bfc), path)
    }

    path <- normalizePath(path, mustWork = FALSE)
  } else {
    path <- NULL
  }

  # Remove the resource from BiocFileCache
  BiocFileCache::bfcremove(bfc, hit$rid[[1]])

  if (!is.null(.amr_bfc_find(bfc, rname))) {
    stop("BV-BRC metadata could not be removed from the BiocFileCache registry.",
         call. = FALSE)
  }

  # Defensive cleanup in case anything remains on disk
  if (!is.null(path)) {
    leftovers <- c(path, paste0(path, ".wal"))

    leftovers <- leftovers[file.exists(leftovers)]

    if (length(leftovers)) {
      unlink(leftovers, force = TRUE)
    }
  }

  invisible(list(registered = TRUE, path = path))
}

# BFC name for a prepared HMMER database
.amr_bfc_hmmer_rname <- function(database, component = NULL) {
  parts <- c("amR_hmmer", database, component)
  parts <- parts[!is.na(parts) & nzchar(parts)]

  paste(parts, collapse = "_")
}

# Stable BFC resource name for the shared BV-BRC CLI metadata DuckDB
.amr_bfc_bvbrc_rname <- function() {
  "amRdata_bvbrc_bacterial_metadata"
}

# Register a prepared HMMER database with BFC
.amr_bfc_register_hmmer <- function(database,
                                    hmm_path,
                                    component = NULL) {
  hmm_path <- normalizePath(hmm_path, mustWork = TRUE)

  pressed <- paste0(
    hmm_path,
    c(".h3m", ".h3i", ".h3f", ".h3p")
  )

  missing <- pressed[!file.exists(pressed)]

  if (length(missing)) {
    stop(
      "Cannot register HMMER database before hmmpress is complete: ",
      paste(missing, collapse = ", ")
    )
  }

  rname <- .amr_bfc_hmmer_rname(
    database = database,
    component = component
  )

  rid <- .amr_bfc_register_local(
    path = hmm_path,
    rname = rname
  )

  invisible(list(
    rid = rid,
    rname = rname,
    hmm = hmm_path,
    pressed = pressed
  ))
}

# List HMMER databases known to the shared amR BFC
.amr_bfc_hmmer_resources <- function() {
  bfc <- .amr_bfc()

  BiocFileCache::bfcquery(
    bfc,
    query = "^amR_hmmer_",
    field = "rname",
    exact = FALSE
  )
}

# Register local files, used for dataset manifests: the manifest remains with data
# while BFC provides for cross-package discovery
.amr_bfc_register_local <- function(path, rname) {
  path <- normalizePath(path, mustWork = TRUE)
  bfc <- .amr_bfc()
  hit <- .amr_bfc_find(bfc, rname)

  if (is.null(hit)) {
    added <- BiocFileCache::bfcadd(
      bfc,
      rname = rname,
      fpath = path,
      rtype = "local",
      action = "asis",
      progress = FALSE
    )

    rid <- names(added)[[1]]

    return(invisible(rid))
  }

  existing_path <- tryCatch(
    normalizePath(hit$rpath[[1]], mustWork = TRUE),
    error = function(e) NA_character_
  )

  if (!identical(existing_path, path)) {
    BiocFileCache::bfcupdate(
      bfc,
      hit$rid[[1]],
      rpath = path,
      rname = rname
    )
  }

  invisible(hit$rid[[1]])
}

# Find datasets that successfully completed prepareGenomes() and still retain
# their working DuckDB for downstream feature processing.
.amr_processing_datasets <- function() {
  bfc <- .amr_bfc()

  resources <- BiocFileCache::bfcquery(
    bfc,
    query = "^amR_dataset_manifest_",
    field = "rname",
    exact = FALSE
  )

  empty_result <- tibble::tibble(
    label = character(),
    dataset_id = character(),
    duckdb_path = character(),
    manifest_path = character(),
    modified = as.POSIXct(character())
  )

  if (!nrow(resources)) {
    return(empty_result)
  }

  candidates <- purrr::map_dfr(
    seq_len(nrow(resources)),
    function(i) {
      manifest_path <- as.character(
        resources$rpath[[i]]
      )

      if (
        is.na(manifest_path) ||
        !nzchar(manifest_path) ||
        !file.exists(manifest_path)
      ) {
        return(NULL)
      }

      manifest <- tryCatch(
        jsonlite::read_json(
          manifest_path,
          simplifyVector = FALSE
        ),
        error = function(e) NULL
      )

      if (is.null(manifest)) {
        return(NULL)
      }

      valid_manifest <- tryCatch(
        {
          .manifest_validate(manifest)
          TRUE
        },
        error = function(e) FALSE
      )

      if (!isTRUE(valid_manifest)) {
        return(NULL)
      }

      prepare_complete <- any(
        purrr::map_lgl(
          manifest$runs %||% list(),
          function(run) {
            if (!identical(
              run$status,
              "success"
            )) {
              return(FALSE)
            }

            any(
              purrr::map_lgl(
                run$stages %||% list(),
                function(stage) {
                  identical(
                    stage$name,
                    "build_genome_file_table"
                  ) &&
                    identical(
                      stage$status,
                      "success"
                    )
                }
              )
            )
          }
        )
      )

      if (!prepare_complete) {
        return(NULL)
      }

      duckdb_path <- as.character(
        manifest$dataset$duckdb %||% ""
      )

      if (
        !nzchar(duckdb_path) ||
        !file.exists(duckdb_path)
      ) {
        return(NULL)
      }

      duckdb_path <- normalizePath(
        duckdb_path,
        mustWork = TRUE
      )

      valid_path <- tryCatch(
        {
          .amr_paths_from_duckdb(
            duckdb_path
          )
          TRUE
        },
        error = function(e) FALSE
      )

      if (!valid_path) {
        return(NULL)
      }

      user_bacs <- unlist(
        manifest$dataset$selection$user_bacs %||% character(),
        use.names = FALSE
      )

      label <- if (length(user_bacs)) {
        paste(
          user_bacs,
          collapse = ", "
        )
      } else {
        basename(
          dirname(
            dirname(duckdb_path)
          )
        )
      }

      tibble::tibble(
        label = label,
        dataset_id = as.character(
          manifest$dataset_id
        ),
        duckdb_path = duckdb_path,
        manifest_path = normalizePath(
          manifest_path,
          mustWork = TRUE
        ),
        modified = file.info(
          manifest_path
        )$mtime
      )
    }
  )

  if (!nrow(candidates)) {
    return(empty_result)
  }

  candidates |>
    dplyr::arrange(
      dplyr::desc(modified)
    ) |>
    dplyr::distinct(
      duckdb_path,
      .keep_all = TRUE
    )
}

# Find amRdata datasets that are still eligible for cleanup, including datasets
# whose mutable work/ directory has already been removed. Cleanup discovery is
# based on the retained registered manifest rather than the working DuckDB.
.amr_cleanup_datasets <- function() {
  bfc <- .amr_bfc()

  resources <- BiocFileCache::bfcquery(
    bfc,
    query = "^amR_dataset_manifest_",
    field = "rname",
    exact = FALSE
  )

  empty_result <- tibble::tibble(
    label = character(),
    dataset_id = character(),
    dataset_path = character(),
    manifest_path = character(),
    modified = as.POSIXct(character())
  )

  if (!nrow(resources)) {
    return(empty_result)
  }

  candidates <- purrr::map_dfr(
    seq_len(nrow(resources)),
    function(i) {
      manifest_path <- as.character(
        resources$rpath[[i]]
      )

      if (
        is.na(manifest_path) ||
        !nzchar(manifest_path) ||
        !file.exists(manifest_path)
      ) {
        return(NULL)
      }

      manifest <- tryCatch(
        jsonlite::read_json(
          manifest_path,
          simplifyVector = FALSE
        ),
        error = function(e) NULL
      )

      if (is.null(manifest)) {
        return(NULL)
      }

      valid_manifest <- tryCatch(
        {
          .manifest_validate(manifest)
          TRUE
        },
        error = function(e) FALSE
      )

      if (
        !isTRUE(valid_manifest) ||
        !identical(
          manifest$manifest_type %||% "",
          "amR_dataset"
        )
      ) {
        return(NULL)
      }

      manifest_path <- normalizePath(
        manifest_path,
        mustWork = TRUE
      )

      orb_dir <- dirname(manifest_path)
      dataset_path <- dirname(orb_dir)
      data_dir <- dirname(dataset_path)

      # Match the directory structure that removeLocalFiles() itself requires.
      if (
        !identical(basename(orb_dir), "orb") ||
        !identical(basename(data_dir), "data") ||
        !dir.exists(dataset_path)
      ) {
        return(NULL)
      }

      dataset_id <- as.character(
        manifest$dataset_id %||% ""
      )

      if (!nzchar(dataset_id)) {
        return(NULL)
      }

      user_bacs <- unlist(
        manifest$dataset$selection$user_bacs %||% character(),
        use.names = FALSE
      )

      label <- if (length(user_bacs)) {
        paste(
          user_bacs,
          collapse = ", "
        )
      } else {
        basename(dataset_path)
      }

      tibble::tibble(
        label = label,
        dataset_id = dataset_id,
        dataset_path = normalizePath(
          dataset_path,
          mustWork = TRUE
        ),
        manifest_path = manifest_path,
        modified = file.info(
          manifest_path
        )$mtime
      )
    }
  )

  if (!nrow(candidates)) {
    return(empty_result)
  }

  candidates |>
    dplyr::arrange(
      dplyr::desc(modified)
    ) |>
    dplyr::distinct(
      dataset_path,
      .keep_all = TRUE
    )
}

# Interactive selector function thing, very cool
.amr_select_processing_dataset <- function() {
  if (!interactive()) {
    stop(
      "`duckdb_path` must be supplied in non-interactive sessions.",
      call. = FALSE
    )
  }

  candidates <- .amr_processing_datasets()

  if (!nrow(candidates)) {
    stop(
      "No prepared amRdata datasets were found.\n",
      "Run prepareGenomes() first.",
      call. = FALSE
    )
  }

  if (nrow(candidates) == 1L) {
    message(
      "Using prepared dataset: ",
      candidates$label[[1]]
    )

    return(
      candidates$duckdb_path[[1]]
    )
  }

  choices <- paste0(
    candidates$label,
    " [",
    candidates$dataset_id,
    "]"
  )

  selection <- utils::menu(
    choices = choices,
    title = "Select an amRdata dataset to process:"
  )

  if (selection == 0L) {
    stop(
      "Dataset selection cancelled.",
      call. = FALSE
    )
  }

  candidates$duckdb_path[[selection]]
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
    parts <- strsplit(line, "[\\t ]", perl = TRUE)[[1]]
    if (length(parts) >= 9) {
      paste(
        c(
          parts[1:8],
          paste(parts[9:length(parts)], collapse = " ")
        ),
        collapse = "\t"
      )
    } else {
      line
    }
  })
  writeLines(out, gff_path, sep = "\n", useBytes = TRUE)
  invisible(TRUE)
}

#' Check BV-BRC data availability for a single taxon
#'
#' Internal worker used by `checkDataAvailability()`. Resolves one taxon to
#' available BV-BRC genomes and summarizes genome metadata, AMR phenotype
#' availability, collection years, genome statistics, and metadata QC.
#'
#' Metadata can be retrieved through the BV-BRC Data API or the legacy CLI/cache
#' workflow.
#'
#' @param user_bac Character scalar. Taxon ID or species name to summarize.
#' @param base_dir Character. Project root. Default: `"."`.
#' @param metadata_method Character. Metadata backend, either `"api"` or
#'   `"cli"`.
#' @param max_checkm_contam Numeric. Maximum allowed CheckM contamination
#'   percentage. Default: `5`.
#' @param min_checkm_complete Numeric. Minimum allowed CheckM completeness
#'   percentage. Default: `95`.
#' @param gc_deviations Numeric or `NULL`. Optional maximum standard deviations
#'   from median GC content.
#' @param length_deviations Numeric or `NULL`. Optional maximum standard
#'   deviations from median genome length.
#' @param cds_deviations Numeric or `NULL`. Optional maximum standard deviations
#'   from median CDS count.
#' @param verbose Logical. Print persistent metadata-query messages.
#'   Default: `TRUE`.
#' @param write_bac_data Logical. If TRUE, allow resolved metadata to be written
#'   to the per-selection `bac_data` table. Default: `FALSE`.
#'
#' @return A one-row tibble summarizing genome, AMR, QC, and collection metadata
#'   availability for the requested taxon.
#' @keywords internal
.checkDataPerTaxon <- function(
    user_bac,
    base_dir = ".",
    metadata_method = c("api", "cli"),
    max_checkm_contam = 5,
    min_checkm_complete = 95,
    gc_deviations = NULL,
    length_deviations = NULL,
    cds_deviations = NULL,
    verbose = TRUE,
    write_bac_data = FALSE
) {
  metadata_method <- match.arg(metadata_method)
  base_dir <- normalizePath(base_dir, mustWork = FALSE)

  # Little cache of internal helper helpers to help the helper
  empty_result <- function() {
    tibble::tibble(
      query = user_bac,
      total_genomes = 0L,
      wgs_genomes = 0L,
      complete_genomes = 0L,
      amr_genomes = 0L,
      amr_records = 0L,
      unique_antibiotics = 0L,
      antibiotics = NA_character_,
      drug_classes = NA_character_,
      checkm_available = 0L,
      qc_pass_genomes = 0L,
      qc_fail_genomes = 0L,
      median_genome_length = NA_real_,
      median_gc_content = NA_real_,
      median_cds = NA_real_,
      median_checkm_completeness = NA_real_,
      median_checkm_contamination = NA_real_,
      collection_year_min = NA_integer_,
      collection_year_max = NA_integer_
    )
  }

  safe_median <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]

    if (!length(x)) {
      return(NA_real_)
    }

    stats::median(x)
  }

  collapse_unique <- function(x) {
    x <- trimws(as.character(x))
    x <- sort(unique(x[!is.na(x) & nzchar(x)]))

    if (!length(x)) {
      return(NA_character_)
    }

    paste(x, collapse = ", ")
  }

  # Resolve this taxon to genome IDs
  genome_ids <- if (identical(metadata_method, "api")) {
    .resolveGenomeIDsApi(
      base_dir = base_dir,
      user_bacs = user_bac,
      verbose = verbose,
      write_bac_data = write_bac_data
    )
  } else {
    # Legacy CLI shenanigans
    bac_input_data <- .retrieveCustomQuery(
      base_dir = base_dir,
      user_bacs = user_bac,
      verbose = verbose
    )

    if (is.null(bac_input_data) || nrow(bac_input_data) == 0L) {
      character(0)
    } else {
      cache_db <- .amr_bfc_bvbrc_path(create = FALSE)

      if (is.null(cache_db) || !file.exists(cache_db)) {
        stop(
          "BV-BRC cache not found in BiocFileCache. ",
          "Run .updateBVBRCdata() first."
        )
      }

      con_cache <- DBI::dbConnect(
        duckdb::duckdb(),
        dbdir = cache_db,
        read_only = TRUE
      )
      on.exit(
        try(
          DBI::dbDisconnect(con_cache),
          silent = TRUE
        ),
        add = TRUE
      )

      taxon_ids <- unique(bac_input_data$genome.taxon_id)
      taxon_sql <- paste(
        DBI::dbQuoteString(con_cache, taxon_ids),
        collapse = ", "
      )

      query <- sprintf(
        paste0(
          "SELECT DISTINCT \"genome.genome_id\" AS genome_id ",
          "FROM bvbrc_bac_data ",
          "WHERE \"genome.taxon_id\" IN (%s)"
        ),
        taxon_sql
      )

      result <- DBI::dbGetQuery(con_cache, query)

      if (nrow(result)) {
        unique(as.character(result$genome_id))
      } else {
        character(0)
      }
    }
  }

  genome_ids <- unique(as.character(genome_ids))
  genome_ids <- genome_ids[
    !is.na(genome_ids) & nzchar(genome_ids)
  ]

  if (!length(genome_ids)) {
    if (isTRUE(verbose)) {
      message("No genomes matched '", user_bac, "'.")
    }
    return(empty_result())
  }

  # Fetch genome metadata
  genome_fields <- paste(
    c(
      "genome_id",
      "genome_name",
      "species",
      "taxon_id",
      "genome_quality",
      "genome_status",
      "collection_year",
      "genome_length",
      "gc_content",
      "cds",
      "checkm_completeness",
      "checkm_contamination"
    ),
    collapse = ","
  )

  genome_data <- if (identical(metadata_method, "api")) {
    .extractGenomeDataApi(
      genome_ids = genome_ids,
      fields = genome_fields,
      verbose = verbose
    )
  } else {
    raw <- .extractGenomeData(
      base_dir = base_dir,
      batch_genome_IDs = genome_ids,
      filter_type = "AMR",
      amr_fields = genome_fields,
      microtrait_fields = genome_fields,
      verbose = verbose
    )

    .parse_bvbrc_tsv(raw)
  }

  genome_data <- tibble::as_tibble(genome_data)

  if (!nrow(genome_data)) {
    if (isTRUE(verbose)) {
      message(
        "No genome metadata were returned for '",
        user_bac,
        "'."
      )
    }
    return(empty_result())
  }

  id_col <- dplyr::case_when(
    "genome.genome_id" %in% names(genome_data) ~
      "genome.genome_id",
    "genome_id" %in% names(genome_data) ~
      "genome_id",
    TRUE ~ NA_character_
  )

  if (is.na(id_col)) {
    stop("Genome metadata did not contain a proper genome ID column.")
  }

  genome_data <- genome_data |>
    dplyr::mutate(
      .genome_id = as.character(.data[[id_col]])
    ) |>
    dplyr::filter(
      !is.na(.data$.genome_id),
      nzchar(.data$.genome_id)
    ) |>
    dplyr::distinct(.data$.genome_id, .keep_all = TRUE)

  # Retrieve AMR phenotype records
  if (isTRUE(verbose)) {
    message("Checking AMR phenotype availability for '", user_bac, "'.")
  }

  amr_data <- if (identical(metadata_method, "api")) {
    .extractAMRtableApi(
      genome_ids = genome_ids,
      abx = "All",
      verbose = verbose
    )
  } else {
    drug_fields <- paste(
      c(
        "genome_id",
        "antibiotic",
        "evidence",
        "laboratory_typing_method",
        "resistant_phenotype"
      ),
      collapse = ","
    )

    raw <- .extractAMRtable(
      base_dir = base_dir,
      batch_genome_IDs = genome_ids,
      abx_filter = "--required antibiotic",
      drug_fields = drug_fields,
      verbose = verbose
    )

    .parse_bvbrc_tsv(raw)
  }

  amr_data <- tibble::as_tibble(amr_data)

  amr_id_col <- dplyr::case_when(
    "genome_drug.genome_id" %in% names(amr_data) ~
      "genome_drug.genome_id",
    "genome_id" %in% names(amr_data) ~
      "genome_id",
    TRUE ~ NA_character_
  )

  amr_ids <- if (!is.na(amr_id_col) && nrow(amr_data)) {
    unique(as.character(amr_data[[amr_id_col]]))
  } else {
    character(0)
  }

  amr_ids <- amr_ids[
    !is.na(amr_ids) & nzchar(amr_ids)
  ]

  antibiotic_col <- dplyr::case_when(
    "genome_drug.antibiotic" %in% names(amr_data) ~
      "genome_drug.antibiotic",
    "antibiotic" %in% names(amr_data) ~
      "antibiotic",
    TRUE ~ NA_character_
  )

  genome_length_col <- dplyr::case_when(
    "genome.genome_length" %in% names(genome_data) ~
      "genome.genome_length",
    "genome_length" %in% names(genome_data) ~
      "genome_length",
    TRUE ~ NA_character_
  )

  gc_content_col <- dplyr::case_when(
    "genome.gc_content" %in% names(genome_data) ~
      "genome.gc_content",
    "gc_content" %in% names(genome_data) ~
      "gc_content",
    TRUE ~ NA_character_
  )

  cds_col <- dplyr::case_when(
    "genome.cds" %in% names(genome_data) ~
      "genome.cds",
    "cds" %in% names(genome_data) ~
      "cds",
    TRUE ~ NA_character_
  )

  # Apply the same metadata QC used by retrieveMetadata()
  qc_out <- .apply_metadata_qc(
    genome_tbl = genome_data,
    max_checkm_contam = max_checkm_contam,
    min_checkm_complete = min_checkm_complete,
    gc_deviations = gc_deviations,
    length_deviations = length_deviations,
    cds_deviations = cds_deviations
  )

  qc_tbl <- tibble::as_tibble(qc_out$qc_tbl)

  total_genomes <- nrow(genome_data)

  wgs_genomes <- if (
    "genome.genome_status" %in% names(genome_data)
  ) {
    sum(
      genome_data$genome.genome_status == "WGS",
      na.rm = TRUE
    )
  } else if ("genome_status" %in% names(genome_data)) {
    sum(genome_data$genome_status == "WGS", na.rm = TRUE)
  } else {
    NA_integer_
  }

  complete_genomes <- if (
    "genome.genome_status" %in% names(genome_data)
  ) {
    sum(
      genome_data$genome.genome_status == "Complete",
      na.rm = TRUE
    )
  } else if ("genome_status" %in% names(genome_data)) {
    sum(
      genome_data$genome_status == "Complete",
      na.rm = TRUE
    )
  } else {
    NA_integer_
  }

  checkm_complete_col <- if (
    "genome.checkm_completeness" %in% names(genome_data)
  ) {
    "genome.checkm_completeness"
  } else if ("checkm_completeness" %in% names(genome_data)) {
    "checkm_completeness"
  } else {
    NA_character_
  }

  checkm_contam_col <- if (
    "genome.checkm_contamination" %in% names(genome_data)
  ) {
    "genome.checkm_contamination"
  } else if ("checkm_contamination" %in% names(genome_data)) {
    "checkm_contamination"
  } else {
    NA_character_
  }

  checkm_available <- if (
    !is.na(checkm_complete_col) &&
    !is.na(checkm_contam_col)
  ) {
    sum(
      !is.na(suppressWarnings(
        as.numeric(genome_data[[checkm_complete_col]])
      )) &
        !is.na(suppressWarnings(
          as.numeric(genome_data[[checkm_contam_col]])
        ))
    )
  } else {
    NA_integer_
  }

  qc_pass <- sum(qc_tbl$qc_keep %in% TRUE, na.rm = TRUE)
  qc_fail <- sum(qc_tbl$qc_keep %in% FALSE, na.rm = TRUE)

  collection_year_col <- if (
    "genome.collection_year" %in% names(genome_data)
  ) {
    "genome.collection_year"
  } else if ("collection_year" %in% names(genome_data)) {
    "collection_year"
  } else {
    NA_character_
  }

  collection_year <- if (!is.na(collection_year_col)) {
    suppressWarnings(
      as.integer(genome_data[[collection_year_col]])
    )
  } else {
    integer(0)
  }

  antibiotics <- if (!is.na(antibiotic_col)) {
    collapse_unique(amr_data[[antibiotic_col]])
  } else {
    NA_character_
  }

  unique_antibiotics <- if (!is.na(antibiotic_col)) {
    x <- trimws(as.character(amr_data[[antibiotic_col]]))
    length(unique(x[!is.na(x) & nzchar(x)]))
  } else {
    0L
  }

observed_drugs <- if (!is.na(antibiotic_col)) {
    x <- trimws(as.character(amr_data[[antibiotic_col]]))
    unique(x[!is.na(x) & nzchar(x)])
  } else {
    character(0)
  }

  drug_classes <- collapse_unique(
    drug_class$drug_class[
      drug_class$drug %in% observed_drugs
    ]
  )

  median_genome_length <- if (!is.na(genome_length_col)) {
    safe_median(genome_data[[genome_length_col]])
  } else {
    NA_real_
  }

  median_gc_content <- if (!is.na(gc_content_col)) {
    safe_median(genome_data[[gc_content_col]])
  } else {
    NA_real_
  }

  median_cds <- if (!is.na(cds_col)) {
    safe_median(genome_data[[cds_col]])
  } else {
    NA_real_
  }

  median_checkm_completeness <- if (!is.na(checkm_complete_col)) {
    safe_median(genome_data[[checkm_complete_col]])
  } else {
    NA_real_
  }

  median_checkm_contamination <- if (!is.na(checkm_contam_col)) {
    safe_median(genome_data[[checkm_contam_col]])
  } else {
    NA_real_
  }

  summary_row <- tibble::tibble(
    query = user_bac,

    total_genomes = as.integer(total_genomes),
    wgs_genomes = as.integer(wgs_genomes),
    complete_genomes = as.integer(complete_genomes),

    amr_genomes = as.integer(
      length(intersect(genome_data$.genome_id, amr_ids))
    ),
    amr_records = as.integer(nrow(amr_data)),
    unique_antibiotics = as.integer(unique_antibiotics),
    antibiotics = antibiotics,
    drug_classes = drug_classes,

    checkm_available = as.integer(checkm_available),
    qc_pass_genomes = as.integer(qc_pass),
    qc_fail_genomes = as.integer(qc_fail),

    median_genome_length = median_genome_length,
    median_gc_content = median_gc_content,
    median_cds = median_cds,
    median_checkm_completeness = median_checkm_completeness,
    median_checkm_contamination = median_checkm_contamination,
    collection_year_min = if (
      length(collection_year) &&
      any(!is.na(collection_year))
    ) {
      min(collection_year, na.rm = TRUE)
    } else {
      NA_integer_
    },
    collection_year_max = if (
      length(collection_year) &&
      any(!is.na(collection_year))
    ) {
      max(collection_year, na.rm = TRUE)
    } else {
      NA_integer_
    }
  )

  if (isTRUE(verbose)) {
    message(
      "Availability summary for '", user_bac, "': ",
      total_genomes, " genomes; ",
      wgs_genomes, " WGS; ",
      complete_genomes, " Complete; ",
      length(intersect(genome_data$.genome_id, amr_ids)),
      " with AMR records; ",
      qc_pass, " pass metadata QC."
    )
  }

  summary_row
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

  manifest_id <- tools::file_path_sans_ext(
    basename(manifest_path)
  )

  manifest <- list(
    schema_version = 1L,
    manifest_type = "amR_dataset",
    manifest_id = manifest_id,
    manifest_created_at = as.character(Sys.time()),
    manifest_updated_at = as.character(Sys.time()),
    dataset_id = dataset_id,
    dataset = list(
      duckdb = duckdb_path,
      selection = selection
    ),
    artifacts = list(),
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

# Registering manifest name and stage
.amr_bfc_manifest_rname <- function(manifest_state) {
  paste0(
    "amR_dataset_manifest_",
    manifest_state$manifest$dataset_id,
    "_",
    manifest_state$manifest$manifest_id
  )
}

# Better status recording for cross-suite hijinx
.manifest_artifact <- function(manifest_state,
                               name,
                               status = "ready",
                               details = list()) {
  if (!inherits(manifest_state, "amr_manifest")) {
    stop("Invalid manifest state.")
  }

  artifact <- c(
    list(
      status = status,
      updated_at = as.character(Sys.time())
    ),
    details
  )

  manifest_state$manifest$artifacts[[name]] <- artifact
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

# Manifest schema validation helper
.manifest_validate <- function(manifest) {
  if (!is.list(manifest)) {
    stop("Manifest must be a list.")
  }

  if (!identical(as.integer(manifest$schema_version), 1L)) {
    stop(
      "Unsupported amR manifest schema version: ",
      manifest$schema_version %||% "missing",
      ". Expected schema version 1."
    )
  }

  if (!identical(manifest$manifest_type, "amR_dataset")) {
    stop("Manifest is not an amR dataset manifest.")
  }

  # We need this stuff
  required <- c(
    "manifest_id",
    "dataset_id",
    "dataset",
    "artifacts",
    "runs"
  )

  missing <- setdiff(required, names(manifest))

  # If we don't have that stuff, hold your horses
  if (length(missing)) {
    stop(
      "Manifest is missing required field(s): ",
      paste(missing, collapse = ", ")
    )
  }

  invisible(TRUE)
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


#' Append a timestamped line to a plain-text progress log
#'
#' A lightweight, human-readable companion to the JSON provenance manifest.
#' The manifest records complete provenance but is rewritten wholesale on
#' every update, which makes it impractical to watch while a long-running
#' pipeline executes. This appends single lines instead, so the file can be
#' tailed (e.g. `tail -f`) to see what stage is currently running.
#'
#' @param log_path Character. Path to the log file. Created if it doesn't exist.
#' @param ... Character fragments pasted together to form the log message.
#'
#' @return Invisibly returns `log_path`.
#' @keywords internal
.log_write <- function(log_path, ...) {
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)

  cat(
    sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)),
    file = log_path,
    append = TRUE
  )

  invisible(log_path)
}

.log_tool_output <- function(log_path, tool, output) {
  if (is.null(log_path) ||
      !length(output)) {
    return(invisible(log_path))
  }

  .log_write(log_path, "----- ", tool, " output -----")

  purrr::walk(as.character(output), ~ .log_write(log_path, .x))

  .log_write(log_path, "----- end ", tool, " output -----")

  invisible(log_path)
}


.log_or_message <- function(log_path = NULL,
                            verbose = TRUE,
                            ...) {
  msg <- paste0(...)

  if (!is.null(log_path)) {
    .log_write(log_path, msg)
  } else if (isTRUE(verbose)) {
    message(msg)
  }

  invisible(msg)
}


#' Find the most recent recorded attempt of a named stage in a manifest run
#'
#' @param run A single run entry from a manifest's `runs` list, or `NULL`.
#' @param name Character. Stage name to look up.
#'
#' @return The matching stage entry (a list), or `NULL` if not found.
#' @keywords internal
.manifest_prior_stage <- function(run, name) {
  if (is.null(run) || !length(run$stages)) {
    return(NULL)
  }

  stage_names <- purrr::map_chr(run$stages, "name")
  idx <- which(stage_names == name)

  if (!length(idx)) {
    return(NULL)
  }

  run$stages[[idx[length(idx)]]]
}


#' Determine which pipeline stages can be safely skipped when resuming
#'
#' Looks at the most recent prior run recorded in a manifest and works out
#' how far into `stage_order` it got before it can be trusted. A stage only
#' counts as done if every earlier stage in `stage_order` also succeeded ---
#' later stages depend on earlier ones' DuckDB writes, so a gap partway
#' through can't be skipped around --- and its recorded output files are
#' still present on disk.
#'
#' @param prev_run A single run entry from a manifest's `runs` list, or `NULL`.
#' @param stage_order Character vector of stage names, in pipeline order.
#'
#' @return Named logical vector (named by `stage_order`) marking which
#'   stages are safe to skip.
#' @keywords internal
.resume_plan <- function(prev_run, stage_order) {
  completed <- stats::setNames(rep(FALSE, length(stage_order)), stage_order)

  for (name in stage_order) {
    stage <- .manifest_prior_stage(prev_run, name)

    if (is.null(stage) || !identical(stage$status, "success")) {
      break
    }

    out_paths <- purrr::map_chr(stage$outputs, "path")

    if (length(out_paths) && !all(file.exists(out_paths))) {
      break
    }

    completed[[name]] <- TRUE
  }

  completed
}

# To distinguish multiple manifests in the same bug directory
.manifest_find_latest <- function(
    dataset_path,
    require_success = TRUE
) {
  dataset_path <- normalizePath(
    dataset_path,
    mustWork = FALSE
  )

  container_dir <- dirname(dataset_path)

  if (!basename(container_dir) %in% c("work", "orb")) {
    stop(
      "Dataset files must be located inside 'work/' or 'orb/'."
    )
  }

  manifest_dir <- file.path(
    dirname(container_dir),
    "orb"
  )

  if (!dir.exists(manifest_dir)) {
    return(NULL)
  }

  manifests <- list.files(
    manifest_dir,
    pattern = "^manifest_.*\\.json$",
    full.names = TRUE
  )

  if (!length(manifests)) {
    return(NULL)
  }

  manifests <- manifests[
    order(
      file.info(manifests)$mtime,
      decreasing = TRUE
    )
  ]

  if (!isTRUE(require_success)) {
    return(manifests[[1]])
  }

  for (path in manifests) {
    manifest <- tryCatch(
      jsonlite::read_json(
        path,
        simplifyVector = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(manifest) || !length(manifest$runs)) {
      next
    }

    if (any(
      purrr::map_lgl(
        manifest$runs,
        ~ identical(.x$status, "success")
      )
    )) {
      return(path)
    }
  }

  NULL
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

  # Is this manifest any good?
  .manifest_validate(manifest)

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
#' features and HMMER annotations registered as virtual relations in the ORB.
#' Feature values are deduplicated and combined into semicolon-separated
#' character fields.
#'
#' @param duckdb_path Character. Path to the processed ORB DuckDB.
#' @param feature_scales Character vector of optional feature types to include.
#'   If NULL, includes `struct` plus all HMMER databases recorded in the
#'   manifest.
#' @param verbose Logical. Print progress messages.
#'
#' @return A tibble with one row per protein-gene dyad.
#'
#' @keywords internal
.exportDyadAnnotations <- function(
    duckdb_path,
    feature_scales = NULL,
    verbose = TRUE
) {
  duckdb_path <- normalizePath(duckdb_path, mustWork = TRUE)

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

  hmmer_stage <- NULL

  for (run in rev(manifest$runs %||% list())) {
    matches <- purrr::keep(
      run$stages %||% list(),
      ~ identical(.x$name, "hmmer") &&
        identical(.x$status, "success")
    )

    if (length(matches)) {
      hmmer_stage <- matches[[length(matches)]]
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

  con <- .amr_connect_dataset_db(
    duckdb_path,
    read_only = TRUE
  )

  on.exit(
    try(
      DBI::dbDisconnect(con),
      silent = TRUE
    ),
    add = TRUE
  )

  available_relations <- DBI::dbListTables(
    con
  )

  required_relations <- c(
    "dyads",
    "dyad_protein_cluster"
  )

  missing_required <- setdiff(
    required_relations,
    available_relations
  )

  if (length(missing_required)) {
    stop(
      "Dyad graph relation(s) not found in ORB: ",
      paste(missing_required, collapse = ", "),
      ". Re-run buildDyadFeatureMap()."
    )
  }

  feature_select <- character()
  feature_joins <- character()

  if (
    "struct" %in% feature_scales &&
    "struct_gene" %in% available_relations
  ) {
    feature_select <- c(
      feature_select,
      "sg.struct"
    )

    feature_joins <- c(
      feature_joins,
      paste0(
        "LEFT JOIN (",
        "SELECT gene, ",
        "string_agg(DISTINCT struct, ';' ORDER BY struct) AS struct ",
        "FROM struct_gene GROUP BY gene",
        ") sg ON d.gene = sg.gene"
      )
    )
  }

  specs <- .dyadHmmerSpecs(
    hmmer_databases
  )

  for (database in intersect(
    hmmer_databases,
    feature_scales
  )) {
    spec <- specs[[database]]

    if (!spec$annotation_view %in% available_relations) {
      if (isTRUE(verbose)) {
        message(
          "Skipping ",
          database,
          ": ",
          spec$annotation_view,
          " was not found in the ORB."
        )
      }
      next
    }

    alias <- paste0(
      "a_",
      spec$key
    )

    annotation_view_sql <- DBI::dbQuoteIdentifier(
      con,
      spec$annotation_view
    )

    output_column_sql <- DBI::dbQuoteIdentifier(
      con,
      spec$output_column
    )

    feature_select <- c(
      feature_select,
      paste0(
        alias,
        ".feature AS ",
        output_column_sql
      )
    )

    feature_joins <- c(
      feature_joins,
      paste0(
        "LEFT JOIN (",
        "SELECT cluster, ",
        "string_agg(DISTINCT feature, ';' ORDER BY feature) AS feature ",
        "FROM ",
        annotation_view_sql,
        " GROUP BY cluster",
        ") ",
        alias,
        " ON pc.cluster = ",
        alias,
        ".cluster"
      )
    )
  }

  select_features <- if (length(feature_select)) {
    paste0(
      ",\n      ",
      paste(
        feature_select,
        collapse = ",\n      "
      )
    )
  } else {
    ""
  }

  join_features <- if (length(feature_joins)) {
    paste0(
      "\n    ",
      paste(
        feature_joins,
        collapse = "\n    "
      )
    )
  } else {
    ""
  }

  result_sql <- paste0(
    "SELECT\n",
    "      d.protein || '|' || d.gene AS dyad,\n",
    "      d.protein,\n",
    "      d.gene",
    select_features,
    "\n    FROM dyads d\n",
    "    LEFT JOIN dyad_protein_cluster pc\n",
    "      ON d.protein = pc.protein",
    join_features
  )

  if (isTRUE(verbose)) {
    message(
      "Building dyad annotation table from virtual ORB relations."
    )
  }

  DBI::dbGetQuery(
    con,
    result_sql
  ) |>
    tibble::as_tibble()
}

.amr_resolve_export_orb <- function(path = NULL) {

  resolve_working_db <- function(duckdb_path) {
    paths <- .amr_paths_from_duckdb(
      duckdb_path
    )

    if (!file.exists(paths$parquet_duckdb)) {
      return(NA_character_)
    }

    normalizePath(
      paths$parquet_duckdb,
      mustWork = TRUE
    )
  }

    if (is.null(path)) {

    bfc <- .amr_bfc()

    resources <- BiocFileCache::bfcquery(
      bfc,
      query = "^amR_dataset_manifest_",
      field = "rname",
      exact = FALSE
    )

    if (!nrow(resources)) {
      stop(
        "No completed amRdata ORBs were found.",
        call. = FALSE
      )
    }

    candidates <- purrr::map_dfr(
      seq_len(nrow(resources)),
      function(i) {

        manifest_path <- as.character(
          resources$rpath[[i]]
        )

        if (
          is.na(manifest_path) ||
          !nzchar(manifest_path) ||
          !file.exists(manifest_path)
        ) {
          return(NULL)
        }

        manifest <- tryCatch(
          jsonlite::read_json(
            manifest_path,
            simplifyVector = FALSE
          ),
          error = function(e) NULL
        )

        if (is.null(manifest)) {
          return(NULL)
        }

        valid_manifest <- tryCatch(
          {
            .manifest_validate(manifest)
            TRUE
          },
          error = function(e) FALSE
        )

        if (!isTRUE(valid_manifest)) {
          return(NULL)
        }

        artifact <- manifest$artifacts$amRml_input %||% NULL

        if (
          is.null(artifact) ||
          !identical(
            artifact$status,
            "ready"
          )
        ) {
          return(NULL)
        }

        orb_dir <- dirname(
          normalizePath(
            manifest_path,
            mustWork = TRUE
          )
        )

        orb_databases <- list.files(
          orb_dir,
          pattern = "_parquet\\.duckdb$",
          full.names = TRUE
        )

        if (length(orb_databases) != 1L) {
          return(NULL)
        }

        orb_path <- normalizePath(
          orb_databases[[1]],
          mustWork = TRUE
        )

        user_bacs <- unlist(
          manifest$dataset$selection$user_bacs %||% character(),
          use.names = FALSE
        )

        label <- if (length(user_bacs)) {
          paste(
            user_bacs,
            collapse = ", "
          )
        } else {
          basename(
            dirname(
              orb_dir
            )
          )
        }

        tibble::tibble(
          label = label,
          dataset_id = as.character(
            manifest$dataset_id
          ),
          orb_path = orb_path,
          manifest_path = normalizePath(
            manifest_path,
            mustWork = TRUE
          ),
          modified = file.info(
            manifest_path
          )$mtime
        )
      }
    )

    if (!nrow(candidates)) {
      stop(
        "No completed amRdata ORBs were found.",
        call. = FALSE
      )
    }

    candidates <- candidates |>
      dplyr::arrange(
        dplyr::desc(modified)
      ) |>
      dplyr::distinct(
        orb_path,
        .keep_all = TRUE
      )

    if (nrow(candidates) == 1L) {
      return(
        candidates$orb_path[[1]]
      )
    }

    if (!interactive()) {
      stop(
        "`duckdb_path` must be supplied when multiple ORBs are available ",
        "in a non-interactive session.",
        call. = FALSE
      )
    }

    choices <- paste0(
      candidates$label,
      " [",
      candidates$dataset_id,
      "]"
    )

    selection <- utils::menu(
      choices = choices,
      title = "Select an amRdata dataset to export:"
    )

    if (selection == 0L) {
      stop(
        "Dataset selection cancelled.",
        call. = FALSE
      )
    }

    return(
      candidates$orb_path[[selection]]
    )
  }

  path <- normalizePath(
    path,
    mustWork = TRUE
  )

  if (dir.exists(path)) {

    orb_databases <- list.files(
      path,
      pattern = "_parquet\\.duckdb$",
      full.names = TRUE
    )

    if (!length(orb_databases)) {
      stop(
        "No ORB DuckDB was found in: ",
        path,
        call. = FALSE
      )
    }

    if (length(orb_databases) > 1L) {
      stop(
        "Multiple ORB DuckDB files were found in: ",
        path,
        call. = FALSE
      )
    }

    return(
      normalizePath(
        orb_databases[[1]],
        mustWork = TRUE
      )
    )
  }

  paths <- .amr_paths_from_dataset_db(
    path
  )

  container <- basename(
    dirname(path)
  )

  if (identical(container, "orb")) {
    return(path)
  }

  if (
    identical(container, "work") &&
    file.exists(paths$parquet_duckdb)
  ) {
    return(
      normalizePath(
        paths$parquet_duckdb,
        mustWork = TRUE
      )
    )
  }

  stop(
    "No completed ORB could be resolved from: ",
    path,
    call. = FALSE
  )
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
#' Streams through an HMM database in bounded line chunks and retains only the
#' `NAME`, `ACC`, and `DESC` fields needed during annotation finalization.
#'
#' @param hmm_file Path to the HMM database file (`.hmm`).
#' @param chunk_lines Number of lines to read per chunk.
#'
#' @returns A tibble with one row per HMM profile.
#'
#' @keywords internal
.parse_hmmer_profiles <- function(hmm_file, chunk_lines = 100000L) {
  input <- file(hmm_file, open = "r")
  on.exit(close(input), add = TRUE)

  names_out <- character()
  accessions_out <- character()
  descriptions_out <- character()

  current_name <- NA_character_
  current_accession <- NA_character_
  current_description <- NA_character_

  flush_profile <- function() {
    if (is.na(current_name)) {
      return(invisible(NULL))
    }

    names_out <<- c(names_out, current_name)
    accessions_out <<- c(accessions_out, current_accession)
    descriptions_out <<- c(descriptions_out, current_description)

    invisible(NULL)
  }

  repeat {
    lines <- readLines(input, n = chunk_lines, warn = FALSE)
    if (!length(lines)) break

    relevant <- lines[grepl("^(NAME|ACC|DESC)\\s+", lines)]
    if (!length(relevant)) next

    for (line in relevant) {
      if (grepl("^NAME\\s+", line)) {
        flush_profile()
        current_name <- sub("^NAME\\s+", "", line)
        current_accession <- NA_character_
        current_description <- NA_character_
      } else if (grepl("^ACC\\s+", line) &&
                 !is.na(current_name) &&
                 is.na(current_accession)) {
        current_accession <- sub("^ACC\\s+", "", line)
      } else if (grepl("^DESC\\s+", line) &&
                 !is.na(current_name) &&
                 is.na(current_description)) {
        current_description <- sub("^DESC\\s+", "", line)
      }
    }
  }

  flush_profile()

  tibble::tibble(
    profile_name = names_out,
    profile_accession = accessions_out,
    profile_description = descriptions_out
  )
}

#' Create or reuse cached HMM profile metadata
#'
#' Stores the small profile lookup table beside the prepared HMM database so
#' dataset finalization does not need to rescan the full HMM file. The sidecar
#' is rebuilt when it is missing or older than the HMM database and is registered
#' as a component of the database in the shared BiocFileCache registry.
#'
#' @param hmm_file Path to the prepared HMM database.
#' @param database Database name used for BiocFileCache registration.
#' @param verbose Logical. Print a message when the sidecar is rebuilt.
#' @param component Character or `NULL`. Optional component name used to keep
#'   profile caches distinct for databases with multiple prepared HMMs.
#'
#' @return A list containing the sidecar path and BiocFileCache identifiers.
#' @keywords internal
.hmmer_profile_cache <- function(hmm_file, database, verbose = FALSE,
                                 component = NULL) {
  hmm_file <- normalizePath(hmm_file, mustWork = TRUE)

  profile_path <- file.path(
    dirname(hmm_file),
    paste0(
      tools::file_path_sans_ext(basename(hmm_file)),
      ".profiles.parquet"
    )
  )

  rebuild <- !file.exists(profile_path)

  if (!rebuild) {
    rebuild <- isTRUE(
      file.info(profile_path)$mtime < file.info(hmm_file)$mtime
    )
  }

  if (rebuild) {
    if (isTRUE(verbose)) {
      cache_parts <- c(database, component)
      cache_parts <- cache_parts[!is.na(cache_parts) & nzchar(cache_parts)]
      message(
        "Caching HMM profile metadata for ",
        paste(cache_parts, collapse = "/")
      )
    }

    profiles <- .parse_hmmer_profiles(hmm_file)

    arrow::write_parquet(
      profiles,
      profile_path,
      compression = "zstd"
    )
  }

  profile_path <- normalizePath(profile_path, mustWork = TRUE)

  rname_parts <- c(component, "profiles")
  rname_parts <- rname_parts[!is.na(rname_parts) & nzchar(rname_parts)]

  rname <- .amr_bfc_hmmer_rname(
    database,
    component = paste(rname_parts, collapse = "_")
  )
  rid <- .amr_bfc_register_local(profile_path, rname)

  invisible(list(
    path = profile_path,
    rid = rid,
    rname = rname
  ))
}

#' Parse HMMER tabular output into a tibble
#'
#' Reads a HMMER `--domtblout` file and returns a tidy tibble with one row per
#' target-query hit. Comment lines are stripped and the free-text description
#' field is reunited from the remaining whitespace-delimited columns.
#'
#' @param file Path to a HMMER `.tbl` output file produced with `--domtblout`.
#'
#' @return A tibble with 23 columns matching the HMMER per-sequence hit table
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
    BiocFileCache::bfccache(.amr_bfc()),
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
