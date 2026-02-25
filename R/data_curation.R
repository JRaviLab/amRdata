#' Helps ensure trailing 0s are retained in genome IDs for proper downloading
#' @keywords internal
.id_checker <- function(x) {
  # Taxon IDs are just numbers, genome IDs have decimals, this tells them apart
  grepl("^[0-9]+$", x)
}

#' Helps tag genomes with their AMR evidence for parsing
#' @keywords internal
.create_amr_tagged_view <- function(con) {
  lab_methods  <- c("Disk diffusion", "MIC", "Broth dilution", "Agar dilution")
  lab_list_sql <- paste(DBI::dbQuoteString(con, lab_methods), collapse = ", ")
  comp_str_sql <- DBI::dbQuoteString(con, "Computational Method")

  DBI::dbExecute(
    con,
    sprintf(
      "CREATE OR REPLACE VIEW amr_phenotype_tagged AS
       SELECT
         *,
         CASE WHEN \"genome_drug.laboratory_typing_method\" IN (%s) THEN 1 ELSE 0 END AS is_lab_row,
         CASE WHEN
           (\"genome_drug.evidence\" = %s)
           OR (COALESCE(\"genome_drug.computational_method\", '') <> '')
           OR (\"genome_drug.laboratory_typing_method\" = 'Computational Prediction')
         THEN 1 ELSE 0 END AS is_comp_row
       FROM amr_phenotype",
      lab_list_sql, comp_str_sql
    )
  )
}

#' Helps appropriately interface with BV-BRC FTPS server, and avoids getting stuck
#' when malformed files can hang an FTPS connection by introducing safeguards
#' @keywords internal
.ftpes_download_one <- function(genomeID, out_dir,
                                connect_timeout = 10L,
                                max_time        = 30L,
                                speed_time      = 30L,    # end if <speed_limit for >speed_time
                                speed_limit     = 2048L,  # B/s
                                min_bytes       = 100L) {

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  exts <- c(".fna", ".PATRIC.faa", ".PATRIC.gff")
  for (ext in exts) {
    dest     <- file.path(out_dir, paste0(genomeID, ext))
    dest_tmp <- paste0(dest, ".tmp")

    # Skip if we already completed these
    if (file.exists(dest) && file.info(dest)$size > min_bytes) next
    if (file.exists(dest) && file.info(dest)$size == 0) try(unlink(dest), silent = TRUE)
    if (file.exists(dest_tmp)) try(unlink(dest_tmp), silent = TRUE)

    url <- sprintf("ftp://ftp.bv-brc.org/genomes/%s/%s%s", genomeID, genomeID, ext)

    args <- c(
      "--fail", "--silent", "--show-error", "-L",
      "--connect-timeout", as.character(connect_timeout),
      "--max-time",        as.character(max_time),
      "--speed-time",      as.character(speed_time),
      "--speed-limit",     as.character(speed_limit),
      "--ftp-ssl",                 # AUTH TLS on port 21 works in testing
      "--ftp-pasv", "--disable-epsv", "--ipv4",
      "--user", "anonymous:",
      "-o", shQuote(dest_tmp), shQuote(url)
    )

    res <- suppressWarnings(system2("curl", args = args, stdout = TRUE, stderr = TRUE))
    status <- attr(res, "status"); if (is.null(status)) status <- 0L

    # Avoiding 0B files cluttering up results -- atomic rename on successful tmp DL
    ok <- (status == 0L && file.exists(dest_tmp) && file.info(dest_tmp)$size > min_bytes)
    if (ok) {
      if (!file.rename(dest_tmp, dest)) {
        # If rename fails for some reason, copy + unlink
        file.copy(dest_tmp, dest, overwrite = TRUE)
        unlink(dest_tmp)
      }
    } else {
      try(unlink(dest_tmp), silent = TRUE)
      return(FALSE)  # Abort early, don't accept partial set
    }
  }

  # Do we have all 3 present?
  .is_complete_set(out_dir, genomeID, min_bytes = min_bytes)
}

#' Helps manage FTPS downloading from BV-BRC, tryng a quick download first, and
#' if that fails, trying a longer timeout 2nd pass at the end in case it was a
#' hiccup. If 2nd pass fails, log and give up on that file.
#' @keywords internal
.ftpes_download_two_pass <- function(genome_ids, out_dir,
                                     workers_first  = 4L,
                                     workers_second = 4L,
                                     log_file       = NULL) {
  genome_ids <- unique(as.character(genome_ids))
  if (!length(genome_ids)) return(character(0))

  if (!is.null(log_file)) {
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("[%s] FTPS run start: %d genomes\n", Sys.time(), length(genome_ids)),
        file = log_file, append = TRUE)
  }

  # Pass 1: 30s per-file cap
  message("FTPS pass 1 (30s timeout)")
  future::plan(future::multisession, workers = max(1, workers_first))
  res1 <- future.apply::future_lapply(
    genome_ids,
    function(gid) {
      ok <- .ftpes_download_one(gid, out_dir,
                                connect_timeout = 10L, max_time = 30L,
                                speed_time = 30L, speed_limit = 2048L)
      list(gid = gid, ok = ok)
    },
    future.seed = TRUE
  )
  ok1 <- vapply(res1, `[[`, logical(1), "ok")
  ok_ids_1 <- genome_ids[ok1]
  fail_ids <- genome_ids[!ok1]
  message(sprintf("Pass 1: ok=%d, fail=%d", length(ok_ids_1), length(fail_ids)))
  if (!is.null(log_file)) {
    cat(sprintf("[%s] Pass1 ok=%d fail=%d\n", Sys.time(), length(ok_ids_1), length(fail_ids)),
        file = log_file, append = TRUE)
  }

  if (!length(fail_ids)) {
    if (!is.null(log_file)) cat(sprintf("[%s] FTPS run end: all OK\n", Sys.time()),
                                file = log_file, append = TRUE)
    return(ok_ids_1)
  }

  # Pass 2: 60s per-file cap where we retry any failures
  message("FTPS pass 2 (60s timeout) for failed genomes")
  future::plan(future::multisession, workers = max(1, workers_second))
  res2 <- future.apply::future_lapply(
    fail_ids,
    function(gid) {
      ok <- .ftpes_download_one(gid, out_dir,
                                connect_timeout = 10L, max_time = 60L,
                                speed_time = 30L, speed_limit = 2048L)
      list(gid = gid, ok = ok)
    },
    future.seed = TRUE
  )
  ok2 <- vapply(res2, `[[`, logical(1), "ok")
  ok_ids_2   <- fail_ids[ok2]
  still_fail <- setdiff(fail_ids, ok_ids_2)
  message(sprintf("Pass 2: ok=%d, still_fail=%d", length(ok_ids_2), length(still_fail)))
  if (!is.null(log_file)) {
    cat(sprintf("[%s] Pass2 ok=%d still_fail=%d\n", Sys.time(), length(ok_ids_2), length(still_fail)),
        file = log_file, append = TRUE)
    if (length(still_fail)) {
      cat("Fail IDs (excluded): ", paste(head(still_fail, 50), collapse = ", "), "\n",
          file = log_file, append = TRUE)
    }
  }

  unique(c(ok_ids_1, ok_ids_2))
}

#' Fetch BV-BRC bacterial genome data
#'
#' Run the BV-BRC CLI (`p3-all-genomes`) inside a Docker image and return results as a tibble.
#'
#' @param image Character scalar. Docker image with BV-BRC CLI preinstalled.
#'   Defaults to `"danylmb/bvbrc:5.3"`.
#' @param verbose Logical. If TRUE, prints informative progress messages. Default: TRUE.
#'
#' @return A tibble with bacterial genome metadata (BV-BRC column names intact).
#' @keywords internal
.fetchBVBRCdata <- function(image = "danylmb/bvbrc:5.3", verbose = TRUE) {
  # Check Docker availability
  docker_path <- Sys.which("docker")
  if (!nzchar(docker_path)) {
    stop("Docker is not available on your PATH, but is required to be installed.
         Please install or load Docker before running this package.")
  }

  if (isTRUE(verbose)) {
    message("Please wait. Fetching BV-BRC bacterial genome metadata.")
  }

  # Construct the BV-BRC CLI command
  cmd <- glue::glue(
    "docker run --rm {image} p3-all-genomes ",
    "--in superkingdom,Bacteria ",
    "--eq genome_quality,Good ",
    "--in genome_status,WGS,Complete ",
    "--attr genome_id,genome_name,genome_quality,genome_status,",
    "taxon_id,species,strain,_version_,",
    "collection_year,state_province,latitude,longitude,",
    "antimicrobial_resistance_evidence,assembly_accession,isolation_country,",
    "isolation_source,disease,host_common_name"
  )

  # Suppress CLI noise
  raw_data <- tryCatch(
    system(cmd, intern = TRUE, ignore.stderr = TRUE),
    error = function(e) {
      stop("BV-BRC CLI call failed.\n", conditionMessage(e))
    }
  )
  # If BV-BRC server is down
  if (length(raw_data) == 0L) {
    stop("BV-BRC returned no data. The service may be unavailable.")
  }

  # Parse and coerce to character to avoid numeric interpretation losing trailing 0s
  df <- utils::read.table(
    text = raw_data, sep = "\t", header = TRUE, fill = TRUE,
    quote = "", check.names = FALSE, comment.char = "", colClasses = "character"
  )
  df <- tibble::as_tibble(df) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ iconv(.x, from = "", to = "UTF-8", sub = "")))

  if (isTRUE(verbose)) {
    message(glue::glue("Retrieved {nrow(df)} rows x {ncol(df)} columns."))
  }

  return(df)
}


# Make sure the BV-BRC metadata live where they're supposed to
.ensure_bvbrc_cache <- function(base_dir = ".",
                                verbose = TRUE,
                                cache_rel = file.path("data", "bvbrc", "bvbrcData.duckdb"),
                                cache_table = "bvbrc_bac_data") {
  base_dir <- normalizePath(base_dir, mustWork = FALSE)
  cache_db <- file.path(base_dir, cache_rel)

  need_build <- !file.exists(cache_db)
  con_cache <- NULL

  if (!need_build) {
    con_cache <- DBI::dbConnect(duckdb::duckdb(), dbdir = cache_db)
    on.exit(try(DBI::dbDisconnect(con_cache, shutdown = TRUE), silent = TRUE), add = TRUE)
    need_build <- !(cache_table %in% DBI::dbListTables(con_cache))
  }

  if (need_build) {
    if (isTRUE(verbose)) message("BV-BRC cache missing or incomplete. Building via .updateBVBRCdata(). Please wait.")
    .updateBVBRCdata(base_dir = base_dir, verbose = verbose)

    if (!is.null(con_cache)) try(DBI::dbDisconnect(con_cache, shutdown = TRUE), silent = TRUE)
    if (!file.exists(cache_db)) stop("After .updateBVBRCdata(), cache DB still missing at: ", cache_db)

    con_cache <- DBI::dbConnect(duckdb::duckdb(), dbdir = cache_db)
    on.exit(try(DBI::dbDisconnect(con_cache, shutdown = TRUE), silent = TRUE), add = TRUE)
    if (!(cache_table %in% DBI::dbListTables(con_cache))) {
      stop("After .updateBVBRCdata(), table '", cache_table, "' still not found in ", cache_db)
    }
  }

  invisible(cache_db)
}

#' Update BV-BRC metadata in DuckDB
#'
#' Fetches bacterial genome metadata from BV-BRC using the BV-BRC CLI and stores
#' it in a DuckDB database under `data/bvbrc/bvbrcData.duckdb` within `base_dir`.
#' If the table exists and is older than `max_age_days`, it refreshes; otherwise,
#' loads the existing table. BV-BRC column names are preserved exactly.
#'
#' @param base_dir Character. Project root. The DuckDB database is created at
#'   `file.path(base_dir, "data", "bvbrc", "bvbrcData.duckdb")`.
#' @param max_age_days Integer. Refresh the table if older than this many days. Default: 30.
#' @param image Character. Docker image used by `.fetchBVBRCdata()`. Default: `"danylmb/bvbrc:5.3"`.
#' @param verbose Logical. If TRUE, prints informative messages. Default: TRUE.
#'
#' @return A tibble containing BV-BRC bacterial genome metadata.
#' @export
.updateBVBRCdata <- function(base_dir = ".",
                             max_age_days = 30L,
                             image = "danylmb/bvbrc:5.3",
                             verbose = TRUE) {
  # base_dir as project root
  base_dir <- normalizePath(base_dir, mustWork = FALSE)
  data_dir <- file.path(base_dir, "data")
  bvbrc_dir <- file.path(data_dir, "bvbrc")
  logs_dir <- file.path(data_dir, "logs")

  dir.create(bvbrc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

  db_path <- file.path(bvbrc_dir, "bvbrcData.duckdb")
  table_name <- "bvbrc_bac_data"
  meta_table <- "__meta"

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path)

  DBI::dbExecute(
    con,
    glue::glue("CREATE TABLE IF NOT EXISTS {meta_table} (
                 table_name TEXT PRIMARY KEY,
                 last_updated TIMESTAMP
               )")
  )

  # Tiny update check helper
  .last_updated <- function(con, table_name) {
    res <- tryCatch(
      DBI::dbGetQuery(
        con,
        glue::glue("SELECT last_updated FROM {meta_table}
                    WHERE table_name = {DBI::dbQuoteString(con, table_name)}")
      ),
      error = function(e) NULL
    )
    if (is.null(res) || nrow(res) == 0L) {
      return(NA)
    }
    as.POSIXct(res$last_updated[[1]], origin = "1970-01-01", tz = "UTC")
  }

  tables <- DBI::dbListTables(con)
  needs_refresh <- TRUE

  if (table_name %in% tables) {
    if (isTRUE(verbose)) message("BV-BRC table exists. Checking time since last update.")
    last_up <- .last_updated(con, table_name)

    if (!is.na(last_up)) {
      age_days <- as.numeric(difftime(Sys.time(), last_up, units = "days"))
      needs_refresh <- isTRUE(age_days > max_age_days)
      if (isTRUE(verbose)) message(paste0(round(age_days, 1), " days since last update."))
    } else {
      file_mtime <- tryCatch(file.info(db_path)$mtime, error = function(e) NA)
      if (!is.na(file_mtime)) {
        age_days <- as.numeric(difftime(Sys.time(), file_mtime, units = "days"))
        needs_refresh <- isTRUE(age_days > max_age_days)
        if (isTRUE(verbose)) message(paste0(round(age_days, 1), " days since last update (file mtime)."))
      } else {
        needs_refresh <- TRUE
      }
    }
  } else {
    if (isTRUE(verbose)) message("BV-BRC metadata table does not exist. Creating.")
    needs_refresh <- TRUE
  }

  if (needs_refresh) {
    if (isTRUE(verbose)) message("Fetching fresh BV-BRC metadata.")
    bvbrc_bacs <- .fetchBVBRCdata(image = image, verbose = verbose)

    if (isTRUE(verbose)) message("Writing table to DuckDB.")
    DBI::dbWriteTable(con, table_name, bvbrc_bacs, overwrite = TRUE)

    DBI::dbExecute(
      con,
      glue::glue("INSERT OR REPLACE INTO {meta_table} (table_name, last_updated)
                  VALUES ({DBI::dbQuoteString(con, table_name)}, CURRENT_TIMESTAMP)")
    )

    if (isTRUE(verbose)) {
      message(glue::glue("Saved {nrow(bvbrc_bacs)} rows x {ncol(bvbrc_bacs)} columns to {db_path} (table '{table_name}')."))
    }
  } else {
    if (isTRUE(verbose)) message("Table is up-to-date. Loading from DuckDB.")
    bvbrc_bacs <- tibble::as_tibble(DBI::dbReadTable(con, table_name))
    if (isTRUE(verbose)) {
      message(glue::glue("Loaded {nrow(bvbrc_bacs)} rows x {ncol(bvbrc_bacs)} columns from existing table."))
    }
  }

  DBI::dbDisconnect(con, shutdown = TRUE)
  invisible(bvbrc_bacs)
}

#' Retrieve BV-BRC records for user-provided bacteria
#'
#' Searches the locally cached BV-BRC bacterial dataset for user-specified inputs.
#' Numeric inputs are treated as taxon IDs; character inputs are matched as
#' case-insensitive substrings against `genome.species`.
#'
#' @param base_dir Character. Project root. BV-BRC cache is expected at
#'   `file.path(base_dir, "data", "bvbrc", "bvbrcData.duckdb")`.
#' @param user_bacs Character vector. Mixed inputs of taxon IDs and/or species strings.
#'
#' @return A tibble with columns `genome.taxon_id` and `genome.species`, or NULL with a message.
.retrieveCustomQuery <- function(base_dir = ".",
                                 user_bacs = c("90371", "Bacillus subtilis")) {
  base_dir <- normalizePath(base_dir, mustWork = FALSE)

  # Ensure global BV-BRC metadata exists/updated at <base_dir>/data/bvbrc/...
  bvbrc_bacs <- .updateBVBRCdata(base_dir = base_dir)

  bac_input_data <- tibble::tibble(
    genome.taxon_id = character(),
    genome.species  = character()
  )

  for (user_bac in user_bacs) {
    if (.id_checker(user_bac)) {
      message("Numeric input detected: ", user_bac)
      if (user_bac %in% bvbrc_bacs$genome.taxon_id) {
        bac_name <- bvbrc_bacs$genome.species[bvbrc_bacs$genome.taxon_id == user_bac]
        bac_df <- tibble::tibble(
          genome.taxon_id = user_bac,
          genome.species  = unique(bac_name)
        )
        bac_input_data <- dplyr::bind_rows(bac_input_data, bac_df)
      } else {
        message("No match in the database for taxon ID: ", user_bac)
      }
    } else {
      message("String input detected: ", user_bac)
      matched <- stringr::str_detect(
        bvbrc_bacs$genome.species,
        stringr::fixed(user_bac, ignore_case = TRUE)
      )
      if (any(matched)) {
        matched_indices <- which(matched)
        bac_df <- tibble::tibble(
          genome.taxon_id = bvbrc_bacs$genome.taxon_id[matched_indices],
          genome.species  = bvbrc_bacs$genome.species[matched_indices]
        ) |> dplyr::distinct()
        bac_input_data <- dplyr::bind_rows(bac_input_data, bac_df)
      } else {
        message("No match in the database for species substring: ", user_bac)
      }
    }
  }

  bac_input_data <- bac_input_data[!is.na(bac_input_data$genome.species), ]

  if (nrow(bac_input_data) > 0) {
    return(bac_input_data)
  } else {
    message("No matches in the database found for the provided inputs.")
    return(NULL)
  }
}


#' Resolve `query value` from `user_bacs` for [.getGenomeIDs()]
#'
#' If query_value is NULL, derive it from user_bacs based on query_type.
#' For species/genome_name: take the first element of user_bacs.
#' For taxon_id: take the first numeric-looking element of user_bacs.
#' If nothing suitable is found, throw a fit and an error.
#' @keywords internal
.resolveQueryValue <- function(query_type, query_value, user_bacs) {
  if (!is.null(query_value) && nzchar(query_value)) {
    return(query_value)
  }
  if (missing(user_bacs) || length(user_bacs) == 0) {
    stop("Provide query_value or user_bacs for the selected query_type.")
  }
  if (query_type %in% c("species", "genome_name")) {
    cand <- user_bacs[1]
    if (is.na(cand) || !nzchar(cand) || !is.character(cand)) {
      stop("Cannot infer query_value for type '", query_type, "'.")
    }
    return(cand)
  }
  if (query_type == "taxon_id") {
    nums <- .id_checker(user_bacs)
    if (!any(nums)) stop("Cannot infer taxon_id from user_bacs. Provide query_value.")
    return(as.character(user_bacs[which(nums)[1]]))
  }

  stop("Unsupported query_type: ", query_type)
}

#' Generate a shortened database name from taxon IDs or species names
#'
#' This function creates an identifier by combining abbreviated species names
#' or taxon IDs. For species names, it uses the first letter of the genus and
#' the first two letters of the species. For single-word names, it appends "sp".
#' For numeric taxon IDs, it prefixes them with "tid_".
#'
#' @param user_bacs Character vector containing mixed inputs of taxon IDs
#' and/or species names.
#'
#' Defaults to `c("90371", "Bacillus subtilis")`.
#'
#' @return A single character string representing the combined shortened name.
#'
#' @examples
#' .generateDBname(c("90371", "Bacillus subtilis"))
#' .generateDBname(c("12345", "Escherichia coli", "Lactobacillus"))
#' @keywords internal
.generateDBname <- function(user_bacs) {
  db_parts <- c()

  for (user_bac in user_bacs) {
    if (.id_checker(user_bac)) {
      db_parts <- c(db_parts, paste0("tid_", user_bac))
    } else {
      parts <- stringr::str_split(user_bac, " ")[[1]]
      if (length(parts) == 1) {
        # If only one word, use "sp" as the second part
        abbreviation <- paste0(stringr::str_sub(parts[1], 1, 1), "sp")
      } else {
        abbreviation <- paste0(stringr::str_sub(parts[1], 1, 1), stringr::str_sub(parts[2], 1, 2))
      }
      db_parts <- c(db_parts, abbreviation)
    }
  }

  db_name <- paste(db_parts, collapse = "_")
  return(db_name)
}

#' Build a DuckDB path for a user-bacs selection
#'
#' Places the per-selection DB at:
#'   <base_dir>/data/<bug_dir>/<abbrev>.duckdb
#' where <bug_dir> is derived from full user_bacs input and <abbrev> from
#' .generateDBname(user_bacs). This function no longer enforces overwrite checks.
#'
#' @param base_dir Character. Project root.
#' @param user_bacs Character vector. The same vector used for DB naming.
#' @param overwrite Logical. Ignored (kept for backward compatibility).
#'
#' @return A list with `db_dir` and `db_path`.
#' @keywords internal
.buildDBpath <- function(base_dir, user_bacs, overwrite = FALSE) {
  base_dir <- normalizePath(base_dir, mustWork = FALSE)
  data_dir <- file.path(base_dir, "data")

  # Directory from full names (order-sensitive by design)
  full_joined <- paste(user_bacs, collapse = "__")
  bug_dirname <- full_joined |>
    stringr::str_replace_all("\\s+", "_") |>
    stringr::str_replace_all("[^A-Za-z0-9._-]", "")

  db_dir <- file.path(data_dir, bug_dirname)
  dir.create(db_dir, recursive = TRUE, showWarnings = FALSE)

  db_file <- paste0(.generateDBname(user_bacs), ".duckdb")
  db_path <- file.path(db_dir, db_file)

  list(db_dir = db_dir, db_path = db_path)
}


#' Retrieve genome IDs from BV-BRC and store them in DuckDB
#'
#' Executes BV-BRC CLI queries in Docker and writes the results to a DuckDB at:
#'   <base_dir>/data/<bug_dir>/<abbrev>.duckdb
#' where <bug_dir> derives from full `user_bacs`, and <abbrev> from `.generateDBname()`.
#' BV-BRC column names are preserved.
#'
#' @param base_dir Character. Project root directory.
#' @param query_type Character. One of "genome_name", "species", or "taxon_id".
#' @param query_value Character or NULL. If NULL, it will be inferred from `user_bacs`
#'   based on `query_type` (first element for species/genome_name, first numeric for taxon_id).
#' @param user_bacs Character vector. Used to construct database location/name and, if needed,
#'   to infer `query_value` when not supplied.
#' @param overwrite Logical. If FALSE and the DuckDB file already exists, abort. Default: FALSE.
#' @param image Character. Docker image containing BV-BRC CLI. Default: "danylmb/bvbrc:5.3".
#' @param verbose Logical. If TRUE, prints messages. Default: TRUE.
#'
#' @return A list with:
#'   - count_result: Integer (count query result)
#'   - duckdbConnection: DBI connection to the DuckDB file
#'   - table_name: "bac_data"
.getGenomeIDs <- function(base_dir = ".",
                          query_type = c("genome_name", "species", "taxon_id"),
                          query_value = NULL,
                          user_bacs,
                          overwrite = FALSE,
                          image = "danylmb/bvbrc:5.3",
                          verbose = TRUE) {
  query_type <- match.arg(query_type)
  query_value <- .resolveQueryValue(query_type, query_value, user_bacs)

  if (isTRUE(verbose)) {
    message("Querying BV-BRC: ", query_type, " == ", query_value)
  }

  # Count
  count_cmd <- paste0(
    "docker run --rm ", image,
    " p3-all-genomes --in ",
    query_type, ",\"", query_value, "\"",
    " --eq genome_quality,Good",
    " --in genome_status,WGS,Complete",
    " --count"
  )
  count_lines <- tryCatch(system(count_cmd, intern = TRUE), error = function(e) character())
  count_result <- suppressWarnings(as.integer(if (length(count_lines) >= 2) count_lines[2] else NA_integer_))
  if (isTRUE(verbose) && !is.na(count_result)) message("Count returned: ", count_result)

  # Details
  data_cmd <- paste0(
    "docker run --rm ", image,
    " p3-all-genomes --in ",
    query_type, ",\"", query_value, "\"",
    " --eq genome_quality,Good",
    " --in genome_status,WGS,Complete",
    " --attr genome_id,genome_name,taxon_id,species,strain"
  )
  data_raw <- tryCatch(system(data_cmd, intern = TRUE), error = function(e) character())
  if (length(data_raw) == 0L) stop("BV-BRC returned no data for: ", query_type, " = ", query_value)

  data_result <- tibble::as_tibble(
    utils::read.table(
      text = data_raw, sep = "\t", header = TRUE, fill = TRUE,
      quote = "", check.names = FALSE, comment.char = "", colClasses = "character"
    )
  ) |>
    dplyr::mutate(
      `genome.genome_id`   = as.character(`genome.genome_id`),
      `genome.genome_name` = as.character(`genome.genome_name`),
      `genome.taxon_id`    = as.character(`genome.taxon_id`),
      `genome.species`     = as.character(`genome.species`),
      `genome.strain`      = as.character(`genome.strain`)
    )

  # Per-bug DB path
  paths <- .buildDBpath(base_dir = base_dir, user_bacs = user_bacs, overwrite = overwrite)
  db_path <- paths$db_path

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path)
  DBI::dbWriteTable(con, "bac_data", data_result, overwrite = TRUE)

  if (isTRUE(verbose)) message("Wrote table 'bac_data' to: ", db_path)

  list(count_result = count_result, duckdbConnection = con, table_name = "bac_data")
}

#' Retrieve genome IDs for each taxon via BV-BRC and DuckDB
#'
#' Resolves user-provided taxa to taxon IDs using the local BV-BRC cache, then
#' returns distinct genome IDs (Good + WGS/Complete). Also writes a 'bac_data'
#' table to the per-selection DuckDB with a cache snapshot of core fields.
#'
#' @return A character vector of distinct `genome.genome_id`, or NULL if none found.
.retrieveQueryIDs <- function(base_dir = ".",
                              user_bacs,
                              overwrite = FALSE,
                              verbose = TRUE) {
  base_dir <- normalizePath(base_dir, mustWork = FALSE)

  if (isTRUE(verbose)) message("Resolving input taxa.")
  bac_input_data <- .retrieveCustomQuery(base_dir = base_dir, user_bacs = user_bacs)
  if (is.null(bac_input_data) || nrow(bac_input_data) == 0) {
    message("No valid input provided or no matches found.")
    return(NULL)
  }

  # Per-selection DB path
  paths   <- .buildDBpath(base_dir = base_dir, user_bacs = user_bacs, overwrite = overwrite)
  db_path <- paths$db_path

  # Got a cache? Use that, it's fast
  cache_db <- file.path(base_dir, "data", "bvbrc", "bvbrcData.duckdb")
  if (!file.exists(cache_db)) {
    stop("BV-BRC cache not found at: ", cache_db, ". Run .updateBVBRCdata() first.")
  }
  con_cache <- DBI::dbConnect(duckdb::duckdb(), dbdir = cache_db, read_only = TRUE)
  on.exit(try(DBI::dbDisconnect(con_cache, shutdown = TRUE), silent = TRUE), add = TRUE)

  taxon_ids <- unique(bac_input_data$genome.taxon_id)
  if (isTRUE(verbose)) message("Querying cache for ", length(taxon_ids), " taxon IDs.")

  taxon_list_sql <- paste(DBI::dbQuoteString(con_cache, taxon_ids), collapse = ", ")
  sql_ids <- sprintf(
    "SELECT DISTINCT
       \"genome.genome_id\"   AS gid,
       \"genome.genome_name\" AS genome_name,
       \"genome.taxon_id\"    AS taxon_id,
       \"genome.species\"     AS species,
       \"genome.strain\"      AS strain
     FROM bvbrc_bac_data
     WHERE \"genome.taxon_id\" IN (%s)
       AND \"genome.genome_quality\" = 'Good'
       AND \"genome.genome_status\"  IN ('WGS','Complete')",
    taxon_list_sql
  )
  cache_rows <- DBI::dbGetQuery(con_cache, sql_ids)
  valid_id_re <- "^[0-9]+\\.[0-9]+$"
  cache_rows <- cache_rows[grepl(valid_id_re, cache_rows$gid), , drop = FALSE]

  if (nrow(cache_rows) == 0L) {
    message("No valid genome IDs found.")
    return(NULL)
  }

  genome_ids <- unique(cache_rows$gid)
  if (isTRUE(verbose)) message("Collected ", length(genome_ids), " distinct genome IDs (cache).")

  # Write 'bac_data' from cache
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)

  bac_data_tbl <- data.frame(
    `genome.genome_id`   = cache_rows$gid,
    `genome.genome_name` = cache_rows$genome_name,
    `genome.taxon_id`    = cache_rows$taxon_id,
    `genome.species`     = cache_rows$species,
    `genome.strain`      = cache_rows$strain,
    stringsAsFactors = FALSE
  )
  DBI::dbWriteTable(con, "bac_data", bac_data_tbl, overwrite = TRUE)
  if (isTRUE(verbose)) message("Wrote table 'bac_data' to: ", db_path)

  genome_ids
}

#' Extract AMR Data Table
#'
#' This function generates and retrieves AMR data for a batch of genome IDs.
#' It uses BV-BRC tools within a Docker container to fetch and process the data based on the input parameters.
#'
#' @param batch_genome_IDs A character vector containing genome IDs to query.
#' @param abx_filter A string specifying the antibiotic filter criteria. This can be in the form
#'   `"--required antibiotic"` or `"--in antibiotic,drug1,drug2"`.
#' @param drug_fields A string specifying the drug fields (attributes) to include in the output.
#'   This corresponds to the attributes retrieved by the `p3-get-genome-drugs` tool.
#' @param path A string representing the file path to a directory where temporary files and data will be stored.
#' @param image Character. Docker image. Default "danylmb/bvbrc:5.3".
#' @param verbose Logical. If TRUE, prints concise messages.
#'
#' @details
#' The function performs the following steps:
#' 1. Ensures the required Docker image (`danylmb/bvbrc:5.3`) is pulled.
#' 2. Creates a genome ID table using `p3-echo` and writes it to a temporary file.
#' 3. Processes the generated file with the `p3-get-genome-drugs` tool to query the AMR data.
#'
#' @return A character vector with the output data from `p3-get-genome-drugs`, typically in tabular format.
#'
#' @note
#' The BV-BRC Docker image (`danylmb/bvbrc:5.3`) must be installed and available on the system.
#'
#' @examples
#' \dontrun{
#' batch_ids <- c("genome1", "genome2", "genome3")
#' abx_filter <- "--required antibiotic"
#' drug_attrs <- "antibiotic,resistance"
#' temp_path <- "/tmp"
#'
#' result <- extractAMRtable(batch_ids, abx_filter, drug_attrs, temp_path)
#' print(result)
#' }
#'
#' @export
.extractAMRtable <- function(base_dir,
                             batch_genome_IDs,
                             abx_filter,
                             drug_fields,
                             image = "danylmb/bvbrc:5.3",
                             verbose = TRUE) {
  base_dir <- normalizePath(base_dir, mustWork = FALSE)
  data_dir <- file.path(base_dir, "data")
  tmp_dir <- file.path(data_dir, "tmp")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(verbose)) {
    message("Preparing AMR query input for ", length(batch_genome_IDs), " genomes.")
  }

  docker_path <- Sys.which("docker")
  if (!nzchar(docker_path)) {
    stop("Docker is not available on your PATH but is required.")
  }

  # Generate genome list with p3-echo
  echo_args <- c(
    "run", "--rm",
    image, "p3-echo",
    "--title=genome_drug.genome_id", paste(batch_genome_IDs, collapse = " ")
  )
  genome_ids_output <- suppressWarnings(
    system2("docker", args = echo_args, stdout = TRUE, stderr = TRUE)
  )

  # Write a temporary file in data/tmp/
  tmp_in <- tempfile(tmpdir = tmp_dir, pattern = "genome_drug_ids_", fileext = ".tsv")
  writeLines(genome_ids_output, con = tmp_in)
  tmp_in_mounted <- file.path("/data", "tmp", basename(tmp_in))

  # Allow abx_filter to be a single string with spaces OR a vector of args
  abx_args <-
    if (length(abx_filter) == 1L) {
      strsplit(abx_filter, "[[:space:]]+", perl = TRUE)[[1]]
    } else {
      abx_filter
    }

  # Query drug data
  drug_args <- c(
    "run", "--rm",
    "-v", paste0(data_dir, ":/data"),
    image, "p3-get-genome-drugs",
    "--input", tmp_in_mounted,
    abx_args,
    "--attr", drug_fields
  )

  if (isTRUE(verbose)) {
    message("Running AMR query.")
  }
  drug_data <- suppressWarnings(system2("docker", args = drug_args, stdout = TRUE, stderr = TRUE))

  # Clean up after yourself
  try(unlink(tmp_in), silent = TRUE)

  return(drug_data)
}

#' This function retrieves metadata for the given genome IDs.
#'
#' @param base_dir Character. Project root.
#' @param batch_genome_IDs Vector of genome IDs.
#' @param filter_type Character. Either "AMR" or "microTraits" (sets data fields).
#' @param amr_fields Character. Attributes (comma-separated) for AMR metadata.
#' @param microtrait_fields Character. Attributes for microtraits.
#' @param image Character.
#' @param verbose Logical.
#'
#' @return A character vector containing the retrieved genome data.
#'
.extractGenomeData <- function(base_dir,
                               batch_genome_IDs,
                               filter_type,
                               amr_fields,
                               microtrait_fields,
                               image = "danylmb/bvbrc:5.3",
                               verbose = TRUE) {
  base_dir <- normalizePath(base_dir, mustWork = FALSE)
  data_dir <- file.path(base_dir, "data")
  tmp_dir <- file.path(data_dir, "tmp")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(verbose)) {
    message("Preparing genome metadata input for ", length(batch_genome_IDs), " genomes.")
  }

  docker_path <- Sys.which("docker")
  if (!nzchar(docker_path)) {
    stop("Docker is not available on your PATH but is required.")
  }

  # Generate genome list with p3-echo
  echo_args <- c(
    "run", "--rm",
    image, "p3-echo",
    "--title=genome.genome_id", paste(batch_genome_IDs, collapse = " ")
  )
  genome_ids_output <- suppressWarnings(
    system2("docker", args = echo_args, stdout = TRUE, stderr = TRUE)
  )

  tmp_in <- tempfile(tmpdir = tmp_dir, pattern = "genome_ids_", fileext = ".tsv")
  writeLines(genome_ids_output, con = tmp_in)
  tmp_in_mounted <- file.path("/data", "tmp", basename(tmp_in))

  # Choose attributes (AMR for this workflow)
  chosen_fields <- if (identical(filter_type, "AMR")) amr_fields else microtrait_fields

  get_args <- c(
    "run", "--rm",
    "-v", paste0(data_dir, ":/data"),
    image, "p3-get-genome-data",
    "--input", tmp_in_mounted,
    "--attr", chosen_fields
  )

  if (isTRUE(verbose)) {
    message("Running genome metadata query.")
  }
  genome_data <- suppressWarnings(system2("docker", args = get_args, stdout = TRUE, stderr = TRUE))

  # Cleaning up
  try(unlink(tmp_in), silent = TRUE)

  return(genome_data)
}

#' Retrieve AMR or Microtrait metadata from BV-BRC and store in DuckDB
#'
#' Queries BV-BRC for AMR or microtrait metadata for genomes corresponding to user inputs.
#' Results are written to a per-selection DuckDB at:
#'   <base_dir>/data/<bug_dir>/<abbrev>.duckdb
#' Tables written:
#'   - amr_phenotype
#'   - genome_data
#'   - metadata (join on genome IDs returned by BV-BRC)
#'
#' @param user_bacs Character vector. Mixed taxon IDs and/or species strings (used for naming).
#' @param filter_type Character. "AMR" or "microTraits". Default "AMR".
#' @param base_dir Character. Project root. Default "results/" in legacy scripts; now default ".".
#' @param abx Character or vector. Antibiotic filter. "All" for all antibiotics, else names.
#' @param overwrite Logical. If FALSE and DuckDB exists already, abort. Default FALSE.
#' @param image Character. Docker image. Default "danylmb/bvbrc:5.3".
#' @param verbose Logical. If TRUE, prints concise messages.
#'
#' @return A list with:
#'   - duckdbConnection: live DBI connection to the created DuckDB
#'   - table_name: "metadata"
retrieveMetadata <- function(user_bacs,
                             filter_type = "AMR",
                             base_dir = ".",
                             abx = "All",
                             overwrite = FALSE,
                             image = "danylmb/bvbrc:5.3",
                             verbose = TRUE) {
  base_dir <- normalizePath(base_dir, mustWork = FALSE)

  if (isTRUE(verbose)) message("Resolving genome IDs for user inputs.")
  genome_ids <- .retrieveQueryIDs(
    base_dir = base_dir, user_bacs = user_bacs,
    overwrite = overwrite, verbose = verbose
  )
  if (length(genome_ids) == 0) {
    message("No genome IDs available for the specified inputs.")
    return(NULL)
  }

  # Desired fields from BV-BRC
  drug_fields <- paste0(
    "antibiotic,computational_method,",
    "evidence,genome_name,id,",
    "laboratory_typing_method,",
    "laboratory_typing_platform,",
    "measurement,measurement_sign,",
    "measurement_unit,measurement_value,",
    "pmid,resistant_phenotype,",
    "source,taxon_id,testing_standard"
  )
  abx_filter <- if (identical(abx, "All")) {
    "--required antibiotic"
  } else {
    paste0("--in antibiotic,", paste(abx, collapse = ","))
  }
  amr_fields <- paste0(
    "assembly_accession,assembly_method,",
    "bioproject_accession,biosample_accession,",
    "body_sample_site,",
    "body_sample_subsite,collection_date,",
    "collection_year,disease,",
    "genbank_accessions,genome_name,",
    "genome_quality,genome_status,",
    "geographic_location,geographic_group,host_age,",
    "host_common_name,host_gender,host_group,",
    "host_health,isolation_country,",
    "isolation_site,isolation_source,",
    "ncbi_project_id,phenotype,publication,",
    "refseq_accessions,",
    "refseq_project_id,sra_accession,species,taxon_id"
  )
  microtrait_fields <- paste0(
    amr_fields, ",genome_length,gram_stain,",
    "habitat,",
    "host_name,",
    "host_scientific_name,",
    "isolation_comments,",
    "lab_host,",
    "latitude,",
    "optimal_temperature,",
    "other_environmental,",
    "oxygen_requirement,",
    "rrna,",
    "salinity,",
    "sporulation,",
    "temperature_range,",
    "trna"
  )

  # Batching downloads and parallel implementation
  batch_size <- 500L
  genome_ids <- as.character(genome_ids)
  genome_batches <- split(genome_ids, ceiling(seq_along(genome_ids) / batch_size))

  n_cores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
  cluster <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cluster), add = TRUE)

  parallel::clusterExport(
    cluster,
    varlist = c(
      ".extractAMRtable", ".extractGenomeData",
      "abx_filter", "drug_fields",
      "filter_type", "amr_fields", "microtrait_fields",
      "base_dir", "image"
    ),
    envir = environment()
  )
  parallel::clusterEvalQ(cluster, {
    library(tibble); library(dplyr)
  })

  # Pull AMR metadata
  if (isTRUE(verbose)) message("Retrieving AMR phenotype data in batches.")
  batch_drug_data <- parallel::parLapply(cluster, genome_batches, function(batch) {
    .extractAMRtable(
      base_dir          = base_dir,
      batch_genome_IDs  = batch,
      abx_filter        = abx_filter,
      drug_fields       = drug_fields,
      image             = image,
      verbose           = FALSE
    )
  })
  combined_drug_data <- unlist(batch_drug_data, use.names = FALSE)
  if (length(combined_drug_data) == 0) {
    message("No drug data returned.")
    return(NULL)
  }
  combined_drug_data_tbl <- tibble::as_tibble(utils::read.table(
    text = combined_drug_data,
    sep = "\t", header = TRUE, fill = TRUE,
    quote = "", check.names = FALSE, comment.char = "", colClasses = "character"
  )) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ iconv(.x, from = "", to = "UTF-8", sub = ""))) |>
    dplyr::mutate(`genome_drug.genome_id` = as.character(`genome_drug.genome_id`))

  # Pull genome metadata
  if (isTRUE(verbose)) message("Retrieving genome metadata in batches.")
  batch_genome_data <- parallel::parLapply(cluster, genome_batches, function(batch) {
    .extractGenomeData(
      base_dir          = base_dir,
      batch_genome_IDs  = batch,
      filter_type       = filter_type,
      amr_fields        = amr_fields,
      microtrait_fields = microtrait_fields,
      image             = image,
      verbose           = FALSE
    )
  })
  combined_genome_data <- unlist(batch_genome_data, use.names = FALSE)
  if (length(combined_genome_data) == 0) {
    message("No genome data returned.")
    return(NULL)
  }
  combined_genome_data_tbl <- tibble::as_tibble(utils::read.table(
    text = combined_genome_data,
    sep = "\t", header = TRUE, fill = TRUE,
    quote = "", check.names = FALSE, comment.char = "", colClasses = "character"
  )) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ iconv(.x, from = "", to = "UTF-8", sub = ""))) |>
    dplyr::mutate(`genome.genome_id` = as.character(`genome.genome_id`))

  # Write & join
  paths   <- .buildDBpath(base_dir = base_dir, user_bacs = user_bacs, overwrite = overwrite)
  db_path <- paths$db_path
  logs_dir <- file.path(base_dir, "data", "logs")
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("[%s] Writing metadata DuckDB: %s\n", Sys.time(), db_path),
      file = file.path(logs_dir, "bvbrc.log"), append = TRUE)

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)

  DBI::dbWriteTable(con, "amr_phenotype", combined_drug_data_tbl, overwrite = TRUE)
  DBI::dbWriteTable(con, "genome_data",  combined_genome_data_tbl, overwrite = TRUE)

  if (isTRUE(verbose)) message("Joining AMR phenotype and genome metadata.")
  DBI::dbExecute(con, '
    CREATE OR REPLACE TABLE metadata AS
    SELECT *
    FROM amr_phenotype
    INNER JOIN genome_data
      ON amr_phenotype."genome_drug.genome_id" = genome_data."genome.genome_id"
  ')

  # Debug summary after writes
  n_targets   <- length(genome_ids)
  n_amr_ids   <- DBI::dbGetQuery(con, 'SELECT COUNT(DISTINCT "genome_drug.genome_id") AS n FROM amr_phenotype')$n
  n_gmeta_ids <- DBI::dbGetQuery(con, 'SELECT COUNT(DISTINCT "genome.genome_id")      AS n FROM genome_data')$n
  if (isTRUE(verbose)) {
    message("Initial summary: targets=", n_targets,
            " | AMR genomes=", n_amr_ids,
            " | genome_data genomes=", n_gmeta_ids)
  }

  # Tagged view
  .create_amr_tagged_view(con)

  # Final debug summary
  n_bac <- if ("bac_data" %in% DBI::dbListTables(con))
    DBI::dbGetQuery(con, 'SELECT COUNT(DISTINCT "genome.genome_id") AS n FROM bac_data')$n else NA_integer_
  n_amr_ids   <- DBI::dbGetQuery(con, 'SELECT COUNT(DISTINCT "genome_drug.genome_id") AS n FROM amr_phenotype')$n
  n_gmeta_ids <- DBI::dbGetQuery(con, 'SELECT COUNT(DISTINCT "genome.genome_id")      AS n FROM genome_data')$n

  ids_zero_amr <- DBI::dbGetQuery(
    con,
    'WITH sel AS (SELECT DISTINCT "genome.genome_id" AS gid FROM genome_data),
          got AS (SELECT DISTINCT "genome_drug.genome_id" AS gid FROM amr_phenotype)
     SELECT sel.gid
     FROM sel LEFT JOIN got USING (gid)
     WHERE got.gid IS NULL
     ORDER BY sel.gid'
  )$gid

  if (isTRUE(verbose)) {
    message("Final summary:")
    message("  targets      : ", length(genome_ids))
    message("  bac_data     : ", n_bac)
    message("  genome_data  : ", n_gmeta_ids)
    message("  amr_phenotype: ", n_amr_ids)
    message("  genomes with 0 AMR rows: ", length(ids_zero_amr))
    if (length(ids_zero_amr)) {
      message("  e.g.: ", paste(utils::head(ids_zero_amr, 10), collapse = ", "))
    }
  }

  list(duckdbConnection = con, table_name = "metadata")
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
  out <- vapply(lines, function(line) {
    if (grepl("^#", line)) {
      return(line)
    }
    parts <- strsplit(line, "[\t ]", perl = TRUE)[[1]]
    if (length(parts) >= 9) {
      paste(c(parts[1:8], paste(parts[9:length(parts)], collapse = " ")), collapse = "\t")
    } else {
      line
    }
  }, character(1))
  writeLines(out, gff_path, sep = "\n", useBytes = TRUE)
  invisible(TRUE)
}


#' Filter genomes by AMR phenotype and metadata, and store results in DuckDB
#'
#' Preferred path: use per-selection DB "metadata" table (from retrieveMetadata()) and
#' apply lab-evidence & genome_quality filters.
#'
#' Fallback path: if "metadata" is missing and fallback_to_bvbrc_cache = TRUE,
#' read BV-BRC cache at <base_dir>/data/bvbrc/bvbrcData.duckdb ("bvbrc_bac_data"),
#' derive genome IDs from user_bacs (taxon IDs or species substring), and
#' write a minimal "filtered" table (without AMR evidence filtering).
#'
#' @param evidence_mode Character. Either "lab_only" (default), "lab_or_comp" (all),
#' "comp_only" (BV-BRC-predicted without lab labels), or "any" (no AMR data required).
#' @return A list with a DuckDB connection and table_name = "filtered"
#' Filter genomes by AMR phenotype and metadata, and store results in DuckDB
#'
#' Preferred path: use per-selection DB "metadata" table (from retrieveMetadata()) and
#' apply evidence & genome_quality filters.
#'
#' Fallback path: if "metadata" is missing and fallback_to_bvbrc_cache = TRUE,
#' read BV-BRC cache at <base_dir>/data/bvbrc/bvbrcData.duckdb ("bvbrc_bac_data"),
#' derive genome IDs from user_bacs (taxon IDs or species substring), and
#' write a minimal "filtered" table (without AMR evidence filtering).
#'
#' @param evidence_mode Character. One of:
#'   "lab_only"   (default) -> only laboratory evidence
#'   "lab_or_comp"          -> laboratory OR computational evidence
#'   "comp_only"            -> only computational evidence
#'   "any"                  -> no AMR required (from genome_data; Good only)
#' @return A list with a DuckDB connection and table_name = "filtered"
.filterGenomes <- function(user_bacs,
                           base_dir = ".",
                           evidence_mode = c("lab_only","lab_or_comp","comp_only","any"),
                           verbose = TRUE,
                           fallback_to_bvbrc_cache = TRUE) {
  evidence_mode <- match.arg(evidence_mode)
  base_dir <- normalizePath(base_dir, mustWork = FALSE)
  paths    <- .buildDBpath(base_dir = base_dir, user_bacs = user_bacs)
  db_path  <- paths$db_path

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path)
  on.exit({ NULL }, add = TRUE)  # keep open for caller

  # The convenient "Metadata Exists" path
  if ("metadata" %in% DBI::dbListTables(con)) {
    if (isTRUE(verbose)) message("Loading metadata for filtering.")

    # If "any", AMR data not needed, just fetch qualifying genome_data
    if (evidence_mode == "any") {
      gd <- DBI::dbReadTable(con, "genome_data")
      if (is.null(gd) || nrow(gd) == 0) {
        DBI::dbDisconnect(con, shutdown = TRUE)
        stop("No data available in 'genome_data'.")
      }
      gd <- tibble::as_tibble(gd) |>
        dplyr::filter(`genome.genome_quality` == "Good") |>
        dplyr::distinct(`genome.genome_id`)
      # Keep BV-BRC column name
      DBI::dbWriteTable(con, "filtered", gd, overwrite = TRUE)
      if (isTRUE(verbose)) {
        message("Post-filter distinct genomes (any): ", nrow(gd))
        message("Wrote table 'filtered' to: ", db_path)
      }
      return(list(duckdbConnection = con, table_name = "filtered"))
    }

    # Otherwise, open up the metadata and start interrogating the AMR data
    md <- DBI::dbReadTable(con, "metadata")
    if (is.null(md) || nrow(md) == 0) {
      DBI::dbDisconnect(con, shutdown = TRUE)
      message("No data available in 'metadata'.")
      return(NULL)
    }

    md <- tibble::as_tibble(md) |>
      dplyr::mutate(
        `genome_drug.evidence` = dplyr::case_when(
          `genome_drug.laboratory_typing_method` %in%
            c("Disk diffusion", "MIC", "Broth dilution", "Agar dilution") ~ "Laboratory Method",
          `genome_drug.laboratory_typing_method` == "Computational Prediction" ~ "Computational Method",
          TRUE ~ `genome_drug.evidence`
        )
      )

    if (evidence_mode == "lab_only") {
      md <- dplyr::filter(md, `genome_drug.evidence` == "Laboratory Method")
    } else if (evidence_mode == "comp_only") {
      md <- dplyr::filter(
        md,
        `genome_drug.evidence` == "Computational Method" |
          (!is.na(`genome_drug.computational_method`) & nzchar(`genome_drug.computational_method`)) |
          (`genome_drug.laboratory_typing_method` == "Computational Prediction")
      )
    } else {
      # Lab_or_comp: Doesn't matter, keep any predictions
      md <- md
    }

    md <- md |>
      dplyr::filter(`genome.genome_quality` == "Good") |>
      dplyr::filter(`genome_drug.resistant_phenotype` %in% c("Resistant","Susceptible","Intermediate")) |>
      dplyr::distinct(`genome.genome_id`)  # keep BV-BRC column name

    DBI::dbWriteTable(con, "filtered", md, overwrite = TRUE)
    if (isTRUE(verbose)) {
      message("Post-filter distinct genomes: ", nrow(md))
      message("Wrote table 'filtered' to: ", db_path)
    }
    return(list(duckdbConnection = con, table_name = "filtered"))
  }

  # No metadata fallback
  if (!isTRUE(fallback_to_bvbrc_cache)) {
    DBI::dbDisconnect(con, shutdown = TRUE)
    stop("No 'metadata' table found in ", db_path, ". Run retrieveMetadata() first.")
  }
  if (isTRUE(verbose)) message("No 'metadata' in per-selection DB. Falling back to BV-BRC cache at data/bvbrc/.")

  cache_db <- file.path(base_dir, "data", "bvbrc", "bvbrcData.duckdb")
  if (!file.exists(cache_db)) {
    DBI::dbDisconnect(con, shutdown = TRUE)
    stop("BV-BRC cache not found at: ", cache_db, ". Run .updateBVBRCdata() first.")
  }

  con_cache <- DBI::dbConnect(duckdb::duckdb(), dbdir = cache_db, read_only = TRUE)
  on.exit(try(DBI::dbDisconnect(con_cache, shutdown = TRUE), silent = TRUE), add = TRUE)

  if (!"bvbrc_bac_data" %in% DBI::dbListTables(con_cache)) {
    DBI::dbDisconnect(con, shutdown = TRUE)
    stop("Table 'bvbrc_bac_data' not found in BV-BRC cache: ", cache_db)
  }

  bv  <- tibble::as_tibble(DBI::dbReadTable(con_cache, "bvbrc_bac_data"))
  sel <- tibble::tibble(`genome.genome_id` = character())

  for (v in user_bacs) {
    if (.id_checker(v)) {
      v_chr   <- as.character(v)
      matches <- bv[bv$`genome.taxon_id` == v_chr, , drop = FALSE]
    } else {
      matches <- bv[stringr::str_detect(
        bv$`genome.species`,
        stringr::fixed(v, ignore_case = TRUE)
      ), , drop = FALSE]
    }
    if (nrow(matches)) {
      sel <- dplyr::bind_rows(sel,
                              tibble::tibble(`genome.genome_id` = as.character(matches$`genome.genome_id`)))
    }
  }
  sel <- dplyr::distinct(sel)

  if (nrow(sel) == 0L) {
    DBI::dbDisconnect(con, shutdown = TRUE)
    stop("No genomes matched user_bacs in BV-BRC cache.")
  }

  DBI::dbWriteTable(con, "filtered", sel, overwrite = TRUE)
  if (isTRUE(verbose)) message("Wrote table 'filtered' to: ", db_path)

  list(duckdbConnection = con, table_name = "filtered")
}

#' Helps check if a complete set exists after DL (.fna + .PATRIC.faa + .PATRIC.gff)
#' @keywords internal
.is_complete_set <- function(dir, genomeID, min_bytes = 100) {
  fna <- file.path(dir, paste0(genomeID, ".fna"))
  faa <- file.path(dir, paste0(genomeID, ".PATRIC.faa"))
  gff <- file.path(dir, paste0(genomeID, ".PATRIC.gff"))
  paths <- c(fna, faa, gff)
  all(file.exists(paths)) &&
    all(vapply(paths, function(x) file.info(x)$size, numeric(1)) > min_bytes)
}

#' Helps collate completed genomes into a set
#' @keywords internal
.list_complete <- function(dir, genome_ids, min_bytes = 100) {
  genome_ids[vapply(genome_ids, .is_complete_set, logical(1), dir = dir, min_bytes = min_bytes)]
}

#' Helps in auditing downloaded files to ensure everything's complete per ID
#' @keywords internal
.audit_gaps <- function(out_dir, ids, min_bytes = 100) {
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  df <- data.frame(
    genome = ids,
    fna = file.exists(file.path(out_dir, paste0(ids, ".fna"))) &
      (file.info(file.path(out_dir, paste0(ids, ".fna")))$size > min_bytes),
    faa = file.exists(file.path(out_dir, paste0(ids, ".PATRIC.faa"))) &
      (file.info(file.path(out_dir, paste0(ids, ".PATRIC.faa")))$size > min_bytes),
    gto = file.exists(file.path(out_dir, paste0(ids, ".gto"))) &
      (file.info(file.path(out_dir, paste0(ids, ".gto")))$size > min_bytes),
    gff = file.exists(file.path(out_dir, paste0(ids, ".PATRIC.gff"))) &
      (file.info(file.path(out_dir, paste0(ids, ".PATRIC.gff")))$size > min_bytes),
    stringsAsFactors = FALSE
  )
  df$complete <- with(df, fna & faa & gff)
  df
}

### BV-BRC CLI downloader [slower by comparison, but does not need FTP server]

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

#' Using p3-dump-genomes in CLI to fetch FASTA and .gto files
#' @keywords internal
.cli_dump_fastas_gto_chunk <- function(image, out_dir, genome_ids, tag, tries = 3L) {
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ids_file <- file.path(out_dir, paste0("ids_", tag, ".txt"))
  writeLines(genome_ids, ids_file)

  shell <- .pick_shell(image)
  mount <- .docker_path(out_dir)

  # Safety against Windows-specific CRLF lines before `p3-dump-genomes`
  sh_cmd <- sprintf(
    'tr -d "\\r" < /out/%s | p3-dump-genomes --outDir /out --fasta --prot --gto -',
    basename(ids_file)
  )

  args <- c("run", "--rm", "-v", paste0(mount, ":/out"), image, shell, "-lc", shQuote(sh_cmd))
  res <- suppressWarnings(system2("docker", args = args, stdout = TRUE, stderr = TRUE))
  st <- attr(res, "status")
  if (is.null(st)) st <- 0L

  if (st != 0L && tries > 1L) {
    Sys.sleep(1)
    return(.cli_dump_fastas_gto_chunk(image, out_dir, genome_ids, tag, tries - 1L))
  }

  # Normalize filenames for this chunk (ensure .fna & PATRIC.faa)
  for (gid in genome_ids) {
    fa <- file.path(out_dir, paste0(gid, ".fa"))
    if (file.exists(fa)) file.rename(fa, file.path(out_dir, paste0(gid, ".fna")))
    fasta <- file.path(out_dir, paste0(gid, ".fasta"))
    if (file.exists(fasta)) file.rename(fasta, file.path(out_dir, paste0(gid, ".fna")))
    faa <- file.path(out_dir, paste0(gid, ".faa"))
    if (file.exists(faa)) file.rename(faa, file.path(out_dir, paste0(gid, ".PATRIC.faa")))
  }

  # Apply sanitizer to FASTA files
  for (gid in genome_ids) {
    fna <- file.path(out_dir, paste0(gid, ".fna"))
    if (file.exists(fna)) .strip_fasta_preamble(fna)
  }

  st == 0L
}

#' Exports GFF files from the downloaded GTO per genome in each chunk
#' @keywords internal
.cli_export_gff_chunk <- function(image, out_dir, genome_ids, tag, tries = 3L) {
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  ids_file <- file.path(out_dir, paste0("ids_", tag, ".txt"))
  writeLines(genome_ids, ids_file)

  exporter <- "/usr/bin/rast-export-genome"
  shell <- .pick_shell(image)
  mount <- .docker_path(out_dir)
  stderr_file <- file.path(out_dir, paste0("gff_chunk_", tag, ".stderr.txt"))

  # GFFs from BV-BRC .gto export are not directly compatible in Panaroo
  # This block reformats the contig IDs and ensures .gff/.fna pairs work together
  sh_cmd <- paste0(
    "set -euo pipefail; ",
    "fail_n=0; : > /out/", basename(stderr_file), "; ",
    'while IFS= read -r b || [ -n "$b" ]; do ',
    '  b=${b%$\'\\r\'}; [ -n "$b" ] || continue; ',
    '  gto="/out/${b}.gto"; gff="/out/${b}.PATRIC.gff"; map="/out/${b}.orig2id.tsv"; ',
    '  if [ ! -s "$gto" ]; then echo "MISSING_GTO $b" >>/out/', basename(stderr_file), "; continue; fi; ",

    # Export GFF
    '  [ -s "$gff" ] || ', exporter, ' -i "$gto" -o "$gff" gff ',
    '     || { echo "EXPORT_FAIL_GFF $b" >>/out/', basename(stderr_file), "; fail_n=$((fail_n+1)); continue; }; ",

    # Build mapping original_id -> id, and default to id if original_id missing
    "  if command -v jq >/dev/null 2>&1; then ",
    '    jq -r \'.contigs[] | [(.original_id // .id), .id] | @tsv\' "$gto" > "$map"; ',
    "  else ",
    '    python3 - "$gto" > "$map" <<\'PY\'\n',
    "import sys, json\n",
    "g = json.load(open(sys.argv[1]))\n",
    'for c in g.get("contigs", []):\n',
    '    o = c.get("original_id") or c.get("id")\n',
    '    i = c.get("id")\n',
    "    if o and i:\n",
    '        print(f"{o}\\t{i}")\n',
    "PY\n",
    "  fi; ",

    # Relabel GFF sequence IDs: original_id -> id
    "  awk 'FNR==NR{m[$1]=$2; next} ",
    "       /^##sequence-region/ { if ($2 in m) {$2=m[$2]} print; next } ",
    "       /^#/ { print; next } ",
    '       { if ($1 in m) $1=m[$1]; print }\' "$map" "$gff" > "${gff}.tmp" && mv "${gff}.tmp" "$gff"; ',
    "done < /out/", basename(ids_file), "; ",
    "exit 0"
  )

  args <- c("run", "--rm", "-v", paste0(mount, ":/out"), image, shell, "-lc", shQuote(sh_cmd))
  res <- suppressWarnings(system2("docker", args = args, stdout = TRUE, stderr = TRUE))
  st <- attr(res, "status")
  if (is.null(st)) st <- 0L

  if (st != 0L && tries > 1L) {
    Sys.sleep(1)
    return(.cli_export_gff_chunk(image, out_dir, genome_ids, tag, tries - 1L))
  }

  # Apply sanitizer for GFFs after they've been extracted
  for (gid in genome_ids) {
    gff <- file.path(out_dir, paste0(gid, ".PATRIC.gff"))
    if (file.exists(gff)) .sanitize_gff(gff)
  }

  TRUE
}

#' Download .fna, .faa, .gff files for filtered BV-BRC genomes
#'
#' Default and fast method="ftp"
#'
#' Alternative path to bypass FTP: method="cli"
#'
#' @param base_dir Project root (results layout preserved).
#' @param user_bacs Input label(s) used to locate per-selection DB path.
#' @param method "ftp" (default) or "cli".
#' @param image Docker image for CLI path (default "danylmb/bvbrc:5.3").
#' @param skip_existing Logical; if TRUE, do not re-download genomes already complete. Default TRUE.
#' @param ftp_workers Parallel workers for FTP path (default 8).
#' @param cli_fasta_workers Parallel chunk containers for FASTA+GTO (default 4).
#' @param cli_gff_workers Parallel chunk containers for GFF export (default 4).
#' @param chunk_size Genomes per chunk container (default 50).
#' @param verbose Verbose messages.
#' @return Character vector of genome IDs with complete file sets on disk.
retrieveGenomes <- function(base_dir = ".",
                            user_bacs,
                            method = c("ftp", "cli"),
                            image = "danylmb/bvbrc:5.3",
                            skip_existing = TRUE,
                            ftp_workers = 4L,
                            cli_fasta_workers = 4L,
                            cli_gff_workers = 4L,
                            chunk_size = 50L,
                            evidence_mode = c("lab_only","lab_or_comp","comp_only","any"),  # NEW
                            verbose = TRUE) {
  method <- match.arg(method)
  evidence_mode <- match.arg(evidence_mode)
  base_dir <- normalizePath(base_dir, mustWork = FALSE)

  # Use 'filtered' if already prepared, or start filtering
  if (isTRUE(verbose)) message("Preparing download set (checking for existing 'filtered').")
  paths <- .buildDBpath(base_dir = base_dir, user_bacs = user_bacs)
  db_path <- paths$db_path
  con0 <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path)
  has_filtered <- "filtered" %in% DBI::dbListTables(con0)
  if (has_filtered) {
    if (isTRUE(verbose)) message("Using existing 'filtered' table (skipping re-filter).")
    con <- con0; tbl <- "filtered"
  } else {
    if (isTRUE(verbose)) message("No 'filtered' table found; filtering now.")
    f_out <- .filterGenomes(base_dir = base_dir,
                            user_bacs = user_bacs,
                            evidence_mode = evidence_mode,
                            verbose = verbose)
    con <- f_out$duckdbConnection
    tbl <- f_out$table_name
  }

  ids <- tibble::as_tibble(DBI::dbReadTable(con, tbl)) |>
    dplyr::distinct(`genome.genome_id`) |>
    dplyr::pull(`genome.genome_id`)

  bug_dir <- dirname(db_path)
  genome_path <- file.path(bug_dir, "genomes")
  logs_dir <- file.path(base_dir, "data", "logs")
  dir.create(genome_path, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(skip_existing)) {
    already <- .list_complete(genome_path, ids)
    if (isTRUE(verbose)) message(length(already), " genomes already completed; skipping.")
    ids <- setdiff(ids, already)
  }

  if (length(ids) == 0L) {
    if (isTRUE(verbose)) message("All genomes already complete.")
    all_complete <- .list_complete(
      genome_path,
      tibble::as_tibble(DBI::dbReadTable(con, tbl)) |>
        dplyr::distinct(`genome.genome_id`) |>
        dplyr::pull(`genome.genome_id`)
    )
    return(all_complete)
  }

  if (identical(method, "ftp")) {
    if (isTRUE(verbose)) message("Downloading by FTPS.")
    ok_ids <- .ftpes_download_two_pass(
      genome_ids     = ids,
      out_dir        = genome_path,
      workers_first  = ftp_workers,
      workers_second = ftp_workers,
      log_file       = file.path(logs_dir, "ftp_download.log")
    )
    if (isTRUE(verbose)) {
      message("Complete file sets for ", length(ok_ids), " genomes total.")
      miss <- setdiff(ids, ok_ids)
      if (length(miss)) message("Excluded (failed after retry): ", length(miss))
    }
    return(ok_ids)
  }

  if (identical(method, "cli")) {
    if (isTRUE(verbose)) {
      message(
        "CLI mode: downloading genomes in parallelized chunks. ",
        "Chunks=", ceiling(length(ids)/chunk_size),
        " | fasta_workers=", cli_fasta_workers,
        " | gff_workers=", cli_gff_workers
      )
    }

    chunks <- split(ids, ceiling(seq_along(ids) / chunk_size))

    # FASTA + GTO in parallel containers
    future::plan(future::multisession, workers = max(1, cli_fasta_workers))
    invisible(future.apply::future_mapply(
      FUN = function(vec, tag) .cli_dump_fastas_gto_chunk(image, genome_path, vec, tag),
      vec = chunks, tag = paste0("fa", seq_along(chunks)),
      SIMPLIFY = TRUE, future.seed = TRUE
    ))

    # GFF export in parallel containers
    future::plan(future::multisession, workers = max(1, cli_gff_workers))
    invisible(future.apply::future_mapply(
      FUN = function(vec, tag) .cli_export_gff_chunk(image, genome_path, vec, tag),
      vec = chunks, tag = paste0("gff", seq_along(chunks)),
      SIMPLIFY = TRUE, future.seed = TRUE
    ))

    ok_ids <- ids[vapply(ids, .is_complete_set, logical(1), dir = genome_path)]
    if (isTRUE(verbose)) {
      message("CLI complete file sets for ", length(ok_ids), " genomes.")
      miss <- setdiff(ids, ok_ids)
      if (length(miss)) message("Excluded (CLI completion failed): ", length(miss))
    }
    return(ok_ids)
  }
}

#' Build a table of local genome file paths and write to DuckDB
#'
#' Scans <base_dir>/data/<bug_dir>/genomes for *.PATRIC.gff, *.fna, *.PATRIC.faa,
#' verifies size > 100 bytes, assembles rows per genome, and writes to the per-bug
#' DuckDB table "files". This writes a Panaroo input file at:
#'   <base_dir>/data/<bug_dir>/<abbrev>.txt
#'
#' @param base_dir Character. Project root.
#' @param user_bacs Character vector. Used to locate per-bug directories and DB.
#' @param verbose Logical. If TRUE, prints messages.
#'
#' @return A list with duckdbConnection and table_name = "files".
genomeList <- function(base_dir = ".",
                       user_bacs,
                       verbose = TRUE) {
  base_dir <- normalizePath(base_dir, mustWork = FALSE)
  paths <- .buildDBpath(base_dir = base_dir, user_bacs = user_bacs)
  db_path <- paths$db_path
  bug_dir <- dirname(db_path)

  genome_path <- file.path(bug_dir, "genomes")
  files_all <- list.files(genome_path, full.names = TRUE)
  files_all <- files_all[file.info(files_all)$size > 100]

  # Separate by type
  gff_files <- files_all[grepl("\\.PATRIC\\.gff$", files_all)]
  fna_files <- files_all[grepl("\\.fna$", files_all)]
  faa_files <- files_all[grepl("\\.PATRIC\\.faa$", files_all)]

  gff_ids <- sub("\\.PATRIC\\.gff$", "", basename(gff_files))
  fna_ids <- sub("\\.fna$", "", basename(fna_files))
  faa_ids <- sub("\\.PATRIC\\.faa$", "", basename(faa_files))

  genome_ids <- unique(c(gff_ids, fna_ids, faa_ids))

  list_of_files <- lapply(genome_ids, function(genomeID) {
    gff_path <- file.path(genome_path, paste0(genomeID, ".PATRIC.gff"))
    fna_path <- file.path(genome_path, paste0(genomeID, ".fna"))
    faa_path <- file.path(genome_path, paste0(genomeID, ".PATRIC.faa"))

    data.frame(
      genome_id = genomeID,
      gff_path = if (file.exists(gff_path) && file.info(gff_path)$size > 100) gff_path else NA,
      fna_path = if (file.exists(fna_path) && file.info(fna_path)$size > 100) fna_path else NA,
      faa_path = if (file.exists(faa_path) && file.info(faa_path)$size > 100) faa_path else NA,
      panaroo_input = if (
        file.exists(gff_path) && file.exists(fna_path) && file.exists(faa_path) &&
        file.info(gff_path)$size > 100 && file.info(fna_path)$size > 100 && file.info(faa_path)$size > 100
      ) {
        paste(gff_path, fna_path)
      } else {
        NA
      },
      stringsAsFactors = FALSE
    )
  })
  list_of_files <- do.call(rbind, list_of_files)
  list_of_files <- tibble::as_tibble(list_of_files) |>
    dplyr::filter(!is.na(panaroo_input))

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path)
  DBI::dbWriteTable(con, "files", list_of_files, overwrite = TRUE)

  # Write Panaroo input next to the DB
  abbrev <- .generateDBname(user_bacs)
  panaroo_txt <- file.path(bug_dir, paste0(abbrev, ".txt"))
  writeLines(na.omit(list_of_files$panaroo_input), con = panaroo_txt)

  if (isTRUE(verbose)) {
    message("Wrote table 'files' and Panaroo input to: ", bug_dir)
  }

  list(duckdbConnection = con, table_name = "files")
}


#' Download and prepare all files for a chosen bacterial species or TaxID
#'
#' This wrapper takes `user_bacs` input (species names and/or taxon IDs),
#' retrieves the corresponding filtered genome set from BV-BRC, downloads all
#' required genome files (.fna, .PATRIC.faa, .PATRIC.gff), and produces the
#' Panaroo input table via `genomeList()`. These outputs are used for the
#' data_processing.R script next.
#'
#' Internally, this runs:
#'   1. `retrieveGenomes()`  – filters BV-BRC metadata, selects genomes, downloads files.
#'   2. `genomeList()`       – scans downloaded files and writes a "files" table
#'                             in the per-selection DuckDB.
#'
#' The per-selection DuckDB is located automatically under:
#'   <base_dir>/data/<bug_dir>/<abbrev>.duckdb
#'
#' @param user_bacs Character vector. Species and/or taxon IDs (e.g.
#'   `c("Shigella flexneri", "623")`).
#' @param base_dir Character. Project root directory. Default `"."`.
#' @param method Character. Download method passed to `retrieveGenomes()`.
#'   `"ftp"` (default) or `"cli"`.
#' @param overwrite Logical. Passed to metadata filtering and DuckDB creation.
#'   Default FALSE.
#' @param evidence_mode Character. Sets what types of AMR evidence is acceptable.
#'    Default `lab_only`. `any` will not require AMR data for downloads. This will
#'    return very large download lists for many species!
#' @param verbose Logical. Print progress messages. Default TRUE.
#'
#' @return A list (the output of `genomeList()`), containing:
#'   - `duckdbConnection`  Active DBI connection to the per-bug DuckDB
#'   - `table_name`        `"files"`
#'
#' @export
prepareGenomes <- function(user_bacs,
                           base_dir = ".",
                           method = c("ftp", "cli"),
                           overwrite = FALSE,
                           evidence_mode = c("lab_only","lab_or_comp","comp_only","any"),
                           verbose = TRUE) {
  method <- match.arg(method)
  evidence_mode <- match.arg(evidence_mode)
  base_dir <- normalizePath(base_dir, mustWork = FALSE)

  .ensure_bvbrc_cache(base_dir = base_dir, verbose = verbose)

  if (isTRUE(verbose)) message("Step 0: Building AMR metadata (retrieveMetadata)")
  invisible(retrieveMetadata(
    user_bacs = user_bacs,
    filter_type = "AMR",
    base_dir = base_dir,
    abx = "All",
    overwrite = overwrite,
    verbose = verbose
  ))

  if (isTRUE(verbose)) message("Step 1: Filtering genomes for download by evidence: ", evidence_mode)
  f_out <- .filterGenomes(
    base_dir      = base_dir,
    user_bacs     = user_bacs,
    evidence_mode = evidence_mode,
    verbose       = verbose,
    fallback_to_bvbrc_cache = FALSE
  )
  if (is.null(f_out)) {
    message("No genomes available after evidence filtering.")
    return(NULL)
  }

  # A little summary of what's left after filtering (or not)
  paths <- .buildDBpath(base_dir = base_dir, user_bacs = user_bacs)
  con   <- DBI::dbConnect(duckdb::duckdb(), dbdir = paths$db_path, read_only = TRUE)
  n_filtered <- DBI::dbGetQuery(con, 'SELECT COUNT(DISTINCT "genome.genome_id") AS n FROM filtered')$n
  n_meta     <- if ("genome_data" %in% DBI::dbListTables(con)) DBI::dbGetQuery(con, 'SELECT COUNT(DISTINCT "genome.genome_id") AS n FROM genome_data')$n else NA
  n_amr      <- if ("amr_phenotype" %in% DBI::dbListTables(con)) DBI::dbGetQuery(con, 'SELECT COUNT(DISTINCT "genome_drug.genome_id") AS n FROM amr_phenotype')$n else NA
  DBI::dbDisconnect(con, shutdown = TRUE)
  if (isTRUE(verbose)) {
    message(sprintf("Evidence filter summary: filtered=%d | genomes with AMR=%s | genomes with genome_data=%s",
                    n_filtered, ifelse(is.na(n_amr), "NA", n_amr), ifelse(is.na(n_meta), "NA", n_meta)))
  }

  if (isTRUE(verbose)) message("Step 2: Downloading genomes from BV-BRC (", method, ")")
  ids <- retrieveGenomes(
    base_dir      = base_dir,
    user_bacs     = user_bacs,
    method        = method,
    skip_existing = !overwrite,
    evidence_mode = evidence_mode,
    verbose       = verbose
  )
  if (length(ids) == 0L) {
    message("No genomes downloaded.")
    return(NULL)
  }

  if (isTRUE(verbose)) message("Step 3: Formatting data into a database for further processing")
  out <- genomeList(
    base_dir  = base_dir,
    user_bacs = user_bacs,
    verbose   = verbose
  )

  if (isTRUE(verbose)) {
    message("Done. Files are ready! Continue with downstream processing with runDataProcessing().")
  }
  out
}
