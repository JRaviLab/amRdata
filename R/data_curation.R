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
    "--attr genome_id,genome_name,taxon_id,species,strain,_version_,",
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
  
  # Parse
  df <- utils::read.table(text = raw_data, sep = "\t", header = TRUE, fill = TRUE,
                          quote = "", check.names = FALSE, comment.char = "")
  df <- tibble::as_tibble(df) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ iconv(.x, from = "", to = "UTF-8", sub = "")))
  
  if (isTRUE(verbose)) {
    message(glue::glue("Retrieved {nrow(df)} rows x {ncol(df)} columns."))
  }
  
  return(df)
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
  logs_dir  <- file.path(data_dir, "logs")
  
  dir.create(bvbrc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir,  recursive = TRUE, showWarnings = FALSE)
  
  db_path    <- file.path(bvbrc_dir, "bvbrcData.duckdb")
  table_name <- "bvbrc_bac_data"
  meta_table <- "__meta"
  
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path)
  
  DBI::dbExecute(
    con,
    glue::glue('CREATE TABLE IF NOT EXISTS {meta_table} (
                 table_name TEXT PRIMARY KEY,
                 last_updated TIMESTAMP
               )')
  )
  
  # Tiny update check helper
  .last_updated <- function(con, table_name) {
    res <- tryCatch(
      DBI::dbGetQuery(
        con,
        glue::glue('SELECT last_updated FROM {meta_table}
                    WHERE table_name = {DBI::dbQuoteString(con, table_name)}')
      ),
      error = function(e) NULL
    )
    if (is.null(res) || nrow(res) == 0L) return(NA)
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
    if (suppressWarnings(!is.na(as.numeric(user_bac)))) {
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



#' Resolve `query value` from `user_bacs` for [getGenomeIDs()]
#'
#' If query_value is NULL, derive it from user_bacs based on query_type.
#' For species/genome_name: take the first element of user_bacs.
#' For taxon_id: take the first numeric-looking element of user_bacs.
#' If nothing suitable is found, throw a fit and an error.
#' @keywords internal
.resolveQueryValue <- function(query_type, query_value, user_bacs) {
  if (!is.null(query_value) && nzchar(query_value)) return(query_value)
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
    nums <- suppressWarnings(!is.na(as.numeric(user_bacs)))
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
#' \dontrun{
#' .generateDBname(c("90371", "Bacillus subtilis"))
#' .generateDBname(c("12345", "Escherichia coli", "Lactobacillus"))
#' }
#'
.generateDBname <- function(user_bacs) {
  db_parts <- c()
  
  for (user_bac in user_bacs) {
    if (suppressWarnings(!is.na(as.numeric(user_bac)))) {
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
  
  db_dir  <- file.path(data_dir, bug_dirname)
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
getGenomeIDs <- function(base_dir = ".",
                         query_type = c("genome_name", "species", "taxon_id"),
                         query_value = NULL,
                         user_bacs,
                         overwrite = FALSE,
                         image = "danylmb/bvbrc:5.3",
                         verbose = TRUE) {
  
  query_type  <- match.arg(query_type)
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
  count_lines  <- tryCatch(system(count_cmd, intern = TRUE), error = function(e) character())
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
      quote = "", check.names = FALSE, comment.char = ""
    )
  ) |>
    dplyr::mutate(
      `genome.genome_id`   = suppressWarnings(as.numeric(`genome.genome_id`)),
      `genome.genome_name` = as.character(`genome.genome_name`),
      `genome.taxon_id`    = suppressWarnings(as.integer(`genome.taxon_id`)),
      `genome.species`     = as.character(`genome.species`),
      `genome.strain`      = as.character(`genome.strain`)
    )
  
  # Per-bug DB path
  paths   <- .buildDBpath(base_dir = base_dir, user_bacs = user_bacs, overwrite = overwrite)
  db_path <- paths$db_path
  
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path)
  DBI::dbWriteTable(con, "bac_data", data_result, overwrite = TRUE)
  
  if (isTRUE(verbose)) message("Wrote table 'bac_data' to: ", db_path)
  
  list(count_result = count_result, duckdbConnection = con, table_name = "bac_data")
}

#' Retrieve genome IDs for each taxon via BV-BRC and DuckDB
#'
#' Resolves user-provided taxa to taxon IDs, queries BV-BRC per unique taxon ID,
#' and returns distinct genome IDs. Uses a per-selection DuckDB located under:
#'   <base_dir>/data/<bug_dir>/<abbrev>.duckdb
#' BV-BRC column names are preserved.
#'
#' @param base_dir Character. Project root directory. Default ".".
#' @param user_bacs Character vector. Mixed inputs of taxon IDs and/or species names.
#' @param overwrite Logical. If FALSE and the DuckDB already exists for this selection, abort. Default: FALSE.
#' @param verbose Logical.
#'
#' @return A numeric vector of distinct `genome.genome_id`, or NULL if none found.
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
  
  # Resolve per-bug DB path (non-enforcing)
  paths   <- .buildDBpath(base_dir = base_dir, user_bacs = user_bacs, overwrite = overwrite)
  db_path <- paths$db_path
  
  bac_data  <- tibble::tibble()
  taxon_ids <- unique(bac_input_data$genome.taxon_id)
  
  if (isTRUE(verbose)) message("Querying BV-BRC for ", length(taxon_ids), " taxon IDs.")
  
  for (i in seq_along(taxon_ids)) {
    tax <- taxon_ids[i]
    if (isTRUE(verbose)) message("Taxon ", i, "/", length(taxon_ids), ": ", tax)
    
    res <- getGenomeIDs(
      base_dir    = base_dir,
      query_type  = "taxon_id",
      query_value = as.character(tax),
      user_bacs   = user_bacs,
      overwrite   = TRUE,   # per-iteration table overwrite is OK
      verbose     = verbose
    )
    
    con <- res$duckdbConnection
    tbl <- res$table_name
    each_bac_data <- tibble::as_tibble(DBI::dbReadTable(con, tbl))
    bac_data <- dplyr::bind_rows(bac_data, each_bac_data)
  }
  
  if (nrow(bac_data) > 0) {
    genome_ids <- bac_data |>
      dplyr::distinct(`genome.genome_id`) |>
      dplyr::pull(`genome.genome_id`)
    genome_ids <- genome_ids[!is.na(genome_ids)]
    
    if (length(genome_ids) > 0) {
      if (isTRUE(verbose)) message("Collected ", length(genome_ids), " distinct genome IDs.")
      return(genome_ids)
    } else {
      message("No valid genome IDs found.")
      return(NULL)
    }
  } else {
    message("No valid genome data found for any input.")
    return(NULL)
  }
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
  tmp_dir  <- file.path(data_dir, "tmp")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (isTRUE(verbose)) {
    message("Preparing AMR query input for ", length(batch_genome_IDs), " genomes.")
  }
  
  docker_path <- Sys.which("docker")
  if (!nzchar(docker_path)) {
    stop("Docker is not available on your PATH but is required.")
  }
  
  # Generate genome list with p3-echo (title must match drug table key)
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
  
  # Query drug data
  drug_args <- c(
    "run", "--rm",
    "-v", paste0(data_dir, ":/data"),
    image, "p3-get-genome-drugs",
    "--input", tmp_in_mounted,
    abx_filter,
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
  tmp_dir  <- file.path(data_dir, "tmp")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (isTRUE(verbose)) {
    message("Preparing genome metadata input for ", length(batch_genome_IDs), " genomes.")
  }
  
  docker_path <- Sys.which("docker")
  
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
  
  # Choose attributes (AMR for this pipeline)
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
#'   - metadata (inner join on genome ID field names as returned by BV-BRC)
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
  genome_ids <- .retrieveQueryIDs(base_dir = base_dir, user_bacs = user_bacs,
                                  overwrite = overwrite, verbose = verbose)
  if (length(genome_ids) == 0) {
    message("No genome IDs available for the specified inputs.")
    return(NULL)
  }
  
  # Fields (unchanged BV-BRC names)
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
  abx_filter <- if (identical(abx, "All")) "--required antibiotic"
  else paste0("--in antibiotic,", paste(abx, collapse = ","))
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
  
  # Batching (as before)
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
  parallel::clusterEvalQ(cluster, { library(tibble); library(dplyr) })
  
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
  if (length(combined_drug_data) == 0) { message("No drug data returned."); return(NULL) }
  
  combined_drug_data_tbl <- tibble::as_tibble(
    utils::read.table(
      text = combined_drug_data,
      sep = "\t", header = TRUE, fill = TRUE,
      quote = "", check.names = FALSE, comment.char = ""
    )
  ) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ iconv(.x, from = "", to = "UTF-8", sub = ""))) |>
    dplyr::mutate(`genome_drug.genome_id` = as.character(`genome_drug.genome_id`))
  
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
  if (length(combined_genome_data) == 0) { message("No genome data returned."); return(NULL) }
  
  combined_genome_data_tbl <- tibble::as_tibble(
    utils::read.table(
      text = combined_genome_data,
      sep = "\t", header = TRUE, fill = TRUE,
      quote = "", check.names = FALSE, comment.char = ""
    )
  ) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ iconv(.x, from = "", to = "UTF-8", sub = ""))) |>
    dplyr::mutate(`genome.genome_id` = as.character(`genome.genome_id`))
  
  # Per-bug DB path (reuse; no enforcement)
  paths   <- .buildDBpath(base_dir = base_dir, user_bacs = user_bacs, overwrite = overwrite)
  db_path <- paths$db_path
  
  logs_dir <- file.path(base_dir, "data", "logs")
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("[%s] Writing metadata DuckDB: %s\n", Sys.time(), db_path),
      file = file.path(logs_dir, "bvbrc.log"), append = TRUE)
  
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path)
  DBI::dbWriteTable(con, "amr_phenotype", combined_drug_data_tbl, overwrite = TRUE)
  DBI::dbWriteTable(con, "genome_data", combined_genome_data_tbl, overwrite = TRUE)
  
  if (isTRUE(verbose)) message("Joining AMR phenotype and genome metadata.")
  initial_metadata <- tibble::as_tibble(DBI::dbGetQuery(
    con,
    'SELECT *
     FROM amr_phenotype
     INNER JOIN genome_data
     ON amr_phenotype."genome_drug.genome_id" = genome_data."genome.genome_id"'
  ))
  DBI::dbWriteTable(con, "metadata", initial_metadata, overwrite = TRUE)
  
  if (isTRUE(verbose)) {
    message("Wrote tables 'amr_phenotype', 'genome_data', and 'metadata' to: ", db_path)
  }
  
  list(duckdbConnection = con, table_name = "metadata")
}


#' Filter genomes by AMR phenotype and metadata, and store results in DuckDB
#'
#' Reads the per-selection DuckDB at:
#'   <base_dir>/data/<bug_dir>/<abbrev>.duckdb
#' Expects a table "metadata" (from retrieveMetadata). Filters to lab-tested evidence,
#' genome_quality == "Good", and resistant_phenotype in {Resistant, Susceptible, Intermediate}.
#'
#' @param base_dir Character. Project root.
#' @param user_bacs Character vector. Used to locate the per-selection DuckDB.
#' @param verbose Logical. If TRUE, prints messages.
#'
#' @return A list with: duckdbConnection and table_name = "filtered"
filterGenomes <- function(base_dir = ".",
                          user_bacs,
                          verbose = TRUE) {
  base_dir <- normalizePath(base_dir, mustWork = FALSE)
  paths    <- .buildDBpath(base_dir = base_dir, user_bacs = user_bacs)
  db_path  <- paths$db_path
  
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path)
  
  if (!"metadata" %in% DBI::dbListTables(con)) {
    DBI::dbDisconnect(con, shutdown = TRUE)
    stop("No 'metadata' table found in ", db_path, ". Run retrieveMetadata() first.")
  }
  
  if (isTRUE(verbose)) message("Loading metadata for filtering.")
  initial_metadata <- DBI::dbReadTable(con, "metadata")
  
  if (is.null(initial_metadata) || nrow(initial_metadata) == 0) {
    DBI::dbDisconnect(con, shutdown = TRUE)
    message("No data available in 'metadata'.")
    return(NULL)
  }
  
  # Map evidence: lab methods to "Laboratory Method"; comp predictions to "Computational Method"
  initial_metadata <- tibble::as_tibble(initial_metadata) |>
    dplyr::mutate(
      `genome_drug.evidence` = dplyr::case_when(
        `genome_drug.laboratory_typing_method` %in%
          c("Disk diffusion", "MIC", "Broth dilution", "Agar dilution") ~ "Laboratory Method",
        `genome_drug.laboratory_typing_method` == "Computational Prediction" ~ "Computational Method",
        TRUE ~ `genome_drug.evidence`
      )
    )
  
  # Filtering
  filtered_metadata <- initial_metadata |>
    dplyr::filter(`genome_drug.evidence` == "Laboratory Method") |>
    dplyr::filter(`genome.genome_quality` == "Good") |>
    dplyr::filter(`genome_drug.resistant_phenotype` %in% c("Resistant", "Susceptible", "Intermediate"))
  
  DBI::dbWriteTable(con, "filtered", filtered_metadata, overwrite = TRUE)
  
  if (nrow(filtered_metadata) == 0) {
    if (isTRUE(verbose)) message("No genomes matched the filtering criteria.")
    return(list(duckdbConnection = con, table_name = "filtered"))
  }
  
  if (isTRUE(verbose)) {
    message("Wrote table 'filtered' to: ", db_path)
  }
  
  list(duckdbConnection = con, table_name = "filtered")
}

#' Download genome files (PATRIC.GFF, FNA, PATRIC.FAA) for filtered BV-BRC genomes
#'
#' Uses filterGenomes() to get filtered metadata, then downloads to:
#'   <base_dir>/data/<bug_dir>/genomes/
#'
#' @param base_dir Character. Project root.
#' @param user_bacs Character vector. Used to locate the DuckDB and output directory.
#' @param verbose Logical. If TRUE, prints messages.
#'
#' @return Character vector of genome IDs with all expected files successfully downloaded.
retrieveGenomes <- function(base_dir = ".",
                            user_bacs,
                            verbose = TRUE) {
  
  base_dir <- normalizePath(base_dir, mustWork = FALSE)
  paths    <- .buildDBpath(base_dir = base_dir, user_bacs = user_bacs)
  db_path  <- paths$db_path
  
  # Ensure filtered table exists
  if (isTRUE(verbose)) message("Filtering genomes before download.")
  filtered_output <- filterGenomes(base_dir = base_dir, user_bacs = user_bacs, verbose = verbose)
  if (is.null(filtered_output)) {
    message("No filtered metadata available.")
    return(character())
  }
  con <- filtered_output$duckdbConnection
  tbl <- filtered_output$table_name
  
  filtered_metadata <- tibble::as_tibble(DBI::dbReadTable(con, tbl))
  filtered_genome_ids <- filtered_metadata |>
    dplyr::distinct(`genome.genome_id`) |>
    dplyr::pull(`genome.genome_id`)
  
  # Destination
  bug_dir     <- dirname(db_path)
  genome_path <- file.path(bug_dir, "genomes")
  dir.create(genome_path, recursive = TRUE, showWarnings = FALSE)
  
  # Simple log
  logs_dir <- file.path(base_dir, "data", "logs")
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(logs_dir, "genome_downloads.log")
  
  if (isTRUE(verbose)) {
    message("Downloading genomes to: ", genome_path)
  }
  cat(sprintf("[%s] Downloading %d genomes to %s\n",
              Sys.time(), length(filtered_genome_ids), genome_path),
      file = log_file, append = TRUE)
  
  # Parallel plan (same as before)
  future::plan(future::multisession, workers = max(1, future::availableCores() - 1))
  
  downloadGenomes <- function(genomeID) {
    files <- c(".PATRIC.gff", ".fna", ".PATRIC.faa")
    success_all <- TRUE
    
    for (ext in files) {
      url  <- paste0("ftps://ftp.bv-brc.org/genomes/", genomeID, "/", genomeID, ext)
      dest <- file.path(genome_path, paste0(genomeID, ext))
      
      success  <- FALSE
      attempts <- 0
      while (!success && attempts < 3) {
        result <- system(paste("wget -q -O", shQuote(dest), shQuote(url)))
        success <- (result == 0 && file.exists(dest) && file.info(dest)$size > 100)
        attempts <- attempts + 1
      }
      
      if (!success) {
        msg <- paste("Failed to download:", url)
        message(msg)
        cat(sprintf("[%s] %s\n", Sys.time(), msg), file = log_file, append = TRUE)
        success_all <- FALSE
      }
    }
    list(genomeID = genomeID, success = success_all)
  }
  
  download_results <- future.apply::future_lapply(filtered_genome_ids, downloadGenomes)
  successful_ids <- purrr::map(download_results, "genomeID")[purrr::map_lgl(download_results, "success")]
  
  if (isTRUE(verbose)) {
    message("Successfully downloaded complete file sets for ", length(successful_ids), " genomes.")
  }
  successful_ids
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
  paths    <- .buildDBpath(base_dir = base_dir, user_bacs = user_bacs)
  db_path  <- paths$db_path
  bug_dir  <- dirname(db_path)
  
  genome_path <- file.path(bug_dir, "genomes")
  files_all   <- list.files(genome_path, full.names = TRUE)
  files_all   <- files_all[file.info(files_all)$size > 100]
  
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
      genome_id   = genomeID,
      gff_path    = if (file.exists(gff_path) && file.info(gff_path)$size > 100) gff_path else NA,
      fna_path    = if (file.exists(fna_path) && file.info(fna_path)$size > 100) fna_path else NA,
      faa_path    = if (file.exists(faa_path) && file.info(faa_path)$size > 100) faa_path else NA,
      panaroo_input = if (
        file.exists(gff_path) && file.exists(fna_path) && file.exists(faa_path) &&
        file.info(gff_path)$size > 100 && file.info(fna_path)$size > 100 && file.info(faa_path)$size > 100
      ) paste(gff_path, fna_path) else NA,
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