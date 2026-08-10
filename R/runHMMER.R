#' Validate if a HMM file has old HMMER3 format and remove them
#'
#' @param hmm_file
#'
#' @returns
#'
#' @keywords internal
#' @examples
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

#' Download and prepare HMMER databases for generating new file types.
#'
#' @param hmmer_db_dir Directory to store HMMER databases
#' @param databases List of databases to prepare (default: c("Pfam", "COG", "AMRFinder"))
#' @param docker_image Docker image containing HMMER (default: "staphb/hmmer")
#' @param hmmer_db_url If the databases contain custom database(s), the url is required to download the database.
#'
#' @returns A list of paths to the database hmm files.
#'
#' @keywords internal
#' @examples
.prepareHmmerDatabases <- function(
    hmmer_db_dir,
    databases = c("Pfam", "COG", "AMRFinder"),
    docker_image = "staphb/hmmer",
    hmmer_db_url = NULL,
    verbose = TRUE
) {

  hmmer_db_dir <- normalizePath(hmmer_db_dir)

  options(timeout = max(3600, getOption("timeout")))

  dbs <- list(

    Pfam = list(
      dir = file.path(hmmer_db_dir, "Pfam"),
      hmm_name = "Pfam-A.hmm",
      url = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
      type = "gz"
    ),

    COG = list(
      dir = file.path(hmmer_db_dir, "COG"),
      hmm_name = "COG_database2024.hmm",
      url = "http://boabio.belozersky.msu.ru/media/COG_database2024.zip",
      type = "zip"
    ),

    AMRFinder = list(
      dir = file.path(hmmer_db_dir, "AMRFinder"),
      hmm_name = NULL,
      url = "https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/latest/NCBIfam-AMRFinder.HMM.tar.gz",
      type = "tar.gz"
    )
  )

  ## Add custom database(s)
  missing_dbs <- setdiff(databases, names(dbs))

  if (length(missing_dbs) > 0) {

    if (is.null(hmmer_db_url)) {
      stop(
        "hmmer_db_url must be supplied when using custom databases"
      )
    }

    get_db_type <- function(url) {

      file <- basename(url)

      if (grepl("\\.(tar\\.gz|tgz)$",
                file,
                ignore.case = TRUE)) {

        return("tar.gz")

      } else if (grepl("\\.zip$",
                       file,
                       ignore.case = TRUE)) {

        return("zip")

      } else if (grepl("\\.gz$",
                       file,
                       ignore.case = TRUE)) {

        return("gz")

      } else {

        stop(
          "Unsupported archive type: ",
          file
        )
      }
    }

    for (db_name in missing_dbs) {

      dbs[[db_name]] <- list(
        dir = file.path(hmmer_db_dir, db_name),
        hmm_name = NULL,
        url = hmmer_db_url,
        type = get_db_type(hmmer_db_url)
      )
    }
  }

  dbs <- dbs[databases]

  db_paths <- list()

  for (db_name in names(dbs)) {

    db <- dbs[[db_name]]

    dir.create(
      db$dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    if(verbose) message("Checking ", db_name)

    hmm_files <- list.files(
      db$dir,
      pattern = "\\.hmm$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )

    if (length(hmm_files) == 0) {

      message("Downloading ", db_name)

      tmp <- tempfile()

      utils::download.file(
        url = db$url,
        destfile = tmp,
        mode = "wb",
        method = "libcurl"
      )

      switch(

        db$type,

        gz = {

          hmm_file <- file.path(
            db$dir,
            db$hmm_name %||% basename(
              sub("\\.gz$",
                  "",
                  basename(db$url),
                  ignore.case = TRUE)
            )
          )

          R.utils::gunzip(
            filename = tmp,
            destname = hmm_file,
            overwrite = TRUE,
            remove = FALSE
          )
        },

        zip = {

          utils::unzip(
            zipfile = tmp,
            exdir = db$dir
          )
        },

        `tar.gz` = {

          utils::untar(
            tarfile = tmp,
            exdir = db$dir
          )
        }
      )

      hmm_files <- list.files(
        db$dir,
        pattern = "\\.hmm$",
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
      )
    }

      hmm_files <- list.files(
      db$dir,
      pattern = "\\.hmm$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )

    if (length(hmm_files) == 0) {

      stop(
        "No .hmm file found for ",
        db_name
      )

    } else if (length(hmm_files) == 1) {

      hmm_file <- hmm_files[1]

    } else {

      hmm_file <- file.path(
        db$dir,
        paste0(db_name, ".hmm")
      )

     source_hmms <- setdiff(
  normalizePath(hmm_files),
  normalizePath(hmm_file, mustWork = FALSE)
)

     # purrr implementation
     valid_hmms <- purrr::map_lgl(source_hmms, .isValidHmmFile)

if (any(!valid_hmms)) {
  bad_files <- basename(hmm_files[!valid_hmms])

  if (isTRUE(verbose)) {
    warning(
      "Ignoring ",
      length(bad_files),
      " invalid HMM file(s):\n",
      paste(bad_files, collapse = "\n"),
      call. = FALSE
    )
  }

  hmm_files <- hmm_files[valid_hmms]
}

if (length(source_hmms) == 0) {

  stop(
    "No valid HMM files found for ",
    db_name
  )
}

      if (!file.exists(hmm_file)) {

        if(verbose) message(
          "Combining ",
          length(source_hmms),
          " HMM files for ",
          db_name
        )

        file.create(hmm_file)

        for (f in sort(source_hmms)) {
          file.append(hmm_file, f)
        }
      }
    }

    pressed_files <- paste0(
      hmm_file,
      c(".h3m", ".h3i", ".h3f", ".h3p")
    )

    if (!all(file.exists(pressed_files))) {

      if(verbose) message(
        "Running hmmpress for ",
        basename(hmm_file)
      )

      output <- system2(
        "docker",
        args = c(
          "run",
          "--rm",
          "-v",
          paste0(dirname(hmm_file), ":/db"),
          docker_image,
          "hmmpress",
          file.path("/db", basename(hmm_file))
        ),
        stdout = TRUE,
        stderr = TRUE
      )

      if (!all(file.exists(pressed_files))) {

        stop(
          "hmmpress failed for ",
          db_name,
          "\n",
          paste(output, collapse = "\n")
        )
      }
    }

    db_paths[[db_name]] <- hmm_file

    message(
      db_name,
      " ready: ",
      hmm_file
    )
  }

  db_paths
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

#' The function to run HMMER with docker
#'
#' @param JOB_NAME
#' @param FASTA
#' @param DB
#' @param Total_proteins
#' @param output_path
#' @param db_paths
#' @param docker_image
#' @param threads
#' @param n_workers
#'
#' @returns
#'
#' @keywords internal
.runHmmerJob <- function(JOB_NAME, FASTA, DB, Total_proteins,
  output_path = NULL, db_paths,
  docker_image = "staphb/hmmer", threads = 8L,
  n_workers = 8L,
  verbose = TRUE
) {
    hmmer_input <- file.path(output_path, FASTA)
    hmmer_output <- file.path(output_path, paste0(JOB_NAME, ".tbl"))

    # database paths
    database_path <- db_paths[[DB]]
    db_host_dir <- dirname(database_path)
    db_filename <- basename(database_path)
    db_cont_dir <- "/opt/hmmer/data"
    db_cont_path <- file.path(db_cont_dir, db_filename)

    # mounts
    mount_host <- output_path
    mount_cont <- "/work"

    threads_per_job <- max(
  1L,
  floor(threads / n_workers)
)

    cmd_args <- c(
      "run", "--rm",
      "-v", paste0(mount_host, ":", mount_cont),
      "-v", paste0(db_host_dir, ":", db_cont_dir),
      docker_image,
      "hmmsearch",
      "--notextw",
      "--cpu", as.character(threads_per_job),
      "-Z", Total_proteins,
      "--domZ", Total_proteins,
      "--domtblout", .to_container(hmmer_output, mount_host, mount_cont),
      db_cont_path,
      .to_container(hmmer_input, mount_host, mount_cont)
    )

    if(verbose) message("Running hmmsearch via Docker...")
    output <- tryCatch(
      {
        system2("docker", args = cmd_args, stdout = TRUE, stderr = TRUE)
      },
      error = function(e) {
        stop("hmmsearch execution failed: ", e$message)
      }
    )

    if (!file.exists(hmmer_output)) {
      stop("hmmsearch failed: output file not found. Check stderr:\n", paste(output, collapse = "\n"))
    }

    if(verbose) message("hmmsearch completed successfully.")

    hmmer_tbl <- .parseHMMEROutput(hmmer_output) |>
      dplyr::select("protein", "query_name")

    hmmer_tbl_filename <- file.path(
      dirname(hmmer_output),
      paste0(tools::file_path_sans_ext(basename(hmmer_output)), ".parquet")
    )

    .write_compressed_parquet(hmmer_tbl, hmmer_tbl_filename)

    hmmer_tbl_filename
  }


#' Wrapper for preparing HMM databases and running HMMER on protein sequences from duckdb and writing them.
#'
#' @param duckdb_path
#' @param output_path
#' @param threads
#' @param hmmer_db_dir
#' @param databases
#' @param docker_image
#' @param num_of_splits
#' @param n_workers
#'
#' @returns
#'
#' @keywords internal
#' @examples
.runHMMER <- function(duckdb_path,
                      output_path,
                      threads = 8L,
                      hmmer_db_dir,
                      databases = c("Pfam", "COG", "AMRFinder"),
                      docker_image = "staphb/hmmer",
                      num_of_splits = 8L,
                      n_workers = 8L,
                      verbose = TRUE
                    ) {
  # Fail fast if Docker is missing
  if (!nzchar(Sys.which("docker"))) {
    stop("Docker is not available on your PATH but is required to run HMMER.")
  }

  # But also check if Docker is on the PATH but isn't running
  docker_ok <- system2(
    "docker",
    "info",
    stdout = FALSE,
    stderr = FALSE
  ) == 0L

  if (!docker_ok) {
    stop(
      "Docker is installed but is not running or cannot be reached. ",
      "Please (re)start Docker Desktop and try again."
    )
  }

  duckdb_path <- .docker_path(duckdb_path)
  if (missing(output_path) || output_path %in% c(".", "results", "results/")) {
    output_path <- dirname(duckdb_path)
  }
  output_path <- .docker_path(output_path)
  if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  prot_seqs <- DBI::dbReadTable(con, "protein_cluster_seq") |>
    tibble::as_tibble()

  # required to define the database size for hmmsearch --Z and --domZ parameters
  Total_proteins <- nrow(prot_seqs)

  if(is.null(hmmer_db_dir)) {
  hmmer_db_dir <- output_path
  }

  # database paths
  if(verbose) message ("Preparing HMM databases")
 db_paths <- .prepareHmmerDatabases(
  hmmer_db_dir = hmmer_db_dir,
  databases = databases,
  docker_image = docker_image,
  verbose = verbose
)

db_paths <- db_paths[databases]

  # clamp splits to the number of sequences available
  chunk_count <- min(as.integer(num_of_splits), nrow(prot_seqs))

  split_fasta <- function(seqs, prefix) {
    records <- paste0(">", seqs$name, "\n", seqs$sequence)
    chunk_size <- ceiling(length(records) / chunk_count)
    chunks <- split(records, ceiling(seq_along(records) / chunk_size))

    purrr::walk2(chunks, seq_along(chunks), function(chunk, i) {
      chunk_path <- file.path(output_path, sprintf("%s_chunk_%02d.fasta", prefix, i))
      readr::write_lines(chunk, chunk_path)
    })
  }

  split_fasta(prot_seqs, "protein")

  job_list <- expand.grid(
    chunk = sprintf("%02d", seq_len(chunk_count)),
    db = databases,
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(
      JOB_NAME = paste0("protein_chunk_", chunk, "_", db),
      FASTA = paste0("protein_chunk_", chunk, ".fasta"),
      DB = db
    ) |>
    dplyr::select(JOB_NAME, FASTA, DB)

  future::plan(
  future::multisession,
  workers = max(1L, n_workers)
)

  if (verbose) message("Running HMMER jobs")
parquet_files <- furrr::future_map_chr(
  seq_len(nrow(job_list)),
  function(i) {

    .runHmmerJob(
      JOB_NAME = job_list$JOB_NAME[i],
      FASTA = job_list$FASTA[i],
      DB = job_list$DB[i],
      Total_proteins = Total_proteins,
      output_path = output_path,
      db_paths = db_paths,
      docker_image = docker_image,
      threads = threads,
      n_workers = n_workers,
      verbose = verbose
    )
  }
)

future::plan(future::sequential)

  parquet_tbl <- tibble::tibble(
  parquet = parquet_files,
  db = job_list$DB
)

 final_parquets <- list()

for (database_name in databases) {

  if(verbose) message("Combining ", database_name)

  db_files <- parquet_tbl |>
    dplyr::filter(
      db == database_name
    ) |>
    dplyr::pull(parquet)

  combined_tbl <- purrr::map(
    db_files,
    arrow::read_parquet
  ) |>
    dplyr::bind_rows() |>
dplyr::left_join(.parse_hmmer_profiles(db_paths[[database_name]]) |>
  dplyr::select(query_name = profile_name, query_accession = profile_accession, description = profile_description),
by = "query_name")

  final_parquet <- file.path(
    output_path,
    paste0(
      "protein_",
      database_name,
      ".parquet"
    )
  )

  .write_compressed_parquet(
    combined_tbl,
    final_parquet
  )

  DBI::dbWriteTable(
    con,
    name = paste0(
      "protein_",
      database_name
    ),
    value = combined_tbl,
    overwrite = TRUE
  )

  final_parquets[[database_name]] <- final_parquet

  message(
    "Created ",
    basename(final_parquet)
  )
}

invisible(final_parquets)

  # purrr::map(parquet_files, arrow::read_parquet) |>
  #   dplyr::bind_rows() |>
  #   .write_compressed_parquet(final_parquet)

  # message("Combined parquet written.")

  # arrow::read_parquet(final_parquet) |>
  #   DBI::dbWriteTable(conn = con, name = tools::file_path_sans_ext(basename(final_parquet)), overwrite = TRUE)
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

#' Map HMMER protein annotations to genome-level count matrix and load into DuckDB
#'
#' Reads a Parquet file of HMMER hits (produced by [.runHMMER()]), joins the
#' annotations to the protein-cluster count matrix already in DuckDB, aggregates
#' counts per genome and annotation, and writes the result both as a Parquet file
#' and as a new table in the DuckDB database.
#'
#' @param annotated_parquet Path to the combined HMMER results Parquet file
#'   (e.g. `"results/Ecoli/protein_COG.parquet"`). The filename stem is used as
#'   the table name in DuckDB.
#' @param duckdb_path Path to the per-selection DuckDB database containing a
#'   `protein_count` table (created by [CDHIT2duckdb()]).
#'
#' @return Invisibly returns the path to the written count Parquet file.
#'
#' @seealso [CDHIT2duckdb()], [runDataProcessing()]
#'
#' @examples
#' \dontrun{
#' proteinAnnotations2Duckdb(
#'   annotated_parquet = "results/Ecoli/protein_COG.parquet",
#'   duckdb_path       = "data/Ecoli/Eco.duckdb"
#' )
#' }
#'
#' @keywords internal
.proteinAnnotations2Duckdb <- function(
    duckdb_path,
    databases = c("Pfam", "COG", "AMRFinder")
) {

  duckdb_path <- .docker_path(duckdb_path)

  con <- DBI::dbConnect(
    duckdb::duckdb(),
    duckdb_path
  )

  on.exit(
    try(
      DBI::dbDisconnect(
        con,
        shutdown = FALSE
      ),
      silent = TRUE
    ),
    add = TRUE
  )

  protein_long <- DBI::dbReadTable(
    con,
    "protein_count"
  ) |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(
      cols = -genome_id,
      names_to = "protein",
      values_to = "count"
    ) |>
    dplyr::filter(count > 0) |>
    dplyr::mutate(
      protein = stringr::str_replace(
        protein,
        "^fig\\.",
        "fig|"
      )
    )

  count_paths <- list()

  for (database in databases) {

    annotation_table <- paste0(
      "protein_",
      database
    )

    if (!DBI::dbExistsTable(con, annotation_table)) {

      warning(
        annotation_table,
        " not found in DuckDB. Skipping."
      )

      next
    }

    message(
      "Processing ",
      annotation_table
    )

    annotation <- DBI::dbReadTable(
      con,
      annotation_table
    ) |>
      tibble::as_tibble()

    genome_annot_matrix <- protein_long |>
      dplyr::inner_join(
        annotation |>
        dplyr::select(
          protein,
          query_name
        ),
        by = "protein"
      ) |>
      dplyr::group_by(
        genome_id,
        query_name
      ) |>
      dplyr::summarise(
        count = sum(count),
        .groups = "drop"
      ) |>
      tidyr::pivot_wider(
        names_from = query_name,
        values_from = count,
        values_fill = 0
      )

    count_table <- paste0(
      annotation_table,
      "_count"
    )

    count_path <- file.path(
      dirname(duckdb_path),
      paste0(
        count_table,
        ".parquet"
      )
    )

    arrow::write_parquet(
      genome_annot_matrix,
      count_path
    )

    DBI::dbWriteTable(
      con,
      count_table,
      genome_annot_matrix,
      overwrite = TRUE
    )

    count_paths[[database]] <- count_path

    message(
      "Created ",
      count_table
    )
  }

  invisible(count_paths)
}

#' Annotate proteins using DefenseFinder + CasFinder HMMs
#' Will add to the duckdb + create the parquet file.
#'
#' @param defense_db_dir Directory used to store downloaded HMMs
#' @param docker_image Docker image containing HMMER
#' @param duckdb_path DuckDB database path
#' @param output_path Output directory
#' @param threads Number of HMMER threads
#'
#' @returns Path to annotation parquet
#' @keywords internal
.defenseHMMER <- function(
    defense_db_dir,
    docker_image = "staphb/hmmer",
    duckdb_path = "inst/extdata/Sfl.duckdb",
    output_path = NULL,
    threads = 8L,
    verbose = TRUE
) {

  if (!nzchar(Sys.which("docker"))) {
    stop("Docker is required.")
  }

  defense_db_dir <- normalizePath(
    defense_db_dir,
    mustWork = FALSE
  )

  if (is.null(output_path)) {
    output_path <- dirname(
      normalizePath(
        duckdb_path,
        mustWork = FALSE
      )
    )
  }

  dir.create(
    defense_db_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    output_path,
    recursive = TRUE,
    showWarnings = FALSE
  )

  ####################################################################
  # download repositories
  ####################################################################

  defense_dir <- file.path(
    defense_db_dir,
    "DefenseFinder"
  )

  cas_dir <- file.path(
    defense_db_dir,
    "CasFinder"
  )

  if (!dir.exists(defense_dir)) {

    if(verbose) message(
      "Downloading DefenseFinder models"
    )

    tmp <- tempfile(fileext = ".zip")

    utils::download.file(
      "https://github.com/mdmparis/defense-finder-models/archive/refs/heads/master.zip",
      tmp,
      mode = "wb",
      method = "libcurl"
    )

    utils::unzip(
      tmp,
      exdir = defense_dir
    )

    unlink(tmp)
  }

  if (!dir.exists(cas_dir)) {

    if(verbose) message(
      "Downloading CasFinder models"
    )

    tmp <- tempfile(fileext = ".zip")

    utils::download.file(
      "https://github.com/macsy-models/CasFinder/archive/refs/heads/main.zip",
      tmp,
      mode = "wb",
      method = "libcurl"
    )

    utils::unzip(
      tmp,
      exdir = cas_dir
    )

    unlink(tmp)
  }

  ####################################################################
  # helper
  ####################################################################

  build_database <- function(
      repo_dir,
      db_name
  ) {

    profile_dirs <- list.dirs(
      repo_dir,
      recursive = TRUE,
      full.names = TRUE
    )

    profile_dirs <- profile_dirs[
      basename(profile_dirs) == "profiles"
    ]

    # moving to purrr implementation
    hmm_files <- profile_dirs |>
      purrr::map(\(x) list.files(x,
                                 pattern = "\\.hmm$",
                                 recursive = TRUE,
                                 full.names = TRUE,
                                 ignore.case = TRUE)) |>
      purrr::flatten_chr() |>
      unique()

    if (length(hmm_files) == 0) {

      stop(
        "No HMM files found for ",
        db_name
      )
    }

    valid_hmms <- purrr::map_lgl(hmm_files, .isValidHmmFile)

    if (any(!valid_hmms)) {
      bad_files <- basename(hmm_files[!valid_hmms])

      if (isTRUE(verbose)) {
        warning(
          "Ignoring ",
          length(bad_files),
          " invalid HMM file(s):\n",
          paste(bad_files, collapse = "\n"),
          call. = FALSE
        )
      }

      hmm_files <- hmm_files[valid_hmms]
    }

if (length(hmm_files) == 0) {
  stop(
    "No valid HMM files found for ",
    db_name
  )
}

    combined_hmm <- file.path(
      repo_dir,
      paste0(
        db_name,
        ".hmm"
      )
    )

    if (file.exists(combined_hmm)) {
      unlink(combined_hmm)
    }

    file.create(combined_hmm)

    for (f in sort(hmm_files)) {

      file.append(
        combined_hmm,
        f
      )
    }

    pressed_files <- paste0(
      combined_hmm,
      c(
        ".h3m",
        ".h3i",
        ".h3f",
        ".h3p"
      )
    )

    if (!all(file.exists(pressed_files))) {

      if(verbose) message(
        "Running hmmpress for ",
        db_name
      )

      output <- system2(
        "docker",
        args = c(
          "run",
          "--rm",
          "-v",
          paste0(
            dirname(combined_hmm),
            ":/db"
          ),
          docker_image,
          "hmmpress",
          file.path(
            "/db",
            basename(combined_hmm)
          )
        ),
        stdout = TRUE,
        stderr = TRUE
      )

      if (!all(file.exists(pressed_files))) {

        stop(
          "hmmpress failed for ",
          db_name,
          "\n",
          paste(output,
                collapse = "\n")
        )
      }
    }

    combined_hmm
  }

  ####################################################################
  # build separate databases
  ####################################################################

  defense_hmm <- build_database(
    defense_dir,
    "DefenseFinder"
  )

  cas_hmm <- build_database(
    cas_dir,
    "CasFinder"
  )

  ####################################################################
  # load proteins
  ####################################################################

  con <- DBI::dbConnect(
    duckdb::duckdb(),
    duckdb_path
  )

  on.exit(
    try(
      DBI::dbDisconnect(
        con,
        shutdown = FALSE
      ),
      silent = TRUE
    ),
    add = TRUE
  )

  prot_seqs <- DBI::dbReadTable(
    con,
    "protein_cluster_seq"
  ) |>
    tibble::as_tibble()

  fasta_file <- file.path(
    output_path,
    "protein_DefenseCas.faa"
  )

   # required to define the database size for hmmsearch --Z and --domZ parameters
  Total_proteins <- nrow(prot_seqs)

  readr::write_lines(
    paste0(
      ">",
      prot_seqs$name,
      "\n",
      prot_seqs$sequence
    ),
    fasta_file
  )

  ####################################################################
  # run hmmsearch separately
  ####################################################################

  databases <- list(
    DefenseFinder = defense_hmm,
    CasFinder = cas_hmm
  )

  all_hits <- list()

  for (db_name in names(databases)) {

    if(verbose) message(
      "Running ",
      db_name
    )

    hmm_file <- databases[[db_name]]

    tbl_file <- file.path(
      output_path,
      paste0(
        "protein_",
        db_name,
        ".tbl"
      )
    )

    output <- system2(
      "docker",
      args = c(
        "run",
        "--rm",
        "-v",
        paste0(output_path, ":/work"),
        "-v",
        paste0(dirname(hmm_file), ":/db"),
        docker_image,
        "hmmsearch",
         "--notextw",
        "--cpu",
        as.character(threads),
         "-Z", Total_proteins,
        "--domZ", Total_proteins,
        "--domtblout",
        file.path(
          "/work",
          basename(tbl_file)
        ),
        file.path(
          "/db",
          basename(hmm_file)
        ),
        "/work/protein_DefenseCas.faa"
      ),
      stdout = TRUE,
      stderr = TRUE
    )

   if (!file.exists(tbl_file)) {

  warning(
    "hmmsearch failed for ",
    db_name,
    "\n",
    paste(output, collapse = "\n")
  )

  next
}

    hits <- .parseHMMEROutput(
      tbl_file
    ) |>
      dplyr::select(
        protein,
        query_name
      ) |>
      dplyr::mutate(
        database = db_name
      )|>
dplyr::left_join(.parse_hmmer_profiles(hmm_file) |>
  dplyr::select(query_name = profile_name, query_accession = profile_accession, description = profile_description),
by = "query_name")

    all_hits[[db_name]] <- hits
  }

  ####################################################################
  # merge at parquet stage
  ####################################################################

  combined_tbl <- dplyr::bind_rows(
    all_hits
  )

  parquet_file <- file.path(
    output_path,
    "protein_DefenseCas.parquet"
  )

  .write_compressed_parquet(
    combined_tbl,
    parquet_file
  )

  DBI::dbWriteTable(
    con,
    "protein_DefenseCas",
    combined_tbl,
    overwrite = TRUE
  )

  message(
    "Created protein_DefenseCas"
  )

  invisible(parquet_file)
}
