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

.runHMMER <- function(duckdb_path,
                      output_path,
                      threads = 8L,
                      database_path,
                      docker_image = "staphb/hmmer",
                      split_jobs = TRUE,
                      num_of_splits = 20L,
                      n_workers = 4L) {
  # Fail fast if Docker is missing
  if (!nzchar(Sys.which("docker"))) {
    stop("Docker is not available on your PATH but is required to run HMMER.")
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

  # derive a clean label from the database filename
  database <- tools::file_path_sans_ext(basename(database_path))

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
    db = database,
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(
      JOB_NAME = paste0("protein_chunk_", chunk, "_", db),
      FASTA = paste0("protein_chunk_", chunk, ".fasta"),
      DB = db
    ) |>
    dplyr::select(JOB_NAME, FASTA, DB)

  .runHmmerJob <- function(JOB_NAME, FASTA, DB) {
    hmmer_input <- file.path(output_path, FASTA)
    hmmer_output <- file.path(output_path, paste0(JOB_NAME, ".tbl"))

    # database paths
    db_host_dir <- dirname(database_path)
    db_filename <- basename(database_path)
    db_cont_dir <- "/opt/hmmer/data"
    db_cont_path <- file.path(db_cont_dir, db_filename)

    # mounts
    mount_host <- output_path
    mount_cont <- "/work"

    cmd_args <- c(
      "run", "--rm",
      "-v", paste0(mount_host, ":", mount_cont),
      "-v", paste0(db_host_dir, ":", db_cont_dir),
      docker_image,
      "hmmscan",
      "--cpu", as.character(threads),
      "--tblout", .to_container(hmmer_output, mount_host, mount_cont),
      db_cont_path,
      .to_container(hmmer_input, mount_host, mount_cont)
    )

    message("Running hmmscan via Docker...")
    output <- tryCatch(
      {
        system2("docker", args = cmd_args, stdout = TRUE, stderr = TRUE)
      },
      error = function(e) {
        stop("hmmscan execution failed: ", e$message)
      }
    )

    if (!file.exists(hmmer_output)) {
      stop("hmmscan failed: output file not found. Check stderr:\n", paste(output, collapse = "\n"))
    }

    message("hmmscan completed successfully.")

    hmmer_tbl <- .parseHMMEROutput(hmmer_output) |>
      dplyr::select("name", "query_name", "description")

    hmmer_tbl_filename <- file.path(
      dirname(hmmer_output),
      paste0(tools::file_path_sans_ext(basename(hmmer_output)), ".parquet")
    )

    .write_compressed_parquet(hmmer_tbl, hmmer_tbl_filename)

    hmmer_tbl_filename
  }

  hmmer_param <- BiocParallel::SnowParam(workers = max(1L, n_workers))
  parquet_files <- BiocParallel::bpmapply(
    FUN = .runHmmerJob,
    JOB_NAME = job_list$JOB_NAME,
    FASTA = job_list$FASTA,
    DB = job_list$DB,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE,
    BPPARAM = hmmer_param
  )

  final_parquet <- file.path(output_path, paste0("protein_", database, ".parquet"))

  purrr::map(parquet_files, arrow::read_parquet) |>
    dplyr::bind_rows() |>
    .write_compressed_parquet(final_parquet)

  message("Combined parquet written.")

  arrow::read_parquet(final_parquet) |>
    DBI::dbWriteTable(conn = con, name = tools::file_path_sans_ext(basename(final_parquet)), overwrite = TRUE)
}


#' Parse HMMER tabular output into a tibble
#'
#' Reads a HMMER `--tblout` file and returns a tidy tibble with one row per
#' target-query hit. Comment lines are stripped and the free-text description
#' field is reunited from the remaining whitespace-delimited columns.
#'
#' @param file Path to a HMMER `.tbl` output file produced with `--tblout`.
#'
#' @return A tibble with 19 columns matching the HMMER per-sequence hit table:
#'   `name`, `accession`, `query_name`, `query_accession`, `sequence_evalue`,
#'   `sequence_score`, `sequence_bias`, `best_evalue`, `best_score`,
#'   `best_bias`, `number_exp`, `number_reg`, `number_clu`, `number_ov`,
#'   `number_env`, `number_dom`, `number_rep`, `number_inc`, `description`.
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
  col_types <- readr::cols(
    name = readr::col_character(),
    accession = readr::col_character(),
    query_name = readr::col_character(),
    query_accession = readr::col_character(),
    sequence_evalue = readr::col_double(),
    sequence_score = readr::col_double(),
    sequence_bias = readr::col_double(),
    best_evalue = readr::col_double(),
    best_score = readr::col_double(),
    best_bias = readr::col_double(),
    number_exp = readr::col_double(),
    number_reg = readr::col_integer(),
    number_clu = readr::col_integer(),
    number_ov = readr::col_integer(),
    number_env = readr::col_integer(),
    number_dom = readr::col_integer(),
    number_rep = readr::col_integer(),
    number_inc = readr::col_character(),
    description = readr::col_character()
  )
  # the line delimiter should always be just "\n", even on Windows
  lines <- readr::read_lines(file, lazy = FALSE, progress = FALSE)

  # drop comment lines
  data_lines <- lines[!grepl("^#", lines)]

  # split: whitespace-separated fields
  split_fields <- strsplit(data_lines, "\\s+", perl = TRUE)

  # count space separated fields
  N <- max(sapply(split_fields, length))

  table <- sub(
    pattern = sprintf("(%s).*", paste0(rep("\\S+", N), collapse = " +")),
    replacement = "\\1",
    x = lines,
    perl = TRUE
  ) |>
    gsub(pattern = "  *", replacement = "\t") |>
    paste0(collapse = "\n") |>
    readr::read_tsv(
      col_names = names(col_types$cols),
      comment = "#",
      na = "-",
      col_types = col_types,
      lazy = FALSE,
      progress = FALSE
    ) |>
    tidyr::unite(description, description:last_col(), sep = " ")
  table$description <- gsub("\t", " ", table$description)

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
#' @export
proteinAnnotations2Duckdb <- function(annotated_parquet, duckdb_path) {
  annotated_parquet <- .docker_path(annotated_parquet)
  duckdb_path <- .docker_path(duckdb_path)

  # derive table name from the annotation filename stem
  database <- tools::file_path_sans_ext(basename(annotated_parquet))

  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  protein_long <- DBI::dbReadTable(con, "protein_count") |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(
      cols = -genome_id,
      names_to = "query_name",
      values_to = "count"
    ) |>
    dplyr::filter(count > 0)

  annotation <- arrow::read_parquet(annotated_parquet)

  genome_annot_matrix <- protein_long |>
    # protein IDs are stored with "." separator in DuckDB but "|" in HMMER output
    dplyr::mutate(query_name = stringr::str_replace(query_name, "^fig\\.", "fig|")) |>
    dplyr::inner_join(
      dplyr::select(annotation, name, query_name),
      by = "query_name"
    ) |>
    dplyr::group_by(genome_id, name) |>
    dplyr::summarise(count = sum(count), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = name,
      values_from = count,
      values_fill = 0
    )

  count_path <- file.path(dirname(duckdb_path), paste0(database, "_count.parquet"))
  arrow::write_parquet(genome_annot_matrix, count_path)

  DBI::dbWriteTable(
    conn = con,
    name = tools::file_path_sans_ext(basename(count_path)),
    value = genome_annot_matrix,
    overwrite = TRUE
  )

  invisible(count_path)
}
