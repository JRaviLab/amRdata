# hmmscan --cpu 32 --tblout "${BUG}_genes_COG.tbl" "$COG_DB" "${BUG}_translated_gene_seqs.fasta"
write_compressed_parquet <- function(df, path) {
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
                      threads = 0,
                      database_path,
                      split_jobs = TRUE,
                      num_of_splits = 20) {
  # Fail fast if Docker is missing
  #  if (!nzchar(Sys.which("docker"))) {
  #    stop("Docker is not available on your PATH but is required to run CD-HIT.")
  #  }

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

  # robust database name
  database <- tools::file_path_sans_ext(basename(database_path))

  # Chunking function
  chunk_count <- num_of_splits

  split_fasta <- function(seqs, prefix) {
    records <- paste0("> ", seqs$name, "\n", seqs$sequence)
    chunk_size <- ceiling(length(records) / chunk_count)
    chunks <- split(records, ceiling(seq_along(records) / chunk_size))

    purrr::walk2(chunks, seq_along(chunks), function(chunk, i) {
      chunk_path <- file.path(output_path, sprintf("%s_chunk_%02d.fasta", prefix, i))
      readr::write_lines(chunk, chunk_path)
    })
  }

  split_fasta(prot_seqs, "protein")

  #  readr::write_lines(paste0("> ", prot_seqs$name, "\n", prot_seqs$sequence),
  #            file.path(output_path, "proteins_for_hmmer.fasta"))

  # Generate job list for ARG and COG
  job_list <- expand.grid(
    chunk = sprintf("%02d", 1:chunk_count),
    db = database,
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(
      JOB_NAME = paste0("protein_chunk_", chunk, "_", db),
      FASTA = paste0("protein_chunk_", chunk, ".fasta"),
      DB = db
    ) |>
    dplyr::select(JOB_NAME, FASTA, DB)

  # readr::write_tsv(job_list, file.path(output_path, paste0("hmmer_jobs_",database,".txt")))

  # number of parallel jobs (NOT threads per hmmscan)
  n_workers <- 4

  # threads per hmmscan
  threads <- 8

  future::plan(future::multisession, workers = n_workers)

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
      "exec",
      "-B", paste0(mount_host, ":", mount_cont),
      "-B", paste0(db_host_dir, ":", db_cont_dir),
      "/scratch/alpine/aghosh5@xsede.org/software/hmmer_latest.sif",
      "hmmscan",
      "--cpu", as.character(threads),
      "--tblout", .to_container(hmmer_output, mount_host, mount_cont),
      db_cont_path,
      .to_container(hmmer_input, mount_host, mount_cont)
    )

    message("Running hmmer via Docker...")
    output <- tryCatch(
      {
        system2("apptainer", args = cmd_args, stdout = TRUE, stderr = TRUE)
      },
      error = function(e) {
        stop("hmmer execution failed: ", e$message)
      }
    )

    if (!file.exists(hmmer_output)) {
      stop("hmmer failed: output file not found. Check stderr:\n", paste(output, collapse = "\n"))
    }

    message("hmmer completed successfully.")

    hmmer_tbl <- .parse_hmmer_output(hmmer_output)

    hmmer_tbl <- hmmer_tbl |>
      dplyr::select("name", "query_name", "description")

    hmmer_tbl_filename <- file.path(dirname(hmmer_output), paste0(tools::file_path_sans_ext(basename(hmmer_output)), ".parquet"))

    hmmer_tbl |>
      write_compressed_parquet(hmmer_tbl_filename)

    return(hmmer_tbl_filename)
  }

  results <- furrr::future_pwalk(
    job_list,
    .runHmmerJob,
    .progress = TRUE
  )
  parquet_files <- results |>
    dplyr::mutate(parquet = paste0(JOB_NAME, ".parquet")) |>
    dplyr::pull(parquet)

  final_parquet <- file.path(output_path, paste0("protein_", database, ".parquet"))

  dataset <- arrow::open_dataset(
    file.path(output_path, parquet_files),
    format = "parquet"
  )

  arrow::write_parquet(
    dataset,
    file.path(output_path, paste0("protein_", database, ".parquet"))
  )

  message("Combined parquet written")

  # arrow::read_parquet("/scratch/alpine/aghosh5@xsede.org/AMR/data/Campylobacter_jejuni/protein_COG_count.parquet") |> DBI::dbWriteTable(conn=con, name="protein_COG_count")

  arrow::read_parquet(final_parquet) |>
    DBI::dbWriteTable(conn = con, name = tools::file_path_sans_ext(basename(final_parquet)), overwrite = TRUE)
}


#' Read a file created as HMMER output
#' modified the rhmmer
#' @param file Filename
#' @return data.frame
#' @export
#'
.parse_hmmer_output <- function(file) {
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
    best_bis = readr::col_double(),
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

  # drop comment lines
  data_lines <- lines[!grepl("^#", lines)]

  # split: whitespace-separated fields
  split_fields <- strsplit(data_lines, "\\s+", perl = TRUE)

  # count space separated fields
  N <- max(sapply(split_fields, length))

  # the line delimiter should always be just "\n", even on Windows
  lines <- readr::read_lines(file, lazy = FALSE, progress = FALSE)

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

proteinAnnotations2Duckdb <- function(annotated_parquet, duckdb_path) {
  annotated_parquet <- .docker_path(annotated_parquet)

  # robust database name
  database <- tools::file_path_sans_ext(basename(annotated_parquet))

  duckdb_path <- .docker_path(duckdb_path)

  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  protein_count <- DBI::dbReadTable(con, "protein_count") |>
    tibble::as_tibble()

  protein_long <- protein_count |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(
      cols = -genome_id,
      names_to = "query_name",
      values_to = "count"
    ) |>
    dplyr::filter(count > 0)

  # protein_COG <- arrow::read_parquet("data/Campylobacter_jejuni/protein_COG.parquet")
  Annotation <- arrow::read_parquet(annotated_parquet)

  protein_long_annot <- protein_long |>
    dplyr::mutate(query_name = stringr::str_replace(query_name, "^fig\\.", "fig|")) |>
    dplyr::inner_join(
      Annotation |>
        dplyr::select(name, query_name),
      by = "query_name"
    )
  genome_annot_counts <- protein_long_annot |>
    dplyr::group_by(genome_id, name) |>
    dplyr::summarise(count = sum(count), .groups = "drop")

  genome_annot_matrix <- genome_annot_counts |>
    tidyr::pivot_wider(
      names_from = name,
      values_from = count,
      values_fill = 0
    )

  arrow::write_parquet(
    genome_annot_matrix,
    file.path(dirname(duckdb_path), paste0(database, "_count", ".parquet"))
  )

  count_path <- file.path(dirname(duckdb_path), paste0(database, "_count", ".parquet"))

  arrow::read_parquet(count_path) |>
    DBI::dbWriteTable(conn = con, name = tools::file_path_sans_ext(basename(count_path)), overwrite = TRUE)
}
