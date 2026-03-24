#' Normalize a host filesystem path for use in Docker
#'
#' Converts Windows and mixed-separator paths to forward slashes
#' and applies `normalizePath()` without requiring the path to exist.
#'
#' @param p Character scalar. A filesystem path on the host OS.
#'
#' @return A normalized path string.
#'
#' @keywords internal
.docker_path <- function(p) gsub("\\\\", "/", normalizePath(p, mustWork = FALSE))

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

# Launch Panaroo to build a pangenome (per batch)
#' processPanaroo()
#'
#' See Panaroo's documentation for details on how the parameters affect your
#' pangenome output: https://gthlab.au/panaroo/#/gettingstarted/params
#'
#' @param batch_input A series of genome IDs for input
#' @param output_path Character scalar. Base directory for Panaroo outputs and temporary files.
#' @param core_threshold Numeric. Core genome threshold for Panaroo (`--core_threshold`). Default `0.90`.
#' @param len_dif_percent Numeric. Length difference percentage (`--len_dif_percent`). Default `0.95`.
#' @param cluster_threshold Numeric. Sequence identity threshold (`--threshold`). Default `0.95`.
#' @param family_seq_identity Numeric. Gene family clustering identity (`-f`). Default `0.5`.
#' @param panaroo_threads_per_job Integer. Number of threads for Panaroo and parallel execution.
#'
#' @returns A list of results for each Panaroo batch in its output directory.
#'
#' @keywords internal
#' @examples NULL
.processPanaroo <- function(batch_input,
                            output_path,
                            core_threshold,
                            len_dif_percent,
                            cluster_threshold,
                            family_seq_identity,
                            panaroo_threads_per_job) {
  output_path <- .docker_path(output_path)
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

  # Fail fast if Docker is missing
  if (!nzchar(Sys.which("docker"))) {
    stop("Docker is not available on your PATH but is required to run Panaroo.")
  }

  # Host mount root = bug directory
  mount_host <- output_path
  mount_cont <- "/work"

  # Write the genome list file (convert each "gff fna" to container-visible paths)
  genome_filepath_host <- tempfile(pattern = "genomeFilepath_", fileext = ".txt", tmpdir = output_path)

  batch_input_cont <- vapply(unlist(batch_input), function(line) {
    parts <- strsplit(line, " +")[[1]]
    parts_cont <- .to_container(parts, host_root = mount_host, container_root = mount_cont)
    paste(parts_cont, collapse = " ")
  }, character(1), USE.NAMES = FALSE)

  # Write with Unix line endings to avoid issues inside Linux container
  con <- file(genome_filepath_host, open = "wb")
  writeLines(batch_input_cont, con = con, sep = "\n", useBytes = TRUE)
  close(con)

  # Create unique output dir by timestamping it
  output_dir_host <- file.path(output_path, paste0("panaroo_out_", format(Sys.time(), "%Y%m%d%H%M%OS4")))
  dir.create(output_dir_host, recursive = TRUE, showWarnings = FALSE)

  # Convert to container-visible paths
  genome_filepath_cont <- .to_container(genome_filepath_host, host_root = mount_host, container_root = mount_cont)
  output_dir_cont <- .to_container(output_dir_host, host_root = mount_host, container_root = mount_cont)

  # Run Panaroo in Docker
  cmd_args <- c(
    "run",
    "--platform", "linux/amd64",
    "--rm",
    "-v", paste0(mount_host, ":", mount_cont),
    "-w", mount_cont,
    "staphb/panaroo:1.5.1",
    "panaroo",
    "-i", genome_filepath_cont,
    "-o", output_dir_cont,
    "--clean-mode", "strict",
    "--merge_paralogs",
    "--remove-invalid-genes",
    "--core_threshold", as.character(core_threshold),
    "--len_dif_percent", as.character(len_dif_percent),
    "--threshold", as.character(cluster_threshold),
    "-f", as.character(family_seq_identity),
    "-t", as.character(panaroo_threads_per_job)
  )

  res <- tryCatch(
    {
      system2("docker", args = cmd_args, stdout = TRUE, stderr = TRUE)
    },
    error = function(e) e
  )

  if (inherits(res, "error")) {
    stop(sprintf("Docker/Panaroo failed to launch: %s", res$message))
  }

  # If Panaroo wrote an error but system2 didn't throw, scan output for clues
  if (length(res) && any(grepl("Traceback|Error|No such file|not found|failed", res, ignore.case = TRUE))) {
    message("Panaroo output:\n", paste(res, collapse = "\n"))
  }

  invisible(res)
}


#' Temporarily set a future plan for parallel execution
#'
#' Sets a `future` plan (sequential or multisession) for the duration of a block
#' and automatically restores the previous plan on exit.
#'
#' @param workers Integer. Number of workers to use; if <= 1, uses sequential mode.
#' @param plan Character. Either `"multisession"` or `"sequential"`.
#'
#' @return Invisibly returns `TRUE` after setting the plan.
#'
#' @keywords internal
.with_future_plan <- function(workers, plan = c("multisession", "sequential")) {
  plan <- match.arg(plan)
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (is.null(workers) || workers <= 1L || identical(plan, "sequential")) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession, workers = workers)
  }
  invisible(TRUE)
}


#' Run Panaroo for Pangenome Analysis in Parallel Batches
#'
#' Executes Panaroo inside a Docker container on genome annotation
#' files prepared by [genomeList()]. The function can optionally split input genomes
#' into batches, runs Panaroo with strict cleaning and clustering options, and
#' returns the results of each batch execution.
#'
#' @param duckdb_path A path to the DuckDB database containing the `"files"` table.
#' @param output_path Character scalar. Base directory for Panaroo outputs and temporary files.
#' @param core_threshold Numeric. Core genome threshold for Panaroo (`--core_threshold`). Default `0.90`.
#' @param len_dif_percent Numeric. Length difference percentage (`--len_dif_percent`). Default `0.95`.
#' @param cluster_threshold Numeric. Sequence identity threshold (`--threshold`). Default `0.95`.
#' @param family_seq_identity Numeric. Gene family clustering identity (`-f`). Default `0.5`.
#' @param threads Integer. Number of threads for Panaroo and parallel execution. Default `8`.
#' @param split_jobs Logical. If TRUE, split into multiple smaller pangenome
#'   generation jobs that can be merged by [.mergePanaroo()]. If FALSE, all isolates in one run.
#'
#' @return A list of results for each Panaroo batch in its output directory.
#'
#' @keywords internal
#' @details
#' - Panaroo uses: `--clean-mode strict`, `--merge_paralogs`, `--remove-invalid-genes`.
#' - Temporary genome file lists are created in `output_path`.
#' - Output directories are named `panaroo_out_<timestamp>` under `output_path`.
#'
.runPanaroo <- function(duckdb_path = "data/{Bug}/{Bug}.duckdb",
                        output_path = "data/{Bug}/",
                        core_threshold = 0.90,
                        len_dif_percent = 0.95,
                        cluster_threshold = 0.95,
                        family_seq_identity = 0.5,
                        threads = 8,
                        split_jobs = FALSE) {
  duckdb_path <- normalizePath(duckdb_path)
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  if (missing(output_path) || output_path %in% c(".", "results", "results/")) {
    output_path <- dirname(duckdb_path)
  }
  output_path <- normalizePath(output_path)

  genome_query_output <- DBI::dbReadTable(con, "files")

  panaroo_input_files <- genome_query_output |>
    dplyr::pull(panaroo_input)

  # Drop true NAs
  panaroo_input_files <- panaroo_input_files[!is.na(panaroo_input_files)]

  split_files <- strsplit(panaroo_input_files, " ")

  # Plan for filtering
  .with_future_plan(workers = threads)
  valid_entries <- furrr::future_map(split_files, function(paths) {
    gff_file <- paths[1]
    if (file.exists(gff_file)) {
      length(readLines(gff_file, n = 5, warn = FALSE)) >= 5
    } else {
      FALSE
    }
  })

  filtered_panaroo_input <- sapply(split_files[unlist(valid_entries)], paste, collapse = " ")

  total_lines <- length(filtered_panaroo_input)
  batch_size <- if (isTRUE(split_jobs)) ceiling(total_lines / 5) else total_lines
  panaroo_batches <- split(filtered_panaroo_input, ceiling(seq_along(filtered_panaroo_input) / batch_size))

  n_jobs <- length(panaroo_batches)
  if (n_jobs == 0L) {
    warning("Panaroo inputs do not exist after filtering. Check your upstream processing.")
    return(invisible(list()))
  }

  # Ensure sum of per-job CPUs does not exceed `threads`
  panaroo_threads_per_job <- max(1L, floor(threads / n_jobs))

  # One worker per batch
  .with_future_plan(workers = n_jobs)
  batch_panaroo_run <- furrr::future_map(
    panaroo_batches,
    ~ .processPanaroo(
      batch_input             = .x,
      output_path             = output_path,
      core_threshold          = core_threshold,
      len_dif_percent         = len_dif_percent,
      cluster_threshold       = cluster_threshold,
      family_seq_identity     = family_seq_identity,
      panaroo_threads_per_job = panaroo_threads_per_job
    ),
    .options = furrr::furrr_options(seed = TRUE)
  )

  invisible(batch_panaroo_run)
}

#' Merge multiple Panaroo batch outputs into a single pangenome result
#'
#' Finds batch output directories under `input_path` that contain `final_graph.gml`,
#' and merges them with `panaroo-merge` inside a Docker container. Output goes to
#' `input_path/merge_output`.
#'
#' @param input_path A directory that contains multiple Panaroo pangenome directories for merging.
#' @param core_threshold Numeric. Core genome threshold for Panaroo (`--core_threshold`). Default `0.90`.
#' @param len_dif_percent Numeric. Length difference percentage (`--len_dif_percent`). Default `0.95`.
#' @param cluster_threshold Numeric. Sequence identity threshold (`--threshold`). Default `0.95`.
#' @param family_seq_identity Numeric. Gene family clustering identity (`-f`). Default `0.5`.
#' @param threads Integer. Number of threads for Panaroo and parallel execution. Default `8`.
#'
#' @returns A a single combined pangenome.
#'
#' @keywords internal
.mergePanaroo <- function(input_path,
                          core_threshold = 0.90,
                          len_dif_percent = 0.95,
                          cluster_threshold = 0.95,
                          family_seq_identity = 0.5,
                          threads = 8) {
  input_path <- .docker_path(input_path)

  # Fail fast if Docker is missing
  if (!nzchar(Sys.which("docker"))) {
    stop("Docker is not available on your PATH but is required to run panaroo-merge.")
  }

  merge_dir <- file.path(input_path, "merge_output")
  dir.create(merge_dir, recursive = TRUE, showWarnings = FALSE)

  all_dirs <- list.dirs(input_path, recursive = FALSE, full.names = TRUE)
  all_dirs <- all_dirs[grepl("^panaroo_out_", basename(all_dirs))]

  valid_dirs <- all_dirs[file.exists(file.path(all_dirs, "final_graph.gml"))]

  if (length(valid_dirs) > 1) {
    mount_host <- input_path
    mount_cont <- "/work"

    # Provide each dir as a separate argv token after "-d"
    dir_args <- as.vector(t(.to_container(valid_dirs, host_root = mount_host, container_root = mount_cont)))

    cmd_args <- c(
      "run",
      "--platform", "linux/amd64",
      "--rm",
      "-v", paste0(mount_host, ":", mount_cont),
      "-w", mount_cont,
      "staphb/panaroo:1.5.1",
      "panaroo-merge",
      "-d", dir_args,
      "-o", file.path(mount_cont, "merge_output"),
      "--merge_paralogs",
      "--core_threshold", as.character(core_threshold),
      "--len_dif_percent", as.character(len_dif_percent),
      "--threshold", as.character(cluster_threshold),
      "-f", as.character(family_seq_identity),
      "-t", as.character(threads)
    )

    system2("docker", args = cmd_args, stdout = TRUE, stderr = TRUE)
  } else {
    stop("No valid Panaroo batch directories found (need >= 2 with final_graph.gml).")
  }
}


#' Load Panaroo gene presence/absence table into DuckDB
#'
#' Reads `gene_presence_absence.csv` and constructs a genome-by-gene count
#' table, writing it into the DuckDB database as `gene_count`.
#'
#' @param panaroo_output_path Path to a Panaroo result directory.
#' @param duckdb_path Path to a DuckDB database file.
#'
#' @return A tibble containing the gene count matrix.
#'
#' @keywords internal
.panaroo2geneTable <- function(panaroo_output_path, duckdb_path) {
  filepath <- file.path(normalizePath(panaroo_output_path), "gene_presence_absence.csv")
  duckdb_path <- normalizePath(duckdb_path)
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  gene_count <- read.table(filepath, sep = ",", header = TRUE, fill = TRUE, quote = "") |>
    tibble::as_tibble() |>
    dplyr::select(-c(Non.unique.Gene.name, Annotation)) |>
    tidyr::pivot_longer(cols = -1) |>
    tidyr::pivot_wider(names_from = Gene, values_from = value) |>
    dplyr::rename("genome_id" = "name") |>
    dplyr::mutate(genome_id = stringr::str_replace_all(genome_id, c("^X" = "", "\\.PATRIC$" = ""))) |>
    dplyr::mutate(across(-genome_id, ~ ifelse(. == "", 0, stringr::str_count(., ";") + 1)))

  DBI::dbWriteTable(con, "gene_count", gene_count, overwrite = TRUE)
  gene_count
}


#' Extract gene names and annotations from Panaroo outputs
#'
#' Reads Panaroo's `gene_presence_absence.csv` to extract gene identifiers
#' and gene annotations, then writes them into the DuckDB table `gene_names`.
#'
#' @inheritParams .panaroo2geneTable
#'
#' @return A tibble with `Gene` and `Annotation` columns.
#'
#' @keywords internal
.panaroo2geneNames <- function(panaroo_output_path, duckdb_path) {
  filepath <- file.path(normalizePath(panaroo_output_path), "gene_presence_absence.csv")
  duckdb_path <- normalizePath(duckdb_path)
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  gene_names <- read.table(filepath, sep = ",", header = TRUE, fill = TRUE, quote = "") |>
    tibble::as_tibble() |>
    dplyr::select(c(Gene, Annotation))

  DBI::dbWriteTable(con, "gene_names", gene_names, overwrite = TRUE)
  gene_names
}


#' Create structural variant presence/absence table from Panaroo outputs
#'
#' Reads `struct_presence_absence.Rtab` and constructs a genome-by-struct
#' presence/absence matrix, writing the result to `gene_struct` in DuckDB.
#'
#' @inheritParams .panaroo2geneTable
#'
#' @return A tibble containing the struct matrix.
#'
#' @keywords internal
.panaroo2StructTable <- function(panaroo_output_path, duckdb_path) {
  struct_filepath <- file.path(normalizePath(panaroo_output_path), "struct_presence_absence.Rtab")
  duckdb_path <- normalizePath(duckdb_path)
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  gene_struct <- read.table(struct_filepath, sep = "\t", header = TRUE, fill = TRUE, quote = "") |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(cols = -1) |>
    tidyr::pivot_wider(names_from = Gene, values_from = value) |>
    dplyr::rename("genome_id" = "name") |>
    dplyr::mutate(genome_id = stringr::str_replace_all(genome_id, c("^X" = "", "\\.PATRIC$" = "")))

  DBI::dbWriteTable(con, "gene_struct", gene_struct, overwrite = TRUE)
  gene_struct
}


#' Import additional Panaroo reference outputs into DuckDB
#'
#' Loads reference sequences and long-format gene–protein mappings from
#' Panaroo outputs and stores them into DuckDB (`gene_ref_seq`, `genome_gene_protein`).
#'
#' @inheritParams .panaroo2geneTable
#'
#' @return Invisibly returns TRUE.
#'
#' @keywords internal
.panaroo2OtherTables <- function(panaroo_output_path, duckdb_path) {
  panaroo_output_path <- normalizePath(panaroo_output_path)
  duckdb_path <- normalizePath(duckdb_path)
  fasta_filepath <- file.path(panaroo_output_path, "pan_genome_reference.fa")
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  gene_fasta <- Biostrings::readDNAStringSet(filepath = fasta_filepath)
  DBI::dbWriteTable(con, "gene_ref_seq",
    tibble::tibble(
      name = names(gene_fasta),
      sequence = as.character(gene_fasta)
    ),
    overwrite = TRUE
  )

  readr::read_csv(file.path(panaroo_output_path, "gene_presence_absence.csv")) |>
    dplyr::select(-`Non-unique Gene name`) |>
    tidyr::pivot_longer(-c("Gene", "Annotation"),
      names_to = "genome_ids",
      values_to = "protein_ids"
    ) |>
    dplyr::mutate(genome_ids = gsub(".PATRIC", "", genome_ids)) |>
    dplyr::select(genome_ids, Gene, protein_ids) |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(protein_ids)) |>
    tidyr::separate_rows(protein_ids, sep = ";") |>
    dplyr::filter(!stringr::str_detect(protein_ids, "_pseudo")) |>
    DBI::dbWriteTable(conn = con, name = "genome_gene_protein", overwrite = TRUE)
}


#' Import all Panaroo-derived outputs into DuckDB
#'
#' Wrapper that loads gene counts, gene names, struct tables, and reference
#' sequence tables from a Panaroo output directory into a DuckDB database.
#'
#' @inheritParams .panaroo2geneTable
#'
#' @return Invisibly returns TRUE.
#'
#' @keywords internal
.panaroo2duckdb <- function(panaroo_output_path, duckdb_path) {
  panaroo_output_path <- normalizePath(panaroo_output_path)
  duckdb_path <- normalizePath(duckdb_path)

  .panaroo2geneTable(panaroo_output_path, duckdb_path)
  .panaroo2geneNames(panaroo_output_path, duckdb_path)
  .panaroo2StructTable(panaroo_output_path, duckdb_path)
  .panaroo2OtherTables(panaroo_output_path, duckdb_path)
  invisible(TRUE)
}


#' Run CD-HIT inside Docker and assemble protein clusters
#'
#' Concatenates `.faa` files, executes CD-HIT in a Docker container,
#' and returns paths to the cluster output files.
#'
#' @param duckdb_path Path to DuckDB containing the `files` table.
#' @param output_path Directory to write concatenated FASTA and CD-HIT results.
#' @param output_prefix String used to prefix CD-HIT output files.
#' @param identity CD-HIT sequence identity threshold (`-c`).
#' @param word_length CD-HIT word size (`-n`).
#' @param threads Integer number of threads.
#' @param memory Integer memory limit (`-M`).
#' @param extra_args Character vector of additional CD-HIT arguments.
#'
#' @return A list containing paths to the concatenated FASTA and cluster FASTA.
#'
#' @keywords internal
.runCDHIT <- function(duckdb_path,
                      output_path,
                      output_prefix = "cdhit_out",
                      identity = 0.9,
                      word_length = 5,
                      threads = 0,
                      memory = 0,
                      extra_args = c("-g", "1")) {
  # Fail fast if Docker is missing
  if (!nzchar(Sys.which("docker"))) {
    stop("Docker is not available on your PATH but is required to run CD-HIT.")
  }

  duckdb_path <- .docker_path(duckdb_path)
  if (missing(output_path) || output_path %in% c(".", "results", "results/")) {
    output_path <- dirname(duckdb_path)
  }
  output_path <- .docker_path(output_path)
  if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  genome_query_output <- DBI::dbReadTable(con, "files")

  cdhit_input_files <- genome_query_output |>
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~ . != "NA")) |>
    dplyr::pull(faa_path)

  if (length(cdhit_input_files) == 0 || !all(file.exists(cdhit_input_files))) {
    stop("Some or all .faa files do not exist.")
  }

  cdhit_input_faa <- file.path(output_path, paste0(output_prefix, "_input.fa"))
  file_conn <- file(cdhit_input_faa, "w")
  for (file in cdhit_input_files) {
    cat(readLines(file), file = file_conn, sep = "\n")
  }
  close(file_conn)

  clustered_faa <- file.path(output_path, paste0(output_prefix, "_proteins"))

  mount_host <- output_path
  mount_cont <- "/work"

  cmd_args <- c(
    "run", "--rm",
    "--platform", "linux/amd64",
    "-v", paste0(mount_host, ":", mount_cont),
    "-w", mount_cont,
    "weizhongli1987/cdhit:4.8.1",
    "cd-hit",
    "-i", .to_container(cdhit_input_faa, mount_host, mount_cont),
    "-o", .to_container(clustered_faa, mount_host, mount_cont),
    "-c", as.character(identity),
    "-n", as.character(word_length),
    "-T", as.character(threads),
    "-M", as.character(memory),
    "-d", "0",
    extra_args
  )

  message("Running cd-hit via Docker...")
  output <- tryCatch(
    {
      system2("docker", args = cmd_args, stdout = TRUE, stderr = TRUE)
    },
    error = function(e) {
      stop("cd-hit execution failed: ", e$message)
    }
  )

  if (!file.exists(clustered_faa)) {
    stop("cd-hit failed: output file not found. Check stderr:\n", paste(output, collapse = "\n"))
  }
  # Ensure .clstr exists (used downstream)
  if (!file.exists(paste0(clustered_faa, ".clstr"))) {
    stop(
      "cd-hit did not produce the expected .clstr file at: ", paste0(clustered_faa, ".clstr"),
      "\nFull output:\n", paste(output, collapse = "\n")
    )
  }

  message("cd-hit completed successfully.")
  list(
    cdhit_input_faa = cdhit_input_faa,
    clustered_faa   = clustered_faa
  )
}

#' Run Panaroo and import pangenome outputs into DuckDB
#'
#' @description
#' `runPanaroo2Duckdb()` executes Panaroo on the genomes registered in a
#' per-selection DuckDB (created earlier by `prepareGenomes()`), optionally in
#' multiple batches, and imports all resulting pangenome tables into the same
#' DuckDB database.
#'
#' It acts as a high-level wrapper around:
#' * **`.runPanaroo()`** — runs Panaroo (single or multi-batch)
#' * **`.mergePanaroo()`** — optionally merges batch outputs
#' * **`.panaroo2duckdb()`** — loads Panaroo results (gene counts, struct variants,
#'   gene names, reference sequences, long tables) into the DuckDB
#'
#' The function determines which Panaroo output directory to use (single-run or merged),
#' verifies that a valid pangenome has been produced, and updates the DuckDB with
#' standardized table names consistent with downstream processing steps.
#'
#' @param duckdb_path Character. Path to the per-selection DuckDB database created by
#'   `prepareGenomes()`. Must contain a `files` table with Panaroo input file paths.
#' @param output_path Character or `NULL`. Directory where Panaroo outputs
#'   (`panaroo_out_*` or merged `merge_output/`) will be written. If `NULL`,
#'   defaults to `dirname(duckdb_path)`.
#'
#' @param core_threshold Numeric. Panaroo `--core_threshold` parameter.
#'   Default: `0.90`.
#' @param len_dif_percent Numeric. Panaroo `--len_dif_percent` parameter.
#'   Default: `0.95`.
#' @param cluster_threshold Numeric. Panaroo global clustering `--threshold`.
#'   Default: `0.95`.
#' @param family_seq_identity Numeric. Panaroo gene family identity `-f`.
#'   Default: `0.5`.
#'
#' @param threads Integer. Total CPU budget to allocate for Panaroo.
#'   If `split_jobs = TRUE`, threads are divided across batches.
#'   Default: `16`.
#'
#' @param split_jobs Logical. If `TRUE`, Panaroo is run in multiple parallel
#'   batches (up to 5, depending on dataset size), and batch outputs are merged
#'   using `.mergePanaroo()`. If `FALSE`, only one Panaroo invocation is run.
#'   Default: `FALSE`.
#'
#' @param verbose Logical. Print status messages during Panaroo execution,
#'   merging, and DuckDB import. Default: `TRUE`.
#'
#' @return
#' Invisibly returns the path to the selected Panaroo output directory
#' (either the single-run output or the merged `merge_output/` directory).
#'
#' @details
#' ### Panaroo Output Discovery
#' After running `.runPanaroo()`, the function scans `output_path` for directories
#' matching `panaroo_out_*` and identifies those containing a `final_graph.gml` file —
#' the minimum requirement for a valid Panaroo run.
#'
#' * If **`split_jobs = TRUE`** and multiple valid outputs are present,
#'   `.mergePanaroo()` is used to combine the outputs.
#' * If **`split_jobs = FALSE`**, the single valid output directory is used directly.
#'
#' ### DuckDB Integration
#' `.panaroo2duckdb()` is then called to import:
#' * gene presence/absence counts (`gene_count`)
#' * gene names (`gene_names`)
#' * structural presence/absence (`gene_struct`)
#' * gene reference FASTA (`gene_ref_seq`)
#' * long-form genome → gene → protein tables
#'
#' These maintain the standardized schema used by downstream feature extraction
#' and modeling steps in `amRdata` and `amRml`.
#'
#' @seealso
#' * `.runPanaroo()` — core Panaroo execution
#' * `.mergePanaroo()` — merge multiple Panaroo batches
#' * `.panaroo2duckdb()` — import Panaroo results into DuckDB
#' * [runDataProcessing()] — full pipeline including CD-HIT & InterProScan
#'
#' @examples
#' \dontrun{
#' # Basic usage:
#' runPanaroo2Duckdb(
#'   duckdb_path = "data/Shigella_flexneri/Sfl.duckdb",
#'   output_path = "data/Shigella_flexneri",
#'   threads     = 12,
#'   split_jobs  = FALSE
#' )
#'
#' # Merging multi-batch pangenomes:
#' runPanaroo2Duckdb(
#'   duckdb_path = "data/Ecoli/Eco.duckdb",
#'   output_path = "data/Ecoli",
#'   split_jobs  = TRUE,
#'   threads     = 24
#' )
#' }
#'
#' @export
runPanaroo2Duckdb <- function(duckdb_path,
                              output_path = NULL,
                              core_threshold = 0.90,
                              len_dif_percent = 0.95,
                              cluster_threshold = 0.95,
                              family_seq_identity = 0.5,
                              threads = 16,
                              split_jobs = FALSE,
                              verbose = TRUE) {
  duckdb_path <- normalizePath(duckdb_path)
  out_dir <- if (is.null(output_path)) dirname(duckdb_path) else normalizePath(output_path)

  if (isTRUE(verbose)) message("Launching Panaroo.")
  .runPanaroo(
    duckdb_path = duckdb_path,
    output_path = out_dir,
    core_threshold = core_threshold,
    len_dif_percent = len_dif_percent,
    cluster_threshold = cluster_threshold,
    family_seq_identity = family_seq_identity,
    threads = threads,
    split_jobs = split_jobs
  )

  # Identify Panaroo outputs that contain a final_graph.gml file
  pan_outs <- list.dirs(out_dir, recursive = FALSE, full.names = TRUE)
  pan_outs <- pan_outs[grepl("^panaroo_out_", basename(pan_outs))]
  valid <- pan_outs[file.exists(file.path(pan_outs, "final_graph.gml"))]

  if (length(valid) == 0L) {
    stop("No valid Panaroo outputs found (no final_graph.gml). Check logs.")
  }

  # If split jobs produced 2+ valid outputs, merge them; else use the single output dir
  target_dir <- NULL
  if (isTRUE(split_jobs) && length(valid) >= 2L) {
    if (isTRUE(verbose)) message("Merging Panaroo batch outputs.")
    .mergePanaroo(
      input_path          = out_dir,
      core_threshold      = core_threshold,
      len_dif_percent     = len_dif_percent,
      cluster_threshold   = cluster_threshold,
      family_seq_identity = family_seq_identity,
      threads             = max(1L, floor(threads / 2))
    )
    target_dir <- file.path(out_dir, "merge_output")
    if (!file.exists(file.path(target_dir, "gene_presence_absence.csv"))) {
      stop("Expected merged Panaroo outputs in merge_output/, but files were not found.")
    }
  } else {
    target_dir <- valid[[1]]
  }

  if (isTRUE(verbose)) message("Writing Panaroo tables to DuckDB.")
  .panaroo2duckdb(panaroo_output_path = target_dir, duckdb_path = duckdb_path)

  invisible(target_dir)
}


#' Parse CD-HIT `.clstr` output into a long-format mapping
#'
#' Reads a CD-HIT `.clstr` file and constructs a mapping of clusters to genome IDs.
#'
#' @param clustered_faa Base path to CD-HIT output (without `.clstr` extension).
#'
#' @return A data.table with columns `cluster` and `genome_id`.
#'
#' @keywords internal
.parseProteinClusters <- function(clustered_faa) {
  clstr <- paste0(clustered_faa, ".clstr")
  if (!file.exists(clstr)) {
    stop(
      "CD-HIT cluster file not found: ", clstr,
      "\nEnsure .runCDHIT() completed successfully and produced the .clstr file."
    )
  }

  lines <- data.table::fread(clstr, sep = "\n", header = FALSE)$V1
  cluster_ids <- grep("^>Cluster", lines)
  cluster_map <- data.table::data.table()

  for (i in seq_along(cluster_ids)) {
    start <- cluster_ids[i] + 1
    end <- if (i < length(cluster_ids)) cluster_ids[i + 1] - 1 else length(lines)
    cluster_lines <- lines[start:end]

    # This finds the reference cluster ID and names the cluster with it
    ref_line <- grep("\\*$", cluster_lines, value = TRUE)
    ref_id <- if (length(ref_line) > 0) {
      stringr::str_extract(ref_line, "fig\\|[0-9]+\\.[0-9]+\\.peg\\.[0-9]+")
    } else {
      paste0("Cluster_", i - 1)
    }

    # Pull genome IDs
    genome_matches <- stringr::str_match(
      cluster_lines,
      "fig\\|([0-9]+\\.[0-9]+)\\.peg\\.[0-9]+"
    )[, 2]
    genome_matches <- genome_matches[!is.na(genome_matches)]

    if (length(genome_matches) > 0) {
      cluster_map <- data.table::rbindlist(list(
        cluster_map,
        data.table::data.table(cluster = ref_id, genome_id = genome_matches)
      ), use.names = TRUE)
    }
  }

  cluster_map
}


#' Build genome-by-protein-cluster count matrix
#'
#' Converts a long-format cluster mapping from `.parseProteinClusters()`
#' into a genome-by-cluster count matrix.
#'
#' @param cluster_map A data.table with `cluster` and `genome_id`.
#'
#' @return A wide-format matrix as a data.frame.
#'
#' @keywords internal
.buildProtMatrices <- function(cluster_map) {
  cluster_map[, count := 1]
  reshape2::dcast(cluster_map, genome_id ~ cluster, value.var = "count", fun.aggregate = sum, fill = 0)
}
# Back-compat wrapper (older external name)
buildMatrices <- function(cluster_map) .buildProtMatrices(cluster_map)


#' Extract per-cluster protein names from CD-HIT cluster FASTA
#'
#' Reads a FASTA file of representative proteins and extracts protein IDs,
#' locus tags, and descriptive names.
#'
#' @param cluster_map Output of `.parseProteinClusters()`.
#' @param cluster_fasta Path to representative FASTA file used by CD-HIT.
#'
#' @return A tibble containing protein metadata.
#'
#' @keywords internal
.clusterNames <- function(cluster_map, cluster_fasta) {
  # Note: cluster_map_unique computed but not used previously—keeping for parity
  cluster_map_unique <- cluster_map |>
    tibble::as_tibble() |>
    dplyr::distinct() |>
    dplyr::group_by(cluster) |>
    dplyr::slice_head(n = 1)

  cdhit_output_faa <- Biostrings::readAAStringSet(cluster_fasta)

  names_faa <- names(cdhit_output_faa) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      proteinID = stringr::str_extract(value, "^fig\\|[0-9]+\\.[0-9]+\\.peg\\.[0-9]+"),
      locus_tag = stringr::str_match(value, "peg\\.[0-9]+\\|([^\\s]+)")[, 2],
      proteinName = stringr::str_trim(stringr::str_match(value, "\\|[^\\s]+\\s+(.*?)\\s+\\[")[, 2])
    ) |>
    dplyr::select(-value)

  names_faa
}

#' Cluster proteins with CD-HIT and write results to DuckDB
#' @export
CDHIT2duckdb <- function(duckdb_path,
                         output_path,
                         output_prefix = "cdhit_out",
                         identity = 0.9,
                         word_length = 5,
                         threads = 0,
                         memory = 0,
                         extra_args = c("-g", "1")) {
  duckdb_path <- normalizePath(duckdb_path)
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  if (missing(output_path) || output_path %in% c(".", "results", "results/")) {
    output_path <- dirname(duckdb_path) # e.g., ./results/<bug>
  }
  output_path <- normalizePath(output_path)

  cdhit_outputs <- .runCDHIT(duckdb_path,
    output_path,
    output_prefix = output_prefix,
    identity = identity,
    word_length = word_length,
    threads = threads,
    memory = memory,
    extra_args = extra_args
  )

  cluster_map <- .parseProteinClusters(cdhit_outputs$clustered_faa)
  cluster_count <- .buildProtMatrices(cluster_map)

  DBI::dbWriteTable(con, "protein_count", cluster_count, overwrite = TRUE)

  cluster_fasta <- cdhit_outputs$cdhit_input_faa
  cluster_name <- .clusterNames(cluster_map, cluster_fasta)
  DBI::dbWriteTable(con, "protein_names", cluster_name, overwrite = TRUE)

  clustered_faa <- Biostrings::readAAStringSet(cdhit_outputs$clustered_faa)
  DBI::dbWriteTable(con, "protein_cluster_seq",
    tibble::tibble(
      name     = names(clustered_faa) |> stringr::str_extract("fig\\|[0-9]+\\.[0-9]+\\.peg\\.[0-9]+"),
      sequence = as.character(clustered_faa)
    ),
    overwrite = TRUE
  )
  invisible(TRUE)
}


#' Check or install InterProScan data bundle
#'
#' Ensures that the InterProScan data directory exists locally, downloading
#' and verifying the appropriate tarball when necessary.
#'
#' @param version InterProScan version string.
#' @param dest_dir Directory where data should be installed.
#' @param docker_image Docker image string for InterProScan.
#' @param platform Character indicating Docker platform (e.g. `"linux/amd64"`).
#' @param curl_bin Path to curl executable.
#' @param verbose Logical; print status messages.
#'
#' @return A list containing `data_dir` and `ready` status.
#'
#' @keywords internal
.checkInterProData <- function(
  version = "5.76-107.0",
  dest_dir = "inst/extdata/interpro",
  docker_image = sprintf("interpro/interproscan:%s", version),
  platform = "linux/amd64",
  curl_bin = "curl",
  verbose = TRUE
) {
  msg <- function(...) if (verbose) message(sprintf(...))

  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  dest_dir <- normalizePath(dest_dir, mustWork = TRUE)

  root_dir <- file.path(dest_dir, sprintf("interproscan-%s", version))
  data_dir <- file.path(root_dir, "data")

  # Simple existence check
  if (dir.exists(data_dir) && length(list.files(data_dir, recursive = TRUE)) > 0) {
    msg("InterProScan data already present at: %s", data_dir)
    return(list(data_dir = normalizePath(data_dir), ready = TRUE))
  }

  # Download bundle if needed
  tar_url <- sprintf(
    "http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/%s/alt/interproscan-data-%s.tar.gz",
    version, version
  )
  md5_url <- paste0(tar_url, ".md5")
  tar_path <- file.path(dest_dir, basename(tar_url))
  md5_path <- paste0(tar_path, ".md5")

  if (!file.exists(tar_path)) {
    msg("Downloading InterProScan data bundle.")
    status_tar <- system2(curl_bin, c("-L", "-o", tar_path, tar_url))
    status_md5 <- system2(curl_bin, c("-L", "-o", md5_path, md5_url))
    if (status_tar != 0 || status_md5 != 0) {
      stop("Failed to download InterProScan data bundle.")
    }
  }

  msg("Verifying MD5 checksum.")
  md5_expected <- sub("\\s+.*$", "", readLines(md5_path)[1])
  md5_actual <- tools::md5sum(tar_path)[[1]]
  if (!identical(tolower(md5_expected), tolower(md5_actual))) {
    stop("MD5 checksum mismatch for InterProScan data bundle.")
  }

  msg("Extracting InterProScan data bundle.")
  utils::untar(tar_path, exdir = dest_dir, tar = "internal")

  msg("Data unpacked successfully.")
  return(list(data_dir = normalizePath(data_dir), ready = TRUE))
}


#' Internal helpers for reading InterProScan TSV outputs
#'
#' Provide standardized column names, types, and a reader wrapper for the
#' InterProScan tab-delimited output format.
#'
#' @param filepath Path to a `.tsv` or `.tsv.gz` InterProScan result file.
#'
#' @return A tibble of parsed InterProScan output.
#'
#' @keywords internal
.getDfIPRColNames <- function() {
  c(
    "AccNum", "SeqMD5Digest", "SLength", "Analysis",
    "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
    "Status", "RunDate", "IPRAcc", "IPRDesc", "placeholder"
  )
}

#' Internal helpers for reading InterProScan TSV outputs
#'
#' Provide standardized column names, types, and a reader wrapper for the
#' InterProScan tab-delimited output format.
#'
#' @param filepath Path to a `.tsv` or `.tsv.gz` InterProScan result file.
#'
#' @return A tibble of parsed InterProScan output.
#'
#' @keywords internal

.getDfIPRColTypes <- function() {
  readr::cols(
    "AccNum"        = readr::col_character(),
    "SeqMD5Digest"  = readr::col_character(),
    "SLength"       = readr::col_integer(),
    "Analysis"      = readr::col_character(),
    "DB.ID"         = readr::col_character(),
    "SignDesc"      = readr::col_character(),
    "StartLoc"      = readr::col_integer(),
    "StopLoc"       = readr::col_integer(),
    "Score"         = readr::col_double(),
    "Status"        = readr::col_character(),
    "RunDate"       = readr::col_character(),
    "IPRAcc"        = readr::col_character(),
    "IPRDesc"       = readr::col_character(),
    "placeholder"   = readr::col_character()
  )
}

#' Internal helpers for reading InterProScan TSV outputs
#'
#' Provide standardized column names, types, and a reader wrapper for the
#' InterProScan tab-delimited output format.
#'
#' @param filepath Path to a `.tsv` or `.tsv.gz` InterProScan result file.
#'
#' @return A tibble of parsed InterProScan output.
#'
#' @keywords internal
.readIPRscanTsv <- function(filepath) {
  readr::read_tsv(filepath,
    col_types = .getDfIPRColTypes(),
    col_names = .getDfIPRColNames()
  )
}


#' Run InterProScan on a sequence chunk inside Docker
#'
#' Executes InterProScan on a subset of protein sequences, writing temporary
#' FASTA and reading back `.tsv` or `.tsv.gz` results.
#'
#' @param chunk A tibble with columns `name` and `sequence`.
#' @param path Working directory used for temporary files.
#' @param ipr_data_path Path to InterProScan data directory.
#' @param out_file_base Output prefix for chunk results.
#' @param appl Character vector of InterProScan applications (e.g. `"Pfam"`).
#' @param chunk_id Integer chunk index.
#' @param threads Number of CPUs for InterProScan container.
#' @param file_format Output format (`"TSV"`).
#' @param docker_image InterProScan Docker image.
#'
#' @return Path to a `.tsv` or `.tsv.gz` InterProScan output file.
#'
#' @keywords internal
.process_chunk <- function(chunk,
                           path,
                           ipr_data_path = "inst/extdata/interpro/data",
                           out_file_base,
                           appl,
                           chunk_id,
                           threads,
                           file_format,
                           docker_image = sprintf("interpro/interproscan:%s", "5.76-107.0")) {
  # Normalize and mount paths
  path <- .docker_path(path)
  bind_data <- .docker_path(ipr_data_path)

  dir.create(file.path(path, "tmp", "iprscan"), recursive = TRUE, showWarnings = FALSE)

  fasta_sequences <- Biostrings::AAStringSet(chunk$sequence)
  names(fasta_sequences) <- chunk$name
  temp_fasta_file <- tempfile(tmpdir = path, fileext = ".fa")
  Biostrings::writeXStringSet(fasta_sequences, temp_fasta_file)

  chunk_out_file_base_host <- file.path(path, sprintf("%s_chunk_%d", out_file_base, chunk_id))
  chunk_out_file_base_cont <- .to_container(chunk_out_file_base_host, path, "/work")

  # Pull image (best-effort)
  try(suppressWarnings(system2("docker", args = c("pull", docker_image))), silent = TRUE)

  appl_str <- paste(appl, collapse = ",")

  cmd_args <- c(
    "run", "--rm",
    "-v", paste0(path, ":", "/work"),
    "-v", paste0(bind_data, ":/opt/interproscan/data"),
    "-w", "/work",
    docker_image,
    "--input", .to_container(temp_fasta_file, path, "/work"),
    "--cpu", as.character(threads),
    "-f", file_format,
    "--appl", appl_str,
    "-b", chunk_out_file_base_cont
  )


  status <- tryCatch(
    {
      system2(
        "docker",
        args = c(
          "run",
          "--rm",
          "--platform", "linux/amd64", # force amd64 for ARM hosts
          "-v", paste0(path, ":", "/work"),
          "-v", paste0(bind_data, ":/opt/interproscan/data"),
          "-w", "/work",
          docker_image,
          "--input", .to_container(temp_fasta_file, path, "/work"),
          "--cpu", as.character(threads),
          "-f", file_format,
          "--appl", appl_str,
          "-b", chunk_out_file_base_cont
        ),
        stdout = TRUE,
        stderr = TRUE
      )
    },
    error = function(e) {
      stop(sprintf("InterProScan execution failed for chunk %d: %s", chunk_id, e$message))
    }
  )

  out_tsv <- paste0(chunk_out_file_base_host, ".tsv")
  out_tsvgz <- paste0(chunk_out_file_base_host, ".tsv.gz")

  if (file.exists(out_tsv)) {
    return(out_tsv)
  } else if (file.exists(out_tsvgz)) {
    return(out_tsvgz)
  } else {
    stop(sprintf(
      "InterProScan produced no output for chunk %d. Checked: %s and %s.\nLast message:\n%s",
      chunk_id, out_tsv, out_tsvgz, paste(status, collapse = "\n")
    ))
  }
}

#' Derive protein domain presence/absence and counts via InterProScan and write to DuckDB
domainFromIPR <- function(duckdb_path,
                          path,
                          out_file_base = "iprscan",
                          appl = c("Pfam"),
                          ipr_version = "5.76-107.0",
                          ipr_dest_dir = "inst/extdata/interpro",
                          ipr_platform = "linux/amd64",
                          auto_prepare_data = TRUE,
                          threads = 8,
                          file_format = "TSV",
                          docker_repo = "interpro/interproscan") {
  duckdb_path <- normalizePath(duckdb_path)
  if (missing(path) || path %in% c(".", "results", "results/")) {
    path <- dirname(duckdb_path)
  }
  path <- normalizePath(path)

  ipr_image <- sprintf("%s:%s", docker_repo, ipr_version)

  # Prepare data if needed
  ipr_info <- if (isTRUE(auto_prepare_data)) {
    .checkInterProData(
      version      = ipr_version,
      dest_dir     = ipr_dest_dir,
      docker_image = ipr_image,
      platform     = ipr_platform,
      verbose      = TRUE
    )
  } else {
    list(
      data_dir = file.path(ipr_dest_dir, sprintf("interproscan-%s", ipr_version), "data"),
      ready = NA
    )
  }
  ipr_data_path <- ipr_info$data_dir

  # Pull image once
  try(suppressWarnings(system2("docker", args = c("pull", ipr_image))), silent = TRUE)

  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  sequences_df <- dplyr::tbl(con, "protein_cluster_seq") |> tibble::as_tibble()
  if (nrow(sequences_df) == 0L) {
    stop("No sequences found in 'protein_cluster_seq'. Please run CDHIT2duckdb() first.")
  }

  # Chunking for parallel (not currently implemented due to memory limits)
  chunks <- list(sequences_df) # Force 1 chunk for RAM limits

  # Forcing 1 container operation for RAM limits
  workers <- 1
  cpu_per_container <- threads

  message(sprintf(
    "InterPro: running in single-container mode with %d CPU(s).",
    cpu_per_container
  ))

  .with_future_plan(workers = 1)

  results <- future.apply::future_lapply(seq_along(chunks), function(i) {
    res <- try(
      .process_chunk(
        chunk         = chunks[[i]],
        path          = path,
        ipr_data_path = ipr_data_path,
        out_file_base = out_file_base,
        appl          = appl,
        chunk_id      = i,
        threads       = cpu_per_container,
        file_format   = file_format,
        docker_image  = ipr_image
      ),
      silent = TRUE
    )
    if (inherits(res, "try-error")) {
      message(sprintf("Chunk %d failed: %s", i, as.character(res)))
      return(NULL)
    }
    res
  })

  # Combine results
  tsvs <- Filter(function(x) !is.null(x) && file.exists(x), results)
  if (length(tsvs) == 0L) {
    stop("InterProScan produced no usable outputs. Check Docker logs above.")
  }

  df_iprscan <- do.call(rbind, lapply(tsvs, .readIPRscanTsv))

  # Load processed tables (unchanged)
  DBI::dbWriteTable(con, "domain_names",
    df_iprscan |>
      dplyr::select(AccNum, DB.ID, SignDesc, IPRAcc, IPRDesc, StartLoc, StopLoc),
    overwrite = TRUE
  )

  df_protein_domain_pa <- df_iprscan |>
    dplyr::select(AccNum, DB.ID, IPRAcc, placeholder) |>
    dplyr::mutate(domain_ID = stringr::str_glue("{DB.ID}_{IPRAcc}")) |>
    dplyr::distinct() |>
    dplyr::mutate(placeholder = stringr::str_replace_all(placeholder, "-", "1")) |>
    tidyr::pivot_wider(
      id_cols = AccNum, names_from = domain_ID, values_from = placeholder,
      values_fill = "0"
    ) |>
    dplyr::group_by(AccNum) |>
    dplyr::summarize(across(everything(), ~ ifelse(any(. == "1"), "1", "0")), .groups = "drop") |>
    dplyr::mutate(across(-AccNum, as.numeric))

  protein_filter <- dplyr::tbl(con, "protein_count") |> tibble::as_tibble()
  accs <- unique(df_protein_domain_pa$AccNum)
  accs_in_matrix <- intersect(accs, colnames(protein_filter))
  if (length(accs_in_matrix) == 0L) {
    stop("No InterPro accessions match protein_count columns.")
  }

  protein_filter <- protein_filter |> dplyr::select(genome_id, dplyr::all_of(accs_in_matrix))
  df_protein_domain_pa <- df_protein_domain_pa |>
    dplyr::filter(AccNum %in% accs_in_matrix) |>
    dplyr::arrange(match(AccNum, accs_in_matrix))

  domain_count <- as.matrix(protein_filter |> dplyr::select(-genome_id)) %*%
    as.matrix(df_protein_domain_pa |> dplyr::select(-AccNum)) |>
    tibble::as_tibble() |>
    dplyr::mutate(genome_id = protein_filter |> dplyr::pull(genome_id)) |>
    dplyr::relocate(genome_id, .before = dplyr::everything())

  DBI::dbWriteTable(conn = con, name = "domain_count", domain_count, overwrite = TRUE)
  invisible(TRUE)
}

# Clean BV-BRC metadata, then save as Parquet files
cleanData <- function(duckdb_path, path, ref_file_path = "data_raw/") {
  duckdb_path <- normalizePath(duckdb_path)
  # If no explicit path is provided (or a generic one), choose results/<bug>/ when
  # the DuckDB lives under data/<bug>/, or else fall back to the DuckDB directory.
  if (missing(path) || path %in% c(".", "results", "results/")) {
    bug_dir <- dirname(duckdb_path)
    mapped_results <- sub(
      paste0(.Platform$file.sep, "data", .Platform$file.sep),
      paste0(.Platform$file.sep, "results", .Platform$file.sep),
      bug_dir,
      fixed = TRUE
    )
    path <- if (!identical(mapped_results, bug_dir)) mapped_results else bug_dir
  }

  path <- normalizePath(path, mustWork = FALSE)
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)

  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)
  ref_file_path <- normalizePath(ref_file_path)

  clean_drug <- readr::read_tsv(file.path(ref_file_path, "clean_drug.tsv"))
  drug_class <- readr::read_tsv(file.path(ref_file_path, "drug_class.tsv"))
  drug_abbr <- readr::read_tsv(file.path(ref_file_path, "drug_abbr.tsv"))
  class_abbr <- readr::read_tsv(file.path(ref_file_path, "class_abbr.tsv"))
  clean_countries <- readr::read_tsv(file.path(ref_file_path, "cleaned_bvbrc_countries.tsv")) |>
    dplyr::select("raw_entry", "clean_name", "short_name") |>
    dplyr::distinct()

  dplyr::tbl(con, "filtered") |>
    tibble::as_tibble() |>
    dplyr::select(
      "genome_drug.genome_id", "genome_drug.antibiotic",
      "genome_drug.genome_name", "genome_drug.laboratory_typing_method",
      "genome_drug.resistant_phenotype", "genome_drug.taxon_id",
      "genome_drug.pmid", "genome.collection_year",
      "genome.isolation_country", "genome.host_common_name",
      "genome.isolation_source", "genome.species"
    ) |>
    dplyr::left_join(clean_drug, by = c("genome_drug.antibiotic" = "original_drug")) |>
    dplyr::filter(!is.na(cleaned_drug)) |>
    dplyr::left_join(drug_class, by = c("cleaned_drug" = "drug")) |>
    dplyr::left_join(drug_abbr, by = c("cleaned_drug" = "drug")) |>
    dplyr::left_join(class_abbr, by = "drug_class") |>
    DBI::dbWriteTable(conn = con, name = "filtered", overwrite = TRUE)

  resistance_summary <- dplyr::tbl(con, "filtered") |>
    tibble::as_tibble() |>
    dplyr::filter(genome_drug.resistant_phenotype == "Resistant") |>
    dplyr::group_by(genome_drug.genome_id) |>
    dplyr::summarise(
      num_resistant_classes = dplyr::n_distinct(drug_class),
      resistant_classes = paste(unique(class_abbr), collapse = "_")
    )

  year_breaks <- seq(1980, 2023, by = 5)
  dplyr::tbl(con, "filtered") |>
    tibble::as_tibble() |>
    dplyr::mutate(genome_drug.antibiotic = cleaned_drug) |>
    dplyr::select(-cleaned_drug) |>
    dplyr::left_join(clean_countries, by = c("genome.isolation_country" = "raw_entry")) |>
    dplyr::rename("cleaned_country" = "clean_name", "country_abbr" = "short_name") |>
    dplyr::mutate(genome.isolation_country = cleaned_country) |>
    dplyr::select(-cleaned_country) |>
    dplyr::left_join(resistance_summary, by = "genome_drug.genome_id") |>
    dplyr::mutate(resistant_classes = dplyr::case_when(
      is.na(resistant_classes) ~ genome_drug.resistant_phenotype,
      TRUE ~ resistant_classes
    )) |>
    dplyr::mutate(num_resistant_classes = dplyr::case_when(
      is.na(num_resistant_classes) ~ 0,
      TRUE ~ num_resistant_classes
    )) |>
    dplyr::mutate(genome.collection_year = as.numeric(genome.collection_year)) |>
    dplyr::mutate(year_bin = cut(genome.collection_year,
      breaks = year_breaks,
      right = FALSE, include.lowest = TRUE,
      labels = paste(year_breaks[-length(year_breaks)],
        year_breaks[-1] - 1,
        sep = "-"
      )
    )) |>
    DBI::dbWriteTable(conn = con, name = "cleaned_metadata", overwrite = TRUE)

  # Parquet output paths
  genes_parquet <- file.path(path, "gene_count.parquet")
  gene_names_parquet <- file.path(path, "gene_names.parquet")
  gene_ref_seq_parquet <- file.path(path, "gene_seqs.parquet")
  genome_gene_protein_parquet <- file.path(path, "genome_gene_protein.parquet")
  struct_parquet <- file.path(path, "struct.parquet")

  proteins_parquet <- file.path(path, "protein_count.parquet")
  domains_parquet <- file.path(path, "domain_count.parquet")

  metadata_parquet <- file.path(path, "metadata.parquet") # cleaned_metadata exported as 'metadata'

  domain_names_parquet <- file.path(path, "domain_names.parquet")
  protein_names_parquet <- file.path(path, "protein_names.parquet")

  protein_cluster_seq_parquet <- file.path(path, "protein_seqs.parquet")

  # Also export AMR/genome/original metadata
  amr_phenotype_parquet <- file.path(path, "amr_phenotype.parquet")
  genome_data_parquet <- file.path(path, "genome_data.parquet")
  original_metadata_parquet <- file.path(path, "original_metadata.parquet")

  writeCompressedParquet <- function(df, path) {
    arrow::write_parquet(
      df,
      path,
      compression = "zstd",
      compression_level = 9,
      use_dictionary = TRUE
    )
  }

  db_name <- duckdb_path |>
    stringr::str_split_i(".duckdb", i = 1) |>
    paste0("_parquet.duckdb")
  con_new <- DBI::dbConnect(duckdb::duckdb(), db_name)
  on.exit(try(DBI::dbDisconnect(con_new, shutdown = FALSE), silent = TRUE), add = TRUE)

  # gene_count -> long parquet + view
  DBI::dbReadTable(con, "gene_count") |>
    tidyr::pivot_longer(-genome_id, names_to = "gene", values_to = "value") |>
    dplyr::filter(!is.na(value) & value != "") |>
    dplyr::mutate(value = as.integer(value)) |>
    writeCompressedParquet(genes_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW gene_count AS SELECT * FROM read_parquet('%s')", genes_parquet))

  # protein_count -> long parquet + view
  DBI::dbReadTable(con, "protein_count") |>
    tidyr::pivot_longer(-genome_id, names_to = "protein", values_to = "value") |>
    dplyr::filter(!is.na(value) & value != "") |>
    dplyr::mutate(value = as.integer(value)) |>
    writeCompressedParquet(proteins_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW protein_count AS SELECT * FROM read_parquet('%s')", proteins_parquet))

  # domain_count -> long parquet + view
  DBI::dbReadTable(con, "domain_count") |>
    tidyr::pivot_longer(-genome_id, names_to = "domain", values_to = "value") |>
    dplyr::filter(!is.na(value) & value != "") |>
    dplyr::mutate(value = as.integer(value)) |>
    writeCompressedParquet(domains_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW domain_count AS SELECT * FROM read_parquet('%s')", domains_parquet))

  # gene_struct -> long parquet + view
  DBI::dbReadTable(con, "gene_struct") |>
    tidyr::pivot_longer(-genome_id, names_to = "struct", values_to = "value") |>
    dplyr::filter(!is.na(value) & value != "") |>
    dplyr::mutate(value = as.integer(value)) |>
    writeCompressedParquet(struct_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW struct AS SELECT * FROM read_parquet('%s')", struct_parquet))

  # cleaned_metadata -> parquet + view (as metadata)
  DBI::dbReadTable(con, "cleaned_metadata") |> writeCompressedParquet(metadata_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW metadata AS SELECT * FROM read_parquet('%s')", metadata_parquet))

  # names/seq tables -> parquet + views
  DBI::dbReadTable(con, "gene_names") |> writeCompressedParquet(gene_names_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW gene_names AS SELECT * FROM read_parquet('%s')", gene_names_parquet))

  DBI::dbReadTable(con, "protein_names") |>
    dplyr::select(-locus_tag) |>
    writeCompressedParquet(protein_names_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW protein_names AS SELECT * FROM read_parquet('%s')", protein_names_parquet))

  DBI::dbReadTable(con, "domain_names") |>
    dplyr::select(-c(IPRAcc, IPRDesc)) |>
    writeCompressedParquet(domain_names_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW domain_names AS SELECT * FROM read_parquet('%s')", domain_names_parquet))

  DBI::dbReadTable(con, "gene_ref_seq") |> writeCompressedParquet(gene_ref_seq_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW gene_seqs AS SELECT * FROM read_parquet('%s')", gene_ref_seq_parquet))

  DBI::dbReadTable(con, "protein_cluster_seq") |> writeCompressedParquet(protein_cluster_seq_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW protein_seqs AS SELECT * FROM read_parquet('%s')", protein_cluster_seq_parquet))

  DBI::dbReadTable(con, "genome_gene_protein") |> writeCompressedParquet(genome_gene_protein_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW genome_gene_protein AS SELECT * FROM read_parquet('%s')", genome_gene_protein_parquet))

  # debug/complete views: amr_phenotype, genome_data, original_metadata
  DBI::dbReadTable(con, "amr_phenotype") |> writeCompressedParquet(amr_phenotype_parquet)
  DBI::dbReadTable(con, "genome_data") |> writeCompressedParquet(genome_data_parquet)
  DBI::dbReadTable(con, "metadata") |> writeCompressedParquet(original_metadata_parquet)

  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW amr_phenotype AS SELECT * FROM read_parquet('%s')", amr_phenotype_parquet))
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW genome_data AS SELECT * FROM read_parquet('%s')", genome_data_parquet))
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW original_metadata AS SELECT * FROM read_parquet('%s')", original_metadata_parquet))

  invisible(TRUE)
}


#' Run the full amRdata processing pipeline (Panaroo → CD-HIT → InterProScan → Parquet)
#'
#' @description
#' `runDataProcessing()` orchestrates the complete feature-extraction pipeline for a
#' BV-BRC selection, starting from a **per-selection DuckDB** (created by
#' [prepareGenomes()] and populated by downstream steps). It:
#' 1. Runs **Panaroo** to build the pangenome and writes gene/struct outputs into DuckDB.
#' 2. Runs **CD-HIT** to cluster proteins and writes protein outputs into DuckDB.
#' 3. Runs **InterProScan** (Pfam) to annotate protein domains and writes domain outputs into DuckDB.
#' 4. **Cleans BV-BRC metadata** (drug names/classes, countries, years) and
#'    exports all feature/metadata tables as **compressed Parquet** files, then creates
#'    a **Parquet-backed DuckDB** with read-only views of those Parquets for downstream ML.
#'
#' The function is a thin controller that delegates each stage to the corresponding
#' internal helpers (Dockerized tools where applicable) and ensures consistent
#' output locations and table schemas across stages.
#'
#' @section Pipeline Steps:
#' \enumerate{
#'   \item **Panaroo** via runPanaroo2Duckdb() → writes:
#'     \itemize{
#'       \item `gene_count` (genome × gene counts)\cr
#'       \item `gene_names`\cr
#'       \item `gene_struct` (structural variants)\cr
#'       \item `gene_ref_seq`, `genome_gene_protein`
#'     }
#'   \item **CD-HIT** via CDHIT2duckdb() (calls internal `.runCDHIT()`) → writes:
#'     \itemize{
#'       \item `protein_count` (genome × protein-cluster counts)\cr
#'       \item `protein_names`\cr
#'       \item `protein_cluster_seq` (representative sequences)
#'     }
#'   \item **InterProScan (Pfam)** via domainFromIPR() → writes:
#'     \itemize{
#'       \item `domain_names`\cr
#'       \item `domain_count` (genome × domain-family matrix)
#'     }
#'   \item **Metadata cleaning + Parquet export** via cleanData() → writes Parquet
#'         files to `output_path`, and builds a **Parquet-backed DuckDB**
#'         (`*_parquet.duckdb`) with views:
#'     \itemize{
#'       \item `gene_count`, `protein_count`, `domain_count`, `struct`\cr
#'       \item `metadata` (cleaned), plus `amr_phenotype`, `genome_data`, `original_metadata`\cr
#'       \item `gene_names`, `protein_names`, `domain_names`\cr
#'       \item `gene_seqs`, `protein_seqs`\cr
#'       \item `genome_gene_protein`
#'     }
#' }
#'
#' @param duckdb_path Character. Path to the **per-selection DuckDB** produced by
#'   [prepareGenomes()] (e.g., `"data/<Bug>/<Abbrev>.duckdb"`). This DB must
#'   already contain at least the tables written by `prepareGenomes()` and subsequent
#'   download steps (e.g., `files`, `filtered`, and metadata tables).
#' @param output_path Character or `NULL`. Base directory for writing Panaroo/CD-HIT/InterProScan
#'   outputs and final Parquet files. If `NULL`, defaults to `dirname(duckdb_path)`.
#'
#' @param threads Integer. Shared concurrency budget used across tools (Panaroo, CD-HIT,
#'   InterProScan). Passed through to each stage as appropriate. Defaults to `16`.
#'
#' @param panaroo_split_jobs Logical. If `TRUE`, Panaroo runs in multiple batches that can be
#'   merged by [.mergePanaroo()]. If `FALSE`, Panaroo runs once on all isolates. Default: `FALSE`.
#' @param panaroo_core_threshold Numeric. Panaroo `--core_threshold`. Default: `0.90`.
#' @param panaroo_len_dif_percent Numeric. Panaroo `--len_dif_percent`. Default: `0.95`.
#' @param panaroo_cluster_threshold Numeric. Panaroo `--threshold`. Default: `0.95`.
#' @param panaroo_family_seq_identity Numeric. Panaroo `-f` (gene family identity). Default: `0.5`.
#'
#' @param cdhit_identity Numeric. CD-HIT `-c` identity threshold. Default: `0.9`.
#' @param cdhit_word_length Integer. CD-HIT `-n` word length. Default: `5`.
#' @param cdhit_memory Integer. CD-HIT `-M` memory limit (MB). Use `0` for unlimited. Default: `0`.
#' @param cdhit_extra_args Character vector. Extra arguments forwarded to `cd-hit`
#'   (e.g., `c("-g","1")`). Default: `c("-g","1")`.
#' @param cdhit_output_prefix Character. Prefix for CD-HIT output files. Default: `"cdhit_out"`.
#'
#' @param ipr_appl Character vector. InterProScan applications to run; typically `c("Pfam")`.
#'   Default: `c("Pfam")`.
#' @param ipr_threads_unused Deprecated/unused. Kept for backward compatibility; ignored.
#' @param ipr_version Character. InterProScan image tag (e.g., `"5.76-107.0"`). Default: `"5.76-107.0"`.
#' @param ipr_dest_dir Character. Local destination for InterProScan data bundle
#'   (used by `.checkInterProData()`). Default: `"inst/extdata/interpro"`.
#' @param ipr_platform Character. Docker platform string for InterProScan containers,
#'   e.g., `"linux/amd64"`. Default: `"linux/amd64"`.
#' @param auto_prepare_data Logical. If `TRUE`, ensure InterProScan data are present
#'   (download/verify if missing). Default: `TRUE`.
#'
#' @param ref_file_path Character. Directory containing reference TSVs used by cleanData()
#'   for metadata harmonization (e.g., `"data_raw/"`). **Required**; defaults to `"data_raw/"`.
#'
#' @param verbose Logical. Print progress messages. Default: `TRUE`.
#'
#' @return
#' Invisibly returns a list with:
#' \itemize{
#'   \item `duckdb_path` – input DuckDB path
#'   \item `panaroo_output` – path to the selected Panaroo output directory used for import
#'   \item `parquet_duckdb_path` – absolute path to the created Parquet-backed DuckDB
#' }
#'
#' @details
#' **Docker & Platform Notes**
#' * All heavy tools (Panaroo, CD-HIT, InterProScan) run inside Docker containers.
#' * On Apple Silicon/ARM hosts, images are forced to `--platform linux/amd64` to ensure compatibility.
#' * Ensure Docker Desktop is running and has sufficient memory/CPUs configured.
#'
#' **Input Requirements**
#' * The `duckdb_path` must reference a per-selection DuckDB that contains:
#'   `files` (paths to `.gff`, `.fna`, `.PATRIC.faa`),
#'   `filtered` (genomes selected for download/filtering), and
#'   BV-BRC metadata tables written by earlier steps.
#'
#' **Outputs & Side Effects**
#' * Writes tool-specific intermediate outputs under `output_path` (e.g., `panaroo_out_*`, CD-HIT files).
#' * Writes Parquet files to `output_path`:
#'   `gene_count.parquet`, `protein_count.parquet`, `domain_count.parquet`, `struct.parquet`,
#'   `gene_names.parquet`, `protein_names.parquet`, `domain_names.parquet`,
#'   `gene_seqs.parquet`, `protein_seqs.parquet`, `genome_gene_protein.parquet`,
#'   `metadata.parquet`, `amr_phenotype.parquet`, `genome_data.parquet`, `original_metadata.parquet`.
#' * Creates a new Parquet-backed DuckDB (`*_parquet.duckdb`) with read-only views pointing to those Parquets.
#'
#' **Threading**
#' * `threads` is a shared budget; each stage uses a portion or all of it.
#' * InterProScan can be memory-intensive; on laptops, single-container mode is used internally.
#'
#' @seealso
#' prepareGenomes(), runPanaroo2Duckdb(), CDHIT2duckdb(), domainFromIPR(), cleanData()
#'
#' @examples
#' \dontrun{
#' # Paths below are illustrative; adapt to your project layout.
#' runDataProcessing(
#'   duckdb_path   = "data/Shigella_flexneri/Sfl.duckdb",
#'   output_path   = "data/Shigella_flexneri",
#'   threads       = 16,
#'   ref_file_path = "data_raw/"
#' )
#'
#' # After completion:
#' #   data/Shigella_flexneri/Sfl_parquet.duckdb
#' # will contain views over the Parquet files for downstream ML.
#' }
#'
#' @export
runDataProcessing <- function(duckdb_path,
                              output_path = NULL,
                              # unified threads for all tools
                              threads = 16,
                              # Panaroo
                              panaroo_split_jobs = FALSE,
                              panaroo_core_threshold = 0.90,
                              panaroo_len_dif_percent = 0.95,
                              panaroo_cluster_threshold = 0.95,
                              panaroo_family_seq_identity = 0.5,
                              # CD-HIT
                              cdhit_identity = 0.9,
                              cdhit_word_length = 5,
                              cdhit_memory = 0,
                              cdhit_extra_args = c("-g", "1"),
                              cdhit_output_prefix = "cdhit_out",
                              # InterPro
                              ipr_appl = c("Pfam"),
                              ipr_threads_unused = NULL,
                              ipr_version = "5.76-107.0",
                              ipr_dest_dir = "inst/extdata/interpro",
                              ipr_platform = "linux/amd64",
                              auto_prepare_data = TRUE,
                              # Metadata cleaning
                              ref_file_path = "data_raw/",
                              verbose = TRUE) {
  duckdb_path <- normalizePath(duckdb_path)
  out_dir <- if (is.null(output_path)) dirname(duckdb_path) else normalizePath(output_path)

  # 1) Panaroo (run + optional merge) -> write Panaroo tables
  pan_dir <- runPanaroo2Duckdb(
    duckdb_path            = duckdb_path,
    output_path            = out_dir,
    core_threshold         = panaroo_core_threshold,
    len_dif_percent        = panaroo_len_dif_percent,
    cluster_threshold      = panaroo_cluster_threshold,
    family_seq_identity    = panaroo_family_seq_identity,
    threads                = threads,
    split_jobs             = panaroo_split_jobs,
    verbose                = verbose
  )

  # 2) CD-HIT -> write `protein` tables
  if (isTRUE(verbose)) message("Running CD-HIT and writing protein tables to DuckDB.")
  CDHIT2duckdb(
    duckdb_path   = duckdb_path,
    output_path   = out_dir,
    output_prefix = cdhit_output_prefix,
    identity      = cdhit_identity,
    word_length   = cdhit_word_length,
    threads       = threads,
    memory        = cdhit_memory,
    extra_args    = cdhit_extra_args
  )

  # 3) InterProScan -> write `domain` tables
  if (isTRUE(verbose)) message("Running InterProScan and writing domain tables to DuckDB.")
  domainFromIPR(
    duckdb_path       = duckdb_path,
    path              = out_dir,
    out_file_base     = "iprscan",
    appl              = ipr_appl,
    ipr_version       = ipr_version,
    ipr_dest_dir      = ipr_dest_dir,
    ipr_platform      = ipr_platform,
    auto_prepare_data = auto_prepare_data,
    threads           = threads,
    file_format       = "TSV"
  )

  # 4) Clean metadata and export Parquet + Parquet-backed DuckDB
  if (missing(ref_file_path) || is.null(ref_file_path)) {
    stop("`ref_file_path` (directory with reference TSVs) must be provided to cleanData().")
  }
  if (isTRUE(verbose)) message("Cleaning metadata and exporting Parquet-backed views.")
  cleanData(duckdb_path = duckdb_path, path = out_dir, ref_file_path = ref_file_path)

  parquet_duckdb_path <- paste0(
    stringr::str_split_i(duckdb_path, ".duckdb", i = 1),
    "_parquet.duckdb"
  )

  if (isTRUE(verbose)) {
    message("\n============================================")
    message("Completed data-processing pipeline successfully.")
    message("Parquet-backed DuckDB created at:")
    message("  ", normalizePath(parquet_duckdb_path))
    message("\nYou can use the amRml package to train machine")
    message("learning models for AMR using this file path.")
    message("For example:")
    message("  runMLmodels(\"", normalizePath(parquet_duckdb_path), "\")")
    message("============================================\n")
  }

  invisible(list(
    duckdb_path = duckdb_path,
    panaroo_output = pan_dir,
    parquet_duckdb_path = normalizePath(parquet_duckdb_path)
  ))
}
