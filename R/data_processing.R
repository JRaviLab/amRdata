#' @importFrom data.table :=
NULL

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
#' @param refind_mode Character. Panaroo's `--refind-mode` (`"off"`, `"default"`, or
#'   `"strict"`). Refinding searches for and recovers gene calls that annotation
#'   tools missed, comparing each candidate against the rest of the pangenome.
#'   Caveat: this search can take substantially longer (and in rare cases fail to
#'   complete within hours) when a genome carries a cluster of CDS with internal
#'   stop codons, which existing upstream genome-quality fields do not flag.
#'   Default `"off"` for now, to avoid that runtime risk; plan to move this back to
#'   `"default"` once a QC step upstream (e.g. in `.apply_metadata_qc()`) can screen
#'   out affected genomes before they reach Panaroo.
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
                            panaroo_threads_per_job,
                            refind_mode = c("off", "default", "strict"),
                            verbose = TRUE) {
  refind_mode <- match.arg(refind_mode)
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  output_path <- .docker_path(output_path)

  # Fail fast if Docker is missing
  if (!nzchar(Sys.which("docker"))) {
    stop("Docker is not available on your PATH but is required to run Panaroo.")
  }

  # Host mount root = bug directory
  mount_host <- output_path
  mount_cont <- "/work"

  # Write the genome list file (convert each "gff fna" to container-visible paths)
  genome_filepath_host <- tempfile(pattern = "genomeFilepath_", fileext = ".txt", tmpdir = output_path)

  batch_input_cont <- purrr::map_chr(unlist(batch_input), function(line) {
    parts <- strsplit(line, " +")[[1]]
    parts_cont <- .to_container(parts, host_root = mount_host, container_root = mount_cont)
    paste(parts_cont, collapse = " ")
  })

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
    "staphb/panaroo:1.7.0",
    "panaroo",
    "-i", genome_filepath_cont,
    "-o", output_dir_cont,
    "--clean-mode", "strict",
    "--merge_paralogs",
    "--remove-invalid-genes",
    "--refind-mode", refind_mode,
    "--core_threshold", as.character(core_threshold),
    "--len_dif_percent", as.character(len_dif_percent),
    "--threshold", as.character(cluster_threshold),
    "-f", as.character(family_seq_identity),
    "-t", as.character(panaroo_threads_per_job)
  )

  # Updating to try controlling Panaroo's print verbosity
  res <- system2("docker", args = cmd_args, stdout = TRUE, stderr = TRUE)
  status <- attr(res, "status")

  if (!is.null(status) && status != 0L) {
    stop(sprintf("Panaroo failed with exit status %s:\n%s", status, paste(res, collapse = "\n")))
  }

  if(isTRUE(verbose)) {
    message("Panaroo output: \n", paste(res, collapse = "\n"))
  }

  if (inherits(res, "error")) {
    stop(sprintf("Docker/Panaroo failed to launch: %s", res$message))
  }

  invisible(res)
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
#' @param refind_mode Character. Panaroo's `--refind-mode` (`"off"`, `"default"`, or
#'   `"strict"`). See [.processPanaroo()] for what refinding does and the runtime
#'   caveat behind the current default. Default `"off"`.
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
                        split_jobs = FALSE,
                        refind_mode = c("off", "default", "strict"),
                        strip_pseudogenes = FALSE,
                        pseudogene_clean_dir = "gff_clean",
                        write_pseudogene_audit = TRUE,
                        verbose = TRUE) {
  refind_mode <- match.arg(refind_mode)
  duckdb_path <- normalizePath(duckdb_path)
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  if (missing(output_path) || output_path %in% c(".", "results", "results/")) {
    output_path <- dirname(duckdb_path)
  }
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  output_path <- normalizePath(output_path)

  genome_query_output <- DBI::dbGetQuery(con, "SELECT * FROM files ORDER BY genome_id")

  panaroo_input_files <- genome_query_output |>
    dplyr::pull(panaroo_input)

  # Drop true NAs
  panaroo_input_files <- panaroo_input_files[!is.na(panaroo_input_files)]

  # Cleaning pseudogene lines out of GFF files
  if (isTRUE(strip_pseudogenes)) {
    cleaned <- .stripPseudogeneGFFs(
      panaroo_input_files = panaroo_input_files,
      output_path = output_path,
      clean_dir = pseudogene_clean_dir
    )

    panaroo_input_files <- cleaned$panaroo_input_files

    if (isTRUE(write_pseudogene_audit)) {
      readr::write_csv(cleaned$audit, file.path(output_path, "panaroo_pseudogene_audit.csv"))
    }

    if (isTRUE(verbose)) {
      audit <- cleaned$audit
      message(sprintf(
        "Pseudogene audit: %d genomes, %d removed pseudogenes, %d features remain.",
        nrow(audit),
        sum(audit$n_pseudogene, na.rm = TRUE),
        sum(audit$n_kept, na.rm = TRUE)
      ))
    }
  }

  split_files <- strsplit(panaroo_input_files, " ")

  valid_entries <- purrr::map_lgl(split_files, function(paths) {
    gff_file <- paths[1]
    if (file.exists(gff_file)) {
      length(readLines(gff_file, n = 5, warn = FALSE)) >= 5
    } else {
      FALSE
    }
  })

  filtered_panaroo_input <- purrr::map_chr(split_files[valid_entries], paste, collapse = " ")

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

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  if (n_jobs <= 1L) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession, workers = n_jobs)
  }

  batch_panaroo_run <- furrr::future_map(
    panaroo_batches,
    ~ .processPanaroo(
      batch_input             = .x,
      output_path             = output_path,
      core_threshold          = core_threshold,
      len_dif_percent         = len_dif_percent,
      cluster_threshold       = cluster_threshold,
      family_seq_identity     = family_seq_identity,
      panaroo_threads_per_job = panaroo_threads_per_job,
      refind_mode             = refind_mode,
      verbose                 = verbose
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
      "staphb/panaroo:1.7.0",
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
    dplyr::mutate(genome_ids = sub("\\.PATRIC\\.\\.\\..*$", "", genome_ids)) |>
    dplyr::select(genome_ids, Gene, protein_ids) |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(protein_ids)) |>
    tidyr::separate_rows(protein_ids, sep = ";") |>
    # dplyr::filter(!stringr::str_detect(protein_ids, "_pseudo")) |>
    dplyr::mutate(protein_ids = gsub("_pseudo", "", protein_ids)) |>
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
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  output_path <- .docker_path(output_path)

  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  genome_query_output <- DBI::dbGetQuery(con, "SELECT * FROM files ORDER BY genome_id")

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
#'   Default: `8`.
#'
#' @param split_jobs Logical. If `TRUE`, Panaroo is run in multiple parallel
#'   batches (up to 5, depending on dataset size), and batch outputs are merged
#'   using `.mergePanaroo()`. If `FALSE`, only one Panaroo invocation is run.
#'   Default: `FALSE`.
#'
#' @param refind_mode Character. Panaroo's `--refind-mode` (`"off"`, `"default"`, or
#'   `"strict"`). See [.processPanaroo()] for what refinding does and the runtime
#'   caveat behind the current default. Default `"off"`.
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
#' * [runDataProcessing()] — full pipeline including CD-HIT & HMMER
#'
#' @examples
#' \dontrun{
#' # Basic usage:
#' runPanaroo2Duckdb(
#'   duckdb_path = "data/Shigella_flexneri/Sfl.duckdb",
#'   output_path = "data/Shigella_flexneri",
#'   threads     = 8,
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
                              threads = 8,
                              split_jobs = FALSE,
                              refind_mode = c("off", "default", "strict"),
                              strip_pseudogenes = FALSE,
                              pseudogene_clean_dir = "gff_clean",
                              write_pseudogene_audit = TRUE,
                              verbose = TRUE) {
  refind_mode <- match.arg(refind_mode)
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
    split_jobs = split_jobs,
    refind_mode = refind_mode,
    strip_pseudogenes = strip_pseudogenes,
    pseudogene_clean_dir = pseudogene_clean_dir,
    write_pseudogene_audit = write_pseudogene_audit,
    verbose = verbose
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
      stringr::str_extract(ref_line, "fig\\|[0-9]+\\.[0-9]+\\.peg(?:sc)?\\.[0-9]+")
    } else {
      paste0("Cluster_", i - 1)
    }

    # Pull genome IDs
    genome_matches <- stringr::str_match(
      cluster_lines,
      "fig\\|([0-9]+\\.[0-9]+)\\.peg(?:sc)?\\.[0-9]+"
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

#' Parse CD-HIT `.clstr` output into a long-format mapping
#'
#' Reads a CD-HIT `.clstr` file and constructs a mapping of clusters to member feature ids.
#'
#' @param clustered_faa Base path to CD-HIT output (without `.clstr` extension).
#'
#' @return A tibble with columns `cluster` and `member`.
#'
#' @keywords internal
.extractMembersInClusters <- function(clustered_faa) {
  clstr <- paste0(clustered_faa, ".clstr")
  if (!file.exists(clstr)) {
    stop(
      "CD-HIT cluster file not found: ", clstr,
      "\nEnsure .runCDHIT() completed successfully and produced the .clstr file."
    )
  }

  lines <- data.table::fread(clstr, sep = "\n", header = FALSE)$V1
  cluster_ids <- grep("^>Cluster", lines)
  cluster_member <- data.table::data.table()

  for (i in seq_along(cluster_ids)) {
    start <- cluster_ids[i] + 1
    end <- if (i < length(cluster_ids)) cluster_ids[i + 1] - 1 else length(lines)
    cluster_lines <- lines[start:end]

    # This finds the reference cluster ID and names the cluster with it
    ref_line <- grep("\\*$", cluster_lines, value = TRUE)
    ref_id <- if (length(ref_line) > 0) {
      stringr::str_extract(ref_line, "fig\\|[0-9]+\\.[0-9]+\\.peg(?:sc)?\\.[0-9]+")
    } else {
      paste0("Cluster_", i - 1)
    }

    # Pull genome IDs
    members <- stringr::str_match(
      cluster_lines,
      "fig\\|([0-9]+\\.[0-9]+)\\.peg(?:sc)?\\.[0-9]+"
    )[, 1]
    members <- members[!is.na(members)]

    if (length(members) > 0) {
      cluster_member <- data.table::rbindlist(list(
        cluster_member,
        data.table::data.table(cluster = ref_id, member = members)
      ), use.names = TRUE)
    }
  }

  tibble::as_tibble(cluster_member)
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
      proteinID = stringr::str_extract(value, "^fig\\|[0-9]+\\.[0-9]+\\.peg(?:sc)?\\.[0-9]+"),
      locus_tag = stringr::str_match(value, "peg(?:sc)?\\.[0-9]+\\|([^\\s]+)")[, 2],
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
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
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
      name     = names(clustered_faa) |> stringr::str_extract("fig\\|[0-9]+\\.[0-9]+\\.peg(?:sc)?\\.[0-9]+"),
      sequence = as.character(clustered_faa)
    ),
    overwrite = TRUE
  )

  cluster_member <- .extractMembersInClusters(cdhit_outputs$clustered_faa)
  DBI::dbWriteTable(con, "protein_members", cluster_member, overwrite = TRUE)

  invisible(TRUE)
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
.prepareHmmerDatabases <- function(
    hmmer_db_dir,
    databases = c("Pfam", "COG", "AMRFinder"),
    docker_image = "staphb/hmmer",
    hmmer_db_url = NULL,
    verbose = TRUE
) {

  hmmer_db_dir <- normalizePath(
    hmmer_db_dir,
    mustWork = FALSE
  )

  dir.create(
    hmmer_db_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

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

  # Add custom database(s)
  missing_dbs <- setdiff(databases, names(dbs))

  if (length(missing_dbs) > 0) {

    if (is.null(hmmer_db_url)) {
      stop(
        "hmmer_db_url must be supplied when using custom databases"
      )
    }

    get_db_type <- function(url) {

      file <- basename(url)

      if (grepl("\\.(tar\\.gz|tgz)$", file, ignore.case = TRUE)) {
        return("tar.gz")
      } else if (grepl("\\.zip$", file, ignore.case = TRUE)) {
        return("zip")
      } else if (grepl("\\.gz$", file, ignore.case = TRUE)) {
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

    if (verbose) {
      message("Checking ", db_name)
    }

    hmm_files <- list.files(
      db$dir,
      pattern = "\\.hmm$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )

    if (length(hmm_files) == 0) {

      if (verbose) {
        message("Downloading ", db_name)
      }

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
              sub(
                "\\.gz$",
                "",
                basename(db$url),
                ignore.case = TRUE
              )
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

      unlink(tmp)

      hmm_files <- list.files(
        db$dir,
        pattern = "\\.hmm$",
        recursive = TRUE,
        full.names = TRUE,
        ignore.case = TRUE
      )
    }

    if (length(hmm_files) == 0) {
      stop(
        "No .hmm file found for ",
        db_name
      )
    }

    if (length(hmm_files) == 1) {

      hmm_file <- hmm_files[[1]]

    } else {

      hmm_file <- file.path(
        db$dir,
        paste0(db_name, ".hmm")
      )

      source_hmms <- setdiff(
        normalizePath(hmm_files),
        normalizePath(hmm_file, mustWork = FALSE)
      )

      valid_hmms <- purrr::map_lgl(
        source_hmms,
        .isValidHmmFile
      )

      if (any(!valid_hmms)) {

        bad_files <- basename(
          source_hmms[!valid_hmms]
        )

        if (isTRUE(verbose)) {
          warning(
            "Ignoring ",
            length(bad_files),
            " invalid HMM file(s):\n",
            paste(bad_files, collapse = "\n"),
            call. = FALSE
          )
        }

        source_hmms <- source_hmms[valid_hmms]
      }

      if (length(source_hmms) == 0) {
        stop(
          "No valid HMM files found for ",
          db_name
        )
      }

      if (!file.exists(hmm_file)) {

        if (verbose) {
          message(
            "Combining ",
            length(source_hmms),
            " HMM files for ",
            db_name
          )
        }

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

      if (verbose) {
        message(
          "Running hmmpress for ",
          basename(hmm_file)
        )
      }

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

    db_paths[[db_name]] <- list(
      hmm = hmm_file,
      source = db$url,
      type = db$type,
      pressed = pressed_files
    )

    if (verbose) {
      message(
        db_name,
        " ready: ",
        hmm_file
      )
    }
  }

  db_paths
}





#' The function to run HMMER with docker
#'
#' @param JOB_NAME protein_chunk id
#' @param FASTA fasta sequences
#' @param DB HMM database
#' @param Total_proteins protein sequence count
#' @param output_path path for saving hmmer outputs
#' @param db_paths path to HMM database 
#' @param docker_image hmmer docker image (ideally from dockerhub)
#' @param threads number of threads
#' @param n_workers number of parallel workers
#'
#' @returns the filename of the parquet file with hmmer output post parsing
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
  database_path <- db_paths[[DB]]$hmm
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

  # Adding an E value cutoff here
  hmmer_tbl <- .parseHMMEROutput(hmmer_output) |>
    dplyr::filter(i_evalue <= 1e-5) |>
    dplyr::select(
      protein,
      query_name,
      query_accession,
      target_description,
      i_evalue,
      domain_score
    )

  hmmer_tbl_filename <- file.path(
    dirname(hmmer_output),
    paste0(tools::file_path_sans_ext(basename(hmmer_output)), ".parquet")
  )

  .write_compressed_parquet(hmmer_tbl, hmmer_tbl_filename)

  hmmer_tbl_filename
}


#' Wrapper for preparing HMM databases and running HMMER on protein sequences from duckdb and writing them.
#'
#' @param duckdb_path path to the duckdb with protein sequences and list
#' @param output_path path where HMMER output will be saved
#' @param threads number of threads
#' @param hmmer_db_dir path to the directory where HMM databases are/will be downloaded
#' @param databases list of HMM databases
#' @param docker_image the docker image of HMMER
#' @param num_of_splits The number of splits of the protein sequence file for parallel processing
#' @param n_workers The number of parallel runs
#'
#' @keywords internal
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

  # Just in case CD-HIT failed to generate sequences somehow
  if (nrow(prot_seqs) == 0L) {
    stop("No sequences found in 'protein_cluster_seq'. Please run CDHIT2duckdb() first.")
  }

  # required to define the database size for hmmsearch --Z and --domZ parameters
  Total_proteins <- nrow(prot_seqs)

  if (is.null(hmmer_db_dir)) {
    hmmer_db_dir <- .defaultHmmerDbDir()
  }

  dir.create(
    hmmer_db_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # database paths
  if(verbose) message ("Preparing HMM databases")
  db_paths <- .prepareHmmerDatabases(
    hmmer_db_dir = hmmer_db_dir,
    databases = databases,
    docker_image = docker_image,
    verbose = verbose
  )

  db_paths <- db_paths[databases]

  # validate split counts before propagating possible hogwash
  num_of_splits <- as.integer(num_of_splits)
  if (is.na(num_of_splits) || num_of_splits < 1L) {
    stop("'num_of_splits' parameter must be a positive integer!")
  }

  # clamp splits to the number of sequences available
  chunk_count <- min(num_of_splits, nrow(prot_seqs))

  split_fasta <- function(seqs, prefix) {
    records <- paste0(">", seqs$name, "\n", seqs$sequence)
    chunk_size <- ceiling(length(records) / chunk_count)
    chunks <- split(records, ceiling(seq_along(records) / chunk_size))

    purrr::walk2(chunks, seq_along(chunks), function(chunk, i) {
      chunk_path <- file.path(output_path, sprintf("%s_chunk_%02d.fasta", prefix, i))
      readr::write_lines(chunk, chunk_path)
    })

    length(chunks)
  }

  actual_chunk_count <- split_fasta(prot_seqs, "protein")

  job_list <- expand.grid(
    chunk = sprintf("%02d", seq_len(actual_chunk_count)),
    db = databases,
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(
      JOB_NAME = paste0("protein_chunk_", chunk, "_", db),
      FASTA = paste0("protein_chunk_", chunk, ".fasta"),
      DB = db
    ) |>
    dplyr::select(JOB_NAME, FASTA, DB)

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

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
      dplyr::left_join(.parse_hmmer_profiles(db_paths[[database_name]]$hmm) |>
                         dplyr::select(query_name = profile_name, description = profile_description),
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

  unlink(
    list.files(
      output_path,
      pattern = "^protein_chunk_.*\\.(fasta|tbl|parquet)$",
      full.names = TRUE
    )
  )

  invisible(list(
    databases = db_paths,
    outputs = final_parquets
  ))

  # purrr::map(parquet_files, arrow::read_parquet) |>
  #   dplyr::bind_rows() |>
  #   .write_compressed_parquet(final_parquet)

  # message("Combined parquet written.")

  # arrow::read_parquet(final_parquet) |>
  #   DBI::dbWriteTable(conn = con, name = tools::file_path_sans_ext(basename(final_parquet)), overwrite = TRUE)
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
#' @keywords internal
.proteinAnnotations2Duckdb <- function(
    duckdb_path,
    databases,
    output_path = dirname(duckdb_path)
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
      tibble::as_tibble() |>
      dplyr::distinct(
        protein,
        query_name
      )

    genome_annot_matrix <- protein_long |>
      dplyr::inner_join(
        annotation |>
          dplyr::select(
            protein,
            query_name
          ),
        by = "protein",
        relationship = "many-to-many"
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
      output_path,
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

    valid_hmms <- purrr::map_lgl(
      hmm_files,
      .isValidHmmFile
    )

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
      stop(
        "hmmsearch failed for ",
        db_name,
        "\n",
        paste(output, collapse = "\n")
      )
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

  unlink(
    c(
      fasta_file,
      file.path(
        output_path,
        paste0("protein_", names(databases), ".tbl")
      )
    )
  )

  invisible(list(
    databases = list(
      DefenseFinder = defense_hmm,
      CasFinder = cas_hmm
    ),
    output = parquet_file
  ))
}


#' Clean BV-BRC metadata, then save as Parquet files
#'
#' @param duckdb_path Path to the **per-selection DuckDB** produced by
#'   [prepareGenomes()] (e.g., `"data/<Bug>/<Abbrev>.duckdb"`). This DB must
#'   already contain the tables written by [prepareGenomes()] and the upstream
#'   genome-processing steps.
#' @param path the path to working directory
#' @param ref_file_path Directory containing reference TSVs used by
#'   [cleanMetaData()] and [cleanData()] for metadata harmonization.
#'   Default: `"data_raw/"`.
#' 
#' @export
cleanMetaData <- function(duckdb_path, path, ref_file_path = "data_raw/") {
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

  # Define lab methods
  lab_methods <- c("Disk diffusion", "MIC", "Broth dilution", "Agar dilution", "Biofosun Gram-positive panels broth dilution",
                  "Vitek_2-P607_card", "cation-adjusted Mueller-Hinton broth", "gradient_diffusion", "kirby-bauer_disc_diffusion")

  dplyr::tbl(con, "filtered") |>
    tibble::as_tibble() |>
    dplyr::select("genome.genome_id") |>
    dplyr::left_join(dplyr::tbl(con, "metadata") |>
      tibble::as_tibble(), by = dplyr::join_by("genome.genome_id" == "genome_drug.genome_id")) |>
    dplyr::select(
      "genome.genome_id", "genome_drug.antibiotic",
      "genome_drug.genome_name", "genome_drug.evidence", "genome_drug.laboratory_typing_method",
      "genome_drug.resistant_phenotype", "genome_drug.taxon_id",
      "genome_drug.pmid", "genome.collection_year",
      "genome.isolation_country", "genome.host_common_name",
      "genome.isolation_source", "genome.species"
    ) |>
    dplyr::mutate(genome_drug.evidence = dplyr::case_when(
      genome_drug.laboratory_typing_method %in% lab_methods ~ "Laboratory Method",
      genome_drug.laboratory_typing_method == "Computational Prediction"  ~ "Computational Method",
      TRUE ~ genome_drug.evidence)) |>
    dplyr::filter(genome_drug.evidence == "Laboratory Method") |>
    dplyr::left_join(clean_drug, by = c("genome_drug.antibiotic" = "original_drug")) |>
    dplyr::filter(!is.na(cleaned_drug)) |>
    dplyr::left_join(drug_class, by = c("cleaned_drug" = "drug")) |>
    dplyr::left_join(drug_abbr, by = c("cleaned_drug" = "drug")) |>
    dplyr::left_join(class_abbr, by = "drug_class") |>
    dplyr::filter(genome_drug.resistant_phenotype %in% c("Resistant", "Susceptible")) |>
    DBI::dbWriteTable(conn = con, name = "filtered_metadata", overwrite = TRUE)

  resistance_summary <- DBI::dbReadTable(con, "filtered_metadata") |>
    dplyr::filter(genome_drug.resistant_phenotype == "Resistant") |>
    dplyr::group_by(genome.genome_id) |>
    dplyr::summarise(
      resistant_classes = paste(sort(unique(class_abbr)), collapse = "_"),
      .groups = "drop"
    ) |>
    dplyr::collect() |>
    dplyr::mutate(
      num_resistant_classes = stringr::str_count(resistant_classes, "_") + 1
    )

  year_breaks <- seq(1980, 2026, by = 5)
  dplyr::tbl(con, "filtered_metadata") |>
    tibble::as_tibble() |>
    dplyr::mutate(genome_drug.antibiotic = cleaned_drug) |>
    dplyr::select(-cleaned_drug) |>
    dplyr::left_join(clean_countries, by = c("genome.isolation_country" = "raw_entry")) |>
    dplyr::rename("cleaned_country" = "clean_name", "country_abbr" = "short_name") |>
    dplyr::mutate(genome.isolation_country = cleaned_country) |>
    dplyr::select(-cleaned_country) |>
    dplyr::left_join(resistance_summary, by = "genome.genome_id") |>
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

  # Parquet output path
  metadata_parquet <- file.path(path, "metadata.parquet") # cleaned_metadata exported as 'metadata'

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

  # Views below reference parquet files by bare filename. Point DuckDB at the
  # parquet directory so schema inference at CREATE VIEW time can resolve them.
  DBI::dbExecute(con_new, sprintf("SET file_search_path='%s'", path))

  # cleaned_metadata -> parquet + view (as metadata)
  DBI::dbReadTable(con, "cleaned_metadata") |> writeCompressedParquet(metadata_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW metadata AS SELECT * FROM read_parquet('%s')", basename(metadata_parquet)))

  # debug/complete views: amr_phenotype, genome_data, original_metadata
  DBI::dbReadTable(con, "amr_phenotype") |> writeCompressedParquet(amr_phenotype_parquet)
  DBI::dbReadTable(con, "genome_data") |> writeCompressedParquet(genome_data_parquet)
  DBI::dbReadTable(con, "metadata") |> writeCompressedParquet(original_metadata_parquet)

  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW amr_phenotype AS SELECT * FROM read_parquet('%s')", basename(amr_phenotype_parquet)))
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW genome_data AS SELECT * FROM read_parquet('%s')", basename(genome_data_parquet)))
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW original_metadata AS SELECT * FROM read_parquet('%s')", basename(original_metadata_parquet)))

  invisible(TRUE)
}

#' Clean feature matrices, then save as Parquet files
#'
#' @param duckdb_path Path to the **per-selection DuckDB** produced by
#'   [prepareGenomes()] (e.g., `"data/<Bug>/<Abbrev>.duckdb"`). This DB must
#'   already contain the tables written by [prepareGenomes()] and the upstream
#'   genome-processing steps.
#' @param path the path to working directory
#'
#' @export
cleanData <- function(duckdb_path, path) {
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

  # Fun new manifest action allows cleanData to find applicable database names
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

  for (run in rev(manifest$runs)) {
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

  hmmer_databases <- unlist(
    hmmer_stage$parameters$databases
  )

  if (!length(hmmer_databases)) {
    stop(
      "HMMER stage in manifest does not contain any databases."
    )
  }

  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  .proteinAnnotations2Duckdb(
    duckdb_path = duckdb_path,
    databases = hmmer_databases,
    output_path = path
  )

  # Parquet output paths
  genes_parquet <- file.path(path, "gene_count.parquet")
  gene_names_parquet <- file.path(path, "gene_names.parquet")
  gene_ref_seq_parquet <- file.path(path, "gene_seqs.parquet")
  genome_gene_protein_parquet <- file.path(path, "genome_gene_protein.parquet")
  struct_parquet <- file.path(path, "struct.parquet")

  proteins_parquet <- file.path(path, "protein_count.parquet")

  protein_names_parquet <- file.path(path, "protein_names.parquet")

  protein_cluster_seq_parquet <- file.path(path, "protein_seqs.parquet")
  protein_cluster_member_parquet <- file.path(path, "protein_members.parquet")

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

  # Views below reference parquet files by bare filename. Point DuckDB at the
  # parquet directory so schema inference at CREATE VIEW time can resolve them.
  DBI::dbExecute(con_new, sprintf("SET file_search_path='%s'", path))

  # gene_count -> long parquet + view
  DBI::dbReadTable(con, "gene_count") |>
    tidyr::pivot_longer(-genome_id, names_to = "gene", values_to = "value") |>
    dplyr::filter(!is.na(value) & value != "") |>
    dplyr::mutate(value = as.integer(value)) |>
    writeCompressedParquet(genes_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW gene_count AS SELECT * FROM read_parquet('%s')", basename(genes_parquet)))

  # protein_count -> long parquet + view
  DBI::dbReadTable(con, "protein_count") |>
    tidyr::pivot_longer(-genome_id, names_to = "protein", values_to = "value") |>
    dplyr::filter(!is.na(value) & value != "") |>
    dplyr::mutate(value = as.integer(value)) |>
    writeCompressedParquet(proteins_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW protein_count AS SELECT * FROM read_parquet('%s')", basename(proteins_parquet)))

  # HMMER annotation counts -> long Parquet + views per database in manifest
  for (database in hmmer_databases) {

    count_table <- paste0(
      "protein_",
      database,
      "_count"
    )

    count_parquet <- file.path(
      path,
      paste0(count_table, ".parquet")
    )

    DBI::dbReadTable(con, count_table) |>
      tidyr::pivot_longer(
        -genome_id,
        names_to = "annotation",
        values_to = "value"
      ) |>
      dplyr::rename(!!database := annotation) |>
      dplyr::filter(!is.na(value) & value != "") |>
      dplyr::mutate(value = as.integer(value)) |>
      writeCompressedParquet(count_parquet)

    DBI::dbExecute(
      con_new,
      sprintf(
        "CREATE OR REPLACE VIEW %s AS SELECT * FROM read_parquet('%s')",
        count_table,
        basename(count_parquet)
      )
    )
  }

  # gene_struct -> long parquet + view
  DBI::dbReadTable(con, "gene_struct") |>
    tidyr::pivot_longer(-genome_id, names_to = "struct", values_to = "value") |>
    dplyr::filter(!is.na(value) & value != "") |>
    dplyr::mutate(value = as.integer(value)) |>
    writeCompressedParquet(struct_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW struct AS SELECT * FROM read_parquet('%s')", basename(struct_parquet)))

  # names/seq tables -> parquet + views
  DBI::dbReadTable(con, "gene_names") |> writeCompressedParquet(gene_names_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW gene_names AS SELECT * FROM read_parquet('%s')", basename(gene_names_parquet)))

  DBI::dbReadTable(con, "protein_names") |>
    dplyr::select(-locus_tag) |>
    writeCompressedParquet(protein_names_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW protein_names AS SELECT * FROM read_parquet('%s')", basename(protein_names_parquet)))

  # Parsing through the different HMMER result Parquets
  for (database in hmmer_databases) {

    annotation_table <- paste0(
      "protein_",
      database
    )

    annotation_parquet <- file.path(
      path,
      paste0(annotation_table, ".parquet")
    )

    DBI::dbReadTable(con, annotation_table) |>
      writeCompressedParquet(annotation_parquet)

    DBI::dbExecute(
      con_new,
      sprintf(
        "CREATE OR REPLACE VIEW %s AS SELECT * FROM read_parquet('%s')",
        annotation_table,
        basename(annotation_parquet)
      )
    )
  }

  DBI::dbReadTable(con, "gene_ref_seq") |> writeCompressedParquet(gene_ref_seq_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW gene_seqs AS SELECT * FROM read_parquet('%s')", basename(gene_ref_seq_parquet)))

  DBI::dbReadTable(con, "protein_cluster_seq") |> writeCompressedParquet(protein_cluster_seq_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW protein_seqs AS SELECT * FROM read_parquet('%s')", basename(protein_cluster_seq_parquet)))

  DBI::dbReadTable(con, "protein_members") |> writeCompressedParquet(protein_cluster_member_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW protein_members AS SELECT * FROM read_parquet('%s')", protein_cluster_member_parquet))

  DBI::dbReadTable(con, "genome_gene_protein") |> writeCompressedParquet(genome_gene_protein_parquet)
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW genome_gene_protein AS SELECT * FROM read_parquet('%s')", basename(genome_gene_protein_parquet)))

  invisible(TRUE)
}


#' Run the full amRdata processing pipeline (Panaroo -> CD-HIT -> HMMER -> Parquet)
#'
#' @description
#' `runDataProcessing()` orchestrates the complete feature-extraction pipeline for a
#' BV-BRC selection, starting from a **per-selection DuckDB** created by
#' [prepareGenomes()] and populated by downstream genome processing steps. It:
#' 1. Runs **Panaroo** to build the pangenome and writes gene/struct outputs into DuckDB.
#' 2. Runs **CD-HIT** to cluster proteins and writes protein outputs into DuckDB.
#' 3. Runs **HMMER** against the requested protein databases and writes annotation
#'    tables into DuckDB.
#' 4. **Cleans BV-BRC metadata** (drug names/classes, countries, years) and
#'    exports feature and metadata tables as compressed Parquet files, then creates
#'    a **Parquet-backed DuckDB** with read-only views for downstream ML.
#'
#' The function is a thin controller that delegates each stage to the corresponding
#' internal helpers (Dockerized tools where applicable) and records processing
#' parameters and provenance in the dataset manifest.
#'
#' @section Pipeline Steps:
#' \enumerate{
#'   \item **Panaroo** via [runPanaroo2Duckdb()] -> writes:
#'     \itemize{
#'       \item `gene_count` (genome x gene counts)\cr
#'       \item `gene_names`\cr
#'       \item `gene_struct` (structural variants)\cr
#'       \item `gene_ref_seq`, `genome_gene_protein`
#'     }
#'   \item **CD-HIT** via [CDHIT2duckdb()] -> writes:
#'     \itemize{
#'       \item `protein_count` (genome x protein-cluster counts)\cr
#'       \item `protein_names`\cr
#'       \item `protein_cluster_seq` (representative sequences)\cr
#'       \item `protein_members`
#'     }
#'   \item **HMMER** via the configured HMMER databases -> writes:
#'     \itemize{
#'       \item `protein_<database>` annotation tables\cr
#'       \item `protein_<database>_count` genome-by-annotation count tables
#'       \item The default databases are `Pfam`, `COG`, `AMRFinder`, and `DefenseCas`.
#'     }
#'   \item **Metadata cleaning + Parquet export** via [cleanData()] -> writes
#'         Parquet files to `output_path`, and builds a **Parquet-backed DuckDB**
#'         (`*_parquet.duckdb`) with views over those Parquets.
#' }
#'
#' @param duckdb_path Character. Path to the **per-selection DuckDB** produced by
#'   [prepareGenomes()] (e.g., `"data/<Bug>/<Abbrev>.duckdb"`). This DB must
#'   already contain the tables written by [prepareGenomes()] and the upstream
#'   genome-processing steps.
#' @param output_path Character or `NULL`. Base directory for writing Panaroo,
#'   CD-HIT, HMMER, and final Parquet outputs. If `NULL`, defaults to
#'   `dirname(duckdb_path)`.
#'
#' @param threads Integer. Shared concurrency budget used across Panaroo, CD-HIT,
#'   and HMMER. Defaults to `8`.
#'
#' @param panaroo_split_jobs Logical. If `TRUE`, Panaroo runs in multiple batches
#'   that can be merged by [.mergePanaroo()]. If `FALSE`, Panaroo runs once on all
#'   isolates. Default: `FALSE`.
#' @param panaroo_core_threshold Numeric. Panaroo `--core_threshold`. Default: `0.90`.
#' @param panaroo_len_dif_percent Numeric. Panaroo `--len_dif_percent`. Default: `0.95`.
#' @param panaroo_cluster_threshold Numeric. Panaroo `--threshold`. Default: `0.95`.
#' @param panaroo_family_seq_identity Numeric. Panaroo `-f` gene family identity.
#'   Default: `0.5`.
#' @param panaroo_refind_mode Character. Panaroo's `--refind-mode` (`"off"`,
#'   `"default"`, or `"strict"`). See [.processPanaroo()] for the runtime caveat
#'   behind refinding. Default: `"off"`.
#' @param panaroo_strip_pseudogenes Logical. If `TRUE`, remove pseudogene feature
#'   records from Panaroo input GFF files before running Panaroo. Default: `FALSE`.
#' @param panaroo_pseudogene_clean_dir Character. Directory name for cleaned GFF
#'   files. Default: `"gff_clean"`.
#' @param panaroo_write_pseudogene_audit Logical. If `TRUE`, write a pseudogene
#'   cleaning audit file. Default: `TRUE`.
#'
#' @param cdhit_identity Numeric. CD-HIT `-c` identity threshold. Default: `0.9`.
#' @param cdhit_word_length Integer. CD-HIT `-n` word length. Default: `5`.
#' @param cdhit_memory Integer. CD-HIT `-M` memory limit in MB. Use `0` for
#'   unlimited. Default: `0`.
#' @param cdhit_extra_args Character vector. Extra arguments forwarded to
#'   `cd-hit`. Default: `c("-g", "1")`.
#' @param cdhit_output_prefix Character. Prefix for CD-HIT output files.
#'   Default: `"cdhit_out"`.
#'
#' @param hmmer_databases Character vector. HMMER annotation databases to run.
#'   Default: `c("Pfam", "COG", "AMRFinder", "DefenseCas")`.
#' @param hmmer_db_dir Character or `NULL`. Directory containing the shared HMMER
#'   database cache. If `NULL`, uses the amRdata user cache.
#' @param hmmer_docker_image Character. Docker image containing HMMER.
#'   Default: `"staphb/hmmer"`.
#' @param hmmer_num_splits Integer. Number of protein-sequence chunks for HMMER.
#'   Default: `8`.
#' @param hmmer_workers Integer. Number of parallel HMMER workers. Default: `8`.
#'
#' @param ref_file_path Character. Directory containing reference TSVs used by
#'   [cleanMetaData()] and [cleanData()] for metadata harmonization.
#'   Default: `"data_raw/"`.
#' @param verbose Logical. Print progress messages. Default: `TRUE`.
#'
#' @return
#' Invisibly returns a list with:
#' \itemize{
#'   \item `duckdb_path` - input DuckDB path
#'   \item `panaroo_output` - path to the selected Panaroo output directory used for import
#'   \item `parquet_duckdb_path` - absolute path to the created Parquet-backed DuckDB
#' }
#'
#' @details
#' **Docker & Platform Notes**
#' * Panaroo, CD-HIT, and HMMER run inside Docker containers.
#' * HMMER databases are stored separately from individual bug directories and
#'   are reused across datasets unless a custom `hmmer_db_dir` is supplied.
#' * Ensure Docker Desktop is running and has sufficient memory and CPU resources.
#'
#' **Input Requirements**
#' * `duckdb_path` must reference a per-selection DuckDB containing the genome
#'   file table, filtered genome selection, and BV-BRC metadata produced by the
#'   upstream curation workflow.
#'
#' **Outputs & Side Effects**
#' * Writes tool-specific intermediate outputs under `output_path`.
#' * Writes feature and metadata Parquet files under `output_path`.
#' * Creates a new Parquet-backed DuckDB (`*_parquet.duckdb`) with read-only views
#'   over the generated Parquet files.
#' * Records processing parameters, software versions, database selections, and
#'   other provenance information in the dataset manifest.
#'
#' **Threading**
#' * `threads` provides the shared CPU budget for the major processing stages.
#' * Panaroo, CD-HIT, and HMMER allocate that budget according to their respective
#'   stage parameters.
#'
#' @seealso
#' [prepareGenomes()], [runPanaroo2Duckdb()], [CDHIT2duckdb()], [cleanMetaData()],
#' [cleanData()]
#'
#' @examples
#' \dontrun{
#' runDataProcessing(
#'   duckdb_path   = "data/Shigella_flexneri/Sfl.duckdb",
#'   output_path   = "data/Shigella_flexneri",
#'   threads       = 8,
#'   ref_file_path = "data_raw/"
#' )
#'
#' # After completion:
#' # data/Shigella_flexneri/Sfl_parquet.duckdb
#' # will contain views over the Parquet files for downstream ML.
#' }
#'
#' @export
runDataProcessing <- function(
    duckdb_path,
    output_path = NULL,
    threads = 8,

    # Panaroo
    panaroo_split_jobs = FALSE,
    panaroo_core_threshold = 0.90,
    panaroo_len_dif_percent = 0.95,
    panaroo_cluster_threshold = 0.95,
    panaroo_family_seq_identity = 0.5,
    panaroo_refind_mode = c("off", "default", "strict"),
    panaroo_strip_pseudogenes = FALSE,
    panaroo_pseudogene_clean_dir = "gff_clean",
    panaroo_write_pseudogene_audit = TRUE,

    # CD-HIT
    cdhit_identity = 0.9,
    cdhit_word_length = 5,
    cdhit_memory = 0,
    cdhit_extra_args = c("-g", "1"),
    cdhit_output_prefix = "cdhit_out",

    # HMMER
    hmmer_databases = c(
      "Pfam",
      "COG",
      "AMRFinder",
      "DefenseCas"
    ),
    hmmer_db_dir = NULL,
    hmmer_docker_image = "staphb/hmmer",
    hmmer_num_splits = 8L,
    hmmer_workers = 8L,

    # Metadata cleaning
    ref_file_path = "data_raw/",
    verbose = TRUE
) {
  panaroo_refind_mode <- match.arg(panaroo_refind_mode)
  duckdb_path <- normalizePath(duckdb_path)
  out_dir <- if (is.null(output_path)) dirname(duckdb_path) else normalizePath(output_path)

  # Find the latest manifest
  manifest_path <- .manifest_find_latest(duckdb_path)

  if (is.null(manifest_path)) {
    stop(
      "No provenance manifest found for: ",
      duckdb_path,
      "\nRun prepareGenomes() first or provide a dataset with an existing manifest."
    )
  }

  # Append a new processing run to the existing manifest
  manifest <- .manifest_resume(
    manifest_path = manifest_path,
    base_dir = dirname(dirname(dirname(duckdb_path))),
    hash_files = FALSE
  )

  run_failed <- TRUE

  on.exit(
    if (run_failed) {
      .manifest_finish(
        manifest,
        status = "failed",
        error = "runDataProcessing() exited before successful completion."
      )
    },
    add = TRUE
  )

  # Record the start of this processing run
  manifest <- .manifest_event(
    manifest,
    message = "Started data-processing run.",
    details = list(
      duckdb_path = duckdb_path,
      output_path = out_dir
    )
  )

  # 1) Panaroo (run + optional merge) -> write Panaroo tables
  if (isTRUE(verbose)) message("Running Panaroo and writing gene & struct tables to DuckDB.")

  # Log!
  manifest <- .manifest_stage(
    manifest,
    name = "panaroo",
    status = "running",
    parameters = list(
      core_threshold = panaroo_core_threshold,
      len_dif_percent = panaroo_len_dif_percent,
      cluster_threshold = panaroo_cluster_threshold,
      family_seq_identity = panaroo_family_seq_identity,
      threads = threads,
      split_jobs = panaroo_split_jobs,
      refind_mode = panaroo_refind_mode,
      strip_pseudogenes = panaroo_strip_pseudogenes
    ),
    inputs = duckdb_path,
    tool = list(
      name = "Panaroo",
      docker_image = "staphb/panaroo:1.7.0"
    )
  )

  pan_dir <- runPanaroo2Duckdb(
    duckdb_path            = duckdb_path,
    output_path            = out_dir,
    core_threshold         = panaroo_core_threshold,
    len_dif_percent        = panaroo_len_dif_percent,
    cluster_threshold      = panaroo_cluster_threshold,
    family_seq_identity    = panaroo_family_seq_identity,
    threads                = threads,
    split_jobs             = panaroo_split_jobs,
    refind_mode            = panaroo_refind_mode,
    strip_pseudogenes      = panaroo_strip_pseudogenes,
    pseudogene_clean_dir   = panaroo_pseudogene_clean_dir,
    write_pseudogene_audit = panaroo_write_pseudogene_audit,
    verbose                = verbose
  )

  manifest <- .manifest_stage(
    manifest,
    name = "panaroo",
    status = "success",
    parameters = list(
      core_threshold = panaroo_core_threshold,
      len_dif_percent = panaroo_len_dif_percent,
      cluster_threshold = panaroo_cluster_threshold,
      family_seq_identity = panaroo_family_seq_identity,
      threads = threads,
      split_jobs = panaroo_split_jobs,
      refind_mode = panaroo_refind_mode,
      strip_pseudogenes = panaroo_strip_pseudogenes
    ),
    inputs = duckdb_path,
    outputs = c(
      pan_dir,
      duckdb_path
    ),
    tool = list(
      name = "Panaroo",
      version = "1.7.0",
      docker_image = "staphb/panaroo:1.7.0"
    )
  )

  # 2) CD-HIT -> write `protein` tables
  if (isTRUE(verbose)) message("Running CD-HIT and writing protein tables to DuckDB.")

  # Log!
  manifest <- .manifest_stage(
    manifest,
    name = "cdhit",
    status = "running",
    parameters = list(
      identity = cdhit_identity,
      word_length = cdhit_word_length,
      memory = cdhit_memory,
      threads = threads,
      extra_args = cdhit_extra_args,
      output_prefix = cdhit_output_prefix
    ),
    inputs = duckdb_path,
    tool = list(
      name = "CD-HIT",
      version = "4.8.1",
      docker_image = "weizhongli1987/cdhit:4.8.1"
    )
  )

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

  manifest <- .manifest_stage(
    manifest,
    name = "cdhit",
    status = "success",
    parameters = list(
      identity = cdhit_identity,
      word_length = cdhit_word_length,
      memory = cdhit_memory,
      threads = threads,
      extra_args = cdhit_extra_args,
      output_prefix = cdhit_output_prefix
    ),
    inputs = duckdb_path,
    outputs = c(
      file.path(out_dir, paste0(cdhit_output_prefix, "_input.fa")),
      file.path(out_dir, paste0(cdhit_output_prefix, "_proteins")),
      file.path(duckdb_path)
    ),
    tool = list(
      name = "CD-HIT",
      version = "4.8.1",
      docker_image = "weizhongli1987/cdhit:4.8.1"
    )
  )

  # 3) HMMER -> write HMM-based match tables for desired databases
  if (isTRUE(verbose)) {
    message(
      "Running HMMER with databases: ",
      paste(hmmer_databases, collapse = ", ")
    )
  }

  hmmer_db_dir <- if (is.null(hmmer_db_dir)) {
    .defaultHmmerDbDir()
  } else {
    normalizePath(hmmer_db_dir, mustWork = FALSE)
  }

  dir.create(
    hmmer_db_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  manifest <- .manifest_stage(
    manifest,
    name = "hmmer",
    status = "running",
    parameters = list(
      databases = hmmer_databases,
      database_dir = hmmer_db_dir,
      docker_image = hmmer_docker_image,
      threads = threads,
      num_of_splits = hmmer_num_splits,
      workers = hmmer_workers
    ),
    inputs = duckdb_path,
    tool = list(
      name = "HMMER",
      version = .hmmer_version(hmmer_docker_image),
      docker_image = hmmer_docker_image
    )
  )

  generic_databases <- intersect(
    hmmer_databases,
    c("Pfam", "COG", "AMRFinder")
  )

  hmmer_result <- NULL
  defense_result <- NULL

  if (length(generic_databases)) {
    hmmer_result <- .runHMMER(
                              duckdb_path = duckdb_path,
                              output_path = out_dir,
                              threads = threads,
                              hmmer_db_dir = hmmer_db_dir,
                              databases = generic_databases,
                              docker_image = hmmer_docker_image,
                              num_of_splits = hmmer_num_splits,
                              n_workers = hmmer_workers,
                              verbose = verbose
                            )
                          }



  if ("DefenseCas" %in% hmmer_databases) {
    defense_result <- .defenseHMMER(
                                    defense_db_dir = if (is.null(hmmer_db_dir)) {
                                      .defaultHmmerDbDir()
                                    } else {
                                      file.path(hmmer_db_dir, "DefenseCas")
                                    },
                                    docker_image = hmmer_docker_image,
                                    duckdb_path = duckdb_path,
                                    output_path = out_dir,
                                    threads = threads,
                                    verbose = verbose
                                  )
  }

  expected_outputs <- file.path(
    out_dir,
    paste0("protein_", hmmer_databases, ".parquet")
  )

  missing_outputs <- expected_outputs[!file.exists(expected_outputs)]

  if (length(missing_outputs)) {
    stop(
      "HMMER did not produce all expected outputs:\n",
      paste(missing_outputs, collapse = "\n")
    )
  }

  con <- DBI::dbConnect(
    duckdb::duckdb(),
    duckdb_path
  )
  on.exit(
    DBI::dbDisconnect(con, shutdown = FALSE),
    add = TRUE
  )

  expected_tables <- paste0("protein_", hmmer_databases)

  missing_tables <- expected_tables[
    !vapply(
      expected_tables,
      DBI::dbExistsTable,
      logical(1),
      conn = con
    )
  ]

  if (length(missing_tables)) {
    stop(
      "HMMER did not produce all expected DuckDB tables:\n",
      paste(missing_tables, collapse = "\n")
    )
  }

  manifest <- .manifest_stage(
    manifest,
    name = "hmmer",
    status = "success",
    parameters = list(
      databases = hmmer_databases,
      database_dir = hmmer_db_dir,
      docker_image = hmmer_docker_image,
      threads = threads,
      num_of_splits = hmmer_num_splits,
      workers = hmmer_workers
    ),
    inputs = duckdb_path,
    outputs = c(
      purrr::map(
        hmmer_databases,
        ~ file.path(out_dir, paste0("protein_", .x, ".parquet"))
      ),
      duckdb_path
    ),
    metrics = list(
      annotation_tables = paste0(
        "protein_",
        hmmer_databases
      ),
      database_provenance = list(
        generic = if (!is.null(hmmer_result)) hmmer_result$databases else NULL,
        DefenseCas = if (!is.null(defense_result)) defense_result$databases else NULL
      )
    ),
    tool = list(
      name = "HMMER",
      docker_image = hmmer_docker_image
    )
  )

  # 4) Clean metadata and export Parquet + Parquet-backed DuckDB
  if (is.null(ref_file_path) || !nzchar(ref_file_path)) {
    stop("`ref_file_path` (directory with reference TSVs) must be provided to cleanData().")
  }
  if (isTRUE(verbose)) message("Cleaning metadata and exporting Parquet-backed views.")
  cleanMetaData(duckdb_path = duckdb_path, path = out_dir, ref_file_path = ref_file_path)
  cleanData(duckdb_path = duckdb_path, path = out_dir)

  parquet_duckdb_path <- paste0(
    stringr::str_split_i(duckdb_path, ".duckdb", i = 1),
    "_parquet.duckdb"
  )

  if (isTRUE(verbose)) {
    message("\n============================================")
    message("Completed data-processing workflow successfully.")
    message("Parquet-backed DuckDB created at:")
    message("  ", normalizePath(parquet_duckdb_path))
    message("\nYou can use the amRml package to train machine")
    message("learning models for AMR using this file path.")
    message("For example:")
    message("  runMLmodels(\"", normalizePath(parquet_duckdb_path), "\")")
    message("============================================\n")
  }

  # Log!
  manifest <- .manifest_stage(
    manifest,
    name = "clean_metadata_and_export",
    status = "success",
    parameters = list(
      reference_path = normalizePath(
        ref_file_path,
        mustWork = FALSE
      )
    ),
    inputs = c(
      duckdb_path,
      ref_file_path
    ),
    outputs = c(
      out_dir,
      parquet_duckdb_path
    ),
    metrics = list(
      parquet_duckdb = parquet_duckdb_path
    )
  )

  run_failed <- FALSE

  .manifest_finish(
    manifest,
    status = "success"
  )

  invisible(list(
    duckdb_path = duckdb_path,
    panaroo_output = pan_dir,
    parquet_duckdb_path = normalizePath(parquet_duckdb_path)
  ))
}

#' Export processed tables from DuckDB database
#'
#' Reads tables from the DuckDB database produced by the `runDataProcessing()` workflow
#' and exports them as CSV, TSV, Parquet, and/or XLSX. This is an optional step
#' that allows users to take their processed data outside our amR workflow for
#' use in their own custom analyses. This is not required to run `amRml`!
#'
#' @param duckdb_path Character. Path to the DuckDB database created by the
#'   workflow (for example, `Sar.duckdb`).
#' @param output_path Character or NULL. Directory for exports. Defaults to
#'   file.path(dirname(duckdb_path), "processed_exports").
#' @param amr_phenotype_mode Character. One of "separate" or "append".
#'   "separate" exports the AMR labels as a separate wide table.
#'   "append" joins those labels onto the main feature tables before export.
#' @param export_formats Character vector. Any of "csv", "tsv", "parquet", "xlsx".
#' @param tables Character vector or NULL. Tables to export. If NULL, exports the
#'   standard processed tables present in the database.
#' @param verbose Logical. If TRUE, prints progress messages.
#'
#' @return Invisibly returns a list containing the export path, table names, and mode.
#' @export
exportProcessedData <- function(duckdb_path,
                                output_path = NULL,
                                amr_phenotype_mode = c("separate", "append"),
                                export_formats = c("csv"),
                                export_sequences = FALSE,
                                tables = NULL,
                                export_tables = TRUE,
                                verbose = TRUE) {
  duckdb_path <- normalizePath(duckdb_path, mustWork = TRUE)

  if (length(amr_phenotype_mode) > 1L) {
    message("`amr_phenotype_mode` not specified; defaulting to 'separate'.")
  }
  amr_phenotype_mode <- match.arg(amr_phenotype_mode)

  export_formats <- unique(tolower(export_formats))
  export_formats[export_formats == "excel"] <- "xlsx"

  allowed_formats <- c("csv", "tsv", "xlsx", "parquet")
  unknown_formats <- setdiff(export_formats, allowed_formats)
  if (length(unknown_formats)) {
    stop("Unsupported export format(s): ", paste(unknown_formats, collapse = ", "))
  }

  if (isTRUE(export_tables) && !length(export_formats)) {
    stop("At least one export format must be supplied when export_tables = TRUE.")
  }

  warn_text_exports <- any(export_formats %in% c("csv", "tsv", "xlsx"))
  if (isTRUE(export_tables) && warn_text_exports) {
    message(
      "\nNote: CSV, TSV, and Excel exports are intended primarily for human readability.\n",
      "However, BV-BRC genome accession are differentiated by trailing zero values.\n",
      "Example: 1282.2280 is a different genome than 1282.228\n",
      "If reading these files into software like Excel, or even re-reading them into R,\n",
      "accession IDs can be read as 'numeric' and trailing zeroes dropped!\n",
      "For programmatic reuse, we suggest using Parquet format, or \n",
      "explicitly import accession ID columns as 'character', not 'numeric'.\n"
    )
  }

  if ("xlsx" %in% export_formats && !requireNamespace("writexl", quietly = TRUE)) {
    stop("Format 'xlsx' was requested but package 'writexl' is not available.")
  }

  if (is.null(output_path)) {
    output_path <- file.path(dirname(duckdb_path), "processed_exports")
  }
  output_path <- normalizePath(output_path, mustWork = FALSE)
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = duckdb_path, read_only = TRUE)
  DBI::dbExecute(
    con,
    sprintf(
      "SET file_search_path='%s'",
      dirname(normalizePath(duckdb_path))
    )
  )
  on.exit(try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)

  available_tables <- DBI::dbListTables(con)
  if (!length(available_tables)) {
    stop("No tables found in DuckDB: ", duckdb_path)
  }

  read_tbl <- function(tbl) {
    tibble::as_tibble(DBI::dbReadTable(con, tbl))
  }

  write_one <- function(df, stem) {
    df <- .preserve_export_id_text(df)

    if ("csv" %in% export_formats) {
      utils::write.table(
        df, file = file.path(output_path, paste0(stem, ".csv")),
        sep = ",", row.names = FALSE, col.names = TRUE,
        quote = TRUE, na = "", qmethod = "double", fileEncoding = "UTF-8"
      )
    }
    if ("tsv" %in% export_formats) {
      utils::write.table(
        df, file = file.path(output_path, paste0(stem, ".tsv")),
        sep = "\t", row.names = FALSE, col.names = TRUE,
        quote = TRUE, na = "", qmethod = "double", fileEncoding = "UTF-8"
      )
    }
    if ("parquet" %in% export_formats) {
      arrow::write_parquet(df, file.path(output_path, paste0(stem, ".parquet")))
    }
    if ("xlsx" %in% export_formats) {
      writexl::write_xlsx(list(data = df), file.path(output_path, paste0(stem, ".xlsx")))
    }
  }

  build_amr_wide <- function() {
    source_tbl <- if ("metadata" %in% available_tables) {
      "metadata"
    } else if ("amr_phenotype" %in% available_tables) {
      "amr_phenotype"
    } else {
      NULL
    }

    if (is.null(source_tbl)) {
      return(NULL)
    }

    md <- read_tbl(source_tbl)

    needed <- c("genome.genome_id", "genome_drug.antibiotic", "genome_drug.resistant_phenotype")
    if (!all(needed %in% names(md))) {
      return(NULL)
    }

    md |>
      dplyr::transmute(
        genome_id = as.character(`genome.genome_id`),
        antibiotic = as.character(`genome_drug.antibiotic`),
        phenotype = as.character(`genome_drug.resistant_phenotype`)
      ) |>
      dplyr::filter(!is.na(genome_id), !is.na(antibiotic), !is.na(phenotype)) |>
      dplyr::distinct() |>
      dplyr::group_by(genome_id, antibiotic) |>
      dplyr::summarise(
        phenotype = paste(sort(unique(phenotype)), collapse = ";"),
        .groups = "drop"
      ) |>
      tidyr::pivot_wider(
        names_from = antibiotic,
        values_from = phenotype,
        values_fill = NA_character_
      ) |>
      dplyr::arrange(genome_id)
  }

  table_specs <- list(
    gene_count = list(source = "gene_count", stem = "gene_count", appendable = TRUE),
    protein_count = list(source = "protein_count", stem = "protein_count", appendable = TRUE),
    domain_count = list(source = "domain_count", stem = "domain_count", appendable = TRUE),
    struct = list(source = "gene_struct", stem = "struct", appendable = TRUE),
    gene_names = list(source = "gene_names", stem = "gene_names", appendable = FALSE),
    protein_names = list(source = "protein_names", stem = "protein_names", appendable = FALSE),
    domain_names = list(source = "domain_names", stem = "domain_names", appendable = FALSE),
    metadata = list(source = "metadata", stem = "metadata", appendable = FALSE),
    genome_data = list(source = "genome_data", stem = "genome_data", appendable = FALSE),
    amr_phenotype_wide = list(source = NULL, stem = "amr_phenotype_wide", appendable = FALSE)
  )

  if (isTRUE(export_sequences)) {
    table_specs$gene_seqs <- list(source = "gene_ref_seq", stem = "gene_seqs", appendable = FALSE)
    table_specs$protein_seqs <- list(source = "protein_cluster_seq", stem = "protein_seqs", appendable = FALSE)
    table_specs$genome_gene_protein <- list(source = "genome_gene_protein", stem = "genome_gene_protein", appendable = FALSE)
  }

  if (is.null(tables)) {
    selected_keys <- c(
      "gene_count",
      "protein_count",
      "domain_count",
      "struct",
      "gene_names",
      "protein_names",
      "domain_names",
      "metadata",
      "genome_data",
      "amr_phenotype_wide"
    )
    if (isTRUE(export_sequences)) {
      selected_keys <- c(selected_keys, "gene_seqs", "protein_seqs", "genome_gene_protein")
    }
  } else {
    selected_keys <- intersect(as.character(tables), names(table_specs))
  }

  if (!length(selected_keys)) {
    stop("No requested tables were found in the DuckDB.")
  }

  # Ensuring we don't accidentally pass a NULL through silently
  phenotype_wide <- build_amr_wide()
  if (!is.null(phenotype_wide)) {
    phenotype_wide <- .preserve_export_id_text(phenotype_wide)
  }
  exported <- character(0)

  for (key in selected_keys) {
    spec <- table_specs[[key]]

    if (key == "amr_phenotype_wide") {
      if (is.null(phenotype_wide)) {
        if (isTRUE(verbose)) message("Skipping amr_phenotype_wide: no AMR source table found.")
        next
      }
      write_one(phenotype_wide, spec$stem)
      exported <- c(exported, spec$stem)
      if (isTRUE(verbose)) message("Exported: ", spec$stem)
      next
    }

    if (is.null(spec$source) || !(spec$source %in% available_tables)) {
      if (isTRUE(verbose)) message("Skipping missing table: ", key)
      next
    }

    df <- .preserve_export_id_text(read_tbl(spec$source))
    out_stem <- spec$stem

    if (identical(amr_phenotype_mode, "append") &&
        isTRUE(spec$appendable) &&
        !is.null(phenotype_wide) &&
        "genome_id" %in% names(df)) {
      df <- dplyr::left_join(df, phenotype_wide, by = "genome_id")
      df <- .preserve_export_id_text(df)
      out_stem <- paste0(spec$stem, "_with_phenotypes")
    }

    write_one(df, out_stem)
    exported <- c(exported, out_stem)

    if (isTRUE(verbose)) message("Exported: ", out_stem)
  }

  invisible(list(
    duckdb_path = duckdb_path,
    output_path = output_path,
    tables = exported,
    amr_phenotype_mode = amr_phenotype_mode,
    export_formats = export_formats,
    export_sequences = isTRUE(export_sequences)
  ))
}
