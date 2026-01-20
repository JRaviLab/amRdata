# Normalize Docker path (reuse from data_curation.R)
.docker_path <- function(p) gsub("\\\\", "/", normalizePath(p, mustWork = FALSE))

# Map host paths under mounted root to container path
.to_container <- function(x, host_root, container_root = "/work") {
  host_root_unix <- .docker_path(host_root)
  x_unix <- .docker_path(x)
  pattern <- paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\\\\\1", host_root_unix))
  sub(pattern, container_root, x_unix)
}

# Launch Panaroo to build a pangenome (per batch)
.processPanaroo <- function(batch_input,
                            output_path,
                            core_threshold,
                            len_dif_percent,
                            cluster_threshold,
                            family_seq_identity,
                            panaroo_threads_per_job) {
  output_path <- .docker_path(output_path)
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

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
  output_dir_cont      <- .to_container(output_dir_host,      host_root = mount_host, container_root = mount_cont)

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
    "--threshold",      as.character(cluster_threshold),
    "-f",               as.character(family_seq_identity),
    "-t",               as.character(panaroo_threads_per_job)
  )

  res <- tryCatch({
    system2("docker", args = cmd_args, stdout = TRUE, stderr = TRUE)
  }, error = function(e) e)

  if (inherits(res, "error")) {
    stop(sprintf("Docker/Panaroo failed to launch: %s", res$message))
  }

  # If Panaroo wrote an error but system2 didn't throw, scan output for clues
  if (length(res) && any(grepl("Traceback|Error|No such file|not found|failed", res, ignore.case = TRUE))) {
    message("Panaroo output:\n", paste(res, collapse = "\n"))
  }

  invisible(res)
}

# Safely set a temporary future plan and restore it on exit
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
#' @param threads Integer. Number of threads for Panaroo and parallel execution. Default `16`.
#' @param split_jobs Logical. If TRUE, split into multiple smaller pangenome
#'   generation jobs that can be merged by [.mergePanaroo()]. If FALSE, all isolates in one run.
#'
#' @return A list of results for each Panaroo batch (stdout/stderr lines).
#'
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
                        threads = 16,
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
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~ . != "NA")) |>
    dplyr::pull(panaroo_input)

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
  batch_size  <- if (isTRUE(split_jobs)) ceiling(total_lines / 5) else total_lines
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
.mergePanaroo <- function(input_path,
                         core_threshold = 0.90,
                         len_dif_percent = 0.95,
                         cluster_threshold = 0.95,
                         family_seq_identity = 0.5,
                         threads = 8){

  input_path <- .docker_path(input_path)

  merge_dir <- file.path(input_path, "merge_output")
  dir.create(merge_dir, recursive = TRUE, showWarnings = FALSE)

  all_dirs <- list.dirs(input_path, recursive = FALSE, full.names = TRUE)
  all_dirs <- all_dirs[grepl("^panaroo_out_", basename(all_dirs))]

  valid_dirs <- all_dirs[file.exists(file.path(all_dirs, "final_graph.gml"))]

  if (length(valid_dirs) > 1) {
    mount_host <- input_path
    mount_cont <- "/work"

    dir_string <- paste(.to_container(valid_dirs, host_root = mount_host, container_root = mount_cont),
                        collapse = " ")

    cmd_args <- c("run",
                  "--platform", "linux/amd64",
                  "--rm",
                  "-v", paste0(mount_host, ":", mount_cont),
                  "-w", mount_cont,
                  "staphb/panaroo:1.5.1",
                  "panaroo-merge",
                  "-d", dir_string,
                  "-o", file.path(mount_cont, "merge_output"),
                  "--merge_paralogs",
                  "--core_threshold", as.character(core_threshold),
                  "--len_dif_percent", as.character(len_dif_percent),
                  "--threshold",      as.character(cluster_threshold),
                  "-f",               as.character(family_seq_identity),
                  "-t",               as.character(threads)
    )

    system2("docker", args = cmd_args, stdout = TRUE, stderr = TRUE)
  } else {
    stop("No valid Panaroo batch directories found (need >= 2 with final_graph.gml).")
  }
}

#' Convert Panaroo gene presence/absence matrix to a per-genome gene count table
.panaroo2geneTable <- function(panaroo_output_path, duckdb_path){
  filepath <- file.path(normalizePath(panaroo_output_path), "gene_presence_absence.csv")
  duckdb_path <- normalizePath(duckdb_path)
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  gene_count <- read.table(filepath, sep=",", header=TRUE, fill=TRUE, quote="") |>
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

#' Extract gene name and annotations from Panaroo gene presence/absence matrix
.panaroo2geneNames <- function(panaroo_output_path, duckdb_path){
  filepath <- file.path(normalizePath(panaroo_output_path), "gene_presence_absence.csv")
  duckdb_path <- normalizePath(duckdb_path)
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  gene_names <- read.table(filepath, sep=",", header=TRUE, fill=TRUE, quote="") |>
    tibble::as_tibble() |>
    dplyr::select(c(Gene, Annotation))

  DBI::dbWriteTable(con, "gene_names", gene_names, overwrite = TRUE)
  gene_names
}

#' Convert Panaroo gene struct presence/absence matrix to a per-genome struct table
.panaroo2StructTable <- function(panaroo_output_path, duckdb_path){
  struct_filepath <- file.path(normalizePath(panaroo_output_path), "struct_presence_absence.Rtab")
  duckdb_path <- normalizePath(duckdb_path)
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  gene_struct <- read.table(struct_filepath, sep="\t", header=TRUE, fill=TRUE, quote="") |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(cols= -1) |>
    tidyr::pivot_wider(names_from = Gene, values_from = value) |>
    dplyr::rename("genome_id" = "name") |>
    dplyr::mutate(genome_id = stringr::str_replace_all(genome_id, c("^X" = "", "\\.PATRIC$" = "")))

  DBI::dbWriteTable(con, "gene_struct", gene_struct, overwrite = TRUE)
  gene_struct
}

#' Write Panaroo reference FASTA and long presence/absence to DuckDB tables
.panaroo2OtherTables <- function(panaroo_output_path, duckdb_path){
  panaroo_output_path <- normalizePath(panaroo_output_path)
  duckdb_path <- normalizePath(duckdb_path)
  fasta_filepath <- file.path(panaroo_output_path, "pan_genome_reference.fa")
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  gene_fasta <- Biostrings::readDNAStringSet(filepath = fasta_filepath)
  DBI::dbWriteTable(con, "gene_ref_seq",
                    tibble::tibble(name = names(gene_fasta),
                                   sequence = as.character(gene_fasta)),
                    overwrite = TRUE)

  readr::read_csv(file.path(panaroo_output_path, "gene_presence_absence.csv")) |>
    dplyr::select(-`Non-unique Gene name`) |>
    tidyr::pivot_longer(-c("Gene", "Annotation"),
                        names_to = "genome_ids",
                        values_to = "protein_ids") |>
    dplyr::mutate(genome_ids = gsub(".PATRIC", "", genome_ids)) |>
    dplyr::select(genome_ids, Gene, protein_ids) |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(protein_ids)) |>
    tidyr::separate_rows(protein_ids, sep = ";") |>
    dplyr::filter(!stringr::str_detect(protein_ids, "_pseudo")) |>
    DBI::dbWriteTable(conn = con, name = "genome_gene_protein", overwrite = TRUE)
}

#' Consolidate Panaroo outputs into DuckDB tables
.panaroo2duckdb <- function(panaroo_output_path, duckdb_path){
  panaroo_output_path <- normalizePath(panaroo_output_path)
  duckdb_path         <- normalizePath(duckdb_path)

  .panaroo2geneTable(panaroo_output_path, duckdb_path)
  .panaroo2geneNames(panaroo_output_path, duckdb_path)
  .panaroo2StructTable(panaroo_output_path, duckdb_path)
  .panaroo2OtherTables(panaroo_output_path, duckdb_path)
  invisible(TRUE)
}

#' Run CD-HIT (via Docker) on concatenated protein FASTA and return output paths
.runCDHIT <- function(duckdb_path,
                      output_path,
                      output_prefix = "cdhit_out",
                      identity = 0.9,
                      word_length = 5,
                      threads = 0,
                      memory = 0,
                      extra_args = c("-g", "1")) {

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
    "-o", .to_container(clustered_faa,   mount_host, mount_cont),
    "-c", as.character(identity),
    "-n", as.character(word_length),
    "-T", as.character(threads),
    "-M", as.character(memory),
    "-d", "0",
    extra_args
  )

  message("Running cd-hit via Docker...")
  output <- tryCatch({
    system2("docker", args = cmd_args, stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    stop("cd-hit execution failed: ", e$message)
  })

  if (!file.exists(clustered_faa)) {
    stop("cd-hit failed: output file not found. Check stderr:\n", paste(output, collapse = "\n"))
  }

  message("cd-hit completed successfully.")
  list(
    cdhit_input_faa = cdhit_input_faa,
    clustered_faa   = clustered_faa
  )
}

.runPanaroo2Duckdb <- function(duckdb_path,
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
    duckdb_path       = duckdb_path,
    output_path       = out_dir,
    core_threshold    = core_threshold,
    len_dif_percent   = len_dif_percent,
    cluster_threshold = cluster_threshold,
    family_seq_identity = family_seq_identity,
    threads           = threads,
    split_jobs        = split_jobs
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

#' Parse CD-HIT protein cluster output
.parseProteinClusters <- function(clustered_faa) {
  lines <- data.table::fread(paste0(clustered_faa, ".clstr"), sep = "\n", header = FALSE)$V1
  cluster_ids <- grep("^>Cluster", lines) # This IDs the distinct clusters
  cluster_map <- data.table::data.table()

  for (i in seq_along(cluster_ids)) {
    start <- cluster_ids[i] + 1
    end   <- if (i < length(cluster_ids)) cluster_ids[i + 1] - 1 else length(lines)
    cluster_lines <- lines[start:end]

    # This finds the reference cluster ID and names the cluster with it
    ref_line <- grep("\\*$", cluster_lines, value = TRUE)
    ref_id <- if (length(ref_line) > 0) {
      stringr::str_extract(ref_line, "fig\\|[0-9]+\\.[0-9]+\\.peg\\.[0-9]+")
    } else {
      paste0("Cluster_", i - 1)
    }

    # Pull genome IDs
    genome_matches <- stringr::str_match(cluster_lines,
                                         "fig\\|([0-9]+\\.[0-9]+)\\.peg\\.[0-9]+")[, 2]
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

#' Build genome-by-cluster count matrix
.buildProtMatrices <- function(cluster_map) {
  cluster_map[, count := 1]
  reshape2::dcast(cluster_map, genome_id ~ cluster, value.var = "count", fun.aggregate = sum, fill = 0)
}
# Back-compat wrapper (older external name)
buildMatrices <- function(cluster_map) .buildProtMatrices(cluster_map)

#' Get cluster names from CD-HIT output
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
      proteinID  = stringr::str_extract(value, "^fig\\|[0-9]+\\.[0-9]+\\.peg\\.[0-9]+"),
      locus_tag  = stringr::str_match(value, "peg\\.[0-9]+\\|([^\\s]+)")[, 2],
      proteinName= stringr::str_trim(stringr::str_match(value, "\\|[^\\s]+\\s+(.*?)\\s+\\[")[, 2])
    ) |>
    dplyr::select(-value)

  names_faa
}
# Back-compat wrapper (older external name)
clusterNames <- function(cluster_map, cluster_fasta) .clusterNames(cluster_map, cluster_fasta)

#' Cluster proteins with CD-HIT and write results to DuckDB
#' @export
CDHIT2duckdb <- function(duckdb_path,
                         output_path,
                         output_prefix = "cdhit_out",
                         identity = 0.9,
                         word_length = 5,
                         threads = 0,
                         memory = 0,
                         extra_args = c("-g", "1")){

  duckdb_path <- normalizePath(duckdb_path)
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  if (missing(output_path) || output_path %in% c(".", "results", "results/")) {
    output_path <- dirname(duckdb_path)  # e.g., ./results/<bug>
  }
  output_path <- normalizePath(output_path)

  cdhit_outputs <- .runCDHIT(duckdb_path,
                             output_path,
                             output_prefix = output_prefix,
                             identity = identity,
                             word_length = word_length,
                             threads = threads,
                             memory = memory,
                             extra_args = extra_args)

  cluster_map   <- .parseProteinClusters(cdhit_outputs$clustered_faa)
  cluster_count <- .buildProtMatrices(cluster_map)

  DBI::dbWriteTable(con, "protein_count", cluster_count, overwrite = TRUE)

  cluster_fasta <- cdhit_outputs$cdhit_input_faa
  cluster_name  <- .clusterNames(cluster_map, cluster_fasta)
  DBI::dbWriteTable(con, "protein_names", cluster_name, overwrite = TRUE)

  clustered_faa <- Biostrings::readAAStringSet(cdhit_outputs$clustered_faa)
  DBI::dbWriteTable(con, "protein_cluster_seq",
                    tibble::tibble(
                      name     = names(clustered_faa) |> stringr::str_extract("fig\\|[0-9]+\\.[0-9]+\\.peg\\.[0-9]+"),
                      sequence = as.character(clustered_faa)
                    ),
                    overwrite = TRUE)
  invisible(TRUE)
}


# Check if InterProScan data are ready for use
.checkInterProData <- function(
    version      = "5.76-107.0",
    dest_dir     = "inst/extdata/interpro",
    docker_image = sprintf("interpro/interproscan:%s", version),
    platform     = "linux/amd64",
    curl_bin     = "curl",
    verbose      = TRUE
) {
  msg <- function(...) if (verbose) message(sprintf(...))

  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  dest_dir <- normalizePath(dest_dir, mustWork = TRUE)

  root_dir <- file.path(dest_dir, sprintf("interproscan-%s", version))
  data_dir <- file.path(root_dir, "data")

  is_indexed <- function() {
    length(Sys.glob(file.path(data_dir, "pfam", "*", "*.h3m"))) > 0
  }

  if (is_indexed()) {
    msg("InterProScan data ready at: %s", data_dir)
    return(list(data_dir = normalizePath(data_dir), ready = TRUE))
  }

  tar_url  <- sprintf("http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/%s/alt/interproscan-data-%s.tar.gz", version, version)
  md5_url  <- paste0(tar_url, ".md5")
  tar_path <- file.path(dest_dir, basename(tar_url))
  md5_path <- paste0(tar_path, ".md5")

  if (!file.exists(tar_path)) {
    msg("Downloading InterProScan data bundle.")
    status_tar  <- system2(curl_bin, c("-L", "-o", tar_path, tar_url))
    status_md5  <- system2(curl_bin, c("-L", "-o", md5_path, md5_url))
    if (status_tar != 0 || status_md5 != 0) stop("Failed to download InterProScan data bundle.")
  }

  msg("Verifying MD5 checksum.")
  md5_expected <- sub("\\s+.*$", "", readLines(md5_path)[1])
  md5_actual   <- tools::md5sum(tar_path)[[1]]
  if (!identical(tolower(md5_expected), tolower(md5_actual)))
    stop("MD5 checksum mismatch for InterProScan data bundle.")

  msg("Extracting InterProScan data bundle.")
  utils::untar(tar_path, exdir = dest_dir, tar = "internal")

  msg("Running InterProScan indexing.")
  docker_args <- c("run", "--rm")
  if (!is.null(platform)) docker_args <- c(docker_args, "--platform", platform)

  bind_data <- gsub("\\\\", "/", normalizePath(data_dir, mustWork = TRUE))

  status_idx <- system2(
    "docker",
    c(
      docker_args,
      "-v", paste0(bind_data, ":/opt/interproscan/data"),
      "-w", "/opt/interproscan",
      docker_image,
      "python3", "setup.py", "-f", "interproscan.properties"
    )
  )

  if (status_idx != 0) {
    warning("InterProScan indexing completed with non-zero exit status.")
  }

  if (is_indexed()) {
    msg("InterProScan Pfam HMMs are indexed and ready: %s", data_dir)
    return(list(data_dir = normalizePath(data_dir), ready = TRUE))
  } else {
    warning("Pfam HMM indices not found. InterProScan may press at runtime, which may be slow the first time it is run.")
    return(list(data_dir = normalizePath(data_dir), ready = FALSE))
  }
}

# Helpers for reading InterPro output
.getDfIPRColNames <- function() {
  c("AccNum", "SeqMD5Digest", "SLength", "Analysis",
    "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
    "Status", "RunDate", "IPRAcc", "IPRDesc", "placeholder")
}
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
.readIPRscanTsv <- function(filepath) {
  readr::read_tsv(filepath,
                  col_types = .getDfIPRColTypes(),
                  col_names = .getDfIPRColNames())
}

# Process a chunk of sequences through InterProScan (Docker)
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
  path      <- .docker_path(path)
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
    "--input",  .to_container(temp_fasta_file, path, "/work"),
    "--cpu", as.character(threads),
    "-f", file_format,
    "--appl", appl_str,
    "-b", chunk_out_file_base_cont
  )

  status <- tryCatch({
    system2("docker", args = cmd_args, stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    stop(sprintf("InterProScan execution failed for chunk %d: %s", chunk_id, e$message))
  })

  out_tsv   <- paste0(chunk_out_file_base_host, ".tsv")
  out_tsvgz <- paste0(chunk_out_file_base_host, ".tsv.gz")

  if (file.exists(out_tsv)) {
    return(out_tsv)
  } else if (file.exists(out_tsvgz)) {
    return(out_tsvgz)
  } else {
    stop(sprintf("InterProScan produced no output for chunk %d. Checked: %s and %s.\nLast message:\n%s",
                 chunk_id, out_tsv, out_tsvgz, paste(status, collapse = "\n")))
  }
}

#' Derive protein domain presence/absence and counts via InterProScan and write to DuckDB
domainFromIPR <- function(duckdb_path,
                          path,
                          out_file_base = "iprscan",
                          appl = c("Pfam"),
                          ipr_version  = "5.76-107.0",
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
  ipr_info <- if (isTRUE(auto_prepare_data)) {
    .checkInterProData(version = ipr_version,
                       dest_dir = ipr_dest_dir,
                       docker_image = ipr_image,
                       platform = ipr_platform,
                       verbose = TRUE)
  } else {
    list(data_dir = file.path(ipr_dest_dir, sprintf("interproscan-%s", ipr_version), "data"),
         ready = NA)
  }
  ipr_data_path <- ipr_info$data_dir

  # Pull the image ONCE (outside the per-chunk function)
  try(suppressWarnings(system2("docker", args = c("pull", ipr_image))), silent = TRUE)

  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  sequences_df <- dplyr::tbl(con, "protein_cluster_seq") |> tibble::as_tibble()
  if (nrow(sequences_df) == 0L) {
    stop("No sequences found in 'protein_cluster_seq'. Please run CDHIT2duckdb() first.")
  }

  # Split (~5000 per chunk) but clamp to at least 1
  chunks <- split(sequences_df, pmax(1, ceiling(seq_along(sequences_df$name) / 5000)))

  # Balance overall concurrency: workers * cpu_per_container <= threads
  total_threads <- max(1L, threads)
  workers <- min(length(chunks), total_threads)
  cpu_per_container <- max(1L, floor(total_threads / workers))

  message(sprintf("InterPro: scheduling %d chunk worker(s) with %d CPU(s) per container (<= %d total).",
                  workers, cpu_per_container, total_threads))

  .with_future_plan(workers = workers)

  results <- future.apply::future_lapply(seq_along(chunks), function(i) {
    res <- try(
      .process_chunk(
        chunk         = chunks[[i]],
        path          = path,
        ipr_data_path = ipr_data_path,
        out_file_base = out_file_base,
        appl          = appl,
        chunk_id      = i,
        threads       = cpu_per_container,  # <= balanced here
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
  if (length(tsvs) == 0L) stop("InterProScan produced no usable outputs. Check Docker logs above.")

  df_iprscan <- do.call(rbind, lapply(tsvs, .readIPRscanTsv))

  DBI::dbWriteTable(con, "domain_names",
                    df_iprscan |>
                      dplyr::select(AccNum, DB.ID, SignDesc, IPRAcc, IPRDesc, StartLoc, StopLoc),
                    overwrite = TRUE)

  df_protein_domain_pa <- df_iprscan |>
    dplyr::select(AccNum, DB.ID, IPRAcc, placeholder) |>
    dplyr::mutate(domain_ID = stringr::str_glue("{DB.ID}_{IPRAcc}")) |>
    dplyr::distinct() |>
    dplyr::mutate(placeholder = stringr::str_replace_all(placeholder, "-", "1")) |>
    tidyr::pivot_wider(id_cols = AccNum, names_from = domain_ID, values_from = placeholder, values_fill = "0") |>
    dplyr::group_by(AccNum) |>
    dplyr::summarize(across(everything(), ~ ifelse(any(. == "1"), "1", "0")), .groups = "drop") |>
    dplyr::mutate(across(-AccNum, as.numeric))

  protein_filter <- dplyr::tbl(con, "protein_count") |>
    tibble::as_tibble()

  accs <- unique(df_protein_domain_pa$AccNum)
  accs_in_matrix <- intersect(accs, colnames(protein_filter))
  if (length(accs_in_matrix) == 0L) {
    stop("No InterPro accessions match protein_count columns.")
  }

  protein_filter <- protein_filter |>
    dplyr::select(genome_id, dplyr::all_of(accs_in_matrix))

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
cleanData <- function(duckdb_path, path, ref_file_path){
  duckdb_path  <- normalizePath(duckdb_path)
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
  drug_abbr  <- readr::read_tsv(file.path(ref_file_path, "drug_abbr.tsv"))
  class_abbr <- readr::read_tsv(file.path(ref_file_path, "class_abbr.tsv"))
  clean_countries <- readr::read_tsv(file.path(ref_file_path, "cleaned_bvbrc_countries.tsv")) |>
    dplyr::select("raw_entry", "clean_name", "short_name")|>
    dplyr::distinct()

  dplyr::tbl(con, "filtered") |>
    tibble::as_tibble() |>
    dplyr::select("genome_drug.genome_id", "genome_drug.antibiotic",
                  "genome_drug.genome_name", "genome_drug.laboratory_typing_method",
                  "genome_drug.resistant_phenotype", "genome_drug.taxon_id",
                  "genome_drug.pmid", "genome.collection_year",
                  "genome.isolation_country", "genome.host_common_name",
                  "genome.isolation_source", "genome.species") |>
    dplyr::left_join(clean_drug, by = c("genome_drug.antibiotic" = "original_drug")) |>
    dplyr::filter(!is.na(cleaned_drug)) |>
    dplyr::left_join(drug_class, by = c("cleaned_drug" = "drug")) |>
    dplyr::left_join(drug_abbr, by = c("cleaned_drug" = "drug")) |>
    dplyr::left_join(class_abbr, by = "drug_class") |>
    DBI::dbWriteTable(conn = con, name = "filtered", overwrite = TRUE)

  resistance_summary <- dplyr::tbl(con, "filtered") |>
    tibble::as_tibble()  |>
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
    dplyr::rename("cleaned_country"="clean_name", "country_abbr"="short_name") |>
    dplyr::mutate(genome.isolation_country = cleaned_country) |>
    dplyr::select(-cleaned_country) |>
    dplyr::left_join(resistance_summary, by = "genome_drug.genome_id") |>
    dplyr::mutate(resistant_classes = dplyr::case_when(
      is.na(resistant_classes) ~ genome_drug.resistant_phenotype,
      TRUE ~ resistant_classes
    )) |>
    dplyr::mutate(num_resistant_classes= dplyr::case_when(
      is.na(num_resistant_classes) ~ 0,
      TRUE ~ num_resistant_classes
    )) |>
    dplyr::mutate(genome.collection_year = as.numeric(genome.collection_year)) |>
    dplyr::mutate(year_bin = cut(genome.collection_year, breaks = year_breaks,
                                 right = FALSE, include.lowest = TRUE,
                                 labels = paste(year_breaks[-length(year_breaks)],
                                                year_breaks[-1] - 1, sep = "-"))) |>
    DBI::dbWriteTable(conn = con, name = "cleaned_metadata", overwrite = TRUE)

  # Parquet output paths
  genes_parquet                  <- file.path(path, "gene_count.parquet")
  gene_names_parquet             <- file.path(path, "gene_names.parquet")
  gene_ref_seq_parquet           <- file.path(path, "gene_seqs.parquet")
  genome_gene_protein_parquet    <- file.path(path, "genome_gene_protein.parquet")
  struct_parquet                 <- file.path(path, "struct.parquet")

  proteins_parquet               <- file.path(path, "protein_count.parquet")
  domains_parquet                <- file.path(path, "domain_count.parquet")

  metadata_parquet               <- file.path(path, "metadata.parquet")  # cleaned_metadata exported as 'metadata'

  domain_names_parquet           <- file.path(path, "domain_names.parquet")
  protein_names_parquet          <- file.path(path, "protein_names.parquet")

  protein_cluster_seq_parquet    <- file.path(path, "protein_seqs.parquet")

  # Also export AMR/genome/original metadata
  amr_phenotype_parquet          <- file.path(path, "amr_phenotype.parquet")
  genome_data_parquet            <- file.path(path, "genome_data.parquet")
  original_metadata_parquet      <- file.path(path, "original_metadata.parquet")

  writeCompressedParquet <- function(df, path) {
    arrow::write_parquet(
      df,
      path,
      compression = "zstd",
      compression_level = 9,
      use_dictionary = TRUE
    )
  }

  db_name <- duckdb_path |> stringr::str_split_i(".duckdb", i = 1) |> paste0("_parquet.duckdb")
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
  DBI::dbReadTable(con, "genome_data")   |> writeCompressedParquet(genome_data_parquet)
  DBI::dbReadTable(con, "metadata")      |> writeCompressedParquet(original_metadata_parquet)

  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW amr_phenotype AS SELECT * FROM read_parquet('%s')", amr_phenotype_parquet))
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW genome_data AS SELECT * FROM read_parquet('%s')", genome_data_parquet))
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW original_metadata AS SELECT * FROM read_parquet('%s')", original_metadata_parquet))

  invisible(TRUE)
}


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
                              cdhit_extra_args = c("-g","1"),
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
  pan_dir <- .runPanaroo2Duckdb(
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
    message("\nYou can use the amR_ml package to train machine")
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

