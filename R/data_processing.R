# Normalize Docker path (reuse from data_curation.R)
.docker_path <- function(p) gsub("\\\\", "/", normalizePath(p, mustWork = FALSE))

# Map host paths under mounted root to container path
.to_container <- function(x, host_root, container_root = "/work") {
  host_root_unix <- .docker_path(host_root)
  x_unix <- .docker_path(x)
  pattern <- paste0("^", gsub("([\\^$.|?*+(){}\\[\\]\\\\])", "\\\\\\1", host_root_unix))
  sub(pattern, container_root, x_unix)
}


# Launch Panaroo to build a pangenome
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

  # Write the genome list file
  genome_filepath_host <- tempfile(pattern = "genomeFilepath_", fileext = ".txt", tmpdir = output_path)

  # Rewrite each "gff fna" pair to container-visible paths
  batch_input_cont <- vapply(unlist(batch_input), function(line) {
    parts <- strsplit(line, " +")[[1]]  # robust to multiple spaces
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

  # Run command
  cmd_args <- c(
    "run",
    "--platform", "linux/amd64",
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


  res <- tryCatch({
    system2("docker", args = cmd_args, stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    e  # return condition for handling below
  })

  if (inherits(res, "error")) {
    stop(sprintf("Docker/Panaroo failed to launch: %s", res$message))
  }

  # If Panaroo wrote an error but system2 didn't throw, scan output for clues
  if (length(res) && any(grepl("Traceback|Error|No such file|not found|failed", res, ignore.case = TRUE))) {
    message("Panaroo output:\n", paste(res, collapse = "\n"))
  }
}

#' Run Panaroo for Pangenome Analysis in Parallel Batches
#'
#' Executes Panaroo inside a Singularity/Apptainer container on genome annotation
#' files prepared by [genomeList()]. The function can optionally split input genomes
#' into batches, runs Panaroo with strict cleaning and clustering options, and
#' returns the results of each batch execution.
#'
#' @param duckdb_path A path to the DuckDB database containing the `"files"` table.
#' @param output_path Character scalar. Base directory for Panaroo outputs and temporary files.
#' @param core_threshold Numeric. Core genome threshold for Panaroo (`--core_threshold`).
#'   Defaults to `0.90`.
#' @param len_dif_percent Numeric. Length difference percentage for clustering
#'   (`--len_dif_percent`). Defaults to `0.95`.
#' @param cluster_threshold Numeric. Sequence identity threshold for clustering
#'   (`--threshold`). Defaults to `0.95`.
#' @param family_seq_identity Numeric. Sequence identity for gene family clustering
#'   (`-f`). Defaults to `0.5`.
#' @param threads Integer. Number of threads for Panaroo and parallel execution.
#'   Defaults to `16`.
#' @param split_jobs Logical. If TRUE, splits into multiple smaller pangenome
#' generation jobs that can be merged by [panarooMerge()]. If FALSE, all isolates
#' will be used in a single pangenome generation step, and merging is not required,
#' but this may take longer for very large numbers of isolates (>5000).
#'
#' @return A list of results from Panaroo batch executions. Each element contains
#'   the stdout/stderr output from the corresponding `system2()` call.
#'
#' @details
#' - Panaroo is run in strict cleaning mode with options:
#'   `--clean-mode strict`, `--merge_paralogs`, `--remove-invalid-genes`.
#' - Temporary genome file lists are created in `output_path` for each batch.
#' - Output directories are named `panaroo_out_<PID>` under `output_path`.
#' - Singularity is invoked with `--bind path:path` to ensure file accessibility.
#'
#' @examples
#' \dontrun{
#' # Run Panaroo on genomes listed in DuckDB
#' res <- runPanaroo(duckdb_path = "results/###.duckdb",
#'   output_path = "results/",
#'   core_threshold = 0.90,
#'   len_dif_percent = 0.95,
#'   cluster_threshold = 0.95,
#'   family_seq_identity = 0.5,
#'   threads = 16,
#'   split_jobs = FALSE)
#' )
runPanaroo <- function(duckdb_path = "results/",
                       output_path = "results/",
                       core_threshold = 0.90,
                       len_dif_percent = 0.95,
                       cluster_threshold = 0.95,
                       family_seq_identity = 0.5,
                       threads = 16,
                       split_jobs = FALSE) {

  # Read genome input table from genomeList()
  duckdb_path <- normalizePath(duckdb_path)
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  # Default Panaroo output to per-bug results directory if generic path supplied
  if (missing(output_path) || output_path %in% c(".", "results", "results/")) {
    output_path <- dirname(duckdb_path)  # e.g., ./results/<bug>
  }
  output_path <- normalizePath(output_path)

  genome_query_output <- DBI::dbReadTable(con, "files")

  # Extract and filter panaroo input paths
  panaroo_input_files <- genome_query_output |>
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~ . != "NA")) |>
    dplyr::pull(panaroo_input)

  split_files <- strsplit(panaroo_input_files, " ")

  # Parallel filtering of valid .gff files
  future::plan(future::multisession, workers = threads)
  valid_entries <- furrr::future_map(split_files, function(paths) {
    gff_file <- paths[1]
    if (file.exists(gff_file)) {
      length(readLines(gff_file, n = 5, warn = FALSE)) >= 5
    } else {
      FALSE
    }
  })

  filtered_panaroo_input <- sapply(split_files[unlist(valid_entries)], paste, collapse = " ")

  # Split into batches for split jobs
  total_lines <- length(filtered_panaroo_input)
  batch_size  <- if (isTRUE(split_jobs)) ceiling(total_lines / 5) else total_lines
  panaroo_batches <- split(filtered_panaroo_input, ceiling(seq_along(filtered_panaroo_input) / batch_size))

  # Per-job thread logistics
  n_jobs <- length(panaroo_batches)
  if (n_jobs == 0L) {
    warning("Panaroo inputs do not exist after filtering. Check your upstream processing.")
    return(invisible(list()))
  }

  # Balance CPUs across jobs
  panaroo_threads_per_job <- max(1L, floor(threads / n_jobs))

  # Reset the plan: one worker per batch
  future::plan(future::multisession, workers = n_jobs)

  # Run one worker per batch and pass all needed params explicitly
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

  return(invisible(batch_panaroo_run))
}

#' Merge multiple Panaroo batch outputs into a single pangenome result
#'
#' Finds Panaroo batch output directories under `input_path` (matching `panaroo_out_*`)
#' that contain a `final_graph.gml`, and merges them using `panaroo-merge`
#' inside a Singularity/Apptainer container. The merged output is written to
#' `input_path/merge_output`.
#'
#' @param input_path Character scalar. Base directory that contains Panaroo batch outputs
#'   named like `panaroo_out_<PID>` (created by `runPanaroo()`) and where the merged
#'   outputs will be written under `merge_output/`.
#' @param core_threshold Numeric. Core genome threshold passed to `--core_threshold`.
#'   Defaults to `0.90`.
#' @param len_dif_percent Numeric. Length difference percentage passed to
#'   `--len_dif_percent`. Defaults to `0.95`.
#' @param cluster_threshold Numeric. Sequence identity threshold for clustering
#'   passed to `--threshold`. Defaults to `0.95`.
#' @param family_seq_identity Numeric. Sequence identity for gene family clustering
#'   passed to `-f`. Defaults to `0.5`.
#' @param threads Integer. Number of threads for `panaroo-merge`. Defaults to `16`.
#'
#' @return
#'   (running `panaroo-merge` and writing merged outputs into `merge_output/`).
#'
#' @details
#' - Only directories containing `final_graph.gml` are considered valid inputs
#'   for merging. At least two valid directories are required.
#'
#' @examples
#' \dontrun{
#' # After running runPanaroo() which creates panaroo_out_* directories:
#' mergePanaroo(
#'   input_path = "results/",
#'   threads = 8
#' )
mergePanaroo <- function(input_path,
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
                  "-v", paste0(mount_host, ":", mount_cont),
                  "-w", mount_cont,
                  "staphb/panaroo:1.5.1",
                  "panaroo-merge",
                  "-d", dir_string,
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
    message("No valid Panaroo batch directories found (need >= 2 with final_graph.gml!")
    quit(status = 1)
  }
}

#' Convert Panaroo gene presence/absence matrix to a per-genome gene count table
#'
#' Reads `gene_presence_absence.csv` (from `runPanaroo()` or `mergePanaroo()` outputs),
#' reshapes it into a genome-by-gene count matrix, and writes the result to DuckDB
#' as a table named `"gene_count"`. Returns the reshaped tibble for immediate use.
#' - Converts presence strings to **counts** by splitting semicolon-delimited entries:
#'   empty -> `0`, `"A;B;C"` -> `3`.
#'
#' @param panaroo_output_path Character scalar. Path to Panaroo output directory;
#' @param duckdb_path Character scalar. Path to the DuckDB database file to use
#'
#' @return A tibble with:
#' \describe{
#'   \item{genome_id}{Character column of normalized genome IDs.}
#'   \item{<gene columns>}{Integer counts per gene (0 if absent, >=1 if present,
#'   counting semicolon-separated entries).}
#' }
#'  (writing table into DuckDB).
#'
#' @examples
#' \dontrun{
#' # Convert Panaroo presence/absence to per-genome counts and store in DuckDB
#' gene_count <- panaroo2geneTable(
#'   panaroo_output_path = "results/merge_output/",
#'   duckdb_path = file.path(".", paste0(generateDBname(c("487000", "Campylobacter coli")), ".duckdb"))
#' )
.panaroo2geneTable <- function(panaroo_output_path, duckdb_path){
  filepath <- file.path(normalizePath(panaroo_output_path), "gene_presence_absence.csv")
  duckdb_path <- normalizePath(duckdb_path)
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)

  gene_count <- read.table(filepath, sep=",", header=T, fill=T, quote="") |>
    tibble::as_tibble() |>
    dplyr::select(-c(Non.unique.Gene.name, Annotation)) |>
    tidyr::pivot_longer(cols= -1) |>
    tidyr::pivot_wider(names_from = Gene,values_from = value) |>
    dplyr::rename("genome_id" = "name") |>
    dplyr::mutate(genome_id = stringr::str_replace_all(genome_id, c("^X" = "", "\\.PATRIC$" = ""))) |>
    dplyr::mutate(across(-genome_id, ~ ifelse(. == "", 0, stringr::str_count(., ";") + 1)))

  DBI::dbWriteTable(con, "gene_count", gene_count, overwrite = T)

  # dbDisconnect(con, shutdown = TRUE)
  return(gene_count)

}

#' Extract gene name and annotations from Panaroo gene presence/absence matrix
#'
#' Reads `gene_presence_absence.csv` (from `runPanaroo()` or `mergePanaroo()` outputs),
#' extract gene names and annotation columns, and writes the result to DuckDB
#' as a table named `"gene_names"`. Returns the tibble for immediate use.
#'
#' @param panaroo_output_path Character scalar. Path to Panaroo output directory;
#' @param duckdb_path Character scalar. Path to the DuckDB database file to open
#'
#' @return A tibble with:
#' \describe{
#'   \item{Gene}{Character column of gene families}
#'   \item{Annotation}{Character column of free text annotation of the genes}
#' }
#' (writing table into DuckDB).
#'
#' @examples
#' \dontrun{
#' # Convert Panaroo presence/absence to per-genome counts and store in DuckDB
#' gene_names <- panaroo2geneNames(
#'   panaroo_output_path = "results/merge_output/",
#'   duckdb_path = file.path(".", paste0(generateDBname(c("487000", "Campylobacter coli")), ".duckdb"))
#' )

.panaroo2geneNames <- function(panaroo_output_path, duckdb_path){
  filepath <- file.path(normalizePath(panaroo_output_path), "gene_presence_absence.csv")
  duckdb_path <- normalizePath(duckdb_path)
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)

  gene_names <- read.table(filepath, sep=",", header=T, fill=T, quote="") |>
    tibble::as_tibble() |>
    dplyr::select(c(Gene, Annotation))

  DBI::dbWriteTable(con, "gene_names", gene_names, overwrite = T)

  return(gene_names)
}

#' Convert Panaroo gene struct presence/absence matrix to a per-genome struct table
#'
#' Reads `struct_presence_absence.Rtab` (from `runPanaroo()` or `mergePanaroo()` outputs),
#' reshapes it into a genome-by-struct matrix, and writes the result to DuckDB
#' as a table named `"gene_struct"`. Returns the reshaped tibble for immediate use.
#'
#' @param panaroo_output_path Character scalar. Path to Panaroo output directory
#' `struct_presence_absence.Rtab`.
#' @param duckdb_path Character scalar. Path to the DuckDB database file to open.
#'
#' @return A tibble with:
#' \describe{
#'   \item{genome_id}{Character column of normalized genome IDs.}
#'   \item{<gene struct columns>}{Integer counts per gene (0 if absent, 1 if present).}
#' }
#'  (writing table into DuckDB).
#'
#' @examples
#' \dontrun{
#' # Convert Panaroo presence/absence to per-genome counts and store in DuckDB
#' gene_struct <- panaroo2StructTable(
#'   panaroo_output_path = "results/merge_output/",
#'   duckdb_path = file.path(".", paste0(generateDBname(c("487000", "Campylobacter coli")), ".duckdb"))
#' )

.panaroo2StructTable <- function(panaroo_output_path, duckdb_path){
  struct_filepath <- file.path(normalizePath(panaroo_output_path), "struct_presence_absence.Rtab")
  duckdb_path <- normalizePath(duckdb_path)
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)

  gene_struct <- read.table(struct_filepath,
                            sep="\t", header=T, fill=T, quote="") |>
    tibble::as_tibble() |>
    tidyr::pivot_longer(cols= -1) |>
    tidyr::pivot_wider(names_from = Gene,values_from = value) |>
    dplyr::rename("genome_id" = "name") |>
    dplyr::mutate(genome_id = stringr::str_replace_all(genome_id, c("^X" = "", "\\.PATRIC$" = "")))

  DBI::dbWriteTable(con, "gene_struct", gene_struct, overwrite = T)

  return(gene_struct)
}

#' Write Panaroo reference FASTA and long presence/absence to DuckDB tables
#'
#' Reads Panaroo outputs produced by [runPanaroo()] or [mergePanaroo()] and
#' writes two tables to DuckDB:
#' \itemize{
#'   \item \code{gene_ref_seq}: Reference gene sequences from
#'         \code{pan_genome_reference.fa} (columns: \code{name}, \code{sequence}).
#'   \item \code{genome_gene_protein}: Long-format mapping of genomes to Panaroo
#'         gene families and their protein IDs, derived from
#'         \code{gene_presence_absence.csv}.
#' }
#'
#' @param panaroo_output_path [character] Base directory containing Panaroo outputs.
#' Normalized via \code{normalizePath()}.
#' @param duckdb_path [character] Path to the DuckDB database file to open
#'   (e.g., created earlier via [generateDBname()]).
#'
#' @return
#'   (writing two tables into DuckDB).
#'
#' @examples
#' \dontrun{
#' # Write Panaroo reference and long presence/absence to DuckDB
#' panaroo2OtherTables(
#'   panaroo_output_path = "/results/merge_output/",
#'   duckdb_path = file.path(".", paste0(generateDBname(c("487000", "Campylobacter coli")), ".duckdb"))
#' )

.panaroo2OtherTables <- function(panaroo_output_path, duckdb_path){
  panaroo_output_path <- normalizePath(panaroo_output_path)
  duckdb_path <- normalizePath(duckdb_path)
  fasta_filepath <- file.path(panaroo_output_path, "pan_genome_reference.fa")
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)

  gene_fasta <- Biostrings::readDNAStringSet(filepath = fasta_filepath)
  DBI::dbWriteTable(con, "gene_ref_seq", tibble::tibble(name = names(gene_fasta),
                                                        sequence = as.character(gene_fasta)),
                    overwrite = T)

  gene_pa_long <- readr::read_csv(file.path(panaroo_output_path, "gene_presence_absence.csv")) |>
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
    DBI::dbWriteTable(conn = con, name = "genome_gene_protein", overwrite=T)

}

#' Consolidate Panaroo outputs into DuckDB tables
#'
#' Runs all Panaroo-to-DuckDB conversion steps in sequence:
#' \itemize{
#'   \item [panaroo2geneTable()]: Converts Panaroo gene presence/absence matrix
#'         into per-genome gene count table.
#'   \item [panaroo2geneNames()]: Writes gene name mapping table.
#'   \item [panaroo2StructTable()]: Writes gene structural annotation table.
#'   \item [panaroo2OtherTables()]: Writes reference FASTA and genome-gene-protein
#'         mapping tables.
#' }
#'
#' @param panaroo_output_path Character scalar. Base directory containing Panaroo outputs
#' @param duckdb_path Character scalar. Path to the DuckDB database file where all
#'   Panaroo tables will be written.
#'
#' @return
#'   (writing multiple tables into DuckDB). After execution, the database will contain:
#'   \itemize{
#'     \item `"gene_count"`: Genome-by-gene count matrix.
#'     \item `"gene_names"`: Panaroo gene name mapping.
#'     \item `"gene_struct"`: Gene struct matrix.
#'     \item `"gene_ref_seq"`: Reference gene sequences.
#'     \item `"genome_gene_protein"`: Genome-to-gene-to-protein mapping.
#'   }
#'
#' @examples
#' \dontrun{
#' # Consolidate Panaroo outputs into DuckDB
#' panaroo2duckdb(
#'   panaroo_output_path = "results/merge_output/",
#'   duckdb_path = file.path(".", paste0(generateDBname(c("487000", "Campylobacter coli")), ".duckdb"))
#' )
panaroo2duckdb <- function(panaroo_output_path, duckdb_path){

  panaroo_output_path <- normalizePath(panaroo_output_path)
  duckdb_path <- normalizePath(duckdb_path)
  gene_count <- .panaroo2geneTable(panaroo_output_path, duckdb_path)
  gene_names <- .panaroo2geneNames(panaroo_output_path, duckdb_path)
  gene_struct <- .panaroo2StructTable(panaroo_output_path, duckdb_path)
  .panaroo2OtherTables(panaroo_output_path, duckdb_path)

}


#' Run CD-HIT (via Docker) on concatenated protein FASTA and return output paths
#'
#' Concatenates protein FASTA files listed in the DuckDB table `"files"` (column
#' `faa_path`), runs CD-HIT in a Docker container, and returns paths to the
#' generated input and clustered output FASTA.
#'
#' The function:
#' - Connects to a DuckDB database and reads the table `"files"`.
#' - Filters rows to ensure no `"NA"` string values and pulls `faa_path`.
#' - Concatenates all `*.faa` files into a single input FASTA:
#'   `<output_path>/<output_prefix>_input.fa`.
#' - Runs CD-HIT via Docker with the specified parameters and writes clustered
#'   proteins to `<output_path>/<output_prefix>_proteins`.
#' - Verifies successful output and returns paths to input/output FASTA files.
#'
#' @param duckdb_path Character scalar. Path to the DuckDB database with a table
#'   named `"files"` containing at least the column `faa_path` with absolute or
#'   relative paths to protein FASTA files.
#' @param output_path Character scalar. Directory where concatenated input and the
#'   clustered output FASTA will be written. Created if it does not exist.
#' @param output_prefix Character scalar. Prefix for output files. Defaults to
#'   `"cdhit_out"`. The function will create:
#'   \itemize{
#'     \item \code{<output_path>/<output_prefix>_input.fa}
#'     \item \code{<output_path>/<output_prefix>_proteins}
#'   }
#' @param identity Numeric. Sequence identity threshold for clustering (CD-HIT `-c`).
#'   Defaults to `0.9`.
#' @param word_length Integer. Word length parameter (CD-HIT `-n`). Defaults to `5`.
#'   Must be compatible with `identity` according to CD-HIT constraints.
#' @param threads Integer. Number of threads (CD-HIT `-T`). `0` uses all available
#'   cores. Defaults to `0`.
#' @param memory Integer. Memory limit in MB (CD-HIT `-M`). `0` means unlimited.
#'   Defaults to `0`.
#' @param extra_args Character vector. Extra CD-HIT flags appended at the end
#'   (e.g., `c("-g","1")` to use the most similar sequences as cluster representatives).
#'
#' @return A list with:
#' \describe{
#'   \item{cdhit_input_faa}{Character scalar. Path to concatenated input FASTA.}
#'   \item{clustered_faa}{Character scalar. Path to clustered output FASTA (CD-HIT `-o`).}
#' }
#'
#' @examples
#' \dontrun{
#' # Run cd-hit using protein FASTAs referenced in a DuckDB "files" table
#' res <- .runCDHIT(
#'   duckdb_path   = "example.duckdb",
#'   output_path    = "trial",
#'   output_prefix = "cdhit",
#'   identity      = 0.9,
#'   word_length   = 5,
#'   threads       = 8,
#'   memory        = 0,
#'   extra_args    = c("-g", "1")
#' )
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


#' Parse CD-HIT protein cluster output
#'
#' Takes a CD-HIT output file (with the -d 0 parameter), finds the clusters,
#' extracts the lines that represent individual sequences in the cluster,
#' and designates which genomes have which clusters.
#'
#' @param clustered_faa Character. Path to the CD-HIT clustered FASTA output file.
#'
#' @return A data.table mapping cluster IDs to genome IDs.
#' @export
#'
#' @examples
#' \dontrun{
#' cluster_map <- .parseProteinClusters("path/to/cdhit_proteins")
#' }
.parseProteinClusters <- function(clustered_faa) {
  lines <- data.table::fread(paste0(clustered_faa, ".clstr"), sep = "\n", header = FALSE)$V1
  cluster_ids <- grep("^>Cluster", lines) # This IDs the distinct clusters
  cluster_map <- data.table::data.table()

  for (i in seq_along(cluster_ids)) {
    start <- cluster_ids[i] + 1
    end <- if (i < length(cluster_ids)) cluster_ids[i + 1] - 1 else length(lines)
    cluster_lines <- lines[start:end]

    # This finds the reference cluster ID and names the cluster with it
    ref_line <- grep("\\*$", cluster_lines, value = TRUE)
    ref_id <- if (length(ref_line) > 0) {
      stringr::str_extract(ref_line,
                           "fig\\|[0-9]+\\.[0-9]+\\.peg\\.[0-9]+")
    } else {
      paste0("Cluster_", i - 1)
    }

    # Pull genome IDs
    genome_matches <- stringr::str_match(cluster_lines,
                                         "fig\\|([0-9]+\\.[0-9]+)\\.peg\\.[0-9]+")[, 2]
    genome_matches <- genome_matches[!is.na(genome_matches)]

    if (length(genome_matches) > 0) {
      cluster_map <- data.table::rbindlist(list(cluster_map,
                                                data.table::data.table(cluster = ref_id,
                                                                       genome_id = genome_matches)), use.names = TRUE)
    }
  }

  return(cluster_map)
}


#' Build genome-by-cluster count matrix
#'
#' Takes a cluster map and builds a count matrix showing which genomes
#' have which protein clusters.
#'
#' @param cluster_map A data.table with columns `cluster` and `genome_id`.
#'
#' @return A data frame with genome IDs as rows and cluster IDs as columns.
#' @export
#'
#' @examples
#' \dontrun{
#' cluster_count <- .buildProtMatrices(cluster_map)
#' }
.buildProtMatrices <- function(cluster_map) {
  cluster_map[, count := 1]
  cluster_count <- reshape2::dcast(cluster_map, genome_id ~ cluster, value.var = "count", fun.aggregate = sum, fill = 0)

  return(cluster_count)
}

#' Get cluster names from CD-HIT output
#'
#' @param cluster_map Data frame mapping sequences to clusters.
#' @param cluster_fasta Path to the CD-HIT output FASTA file.
#'
#' @return A data frame with cluster names.
#' @export
#'
#' @examples
#' \dontrun{
#' .clusterNames(cluster_map, "path/to/clusters.fasta")
#' }
.clusterNames <- function(cluster_map, cluster_fasta) {

  cluster_map_unique <- cluster_map |>
    tibble::as_tibble() |>
    dplyr::distinct() |>
    dplyr::group_by(cluster) |>
    dplyr::slice_head(n=1)

  cdhit_output_faa <- Biostrings::readAAStringSet(cluster_fasta)

  names_faa <- names(cdhit_output_faa) |>
    tibble::as_tibble()

  names_faa <- names_faa |>
    dplyr::mutate(proteinID = stringr::str_extract(value,
                                                   "^fig\\|[0-9]+\\.[0-9]+\\.peg\\.[0-9]+")) |>
    dplyr::mutate(locus_tag = stringr::str_match(value,
                                                 "peg\\.[0-9]+\\|([^\\s]+)")[, 2]) |>
    dplyr::mutate(proteinName = stringr::str_trim(stringr::str_match(value,
                                                                     "\\|[^\\s]+\\s+(.*?)\\s+\\[")[, 2])) |>
    dplyr::select(-value)

  return(names_faa)
}

#' Cluster proteins with CD-HIT and write results to DuckDB
#'
#' Runs CD-HIT inside a Docker container to cluster protein sequences,
#' parses the resulting clusters, builds per-genome protein count matrices, and writes
#' multiple tables into a DuckDB database:
#' \itemize{
#'   \item \code{protein_count}: genome-by-protein-cluster count matrix.
#'   \item \code{protein_names}: mapping of cluster IDs to representative protein names.
#'   \item \code{protein_cluster_seq}: clustered protein sequences (AA) with IDs.
#' }
#'
#' @param duckdb_path Character scalar. Path to the DuckDB database file.
#' @param output_path Character scalar. Directory where CD-HIT outputs will be written.
#' @param output_prefix Character scalar. Prefix for CD-HIT output files. Defaults to
#'   \code{"cdhit_out"}.
#' @param identity Numeric. Sequence identity threshold for clustering (\code{-c}).
#'   Defaults to \code{0.9}.
#' @param word_length Integer. Word length parameter (\code{-n}). Defaults to \code{5}.
#'   Must be compatible with \code{identity} per CD-HIT rules.
#' @param threads Integer. Number of threads for CD-HIT (\code{-T}). \code{0} uses all
#'   available cores. Defaults to \code{0}.
#' @param memory Integer. Memory limit in MB (\code{-M}). \code{0} means unlimited.
#'   Defaults to \code{0}.
#' @param extra_args Character vector. Additional CD-HIT flags (e.g., \code{"-g", "1"}).
#'   Defaults to \code{c("-g", "1")}.
#'
#' @return
#'   writing the following DuckDB tables:
#'   \itemize{
#'     \item \code{"protein_count"}
#'     \item \code{"protein_names"}
#'     \item \code{"protein_cluster_seq"}
#'   }
#'
#' @examples
#' \dontrun{
#' # Cluster proteins and write results to DuckDB
#' CDHIT2duckdb(
#'   duckdb_path = "example.duckdb",
#'   output_path = "trial",
#'   identity = 0.9,
#'   word_length = 5,
#'   threads = 8
#' )
#' }
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
  aa <- Biostrings::readAAStringSet(cdhit_outputs$clustered_faa)

  protein_cluster_seq <- tibble::tibble(
    name     = names(aa),
    sequence = as.character(aa)
  )

  DBI::dbWriteTable(con, "protein_names", cluster_name, overwrite = TRUE)

  clustered_faa <- Biostrings::readAAStringSet(cdhit_outputs$clustered_faa)

  DBI::dbWriteTable(con, "protein_cluster_seq",
                    tibble::tibble(
                      name    = names(clustered_faa) |> stringr::str_extract("fig\\|[0-9]+\\.[0-9]+\\.peg\\.[0-9]+"),
                      sequence = as.character(clustered_faa)
                    ),
                    overwrite = TRUE)

  DBI::dbDisconnect(con)
}

#' Check if InterProScan data are ready for use
#'
#' @param version Character. InterProScan 5 version tag, e.g. "5.76-107.0".
#' @param dest_dir Character. Target directory for IPR data.
#' @param docker_image Character. Docker image tag to use.
#' @param platform Character or NULL.
#' @param curl_bin Character. System curl binary to use; defaults to "curl".
#' @param verbose Logical. Chatty progress updates.
#'
#' @return List of:
#'   - data_dir: absolute path to <dest_dir>/interproscan-<version>/data
#'   - ready: TRUE/FALSE indicating data presence and (basic) index check
#' @examples
#' \dontrun{
#'   res <- checkInterProData(
#'     version = "5.76-107.0",
#'     dest_dir = "results/interpro",
#'     platform = "linux/amd64"
#'   )
#'   # In .process_chunk, we can mount res$data_dir -> /opt/interproscan/data
#' }
checkInterProData <- function(
    version      = "5.76-107.0",
    dest_dir     = "inst/extdata/interpro",
    docker_image = sprintf("interpro/interproscan:%s", version),
    platform     = "linux/amd64",
    curl_bin     = "curl",
    verbose      = TRUE
) {
  msg <- function(...) if (verbose) message(sprintf(...))

  # Create base directory
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  dest_dir <- normalizePath(dest_dir, mustWork = TRUE)

  # Layout under dest_dir
  root_dir <- file.path(dest_dir, sprintf("interproscan-%s", version))
  data_dir <- file.path(root_dir, "data")

  # Is Pfam already pressed?
  is_indexed <- function() {
    length(Sys.glob(file.path(data_dir, "pfam", "*", "*.h3m"))) > 0
  }

  if (is_indexed()) {
    msg("InterProScan data ready at: %s", data_dir)
    return(list(data_dir = normalizePath(data_dir), ready = TRUE))
  }

  # Download the InterPro data bundle
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

  # Verify MD5 like InterPro asks nicely for us to do
  msg("Verifying MD5 checksum.")
  md5_expected <- sub("\\s+.*$", "", readLines(md5_path)[1])
  md5_actual   <- tools::md5sum(tar_path)[[1]]
  if (!identical(tolower(md5_expected), tolower(md5_actual)))
    stop("MD5 checksum mismatch for InterProScan data bundle.")

  # Extract using base R
  msg("Extracting InterProScan data bundle.")
  utils::untar(tar_path, exdir = dest_dir, tar = "internal")

  # Run InterProScan's pressing step
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

  # Final check
  if (is_indexed()) {
    msg("InterProScan Pfam HMMs are indexed and ready: %s", data_dir)
    return(list(data_dir = normalizePath(data_dir), ready = TRUE))
  } else {
    warning("Pfam HMM indices not found. InterProScan may press at runtime, which may be slow the first time it is run.")
    return(list(data_dir = normalizePath(data_dir), ready = FALSE))
  }
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
.getDfIPRColNames <- function() {
  column_names <- c(
    "AccNum", "SeqMD5Digest", "SLength", "Analysis",
    "DB.ID", "SignDesc", "StartLoc", "StopLoc", "Score",
    "Status", "RunDate", "IPRAcc", "IPRDesc", "placeholder"
  )
  return(column_names)
}

#' construct column types for reading interproscan output TSVs
#' (based upon the global variable written in
#' molevol_scripts/R/colnames_molevol.R)
#' @return [collector] a named vector of type expectations
#' for interproscan columns
#'
.getDfIPRColTypes <- function() {
  column_types <- readr::cols(
    "AccNum" = readr::col_character(),
    "SeqMD5Digest" = readr::col_character(),
    "SLength" = readr::col_integer(),
    "Analysis" = readr::col_character(),
    "DB.ID" = readr::col_character(),
    "SignDesc" = readr::col_character(),
    "StartLoc" = readr::col_integer(),
    "StopLoc" = readr::col_integer(),
    "Score" = readr::col_double(),
    "Status" = readr::col_character(),
    "RunDate" = readr::col_character(),
    "IPRAcc" = readr::col_character(),
    "IPRDesc" = readr::col_character(),
    "placeholder" = readr::col_character()
  )
  return(column_types)
}

#' Read an interproscan output TSV with standardized
#' column names and types
#' @param filepath [chr] path to interproscan output TSV
#' @return [tbl_df] interproscan output table
.readIPRscanTsv <- function(filepath) {
  df_ipr <- readr::read_tsv(filepath,
                            col_types = .getDfIPRColTypes(),
                            col_names = .getDfIPRColNames()
  )
  return(df_ipr)
}

# Function to process a chunk of sequences
#' Process a chunk of sequences through InterProScan
#'
#' @param chunk Character vector of sequences to process.
#' @param path Character. Path to the working directory.
#' @param ipr_data_path Character. Path to InterProScan data directory.
#' @param out_file_base Character. Base name for output files.
#' @param appl Character. InterProScan applications to run.
#' @param chunk_id Integer. Identifier for this chunk.
#' @param threads Integer. Number of threads to use.
#' @param file_format Character. Output file format.
#'
#' @return Data frame with InterProScan results.
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' .process_chunk(chunk, path, ipr_data_path, out_file_base, appl, chunk_id, threads, file_format)
#' }
.process_chunk <- function(chunk,
                           path,
                           ipr_data_path = "inst/extdata/interpro/data",
                           out_file_base,
                           appl,
                           chunk_id,
                           threads,
                           file_format,
                           docker_image = sprintf("interpro/interproscan:%s", "5.76-107.0")) {
  # Normalize and quote paths
  path      <- .docker_path(path)
  bind_data <- .docker_path(ipr_data_path)

  dir.create(file.path(path, "tmp", "iprscan"), recursive = TRUE, showWarnings = FALSE)

  fasta_sequences <- Biostrings::AAStringSet(chunk$sequence)
  names(fasta_sequences) <- chunk$name
  temp_fasta_file <- tempfile(tmpdir = path, fileext = ".fa")
  Biostrings::writeXStringSet(fasta_sequences, temp_fasta_file)

  chunk_out_file_base_host <- file.path(path, sprintf("%s_chunk_%d", out_file_base, chunk_id))
  chunk_out_file_base_cont <- .to_container(chunk_out_file_base_host, path, "/work")

  # pull the image to run
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

  # Determine output file path
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
#'
#' @param duckdb_path   Path to the DuckDB database file to open.
#' @param path          Working directory for chunk files and InterProScan outputs. Defaults to dirname(duckdb).
#' @param out_file_base Base prefix for InterProScan outputs per chunk. Default "iprscan".
#' @param appl          Vector of InterProScan applications, e.g. c("Pfam").
#' @param ipr_version   InterProScan version tag, e.g. "5.76-107.0".
#' @param ipr_dest_dir  Host directory where InterPro data live or will be created.
#' @param ipr_platform  Optional Docker platform string (e.g., "linux/amd64") for indexing step.
#' @param auto_prepare_data Logical. If TRUE, ensures the InterPro data are present & indexed (download if missing).
#' @param threads       Number of cores for parallel chunks and IPR --cpu per chunk.
#' @param file_format   InterProScan output format, default "TSV".
#' @param docker_repo   InterProScan Docker repo (prefix), default "interpro/interproscan".
#'
#' @return Writes domain tables into DuckDB; returns invisibly the domain_count write result.
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
    path <- dirname(duckdb_path)  # e.g., ./results/<bug>
  }
  path <- normalizePath(path)

  # Ensure InterPro data are present and indexed

  if (is.null(ipr_dest_dir)) {
    ipr_dest_dir <- file.path(path, "interpro")
  }

  ipr_image <- sprintf("%s:%s", docker_repo, ipr_version)

  ipr_info <- checkInterProData(
    version      = ipr_version,
    dest_dir     = ipr_dest_dir,
    docker_image = ipr_image,
    platform     = ipr_platform,
    verbose      = TRUE
  )

  ipr_data_path <- ipr_info$data_dir

#  ipr_data_path <- .docker_path(ipr_data_path)  # for -v mount

  # Load protein sequences to scan
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE), add = TRUE)

  sequences_df <- dplyr::tbl(con, "protein_cluster_seq") |> tibble::as_tibble()
  if (nrow(sequences_df) == 0L) {
    stop("No sequences found in 'protein_cluster_seq'. Please run CDHIT2duckdb() first.")
  }

  # Split sequences into chunks
  chunks <- split(sequences_df, ceiling(seq_along(sequences_df$name) / 5000))

  # OS-agnostic parallel implementation
  future::plan(future::multisession, workers = threads)
  results <- future.apply::future_lapply(seq_along(chunks), function(i) {
    res <- try(
      .process_chunk(
        chunk         = chunks[[i]],
        path          = path,
        ipr_data_path = ipr_data_path,
        out_file_base = out_file_base,
        appl          = appl,
        chunk_id      = i,
        threads       = threads,
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
  df_iprscan <- do.call(rbind, lapply(results, function(file) {
    if (!is.null(file) && file.exists(file)) .readIPRscanTsv(file)
  }))

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
    dplyr::summarize(across(everything(), ~ ifelse(any(. == "1"), "1", "0"))) |>
    dplyr::mutate(across(-AccNum, as.numeric))

  protein_filter <- dplyr::tbl(con, "protein_count") |>
    tibble::as_tibble() |>
    dplyr::select(genome_id, df_protein_domain_pa |> dplyr::pull(AccNum))

  domain_count <- as.matrix(protein_filter |> dplyr::select(-genome_id)) %*%
    as.matrix(df_protein_domain_pa |> dplyr::select(-AccNum)) |>
    tibble::as_tibble() |>
    dplyr::mutate(genome_id = protein_filter |> dplyr::pull(genome_id)) |>
    dplyr::relocate(genome_id, .before = dplyr::everything()) |>
    DBI::dbWriteTable(conn = con, name = "domain_count", overwrite = TRUE)

  invisible(domain_count)
}

#' Clean BV-BRC metadata, then save as Parquet files
#'
#' Cleans the AMR metadata in the DuckDB table `"filtered"` using reference
#' dictionaries (drug names, classes, abbreviations, countries), produces a
#' `"cleaned_metadata"` table, and then exports feature matrices and metadata
#' to compressed Parquet files. A new DuckDB database (suffix `_parquet.duckdb`)
#' is created with **views** that read directly from these Parquet files for
#' fast, columnar analytics.
#'
#' The function:
#' - Reads reference TSVs (`clean_drug.tsv`, `drug_class.tsv`, `drug_abbr.tsv`,
#'   `class_abbr.tsv`, `cleaned_bvbrc_countries.tsv`) from `ref_file_path`.
#' - Cleans `"filtered"`: normalizes antibiotic names, joins class & abbreviations,
#'   and rewrites `"filtered"` in the original DuckDB.
#' - Builds a per-genome **resistance summary** (distinct classes & concatenated
#'   class abbreviations) and writes a `"cleaned_metadata"` table.
#' - Exports feature tables (`gene_count`, `protein_count`, `domain_count`,
#'   `gene_struct`) as **long** Parquet files with Zstandard compression.
#' - Exports lookup/sequence tables (`metadata`, `gene_names`, `protein_names`,
#'   `domain_names`, `gene_ref_seq`, `protein_cluster_seq`, `genome_gene_protein`)
#'   as Parquet.
#'
#' @param duckdb_path Character scalar. Path to the source DuckDB file that already
#'   contains tables: `"filtered"`, `"gene_count"`, `"protein_count"`, `"domain_count"`,
#'   `"gene_struct"`, `"gene_names"`, `"protein_names"`, `"domain_names"`,
#'   `"gene_ref_seq"`, `"protein_cluster_seq"`, `"genome_gene_protein"`.
#' @param path Character scalar. Output directory for Parquet files and the new
#'   DuckDB (`*_parquet.duckdb`). Normalized via `normalizePath()`.
#' @param ref_file_path Character scalar. Directory containing reference TSV files:
#'   \itemize{
#'     \item \code{clean_drug.tsv} (columns should include \code{original_drug}, \code{cleaned_drug})
#'     \item \code{drug_class.tsv} (e.g., \code{drug}, \code{drug_class})
#'     \item \code{drug_abbr.tsv} (e.g., \code{drug}, \code{drug_abbr})
#'     \item \code{class_abbr.tsv} (e.g., \code{drug_class}, \code{class_abbr})
#'     \item \code{cleaned_bvbrc_countries.tsv} (at least \code{raw_entry}, \code{clean_name}, \code{short_name})
#'   }
#'
#' @return
#'   rewriting `"filtered"`, writing `"cleaned_metadata"`, exporting Parquet files,
#'   and creating Parquet-backed views in a new DuckDB (`*_parquet.duckdb`).
#'
#' @examples
#' \dontrun{
#' # Clean and export to Parquet-backed views
#' cleanData(
#'   duckdb_path    = "trial/example.duckdb",
#'   path          = "trial,
#'   ref_file_path = "/path/to/reference-tsvs"
#' )

#' Clean BV-BRC metadata, then save as Parquet files
#'
#' Cleans the AMR metadata in the DuckDB table `"filtered"` using reference
#' dictionaries (drug names, classes, abbreviations, countries), produces a
#' `"cleaned_metadata"` table, and then exports feature matrices and metadata
#' to compressed Parquet files. A new DuckDB database (suffix `_parquet.duckdb`)
#' is created with **views** that read directly from these Parquet files for
#' fast, columnar analytics.
#'
#' The function:
#' - Reads reference TSVs (`clean_drug.tsv`, `drug_class.tsv`, `drug_abbr.tsv`,
#'   `class_abbr.tsv`, `cleaned_bvbrc_countries.tsv`) from `ref_file_path`.
#' - Cleans `"filtered"`: normalizes antibiotic names, joins class & abbreviations,
#'   and rewrites `"filtered"` in the original DuckDB.
#' - Builds a per-genome **resistance summary** (distinct classes & concatenated
#'   class abbreviations) and writes a `"cleaned_metadata"` table.
#' - Exports feature tables (`gene_count`, `protein_count`, `domain_count`,
#'   `gene_struct`) as **long** Parquet files with Zstandard compression.
#' - Exports lookup/sequence tables (`metadata`, `gene_names`, `protein_names`,
#'   `domain_names`, `gene_ref_seq`, `protein_cluster_seq`, `genome_gene_protein`)
#'   as Parquet.
#' - Also exports `amr_phenotype`, `genome_data`, and the original joined `metadata`
#'   (as `original_metadata`) so **no information is lost** in the Parquet-backed DB.
#'
#' @param duckdb_path Character scalar. Path to the source DuckDB file that already
#'   contains tables: `"filtered"`, `"gene_count"`, `"protein_count"`, `"domain_count"`,
#'   `"gene_struct"`, `"gene_names"`, `"protein_names"`, `"domain_names"`,
#'   `"gene_ref_seq"`, `"protein_cluster_seq"`, `"genome_gene_protein"`.
#' @param path Character scalar. Output directory for Parquet files and the new
#'   DuckDB (`*_parquet.duckdb`). If missing or set to ".", "results", or "results/",
#'   it defaults to `dirname(duckdb_path)` (e.g., `./results/<bug>`). Normalized via `normalizePath()`.
#' @param ref_file_path Character scalar. Directory containing reference TSV files:
#'   \itemize{
#'     \item \code{clean_drug.tsv} (columns should include \code{original_drug}, \code{cleaned_drug})
#'     \item \code{drug_class.tsv} (e.g., \code{drug}, \code{drug_class})
#'     \item \code{drug_abbr.tsv} (e.g., \code{drug}, \code{drug_abbr})
#'     \item \code{class_abbr.tsv} (e.g., \code{drug_class}, \code{class_abbr})
#'     \item \code{cleaned_bvbrc_countries.tsv} (at least \code{raw_entry}, \code{clean_name}, \code{short_name})
#'   }
#'
#' @return
#'   rewriting `"filtered"`, writing `"cleaned_metadata"`, exporting Parquet files,
#'   and creating Parquet-backed views in a new DuckDB (`*_parquet.duckdb`).
#'
#' @examples
#' \dontrun{
#' # Clean and export to Parquet-backed views
#' cleanData(
#'   duckdb_path    = "./results/Staphylococcus_epidermidis/Sep.duckdb",
#'   path           = ".",   # defaults to ./results/Staphylococcus_epidermidis
#'   ref_file_path  = "/path/to/reference-tsvs"
#' )
cleanData <- function(duckdb_path, path, ref_file_path){
  duckdb_path  <- normalizePath(duckdb_path)
  # If no explicit path is provided (or a generic one), choose results/<bug>/ when
  # the DuckDB lives under data/<bug>/, or else fall back to the DuckDB directory.
  if (missing(path) || path %in% c(".", "results", "results/")) {
    bug_dir <- dirname(duckdb_path)                      # e.g., ./data/<bug>
    # map .../data/<bug>  ->  .../results/<bug>
    mapped_results <- sub(
      paste0(.Platform$file.sep, "data", .Platform$file.sep),
      paste0(.Platform$file.sep, "results", .Platform$file.sep),
      bug_dir,
      fixed = TRUE
    )
    # If mapping changed anything, use the results path; otherwise fallback to bug_dir
    if (!identical(mapped_results, bug_dir)) {
      path <- mapped_results
    } else {
      path <- bug_dir
    }
  }

  path <- normalizePath(path, mustWork = FALSE)
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)

  con          <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
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

  # create the resistance summary
  resistance_summary <- dplyr::tbl(con, "filtered") |>
    tibble::as_tibble()  |>
    dplyr::filter(genome_drug.resistant_phenotype == "Resistant") |>
    dplyr::group_by(genome_drug.genome_id) |>
    dplyr::summarise(
      num_resistant_classes = dplyr::n_distinct(drug_class),
      resistant_classes = paste(unique(class_abbr), collapse = "_")
    )

  year_breaks <- seq(1980, 2023, by = 5)
  filtered <- dplyr::tbl(con, "filtered") |>
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

  # Also export AMR/genome/original metadata so parquet DB is complete
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

  db_name  <- duckdb_path |> stringr::str_split_i(".duckdb", i = 1) |>
    paste0("_parquet.duckdb")
  con_new <- DBI::dbConnect(duckdb::duckdb(), db_name)

  # Read and transform gene data
  DBI::dbReadTable(con, "gene_count") |>
    tidyr::pivot_longer(-genome_id, names_to = "gene", values_to = "value") |>
    dplyr::filter(!is.na(value) & value != "") |>
    dplyr::mutate(value = as.integer(value)) |>
    writeCompressedParquet(genes_parquet)

  # Create a view that reads directly from the Parquet file
  DBI::dbExecute(con_new, sprintf("
  CREATE OR REPLACE VIEW gene_count AS
  SELECT * FROM read_parquet('%s')
", genes_parquet))

  # Protein
  DBI::dbReadTable(con, "protein_count") |>
    tidyr::pivot_longer(-genome_id, names_to = "protein", values_to = "value") |>
    dplyr::filter(!is.na(value) & value != "") |>
    dplyr::mutate(value = as.integer(value)) |>
    writeCompressedParquet(proteins_parquet)

  DBI::dbExecute(con_new, sprintf("
  CREATE OR REPLACE VIEW protein_count AS
  SELECT * FROM read_parquet('%s')
", proteins_parquet))

  # Domain features
  DBI::dbReadTable(con, "domain_count") |>
    tidyr::pivot_longer(-genome_id, names_to = "domain", values_to = "value") |>
    dplyr::filter(!is.na(value) & value != "") |>
    dplyr::mutate(value = as.integer(value)) |>
    writeCompressedParquet(domains_parquet)

  DBI::dbExecute(con_new, sprintf("
  CREATE OR REPLACE VIEW domain_count AS
  SELECT * FROM read_parquet('%s')
", domains_parquet))

  # Struct features
  DBI::dbReadTable(con, "gene_struct") |>
    tidyr::pivot_longer(-genome_id, names_to = "struct", values_to = "value") |>
    dplyr::filter(!is.na(value) & value != "") |>
    dplyr::mutate(value = as.integer(value)) |>
    writeCompressedParquet(struct_parquet)

  DBI::dbExecute(con_new, sprintf("
  CREATE OR REPLACE VIEW struct AS
  SELECT * FROM read_parquet('%s')
", struct_parquet))

  # Metadata (export cleaned as 'metadata')
  DBI::dbReadTable(con, "cleaned_metadata") |>
    writeCompressedParquet(metadata_parquet)

  DBI::dbExecute(con_new, sprintf("
  CREATE OR REPLACE VIEW metadata AS
  SELECT * FROM read_parquet('%s')
", metadata_parquet))

  # gene names
  DBI::dbReadTable(con, "gene_names") |>
    writeCompressedParquet(gene_names_parquet)

  DBI::dbExecute(con_new, sprintf("
  CREATE OR REPLACE VIEW gene_names AS
  SELECT * FROM read_parquet('%s')
", gene_names_parquet))

  # protein names
  DBI::dbReadTable(con, "protein_names") |>
    dplyr::select(-locus_tag) |>
    writeCompressedParquet(protein_names_parquet)

  DBI::dbExecute(con_new, sprintf("
  CREATE OR REPLACE VIEW protein_names AS
  SELECT * FROM read_parquet('%s')
", protein_names_parquet))

  # domain names
  DBI::dbReadTable(con, "domain_names") |>
    dplyr::select(-c(IPRAcc, IPRDesc)) |>
    writeCompressedParquet(domain_names_parquet)

  DBI::dbExecute(con_new, sprintf("
  CREATE OR REPLACE VIEW domain_names AS
  SELECT * FROM read_parquet('%s')
", domain_names_parquet))

  # gene ref seq
  DBI::dbReadTable(con, "gene_ref_seq") |>
    writeCompressedParquet(gene_ref_seq_parquet)

  DBI::dbExecute(con_new, sprintf("
  CREATE OR REPLACE VIEW gene_seqs AS
  SELECT * FROM read_parquet('%s')
", gene_ref_seq_parquet))

  # protein cluster seq
  DBI::dbReadTable(con, "protein_cluster_seq") |>
    writeCompressedParquet(protein_cluster_seq_parquet)

  DBI::dbExecute(con_new, sprintf("
  CREATE OR REPLACE VIEW protein_seqs AS
  SELECT * FROM read_parquet('%s')
", protein_cluster_seq_parquet))

  # genome gene protein mapping
  DBI::dbReadTable(con, "genome_gene_protein") |>
    writeCompressedParquet(genome_gene_protein_parquet)

  DBI::dbExecute(con_new, sprintf("
  CREATE OR REPLACE VIEW genome_gene_protein AS
  SELECT * FROM read_parquet('%s')
", genome_gene_protein_parquet))

  ### This is for debugging and is not currently used downstream ###
  DBI::dbReadTable(con, "amr_phenotype") |> writeCompressedParquet(amr_phenotype_parquet)
  DBI::dbReadTable(con, "genome_data")   |> writeCompressedParquet(genome_data_parquet)
  DBI::dbReadTable(con, "metadata")      |> writeCompressedParquet(original_metadata_parquet)

  DBI::dbExecute(con_new, sprintf("
  CREATE OR REPLACE VIEW amr_phenotype AS
  SELECT * FROM read_parquet('%s')
", amr_phenotype_parquet))

  DBI::dbExecute(con_new, sprintf("
  CREATE OR REPLACE VIEW genome_data AS
  SELECT * FROM read_parquet('%s')
", genome_data_parquet))

  DBI::dbExecute(con_new, sprintf("
  CREATE OR REPLACE VIEW original_metadata AS
  SELECT * FROM read_parquet('%s')
", original_metadata_parquet))
}
