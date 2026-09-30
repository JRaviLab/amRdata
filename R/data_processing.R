#' @importFrom data.table :=
NULL

#' Clear cached HMMER databases
#'
#' Removes user-specified HMMER databases from the shared amRdata BFC registry
#' and deletes their local databases. Databases will be downloaded and prepared
#' again the next time they are requested. This can help resolve corrupt database
#' issues that may arise from time to time, especially on certain environments
#' with unstable network connections.
#'
#' If specific `databases` are not supplied in an interactive R session, a menu
#' allows the user to select a specific database, or remove all databases.
#'
#' @param databases Character vector of HMMER databases to remove.
#'   Supported values are `"Pfam"`, `"COG"`, `"AMRFinder"`, and `"DefenseCas"`.
#'   If `NULL` in an interactive session, the user is prompted to choose.
#' @param verbose Logical. Print information about removed databases.
#'   Default: `TRUE`.
#'
#' @return Invisibly returns the names of databases removed.
#'
#' @export
clearHMMERdatabases <- function(
    databases = NULL,
    verbose = TRUE
) {
  supported <- c(
    "Pfam",
    "COG",
    "AMRFinder",
    "DefenseCas"
  )

  # Interactive selection if no database supplied
  if (is.null(databases)) {
    if (!interactive()) {
      stop(
        "Values for `databases` must be supplied in non-interactive sessions.",
        call. = FALSE
      )
    }

    selection <- utils::menu(
      choices = c(supported, "All"),
      title = "Which HMMER database would you like to remove?"
    )

    # utils::menu() returns a 0 value when cancelled. Reassure user that
    # no damage was done to their precious databases
    if (selection == 0L) {
      if (isTRUE(verbose)) message("No HMMER databases removed.")

      return(invisible(character(0)))
    }

    # Final numbered option is "All" if you want the nuclear option
    if (selection == length(supported) + 1L) {
      databases <- supported
    } else {
      databases <- supported[[selection]]
    }
  }

  databases <- unique(as.character(databases))

  if (!length(databases)) {
    stop("At least one HMMER database must be specified.", call. = FALSE)
  }

  # For when you either have a typo or forget what databases there are
  unknown <- setdiff(databases, supported)

  if (length(unknown)) {
    stop(
      "Unknown HMMER database(s): ",
      paste(unknown, collapse = ", "),
      ". Supported databases are: ",
      paste(supported, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  # What's in the file cache already?
  bfc <- .amr_bfc()
  resources <- .amr_bfc_hmmer_resources()
  hmmer_dir <- .defaultHmmerDbDir()

  for (db in databases) {
    prefix <- .amr_bfc_hmmer_rname(db)

    # Match the database itself and any registered components
    hits <- resources[resources$rname == prefix | startsWith(resources$rname, paste0(prefix, "_")),, drop = FALSE]

    if (nrow(hits)) {
      BiocFileCache::bfcremove(bfc, hits$rid)
    }

    # Purge the whole directory so partially downloaded/extracted files,
    # combined HMMs, and hmmpress fluff cannot persist in vile ways
    db_dir <- file.path(hmmer_dir, db)

    if (dir.exists(db_dir)) {
      unlink(db_dir, recursive = TRUE, force = TRUE)
    }

    if (isTRUE(verbose)) message("Cleared HMMER database: ", db)
  }

  invisible(databases)
}

#' Run one Panaroo pangenome job
#'
#' Runs Panaroo inside Docker for one batch of prepared genome GFF/FASTA inputs.
#' Host paths are translated to container-visible paths and Panaroo output is
#' written to a uniquely named directory under `output_path`.
#'
#' Panaroo is run with strict cleaning, paralog merging, and invalid-gene
#' removal. Refinding behavior is controlled by `refind_mode`.
#'
#' See Panaroo's documentation for details on how the parameters affect your
#' pangenome output: https://gthlab.au/panaroo/#/gettingstarted/params
#'
#' @param batch_input A series of genome IDs for input
#' @param output_path Character scalar. Base directory for Panaroo outputs and temporary files.
#' @param mount_root Character. Host dataset directory mounted into the Docker
#'   container.
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
#' @param log_path Character or `NULL`. Optional path for Panaroo command output
#'   logging.
#'
#' @returns A list of results for each Panaroo batch in its output directory.
#'
#' @keywords internal
#' @examples NULL
.processPanaroo <- function(batch_input,
                            output_path,
                            mount_root,
                            core_threshold,
                            len_dif_percent,
                            cluster_threshold,
                            family_seq_identity,
                            panaroo_threads_per_job,
                            refind_mode = c("off", "default", "strict"),
                            log_path = NULL) {
  refind_mode <- match.arg(refind_mode)
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  output_path <- .docker_path(output_path)
  mount_root <- .docker_path(mount_root)

  # Fail fast if Docker is missing
  if (!nzchar(Sys.which("docker"))) {
    stop("Docker is not available on your PATH but is required to run Panaroo.")
  }

  mount_cont <- "/work"

  # Write the genome list file (convert each "gff fna" to container-visible paths)
  genome_filepath_host <- tempfile(pattern = "genomeFilepath_", fileext = ".txt", tmpdir = output_path)

  batch_input_cont <- purrr::map_chr(unlist(batch_input), function(line) {
    parts <- strsplit(line, " +")[[1]]
    parts_cont <- .to_container(parts, host_root = mount_root, container_root = mount_cont)
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
  genome_filepath_cont <- .to_container(genome_filepath_host, host_root = mount_root, container_root = mount_cont)
  output_dir_cont <- .to_container(output_dir_host, host_root = mount_root, container_root = mount_cont)

  # Run Panaroo in Docker
  cmd_args <- c(
    "run",
    "--platform", "linux/amd64",
    "--rm",
    "-v", paste0(mount_root, ":", mount_cont),
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
  res <- suppressWarnings(
  system2(
    "docker",
    args = cmd_args,
    stdout = TRUE,
    stderr = TRUE
  )
)

status <- attr(
  res,
  "status"
)

.log_tool_output(
  log_path,
  paste0(
    "Panaroo: ",
    basename(output_dir_host)
  ),
  res
)

if (!is.null(status) && status != 0L) {
  stop(
    "Panaroo failed with exit status ",
    status,
    ":\n",
    paste(
      utils::tail(
        res,
        40L
      ),
      collapse = "\n"
    ),
    "\n\nFull Panaroo output was written to the processing log.",
    call. = FALSE
  )
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
.runPanaroo <- function(duckdb_path,
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
                        verbose = TRUE,
                        log_path = NULL) {
  refind_mode <- match.arg(refind_mode)
  threads <- .resolve_workers(requested = threads)
  duckdb_path <- normalizePath(duckdb_path, mustWork = TRUE)
  paths <- .amr_paths_from_duckdb(duckdb_path)

  if (is.null(output_path)) {
    output_path <- paths$panaroo
  }
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

  output_path <- normalizePath(output_path, mustWork = TRUE)

  # Avoiding keeping DuckDB connections live when playing with futures
  genome_query_output <- local({con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)

    on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

    DBI::dbGetQuery(con, "SELECT * FROM files ORDER BY genome_id")
  })

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

    audit <- cleaned$audit

    .log_or_message(
      log_path,
      verbose,
      sprintf(
        "Pseudogene audit: %d genomes, %d removed pseudogenes, %d features remain.",
        nrow(audit),
        sum(audit$n_pseudogene, na.rm = TRUE),
        sum(audit$n_kept, na.rm = TRUE)
      ))
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

  # Never allow more simultaneous Panaroo jobs than total CPUs -- that's bad
  n_parallel_jobs <- .resolve_workers(requested = min(n_jobs, threads),
    n_tasks = n_jobs, warn = FALSE)

  # Divide the total CPU budget across simultaneous Panaroo jobs
  panaroo_threads_per_job <- max(1L, floor(threads / n_parallel_jobs))

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  .amr_set_future_plan(n_parallel_jobs)

  batch_panaroo_run <- furrr::future_map(
    panaroo_batches,
    ~ .processPanaroo(
      batch_input             = .x,
      output_path             = output_path,
      mount_root              = paths$root,
      core_threshold          = core_threshold,
      len_dif_percent         = len_dif_percent,
      cluster_threshold       = cluster_threshold,
      family_seq_identity     = family_seq_identity,
      panaroo_threads_per_job = panaroo_threads_per_job,
      refind_mode             = refind_mode,
      log_path                = log_path
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
                          threads = 8,
                          log_path = NULL) {
  input_path <- .docker_path(input_path)

  # Fail fast if Docker is missing
  if (!nzchar(Sys.which("docker"))) {
    stop("Docker is not available on your PATH but is required to run panaroo-merge.")
  }

  threads <- .resolve_workers(requested = threads)
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

    output <- system2("docker", args = cmd_args, stdout = TRUE, stderr = TRUE)

    .log_tool_output(log_path, "Panaroo merge", output)

    status <- attr(output, "status")

    if (!is.null(status) && status != 0L) {
      stop("Panaroo merge failed with exit status ",
           status,
           ":\n",
           paste(output, collapse = "\n"))
    }
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
  on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

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
  on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

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
  on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

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
  on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

  gene_fasta <- Biostrings::readDNAStringSet(filepath = fasta_filepath)
  DBI::dbWriteTable(con, "gene_ref_seq",
    tibble::tibble(
      name = names(gene_fasta),
      sequence = as.character(gene_fasta)
    ),
    overwrite = TRUE
  )
  # col_types = FALSE reduces console spam
  readr::read_csv(file.path(panaroo_output_path, "gene_presence_absence.csv"), show_col_types = FALSE) |>
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
    dplyr::mutate(protein_ids = gsub("_len", "", protein_ids)) |>
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
                      output_path = NULL,
                      output_prefix = "cdhit_out",
                      identity = 0.9,
                      word_length = 5,
                      threads = 0,
                      memory = 0,
                      extra_args = c("-g", "1"),
                      verbose = TRUE,
                      log_path = NULL) {
  # Fail fast if Docker is missing
  if (!nzchar(Sys.which("docker"))) {
    stop("Docker is not available on your PATH but is required to run CD-HIT.")
  }

  # CD-HIT sets threads = 0 to mean all CPUs, we limit to CPUs available to this R session
  if (identical(threads, 0L) || identical(threads, 0)) {
    threads <- .resolve_workers(requested = NULL)
  } else {
    threads <- .resolve_workers(requested = threads)
  }
  duckdb_path <- normalizePath(duckdb_path, mustWork = TRUE)
  paths <- .amr_paths_from_duckdb(duckdb_path)
  if (is.null(output_path)) {
    output_path <- paths$cdhit
  }
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  duckdb_path <- .docker_path(duckdb_path)
  output_path <- .docker_path(output_path)

  genome_query_output <- local({
    con <- DBI::dbConnect(duckdb::duckdb(),duckdb_path)
    on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

    DBI::dbGetQuery(con, "SELECT * FROM files ORDER BY genome_id")
  })

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

  .log_or_message(log_path, verbose, "Running CD-HIT via Docker...")
  output <- tryCatch({
    system2("docker",
            args = cmd_args,
            stdout = TRUE,
            stderr = TRUE)
  }, error = function(e) {
    stop("CD-HIT execution failed: ", e$message)
  })
  .log_tool_output(log_path, "CD-HIT", output)

  if (!file.exists(clustered_faa)) {
    stop("CD-HIT failed: output file not found. Check stderr:\n", paste(output, collapse = "\n"))
  }
  # Ensure .clstr exists (used downstream)
  if (!file.exists(paste0(clustered_faa, ".clstr"))) {
    stop(
      "CD-HIT did not produce the expected .clstr file at: ", paste0(clustered_faa, ".clstr"),
      "\nFull output:\n", paste(output, collapse = "\n")
    )
  }

  if (isTRUE(verbose)) {
    .log_or_message(log_path, verbose, "CD-HIT completed successfully.")
  }
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
                              verbose = TRUE,
                              log_path = NULL) {
  refind_mode <- match.arg(refind_mode)
  duckdb_path <- normalizePath(duckdb_path, mustWork = TRUE)
  paths <- .amr_paths_from_duckdb(duckdb_path)
  out_dir <- if (is.null(output_path)) {
    paths$panaroo
  } else {
    normalizePath(output_path, mustWork = FALSE)
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  .log_or_message(log_path, verbose, "Launching Panaroo.")

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
    verbose = verbose,
    log_path = log_path
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
    .log_or_message(log_path, verbose, "Merging Panaroo batch outputs.")
    .mergePanaroo(
      input_path          = out_dir,
      core_threshold      = core_threshold,
      len_dif_percent     = len_dif_percent,
      cluster_threshold   = cluster_threshold,
      family_seq_identity = family_seq_identity,
      threads             = max(1L, floor(threads / 2)),
      log_path            = log_path
    )
    target_dir <- file.path(out_dir, "merge_output")
    if (!file.exists(file.path(target_dir, "gene_presence_absence.csv"))) {
      stop("Expected merged Panaroo outputs in merge_output/, but files were not found.")
    }
  } else {
    target_dir <- valid[[1]]
  }

  .log_or_message(log_path, verbose, "Writing Panaroo tables to DuckDB.")
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
#' @param cluster_fasta Path to representative FASTA file used by CD-HIT.
#'
#' @return A tibble containing protein metadata.
#'
#' @keywords internal
.clusterNames <- function(cluster_fasta) {
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
                         output_path = NULL,
                         output_prefix = "cdhit_out",
                         identity = 0.9,
                         word_length = 5,
                         threads = 0,
                         memory = 0,
                         extra_args = c("-g", "1"),
                         verbose = TRUE,
                         log_path = NULL) {
  duckdb_path <- normalizePath(duckdb_path, mustWork = TRUE)

  paths <- .amr_paths_from_duckdb(duckdb_path)

  if (is.null(output_path)) {
    output_path <- paths$cdhit
  }

  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

  output_path <- normalizePath(output_path, mustWork = TRUE)

  cdhit_outputs <- .runCDHIT(
    duckdb_path,
    output_path,
    output_prefix = output_prefix,
    identity = identity,
    word_length = word_length,
    threads = threads,
    memory = memory,
    extra_args = extra_args,
    verbose = verbose,
    log_path = log_path
  )

  cluster_map <- .parseProteinClusters(cdhit_outputs$clustered_faa)
  cluster_count <- .buildProtMatrices(cluster_map)

  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

  DBI::dbWriteTable(con, "protein_count", cluster_count, overwrite = TRUE)

  cluster_fasta <- cdhit_outputs$cdhit_input_faa
  cluster_name <- .clusterNames(cluster_fasta)
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
#' @param hmmer_db_dir Character. Directory where HMMER databases are cached.
#' @param databases Character vector of database names to prepare. Supported
#'   built-in databases are `Pfam`, `COG`, and `AMRFinder`; `DefenseCas` is
#'   prepared separately by the DefenseFinder/CasFinder workflow.
#' @param docker_image Character. Docker image containing HMMER tools used to
#'   press the prepared databases. Default: `"staphb/hmmer"`.
#' @param hmmer_db_url NON-FUNCTIONAL. Character or `NULL`. URL used to download
#'   a custom HMMER database when `databases` contains names not covered by the
#'   built-in database definitions. This function is not currently active!
#' @param verbose Logical. Print status messages while checking, downloading,
#'   combining, and pressing databases. Default: `TRUE`.
#' @param log_path Character or `NULL`. Optional path for database preparation
#'   and external-tool logging.
#' @param progress Logical. Show download progress when HMM database archives
#'   must be retrieved. Default: `TRUE`.
#'
#' @return A named list containing the prepared HMM path and associated files
#'   for each requested database..
#'
#' @keywords internal
.prepareHmmerDatabases <- function(
    hmmer_db_dir,
    databases = c("Pfam", "COG", "AMRFinder"),
    docker_image = "staphb/hmmer",
    hmmer_db_url = NULL,
    verbose = TRUE,
    log_path = NULL,
    progress = TRUE
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

      .amr_download_file(db$url, tmp, progress = progress)

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

      output <- .amr_progress_step(
          paste0("Preparing ", db_name, " HMM database"),
          system2(
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
          ),
          progress = progress,
          verbose = verbose,
          log_path = log_path
        )

      .log_tool_output(log_path, "HMMER", output)

      if (!all(file.exists(pressed_files))) {
        stop(
          "hmmpress failed for ",
          db_name,
          "\n",
          paste(output, collapse = "\n")
        )
      }
    }

    # Registering HMMER databases in BiocFileCache for later use
    bfc_resource <- .amr_bfc_register_hmmer(
      database = db_name,
      hmm_path = hmm_file
    )

    db_paths[[db_name]] <- list(
      hmm = hmm_file,
      source = db$url,
      type = db$type,
      pressed = pressed_files,
      bfc_rid = bfc_resource$rid,
      bfc_rname = bfc_resource$rname
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





#' Run one HMMER chunk job
#'
#' Runs `hmmsearch` for one protein FASTA chunk against one prepared HMM
#' database. The available CPU budget is divided across concurrently running
#' HMMER workers. Raw HMMER output is parsed and written to a Parquet file for
#' later combination by [.runHMMER()].
#'
#' @param JOB_NAME Character. Unique name for this HMMER job, used to name
#'   intermediate and output files.
#' @param FASTA Character. Filename of the protein FASTA chunk to search.
#' @param DB Character. Name of the prepared HMM database to search.
#' @param total_proteins Integer. Total number of representative proteins in the
#'   dataset, used to set HMMER's `-Z` and `--domZ` search-space values.
#' @param output_path Character. Directory containing the FASTA input and where
#'   HMMER intermediate and Parquet outputs should be written.
#' @param db_paths List. Prepared HMM database metadata indexed by database name.
#' @param docker_image Character. Docker image containing HMMER.
#'   Default: `"staphb/hmmer"`.
#' @param threads Integer. Total CPU budget available to HMMER. Default: `8`.
#' @param n_workers Integer. Number of HMMER jobs running concurrently. Used to
#'   divide the CPU budget among jobs. Default: `8`.
#' @param log_path Character or `NULL`. Optional path for HMMER command output
#'   logging.
#'
#' @return Character path to the parsed Parquet output for this HMMER job.
#' @keywords internal
.runHmmerJob <- function(JOB_NAME, FASTA, DB, total_proteins,
                         output_path = NULL, db_paths,
                         docker_image = "staphb/hmmer", threads = 8L,
                         n_workers = 8L,
                         log_path = NULL
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
    "-Z", total_proteins,
    "--domZ", total_proteins,
    "--domtblout", .to_container(hmmer_output, mount_host, mount_cont),
    db_cont_path,
    .to_container(hmmer_input, mount_host, mount_cont)
  )

  stderr_file <- tempfile(pattern = paste0(JOB_NAME, "_"), fileext = ".stderr")

  on.exit(unlink(stderr_file, force = TRUE), add = TRUE)

  status <- tryCatch(
    system2(
      "docker",
      args = cmd_args,
      stdout = FALSE,
      stderr = stderr_file
    ),
    error = function(e) {
      stop("hmmsearch execution failed: ", e$message)
    }
  )

  if (!identical(status, 0L) ||
      !file.exists(hmmer_output)) {
    diagnostics <- if (file.exists(stderr_file)) {
      readLines(stderr_file, warn = FALSE)
    } else {
      character()
    }

    if (length(diagnostics)) {
      .log_tool_output(log_path, paste0("HMMER failure: ", JOB_NAME), diagnostics)
    }

    stop("hmmsearch failed for ",
         JOB_NAME,
         ". See the processing log for diagnostics.",
         call. = FALSE)
  }

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


#' Run HMMER against prepared protein reference databases
#'
#' Splits representative protein-cluster sequences into FASTA chunks and runs
#' HMMER searches against each requested reference database. Chunks for each
#' database are processed in parallel, then combined into database-specific
#' Parquet files and DuckDB tables.
#'
#' HMM databases are prepared automatically when needed. Progress reports the
#' number of completed chunks for the database currently being processed.
#'
#' @param duckdb_path Character. Path to the DuckDB database containing
#'   `protein_cluster_seq`, which provides the protein sequences to analyze.
#' @param output_path Character. Directory for HMMER intermediate and final
#'   Parquet outputs.
#' @param threads Integer. Total CPU budget used by HMMER jobs. Default: `8`.
#' @param hmmer_db_dir Character. Directory containing the prepared HMMER
#'   databases. If `NULL`, the default `amRdata` HMMER database cache is used.
#' @param databases Character vector of HMMER databases to run.
#' @param docker_image Character. Docker image containing HMMER. Default:
#'   `"staphb/hmmer"`.
#' @param num_of_splits Integer. Number of chunks into which the protein
#'   sequences should be divided. Must be a positive integer. The requested
#'   value is automatically reduced when fewer protein sequences are available.
#'   Default: `8`.
#' @param n_workers Integer. Maximum number of HMMER chunk jobs to run in
#'   parallel. Automatically reduced when fewer chunks are available.
#'   Default: `8`.
#' @param verbose Logical. Print persistent HMMER status and diagnostic
#'   messages. Default: `TRUE`.
#' @param log_path Character or `NULL`. Optional path for HMMER and external
#'   tool output logging.
#' @param progress Logical. Show temporary completed-chunk progress.
#'   Default: `TRUE`.
#' @return Invisibly returns a list containing the prepared database paths and
#'   final HMMER output Parquet files.
#'
#' @keywords internal
.runHMMER <- function(duckdb_path,
                      output_path = NULL,
                      threads = 8L,
                      hmmer_db_dir = NULL,
                      databases = c("Pfam", "COG", "AMRFinder"),
                      docker_image = "staphb/hmmer",
                      num_of_splits = 8L,
                      n_workers = 8L,
                      verbose = TRUE,
                      log_path = NULL,
                      progress = TRUE)
  {
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

  threads <- .resolve_workers(requested = threads)
  duckdb_path <- normalizePath(duckdb_path, mustWork = TRUE)
  paths <- .amr_paths_from_duckdb(duckdb_path)

  if (is.null(output_path)) {
    output_path <- paths$hmmer
  }

  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  duckdb_path <- .docker_path(duckdb_path)
  output_path <- .docker_path(output_path)

  prot_seqs <- local({
    con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
    on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)
    DBI::dbReadTable(con, "protein_cluster_seq") |>
      tibble::as_tibble()
  })

  # Just in case CD-HIT failed to generate sequences somehow
  if (nrow(prot_seqs) == 0L) {
    stop("No sequences found in 'protein_cluster_seq'. Please run CDHIT2duckdb() first.")
  }

  # required to define the database size for hmmsearch --Z and --domZ parameters
  total_proteins <- nrow(prot_seqs)

  if (is.null(hmmer_db_dir)) {
    hmmer_db_dir <- .defaultHmmerDbDir()
  }

  dir.create(
    hmmer_db_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # database paths
  .log_or_message(log_path, verbose, "Preparing HMM databases")
  db_paths <- .prepareHmmerDatabases(
    hmmer_db_dir = hmmer_db_dir,
    databases = databases,
    docker_image = docker_image,
    verbose = verbose,
    log_path = log_path,
    progress = progress
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
    dplyr::select(JOB_NAME, FASTA, DB, chunk)

  n_workers <- .resolve_workers(requested = n_workers, n_tasks = min(actual_chunk_count, threads))

  parquet_files <- local({
    old_plan <- future::plan()

    on.exit(
      .amr_progress_step(
        "Closing HMMER workers",
        future::plan(old_plan),
        progress = progress,
        verbose = verbose,
        log_path = log_path
      ),
      add = TRUE
    )

    .amr_set_future_plan(n_workers)

  .log_or_message(log_path, verbose, "Running ", nrow(job_list), " HMMER jobs.")

  purrr::map(seq_along(databases), function(db_i) {
    db <- databases[[db_i]]
    db_rows <- which(job_list$DB == db)

    progress_message <- sprintf("HMMER %d/%d: %s chunk", db_i, length(databases), db)

    .amr_with_progress({
      p <- .amr_progressor(
        steps = length(db_rows),
        progress = progress,
        label = db,
        message = progress_message
      )

      p(amount = 0, message = progress_message)

      furrr::future_map_chr(db_rows, function(i) {
        result <- .runHmmerJob(
          JOB_NAME = job_list$JOB_NAME[[i]],
          FASTA = job_list$FASTA[[i]],
          DB = db,
          total_proteins = total_proteins,
          output_path = output_path,
          db_paths = db_paths,
          docker_image = docker_image,
          threads = threads,
          n_workers = n_workers,
          log_path = log_path
        )

        p(message = progress_message)

        result
      }, .options = furrr::furrr_options(seed = TRUE))
    }, progress = progress, type = "steps")}) |>
    unlist(use.names = FALSE)
  })

  parquet_tbl <- tibble::tibble(parquet = parquet_files, db = job_list$DB)

  final_parquets <- .amr_with_progress(
    {
      p <- .amr_progressor(steps = length(databases),
                           progress = progress,
                           message = "Finalizing HMMER annotations")

      con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)

      on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

      purrr::set_names(databases) |>
        purrr::map(function(database_name) {

          progress_message <- paste0("Finalizing HMMER annotations: ", database_name)

          p(amount = 0, message = progress_message)

          .log_or_message(log_path, verbose, "Combining ", database_name)

          db_files <- parquet_tbl |>
            dplyr::filter(db == database_name) |>
            dplyr::pull(parquet)

          combined_tbl <- db_files |>
            purrr::map(arrow::read_parquet) |>
            dplyr::bind_rows() |>
            dplyr::left_join(
              .parse_hmmer_profiles(db_paths[[database_name]]$hmm) |>
                dplyr::select(query_name = profile_name, description = profile_description),
              by = "query_name"
            )

          final_parquet <- file.path(output_path, paste0("protein_", database_name, ".parquet"))

          .write_compressed_parquet(combined_tbl, final_parquet)

          DBI::dbWriteTable(
            con,
            name = paste0("protein_", database_name),
            value = combined_tbl,
            overwrite = TRUE
          )

          .log_or_message(log_path, verbose, "Created ", basename(final_parquet))

          p(message = paste0("Finalized ", database_name))

          final_parquet
        })
    }, progress = progress, type = "steps")

  unlink(
    list.files(output_path, pattern = "^protein_chunk_.*\\.(fasta|tbl|parquet)$", full.names = TRUE)
  )

  invisible(list(databases = db_paths, outputs = final_parquets))
}

#' Map HMMER protein annotations to genome-level count matrix and load into DuckDB
#'
#' Reads a Parquet file of HMMER hits (produced by [.runHMMER()]), joins the
#' annotations to the protein-cluster count matrix already in DuckDB, aggregates
#' counts per genome and annotation, and writes the result both as a Parquet file
#' and as a new table in the DuckDB database.
#'
#' @param duckdb_path Character. Path to the per-selection DuckDB database
#'   containing the `protein_count` table created by [CDHIT2duckdb()].
#' @param databases Character vector of HMMER database names to process.
#'   Each database must correspond to a `protein_<database>` annotation table
#'   already present in the DuckDB.
#' @param output_path Character. Directory where the genome-by-annotation
#'   Parquet files will be written. Defaults to `dirname(duckdb_path)`.
#'
#' @return Invisibly returns the path to the written count Parquet file.
#'
#' @seealso [CDHIT2duckdb()], [runDataProcessing()]
#'
#' @keywords internal
.proteinAnnotations2Duckdb <- function(
    duckdb_path,
    databases,
    output_path = NULL,
    verbose = TRUE
) {
  duckdb_path <- normalizePath(duckdb_path, mustWork = TRUE)
  paths <- .amr_paths_from_duckdb(duckdb_path)
  if (is.null(output_path)) {
    output_path <- paths$orb
  }

  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  output_path <- normalizePath(output_path, mustWork = TRUE)
  duckdb_path <- .docker_path(duckdb_path)
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

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

    if (isTRUE(verbose)) {
      message(
        "Processing ",
        annotation_table
      )
    }

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

    if (isTRUE(verbose)) {
      message(
        "Created ",
        count_table
      )
    }
  }

  invisible(count_paths)
}

#' Prepare DefenseFinder and CasFinder HMMs and annotate proteins
#'
#' Downloads and prepares the DefenseFinder and CasFinder model collections,
#' combines their HMM profiles, runs HMMER against the representative protein
#' sequences in the selected dataset, and writes the resulting annotations for
#' downstream processing.
#'
#' Prepared model files are reused when already present.
#'
#' @param defense_db_dir Character. Directory used to cache DefenseFinder and
#'   CasFinder model files.
#' @param docker_image Character. Docker image containing HMMER.
#'   Default: `"staphb/hmmer"`.
#' @param duckdb_path Character. Path to the working dataset DuckDB.
#' @param output_path Character or `NULL`. Directory for HMMER outputs. If
#'   `NULL`, the dataset HMMER directory is used.
#' @param threads Integer. CPU budget available to HMMER. Default: `8`.
#' @param verbose Logical. Print persistent preparation and execution messages.
#'   Default: `TRUE`.
#' @param log_path Character or `NULL`. Optional path for external-tool logging.
#' @param progress Logical. Show download progress when DefenseFinder or
#'   CasFinder model archives must be retrieved. Default: `TRUE`.
#'
#' @return Invisibly returns the path to the DefenseCas annotation Parquet file.
#' @keywords internal
.defenseHMMER <- function(
    defense_db_dir,
    docker_image = "staphb/hmmer",
    duckdb_path,
    output_path = NULL,
    threads = 8L,
    verbose = TRUE,
    log_path = NULL,
    progress = TRUE
) {

  if (!nzchar(Sys.which("docker"))) {
    stop("Docker is required.")
  }

  threads <- .resolve_workers(requested = threads)
  duckdb_path <- normalizePath(duckdb_path, mustWork = TRUE)
  paths <- .amr_paths_from_duckdb(duckdb_path)
  defense_db_dir <- normalizePath(defense_db_dir, mustWork = FALSE)
  if (is.null(output_path)) {
    output_path <- paths$hmmer
  }
  dir.create(defense_db_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  output_path <- normalizePath(output_path, mustWork = TRUE)

  ####################################################################
  # download repositories
  ####################################################################

  defense_dir <- file.path(defense_db_dir, "DefenseFinder")

  cas_dir <- file.path(defense_db_dir, "CasFinder")

  if (!dir.exists(defense_dir)) {

    .log_or_message(log_path, verbose, "Downloading DefenseFinder models")

    tmp <- tempfile(fileext = ".zip")

    .amr_download_file(
      url = "https://github.com/mdmparis/defense-finder-models/archive/refs/heads/master.zip",
      destfile = tmp,
      progress = progress
    )

    utils::unzip(
      tmp,
      exdir = defense_dir
    )

    unlink(tmp)
  }

  if (!dir.exists(cas_dir)) {

    .log_or_message(log_path, verbose, "Downloading CasFinder models")

    tmp <- tempfile(fileext = ".zip")

    .amr_download_file(url = "https://github.com/macsy-models/CasFinder/archive/refs/heads/main.zip",
      tmp,
      progress = progress)

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

      .log_or_message(log_path, verbose, "Running hmmpress for ", db_name)

      stderr_file <- tempfile(
        pattern = paste0(
          "hmmpress_",
          db_name,
          "_"
        ),
        fileext = ".stderr"
      )

      on.exit(
        unlink(
          stderr_file,
          force = TRUE
        ),
        add = TRUE
      )

      status <- tryCatch(
        system2(
          "docker",
          args = c(
            "run",
            "--rm",
            "-v",
            paste0(dirname(combined_hmm), ":/db"),
            docker_image,
            "hmmpress",
            file.path("/db", basename(combined_hmm))
          ),
          stdout = FALSE,
          stderr = stderr_file
        ),
        error = function(e) {
          stop(
            "hmmpress execution failed for ",
            db_name,
            ": ",
            e$message
          )
        }
      )

      if (
        !identical(status, 0L) ||
        !all(file.exists(pressed_files))
      ) {

        diagnostics <- if (file.exists(stderr_file)) {
          readLines(
            stderr_file,
            warn = FALSE
          )
        } else {
          character()
        }

        if (length(diagnostics)) {
          .log_tool_output(
            log_path,
            paste0(
              "HMMER failure: ",
              db_name
            ),
            diagnostics
          )
        }

        stop(
          "hmmpress failed for ",
          db_name,
          ". See the processing log for diagnostics.",
          call. = FALSE
        )
      }
    }

    combined_hmm
  }

  ####################################################################
  # build separate databases
  ####################################################################

  defense_hmm <- .amr_progress_step(
  "Preparing DefenseFinder HMM database",
  build_database(
    defense_dir,
    "DefenseFinder"
  ),
  progress = progress,
  verbose = verbose,
  log_path = log_path
)

cas_hmm <- .amr_progress_step(
  "Preparing CasFinder HMM database",
  build_database(
    cas_dir,
    "CasFinder"
  ),
  progress = progress,
  verbose = verbose,
  log_path = log_path
)

  defense_bfc <- .amr_bfc_register_hmmer(
    database = "DefenseCas",
    component = "DefenseFinder",
    hmm_path = defense_hmm
  )

  cas_bfc <- .amr_bfc_register_hmmer(
    database = "DefenseCas",
    component = "CasFinder",
    hmm_path = cas_hmm
  )

  ####################################################################
  # load proteins
  ####################################################################

  prot_seqs <- local({
    con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
    on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

    DBI::dbReadTable(con, "protein_cluster_seq") |>
      tibble::as_tibble()
    }
    )

  fasta_file <- file.path(output_path, "protein_DefenseCas.faa")

  # required to define the database size for hmmsearch --Z and --domZ parameters
  total_proteins <- nrow(prot_seqs)

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

  combined_tbl <- purrr::imap_dfr(
    databases, function(hmm_file, db_name) {

      .log_or_message(log_path, verbose, "Running ", db_name)

      tbl_file <- file.path(output_path, paste0("protein_", db_name, ".tbl")
      )

      stderr_file <- tempfile(
        pattern = paste0(
          "hmmer_",
          db_name,
          "_"
        ),
        fileext = ".stderr"
      )

      status <- .amr_progress_step(
                  paste0("Running ", db_name, " HMMER"),
                  tryCatch(
                    system2(
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
                        "-Z",
                        total_proteins,
                        "--domZ",
                        total_proteins,
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
                      stdout = FALSE,
                      stderr = stderr_file
                    ),
                    error = function(e) {
                      stop(
                        "hmmsearch execution failed for ",
                        db_name,
                        ": ",
                        e$message
                      )
                    }
                  ),
                  progress = progress,
                  verbose = verbose,
                  log_path = log_path
                )

      if (
        !identical(status, 0L) ||
        !file.exists(tbl_file)
      ) {
        diagnostics <- if (file.exists(stderr_file)) {
          readLines(
            stderr_file,
            warn = FALSE
          )
        } else {
          character()
        }

        if (length(diagnostics)) {
          .log_tool_output(
            log_path,
            paste0(
              "HMMER failure: ",
              db_name
            ),
            diagnostics
          )
        }

        unlink(
          stderr_file,
          force = TRUE
        )

        stop(
          "hmmsearch failed for ",
          db_name,
          ". See the processing log for diagnostics.",
          call. = FALSE
        )
      }

      unlink(
        stderr_file,
        force = TRUE
      )

      .parseHMMEROutput(tbl_file) |>
        dplyr::select(
          protein,
          query_name
        ) |>
        dplyr::mutate(
          database = db_name
        ) |>
        dplyr::left_join(
          .parse_hmmer_profiles(hmm_file) |>
            dplyr::select(
              query_name = profile_name,
              query_accession = profile_accession,
              description = profile_description
            ),
          by = "query_name"
        )
    }
  )

  parquet_file <- file.path(output_path, "protein_DefenseCas.parquet")

  .write_compressed_parquet(combined_tbl, parquet_file)

  local({
    con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)

    on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

    DBI::dbWriteTable(con, "protein_DefenseCas", combined_tbl, overwrite = TRUE)
  })

  .log_or_message(log_path, verbose, "Created protein_DefenseCas")


  unlink(
    c(fasta_file, file.path(output_path, paste0("protein_", names(databases), ".tbl")))
  )

  invisible(list(
    databases = list(
      DefenseFinder = list(
        hmm = defense_hmm,
        bfc_rid = defense_bfc$rid,
        bfc_rname = defense_bfc$rname
      ),
      CasFinder = list(
        hmm = cas_hmm,
        bfc_rid = cas_bfc$rid,
        bfc_rname = cas_bfc$rname
      )
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
#'
#' @export
cleanMetaData <- function(duckdb_path, path = NULL) {
  duckdb_path <- normalizePath(duckdb_path, mustWork = TRUE)
  paths <- .amr_paths_from_duckdb(duckdb_path)
  if (is.null(path)) {
    path <- paths$orb
  }

  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path <- normalizePath(path, mustWork = TRUE)

  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)
  clean_countries <- cleaned_bvbrc_countries |>
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

  # What year is it?!
  current_year <- as.integer(format(Sys.Date(), "%Y"))
  last_year_break <- (floor(current_year / 5) + 1L) * 5L
  year_breaks <- seq(1980, last_year_break, by = 5)

  # 2025-era code doesn't handle 2026 very well, so updating above
  # year_breaks <- seq(1980, 2026, by = 5)
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

  # Final metadata Parquets
  metadata_parquet <- file.path(path, "metadata.parquet")
  amr_phenotype_parquet <- file.path(path, "amr_phenotype.parquet")
  genome_data_parquet <- file.path(path, "genome_data.parquet")
  original_metadata_parquet <- file.path(path, "original_metadata.parquet")
  metadata_qc_parquet <- file.path(path, "metadata_qc.parquet")
  metadata_qc_rejections_parquet <- file.path(path, "metadata_qc_rejections.parquet")
  selected_genomes_parquet <- file.path(path,"selected_genomes.parquet")

  writeCompressedParquet <- function(df, path) {
    arrow::write_parquet(
      df,
      path,
      compression = "zstd",
      compression_level = 9,
      use_dictionary = TRUE
    )
  }

  db_name <- file.path(path, paste0(tools::file_path_sans_ext(basename(duckdb_path)), "_parquet.duckdb"))

  con_new <- .amr_connect_dataset_db(db_name)
  on.exit(try(DBI::dbDisconnect(con_new), silent = TRUE), add = TRUE)

  # Cleaned analysis metadata
  DBI::dbReadTable(con,"cleaned_metadata") |>
    writeCompressedParquet(metadata_parquet)

  DBI::dbExecute(con_new,sprintf(
      "CREATE OR REPLACE VIEW metadata AS SELECT * FROM read_parquet('%s')",
      basename(metadata_parquet)
    )
  )

  # Original metadata tables retained in the ORB
  DBI::dbReadTable(con, "amr_phenotype") |>
    writeCompressedParquet(amr_phenotype_parquet)
  DBI::dbReadTable(con, "genome_data") |>
    writeCompressedParquet(genome_data_parquet)
  DBI::dbReadTable(con, "metadata") |>
    writeCompressedParquet(original_metadata_parquet)
  DBI::dbExecute(con_new, sprintf(
      "CREATE OR REPLACE VIEW amr_phenotype AS SELECT * FROM read_parquet('%s')",
      basename(amr_phenotype_parquet)
    )
  )

  DBI::dbExecute(con_new, sprintf(
      "CREATE OR REPLACE VIEW genome_data AS SELECT * FROM read_parquet('%s')",
      basename(genome_data_parquet)
    )
  )
  DBI::dbExecute(con_new, sprintf(
      "CREATE OR REPLACE VIEW original_metadata AS SELECT * FROM read_parquet('%s')",
      basename(original_metadata_parquet)
    )
  )

  # Metadata QC audit
  DBI::dbReadTable(con, "metadata_qc") |>
    writeCompressedParquet(metadata_qc_parquet)
  DBI::dbReadTable(con, "metadata_qc_rejections") |>
    writeCompressedParquet(metadata_qc_rejections_parquet)
  DBI::dbExecute(con_new, sprintf(
      "CREATE OR REPLACE VIEW metadata_qc AS SELECT * FROM read_parquet('%s')",
      basename(metadata_qc_parquet)
    )
  )
  DBI::dbExecute(con_new, sprintf(
      paste0(
        "CREATE OR REPLACE VIEW metadata_qc_rejections ",
        "AS SELECT * FROM read_parquet('%s')"
      ),
      basename(metadata_qc_rejections_parquet)
    )
  )

  # Snapshot of selected genomes
  selected_genomes <- DBI::dbReadTable(con, "filtered") |>
    tibble::as_tibble()

  if (!"genome.genome_id" %in% names(selected_genomes)) {
    stop("Table 'filtered' does not contain 'genome.genome_id'.")
  }

  selected_genomes <- selected_genomes |>
    dplyr::transmute(genome_id = as.character(.data[["genome.genome_id"]])) |>
    dplyr::filter(!is.na(genome_id), nzchar(genome_id)) |>
    dplyr::distinct() |>
    dplyr::arrange(genome_id)

  writeCompressedParquet(selected_genomes, selected_genomes_parquet)

  DBI::dbExecute(con_new,sprintf(
      "CREATE OR REPLACE VIEW selected_genomes AS SELECT * FROM read_parquet('%s')",
      basename(selected_genomes_parquet)
    )
  )

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
cleanData <- function(duckdb_path, path = NULL, verbose = TRUE) {
  duckdb_path <- normalizePath(duckdb_path, mustWork = TRUE)

  paths <- .amr_paths_from_duckdb(duckdb_path)

  if (is.null(path)) {
    path <- paths$orb
  }

  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  path <- normalizePath(path, mustWork = TRUE)

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
    stop("No successful HMMER stage found in manifest: ", manifest_path)
  }

  hmmer_databases <- unlist(
    hmmer_stage$parameters$databases
  )

  if (!length(hmmer_databases)) {
    stop(
      "HMMER stage in manifest does not contain any databases."
    )
  }


  .proteinAnnotations2Duckdb(
    duckdb_path = duckdb_path,
    databases = hmmer_databases,
    output_path = path,
    verbose = verbose
  )

  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

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

  db_name <- file.path(
    path,
    paste0(
      tools::file_path_sans_ext(
        basename(duckdb_path)
      ),
      "_parquet.duckdb"
    )
  )

  con_new <- .amr_connect_dataset_db(db_name)
  on.exit(try(DBI::dbDisconnect(con_new), silent = TRUE), add = TRUE)

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
  DBI::dbExecute(con_new, sprintf("CREATE OR REPLACE VIEW protein_members AS SELECT * FROM read_parquet('%s')", basename(protein_cluster_member_parquet)))

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
#'
#'   \item **Dyad feature mapping** via [buildDyadFeatureMap()] -> writes the
#'         protein-gene dyad feature network used for downstream graph analysis.
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
#' @param resume Logical. If `TRUE`, checks the dataset manifest for the most
#'   recent prior `runDataProcessing()` attempt and skips Panaroo, CD-HIT,
#'   and/or HMMER if they already completed successfully and their recorded
#'   output files are still present on disk. Only a contiguous run of
#'   successes from the start of the pipeline is honored (e.g. if CD-HIT
#'   failed, HMMER is always re-run even if it previously succeeded, since it
#'   depends on CD-HIT's output). Metadata cleaning and Parquet export always
#'   run, since they're fast and idempotent. Default: `FALSE`.
#'
#' @param export_tabular_data Logical. If TRUE, automatically call
#'   [exportProcessedData()] after the processing pipeline completes to create
#'   the default human-readable exports. Default: `FALSE`.
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
#' @param hmmer_db_dir Character. Directory containing the prepared HMMER
#'   databases. If `NULL`, the default BiocFileCache-managed HMMER directory
#'   is used.
#' @param hmmer_docker_image Character. Docker image containing HMMER.
#'   Default: `"staphb/hmmer"`.
#' @param hmmer_num_splits Integer. Number of protein-sequence chunks for HMMER.
#'   Default: `8`.
#' @param hmmer_workers Integer. Number of parallel HMMER workers. Default: `8`.
#'
#' @param verbose Logical. Print persistent status and diagnostic messages.
#'   Routine workflow output is suppressed when `FALSE`; warnings and errors are
#'   still reported. Default: `FALSE`.
#' @param progress Logical. Show temporary progress for long-running operations.
#'   This is independent of `verbose`. Defaults to `TRUE`, unless overridden by
#'   the `amRdata.progress` option.
#'
#' @return
#' Invisibly returns a list with:
#' \itemize{
#'   \item `duckdb_path` - input DuckDB path
#'   \item `panaroo_output` - path to the selected Panaroo output directory used for import
#'   \item `parquet_duckdb_path` - absolute path to the created Parquet-backed DuckDB
#'   \item `log_path` - path to a plain-text progress log, appended to as each
#'     stage starts and finishes. Tail it (e.g. `tail -f <log_path>`) to watch
#'     progress without parsing the full JSON provenance manifest.
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
#' * Creates the protein-gene dyad feature map for downstream graph analysis.
#' * If `export_tabular_data = TRUE`, exports features in human-readable CSVs.
#'   Additional export options are available by calling `exportProcessedData()`
#'   after `runDataProcessing()` completes.
#'
#' **Threading**
#' * `threads` provides the shared CPU budget for the major processing stages.
#' * Panaroo, CD-HIT, and HMMER allocate that budget according to their respective
#'   stage parameters.
#'
#' **Resuming a Failed Run**
#' * Set `resume = TRUE` to avoid re-running stages that already completed
#'   successfully in the most recent prior attempt (per the dataset manifest).
#' * Resuming is driven entirely by the manifest plus a check that the
#'   previously recorded output files still exist --- it does not re-validate
#'   the contents of those outputs, so don't rely on it if files may have
#'   been altered or deleted since the failed run.
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
#'   threads       = 8
#' )
#'
#' # After completion:
#' # data/Shigella_flexneri/Sfl_parquet.duckdb
#' # will contain views over the Parquet files for downstream ML.
#' }
#'
#' @export
runDataProcessing <- function(
    duckdb_path = NULL,
    threads = 8,
    resume = FALSE,
    export_tabular_data = FALSE,


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
    verbose = FALSE,
    progress = getOption("amRdata.progress", TRUE)
) {
  panaroo_refind_mode <- match.arg(panaroo_refind_mode)

  # User supplied nothing? Go fetch what's eligible to run
  if (is.null(duckdb_path)) {
    duckdb_path <- .amr_select_processing_dataset()
  }

  requested_threads <- threads
  threads <- .resolve_workers(requested = threads)
  duckdb_path <- normalizePath(duckdb_path, mustWork = TRUE)
  paths <- .amr_paths_from_duckdb(duckdb_path)

  purrr::walk(
    c(
      paths$panaroo,
      paths$cdhit,
      paths$hmmer,
      paths$orb
    ),
    ~ dir.create(
      .x,
      recursive = TRUE,
      showWarnings = FALSE
    )
  )

  log_path <- paths$processing_log

  if (file.exists(log_path)) {
    unlink(log_path, force = TRUE)
  }

  .log_write(log_path, "Started runDataProcessing().")

  if (isTRUE(verbose)) {
    message(
      "Individual tool output will be written to:\n",
      "  ",
      normalizePath(log_path, mustWork = FALSE)
    )
  }



  status <- function(...) {
    .amr_status(
      ...,
      verbose = verbose,
      log_path = log_path
    )
  }

  status("Extracting features...")

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
    base_dir = dirname(dirname(paths$root)),
    hash_files = FALSE
  )

  # Work out which stages can be skipped, based on the most recent prior run
  stage_order <- c("panaroo", "cdhit", "hmmer")
  prev_run <- if (isTRUE(resume) && manifest$run_index > 1L) {
    manifest$manifest$runs[[manifest$run_index - 1L]]
  } else {
    NULL
  }
  resume_completed <- .resume_plan(prev_run, stage_order)

  if (any(resume_completed)) {
    status(
      "Resuming: skipping already-completed stage(s): ",
      paste(stage_order[resume_completed], collapse = ", ")
    )
  }

  run_failed <- TRUE

  on.exit(
    if (run_failed) {
      status("FAILED: runDataProcessing() exited before successful completion.")
      .manifest_finish(
        manifest,
        status = "failed",
        error = "runDataProcessing() exited before successful completion."
      )
    },
    add = TRUE
  )

  processing_run_id <-
    manifest$manifest$runs[[manifest$run_index]]$run_id

  manifest <- .manifest_artifact(
    manifest,
    name = "amRml_input",
    status = "building",
    details = list(
      producer = "amRdata",
      producer_run_id = processing_run_id
    )
  )

  # Record the start of this processing run
  manifest <- .manifest_event(
    manifest,
    message = "Started data-processing run.",
    details = list(
      duckdb_path = duckdb_path,
      panaroo_dir = paths$panaroo,
      cdhit_dir = paths$cdhit,
      hmmer_dir = paths$hmmer,
      orb_dir = paths$orb
    )
  )

  # 1) Panaroo (run + optional merge) -> write Panaroo tables
  if (resume_completed[["panaroo"]]) {
    status("Skipping Panaroo (resume): reusing output from the previous run.")

    pan_dir <- .manifest_prior_stage(prev_run, "panaroo")$outputs[[1]]$path

    manifest <- .manifest_stage(
      manifest,
      name = "panaroo",
      status = "skipped",
      inputs = duckdb_path,
      outputs = c(
        pan_dir,
        duckdb_path
      ),
      tool = list(
        name = "Panaroo",
        docker_image = "staphb/panaroo:1.7.0"
      ),
      message = "Resumed: reused successful output from a previous run."
    )
  } else {
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
        threads_requested = requested_threads,
        threads_used = threads,
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

    pan_dir <- .amr_progress_step(
      "Running Panaroo",
      runPanaroo2Duckdb(
        duckdb_path            = duckdb_path,
        output_path            = paths$panaroo,
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
        verbose                = verbose,
        log_path               = log_path
      ),
      progress = progress,
      verbose = verbose,
      log_path = log_path
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
  }

  # 2) CD-HIT -> write `protein` tables
  if (resume_completed[["cdhit"]]) {
    status("Skipping CD-HIT (resume): reusing output from the previous run.")

    manifest <- .manifest_stage(
      manifest,
      name = "cdhit",
      status = "skipped",
      inputs = duckdb_path,
      outputs = c(
        file.path(paths$cdhit, paste0(cdhit_output_prefix, "_input.fa")),
        file.path(paths$cdhit, paste0(cdhit_output_prefix, "_proteins")),
        file.path(duckdb_path)
      ),
      tool = list(
        name = "CD-HIT",
        version = "4.8.1",
        docker_image = "weizhongli1987/cdhit:4.8.1"
      ),
      message = "Resumed: reused successful output from a previous run."
    )
  } else {
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

    .amr_progress_step(
      "Running CD-HIT",
      CDHIT2duckdb(
        duckdb_path   = duckdb_path,
        output_path   = paths$cdhit,
        output_prefix = cdhit_output_prefix,
        identity      = cdhit_identity,
        word_length   = cdhit_word_length,
        threads       = threads,
        memory        = cdhit_memory,
        extra_args    = cdhit_extra_args,
        verbose       = verbose,
        log_path      = log_path
      ),
      progress = progress,
      verbose = verbose,
      log_path = log_path
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
        file.path(paths$cdhit, paste0(cdhit_output_prefix, "_input.fa")),
        file.path(paths$cdhit, paste0(cdhit_output_prefix, "_proteins")),
        duckdb_path
      ),
      tool = list(
        name = "CD-HIT",
        version = "4.8.1",
        docker_image = "weizhongli1987/cdhit:4.8.1"
      )
    )
  }

  # 3) HMMER -> write HMM-based match tables for desired databases
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

  hmmer_result <- NULL
  defense_result <- NULL

  if (resume_completed[["hmmer"]]) {
    status("Skipping HMMER (resume): reusing output from the previous run.")

    manifest <- .manifest_stage(
      manifest,
      name = "hmmer",
      status = "skipped",
      inputs = duckdb_path,
      tool = list(
        name = "HMMER",
        docker_image = hmmer_docker_image
      ),
      message = "Resumed: reused successful output from a previous run."
    )
  } else {
    status("Running HMMER with databases: ", paste(hmmer_databases, collapse = ", "))

    manifest <- .manifest_stage(
      manifest,
      name = "hmmer",
      status = "running",
      parameters = list(
        databases = hmmer_databases,
        database_dir = hmmer_db_dir,
        database_cache = "BiocFileCache",
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

    if (length(generic_databases)) {
      hmmer_result <- .runHMMER(
        duckdb_path = duckdb_path,
        output_path = paths$hmmer,
        threads = threads,
        hmmer_db_dir = hmmer_db_dir,
        databases = generic_databases,
        docker_image = hmmer_docker_image,
        num_of_splits = hmmer_num_splits,
        n_workers = hmmer_workers,
        verbose = verbose,
        log_path = log_path,
        progress = progress
      )
    }

    if ("DefenseCas" %in% hmmer_databases) {
  defense_result <- .defenseHMMER(
    defense_db_dir = file.path(
      hmmer_db_dir,
      "DefenseCas"
    ),
    docker_image = hmmer_docker_image,
    duckdb_path = duckdb_path,
    output_path = paths$hmmer,
    threads = threads,
    verbose = verbose,
    log_path = log_path,
    progress = progress
  )
}
  }

  # Verify expected outputs regardless of whether HMMER just ran or was skipped
  expected_outputs <- file.path(
    paths$hmmer,
    paste0(
      "protein_",
      hmmer_databases,
      ".parquet"
    )
  )

  missing_outputs <- expected_outputs[!file.exists(expected_outputs)]

  if (length(missing_outputs)) {
    stop(
      "HMMER did not produce all expected outputs:\n",
      paste(missing_outputs, collapse = "\n")
    )
  }

  expected_tables <- paste0("protein_", hmmer_databases)

  missing_tables <- local({
    con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)

    on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

    expected_tables[
      !vapply(
        expected_tables,
        function(tbl) {
          DBI::dbExistsTable(con, tbl)
        },
        logical(1)
      )
    ]
  })

  if (length(missing_tables)) {
    stop(
      "HMMER did not produce all expected DuckDB tables:\n",
      paste(missing_tables, collapse = "\n")
    )
  }

  if (!resume_completed[["hmmer"]]) {
    manifest <- .manifest_stage(
      manifest,
      name = "hmmer",
      status = "success",
      parameters = list(
        databases = hmmer_databases,
        database_dir = hmmer_db_dir,
        database_cache = "BiocFileCache",
        docker_image = hmmer_docker_image,
        threads = threads,
        num_of_splits = hmmer_num_splits,
        workers = hmmer_workers
      ),
      inputs = duckdb_path,
      outputs = c(
        file.path(
          paths$hmmer,
          paste0(
            "protein_",
            hmmer_databases,
            ".parquet"
          )
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
    status("Finished HMMER.")
  }

  # 4) Clean metadata and export Parquet + Parquet-backed DuckDB
  .amr_progress_step("Building resource bundle", {
      cleanMetaData(duckdb_path = duckdb_path, path = paths$orb)

      cleanData(
        duckdb_path = duckdb_path,
        path = paths$orb,
        verbose = verbose
      )
    },
    progress = progress,
    verbose = verbose,
    log_path = log_path
  )

  parquet_duckdb_path <- paths$parquet_duckdb

  dyad_parquet <- .amr_progress_step(
    "Building protein-gene dyad map",
    buildDyadFeatureMap(
      duckdb_path = duckdb_path,
      output_path = paths$orb
    ),
    progress = progress,
    verbose = verbose,
    log_path = log_path
  )

  local({
    con_orb <- .amr_connect_dataset_db(
      paths$parquet_duckdb
    )

    on.exit(
      try(
        DBI::dbDisconnect(con_orb),
        silent = TRUE
      ),
      add = TRUE
    )

    DBI::dbExecute(
      con_orb,
      sprintf(
        "CREATE OR REPLACE VIEW dyad_feature AS SELECT * FROM read_parquet('%s')",
        basename(dyad_parquet)
      )
    )
  })

  parquet_duckdb_path <- normalizePath(
    parquet_duckdb_path,
    mustWork = TRUE
  )

  if (isTRUE(export_tabular_data)) {
    .amr_progress_step(
      "Exporting processed data",
      exportProcessedData(
        duckdb_path = paths$parquet_duckdb,
        output_path = paths$exports,
        verbose = verbose
      ),
      progress = progress,
      verbose = verbose,
      log_path = log_path
    )
  }

  expected_orb_contents <- c(
    "metadata",
    "amr_phenotype",
    "genome_data",
    "original_metadata",
    "metadata_qc",
    "metadata_qc_rejections",
    "selected_genomes",

    "gene_count",
    "gene_names",
    "gene_seqs",
    "genome_gene_protein",
    "struct",

    "protein_count",
    "protein_names",
    "protein_seqs",
    "protein_members",

    paste0(
      "protein_",
      hmmer_databases
    ),
    paste0(
      "protein_",
      hmmer_databases,
      "_count"
    ),

    "dyad_feature"
  )

  local({
    con_orb <- .amr_connect_dataset_db(
      parquet_duckdb_path,
      read_only = TRUE
    )

    on.exit(
      try(
        DBI::dbDisconnect(con_orb),
        silent = TRUE
      ),
      add = TRUE
    )

    available_relations <- DBI::dbListTables(
      con_orb
    )

    missing_relations <- setdiff(
      expected_orb_contents,
      available_relations
    )

    if (length(missing_relations)) {
      stop(
        "ORB validation failed. Missing expected relation(s):\n",
        paste(
          missing_relations,
          collapse = "\n"
        )
      )
    }

    purrr::walk(
      expected_orb_contents,
      function(relation) {
        relation_sql <- DBI::dbQuoteIdentifier(
          con_orb,
          relation
        )

        DBI::dbGetQuery(
          con_orb,
          paste0(
            "SELECT * FROM ",
            relation_sql,
            " LIMIT 0"
          )
        )
      }
    )
  })

  # Final Parquets generated by this run
  parquet_files <- c(
    # Metadata and QC outputs
    file.path(paths$orb, "metadata.parquet"),
    file.path(paths$orb, "amr_phenotype.parquet"),
    file.path(paths$orb, "genome_data.parquet"),
    file.path(paths$orb, "original_metadata.parquet"),
    file.path(paths$orb, "metadata_qc.parquet"),
    file.path(paths$orb, "metadata_qc_rejections.parquet"),
    file.path(paths$orb, "selected_genomes.parquet"),

    # Core feature outputs
    file.path(paths$orb, "gene_count.parquet"),
    file.path(paths$orb, "gene_names.parquet"),
    file.path(paths$orb, "gene_seqs.parquet"),
    file.path(paths$orb, "genome_gene_protein.parquet"),
    file.path(paths$orb, "struct.parquet"),
    file.path(paths$orb, "protein_count.parquet"),
    file.path(paths$orb, "protein_names.parquet"),
    file.path(paths$orb, "protein_seqs.parquet"),
    file.path(paths$orb, "protein_members.parquet"),

    # HMMER outputs
    file.path(paths$orb, paste0("protein_", hmmer_databases, ".parquet")),
    file.path(
      paths$orb,
      paste0("protein_", hmmer_databases, "_count.parquet")
    ),

    # Dyad feature map
    file.path(paths$orb, "dyad_feature.parquet")
  )

  missing_parquet_files <- parquet_files[!file.exists(parquet_files)]

  if (length(missing_parquet_files)) {
    stop(
      "ORB build is incomplete. Missing expected file(s):\n",
      paste(missing_parquet_files, collapse = "\n")
    )
  }

  if (!file.exists(parquet_duckdb_path)) {
    stop(
      "ORB build is incomplete. Parquet-backed DuckDB was not created: ",
      parquet_duckdb_path
    )
  }

  parquet_files <- normalizePath(parquet_files, mustWork = TRUE)

  parquet_duckdb_path <- normalizePath(parquet_duckdb_path, mustWork = TRUE)

  # Log!
  manifest <- .manifest_stage(
    manifest,
    name = "clean_metadata_and_export",
    status = "success",
    parameters = list(
      reference_data = c(
        "clean_drug",
        "drug_class",
        "drug_abbr",
        "class_abbr",
        "cleaned_bvbrc_countries"
      )
    ),
    inputs = duckdb_path,
    outputs = c(parquet_files, parquet_duckdb_path),
    metrics = list(parquet_duckdb = parquet_duckdb_path, parquet_files = parquet_files)
  )

  # Labeling that this run is ready for amRml in the next package
  manifest <- .manifest_artifact(
    manifest,
    name = "amRml_input",
    status = "ready",
    details = list(
      producer = "amRdata",
      producer_run_id = processing_run_id,
      path_mode = "relative_to_manifest",
      directory = ".",
      parquet_duckdb = basename(
        paths$parquet_duckdb
      ),
      metadata_parquet = "metadata.parquet",
      metadata_qc_parquet = "metadata_qc.parquet",
      metadata_qc_rejections_parquet =
        "metadata_qc_rejections.parquet",
      selected_genomes_parquet =
        "selected_genomes.parquet"
    )
  )

  manifest <- .manifest_finish(manifest, status = "success")

  status("Features ready.")
  message("Data processing complete.")

  run_failed <- FALSE

  invisible(
    list(
      duckdb_path = duckdb_path,
      panaroo_output = pan_dir,
      parquet_duckdb_path = normalizePath(parquet_duckdb_path),
      log_path = log_path
    )
  )
}



#' Export processed data from an amRdata ORB
#'
#' Exports processed relations from the resource bundle (ORB) produced by
#' [runDataProcessing()] as CSV, TSV, Parquet, and/or XLSX. This is an optional
#' step that allows users to take their processed data outside our amR workflow
#' for use in their own custom analyses. This is not required to run `amRml`!
#'
#' If `duckdb_path` is `NULL`, available ORBs are discovered automatically. A
#' single ORB is selected automatically; when multiple ORBs are available in an
#' interactive session, a selection menu is shown.
#'
#' @param duckdb_path Character or NULL. Path to an amRdata ORB directory,
#'   ORB DuckDB, or working DuckDB associated with a completed ORB. If NULL,
#'   available ORBs are discovered automatically.
#' @param output_path Character or NULL. Directory for exports. Defaults to the
#'   dataset `exports/` directory.
#' @param amr_phenotype_mode Character. One of "separate" or "append".
#'   "separate" exports the AMR labels as a separate wide table.
#'   "append" also joins those labels onto genome-level feature count tables.
#' @param export_formats Character vector. Any of "csv", "tsv", "parquet", "xlsx".
#' @param export_sequences Logical. If TRUE, also exports gene and protein
#'   sequence and mapping relations detected in the ORB. Default FALSE.
#' @param export_dyads Logical. If TRUE, exports the optional dyad annotation
#'   table. Each row represents a protein-gene dyad with semicolon-separated
#'   mapped feature values. The export uses `export_formats`. Default FALSE.
#' @param tables Character vector or NULL. Relations to export. If NULL,
#'   exportable relations are discovered automatically from the ORB.
#' @param export_tables Logical. If TRUE, write the selected tables to disk.
#'   Default TRUE.
#' @param verbose Logical. If TRUE, prints progress messages.
#'
#' @return Invisibly returns a list containing the ORB path, export path,
#'   exported table names, and export settings.
#' @export
exportProcessedData <- function(duckdb_path = NULL,
                                output_path = NULL,
                                amr_phenotype_mode = c("separate", "append"),
                                export_formats = c("csv"),
                                export_sequences = FALSE,
                                export_dyads = FALSE,
                                tables = NULL,
                                export_tables = TRUE,
                                verbose = TRUE) {

  requested_duckdb_path <- duckdb_path

  duckdb_path <- .amr_resolve_export_orb(
    duckdb_path
  )

  paths <- .amr_paths_from_dataset_db(
    duckdb_path
  )

  if (isTRUE(verbose)) {
    message(
      "Using processed ORB: ",
      duckdb_path
    )
  }

  if (length(amr_phenotype_mode) > 1L &&
      isTRUE(verbose)) {
    message(
      "`amr_phenotype_mode` not specified; defaulting to 'separate'."
    )
  }

  amr_phenotype_mode <- match.arg(
    amr_phenotype_mode
  )

  export_formats <- unique(
    tolower(export_formats)
  )

  export_formats[
    export_formats == "excel"
  ] <- "xlsx"

  allowed_formats <- c(
    "csv",
    "tsv",
    "xlsx",
    "parquet"
  )

  unknown_formats <- setdiff(
    export_formats,
    allowed_formats
  )

  if (length(unknown_formats)) {
    stop(
      "Unsupported export format(s): ",
      paste(
        unknown_formats,
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  if (
    isTRUE(export_tables) &&
    !length(export_formats)
  ) {
    stop(
      "At least one export format must be supplied when export_tables = TRUE.",
      call. = FALSE
    )
  }

  warn_text_exports <- any(
    export_formats %in% c(
      "csv",
      "tsv",
      "xlsx"
    )
  )

  if (
    isTRUE(export_tables) &&
    warn_text_exports &&
    isTRUE(verbose)
  ) {

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

  if (
    "xlsx" %in% export_formats &&
    !requireNamespace(
      "writexl",
      quietly = TRUE
    )
  ) {
    stop(
      "Format 'xlsx' was requested but package 'writexl' is not available.",
      call. = FALSE
    )
  }

  if (is.null(output_path)) {
    output_path <- paths$exports
  }

  dir.create(
    output_path,
    recursive = TRUE,
    showWarnings = FALSE
  )

  output_path <- normalizePath(
    output_path,
    mustWork = TRUE
  )

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

  available_tables <- sort(
    DBI::dbListTables(con)
  )

  if (!length(available_tables)) {
    stop(
      "No relations found in ORB: ",
      duckdb_path,
      call. = FALSE
    )
  }

  read_tbl <- function(tbl) {
    tibble::as_tibble(
      DBI::dbReadTable(
        con,
        tbl
      )
    )
  }

  write_one <- function(df, stem) {

    df <- .preserve_export_id_text(
      df
    )

    if ("csv" %in% export_formats) {
      utils::write.table(
        df,
        file = file.path(
          output_path,
          paste0(
            stem,
            ".csv"
          )
        ),
        sep = ",",
        row.names = FALSE,
        col.names = TRUE,
        quote = TRUE,
        na = "",
        qmethod = "double",
        fileEncoding = "UTF-8"
      )
    }

    if ("tsv" %in% export_formats) {
      utils::write.table(
        df,
        file = file.path(
          output_path,
          paste0(
            stem,
            ".tsv"
          )
        ),
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = TRUE,
        na = "",
        qmethod = "double",
        fileEncoding = "UTF-8"
      )
    }

    if ("parquet" %in% export_formats) {
      arrow::write_parquet(
        df,
        file.path(
          output_path,
          paste0(
            stem,
            ".parquet"
          )
        )
      )
    }

    if ("xlsx" %in% export_formats) {
      writexl::write_xlsx(
        list(
          data = df
        ),
        file.path(
          output_path,
          paste0(
            stem,
            ".xlsx"
          )
        )
      )
    }
  }

  # ID HMMER outputs from ORB
  hmmer_count_tables <- grep(
    "^protein_.+_count$",
    available_tables,
    value = TRUE
  )

  hmmer_databases <- sub(
    "^protein_(.+)_count$",
    "\\1",
    hmmer_count_tables
  )

  hmmer_databases <- hmmer_databases[
    paste0(
      "protein_",
      hmmer_databases
    ) %in% available_tables
  ]

  hmmer_databases <- unique(
    hmmer_databases
  )

  if (
    length(hmmer_databases) &&
    isTRUE(verbose)
  ) {
    message(
      "HMMER outputs detected in ORB: ",
      paste(
        hmmer_databases,
        collapse = ", "
      )
    )
  }

  # AMR phenotype table
  build_amr_wide <- function() {

    source_tbl <- if (
      "amr_phenotype" %in% available_tables
    ) {
      "amr_phenotype"
    } else if (
      "metadata" %in% available_tables
    ) {
      "metadata"
    } else {
      NULL
    }

    if (is.null(source_tbl)) {
      return(NULL)
    }

    md <- read_tbl(
      source_tbl
    )

    genome_col <- c(
      "genome_drug.genome_id",
      "genome.genome_id",
      "genome_id"
    )

    genome_col <- genome_col[
      genome_col %in% names(md)
    ]

    antibiotic_col <- c(
      "genome_drug.antibiotic",
      "antibiotic"
    )

    antibiotic_col <- antibiotic_col[
      antibiotic_col %in% names(md)
    ]

    phenotype_col <- c(
      "genome_drug.resistant_phenotype",
      "phenotype"
    )

    phenotype_col <- phenotype_col[
      phenotype_col %in% names(md)
    ]

    if (
      !length(genome_col) ||
      !length(antibiotic_col) ||
      !length(phenotype_col)
    ) {
      return(NULL)
    }

    genome_col <- genome_col[[1]]
    antibiotic_col <- antibiotic_col[[1]]
    phenotype_col <- phenotype_col[[1]]

    md |>
      dplyr::transmute(
        genome_id = as.character(
          .data[[genome_col]]
        ),
        antibiotic = as.character(
          .data[[antibiotic_col]]
        ),
        phenotype = as.character(
          .data[[phenotype_col]]
        )
      ) |>
      dplyr::filter(
        !is.na(genome_id),
        nzchar(genome_id),
        !is.na(antibiotic),
        nzchar(antibiotic),
        !is.na(phenotype),
        nzchar(phenotype)
      ) |>
      dplyr::distinct() |>
      dplyr::group_by(
        genome_id,
        antibiotic
      ) |>
      dplyr::summarise(
        phenotype = paste(
          sort(
            unique(
              phenotype
            )
          ),
          collapse = ";"
        ),
        .groups = "drop"
      ) |>
      tidyr::pivot_wider(
        names_from = antibiotic,
        values_from = phenotype,
        values_fill = NA_character_
      ) |>
      dplyr::arrange(
        genome_id
      )
  }

  phenotype_wide <- build_amr_wide()

  if (!is.null(phenotype_wide)) {
    phenotype_wide <- .preserve_export_id_text(
      phenotype_wide
    )
  }

  # Determine the exports by beholding the ORB
  sequence_tables <- c(
    "gene_seqs",
    "protein_seqs",
    "protein_members",
    "genome_gene_protein"
  )

  default_exclusions <- c(
    "amr_phenotype",
    sequence_tables,
    "dyad_feature"
  )

  if (is.null(tables)) {

    selected_tables <- setdiff(
      available_tables,
      default_exclusions
    )

    if (isTRUE(export_sequences)) {
      selected_tables <- unique(c(
        selected_tables,
        intersect(
          sequence_tables,
          available_tables
        )
      ))
    }

    selected_tables <- unique(c(
      selected_tables,
      "amr_phenotype_wide"
    ))

  } else {

    requested_tables <- unique(
      as.character(tables)
    )

    selectable_tables <- c(
      available_tables,
      "amr_phenotype_wide"
    )

    missing_tables <- setdiff(
      requested_tables,
      selectable_tables
    )

    if (
      length(missing_tables) &&
      isTRUE(verbose)
    ) {
      message(
        "Requested relation(s) not present in ORB: ",
        paste(
          missing_tables,
          collapse = ", "
        )
      )
    }

    selected_tables <- intersect(
      requested_tables,
      selectable_tables
    )
  }

  if (!length(selected_tables)) {
    stop(
      "No requested relations were found in the ORB.",
      call. = FALSE
    )
  }

  exported <- character()

  # Optional dyad annotations
  if (isTRUE(export_dyads)) {

    dyad_tbl <- .exportDyadAnnotations(
      duckdb_path = duckdb_path,
      verbose = verbose
    )

    if (isTRUE(export_tables)) {

      write_one(
        dyad_tbl,
        "dyad_annotations"
      )

      exported <- c(
        exported,
        "dyad_annotations"
      )

      if (isTRUE(verbose)) {
        message(
          "Exported: dyad_annotations"
        )
      }
    }
  }

  for (table in selected_tables) {

    if (identical(
      table,
      "amr_phenotype_wide"
    )) {

      if (is.null(phenotype_wide)) {

        if (isTRUE(verbose)) {
          message(
            "Skipping amr_phenotype_wide: no AMR source relation found."
          )
        }

        next
      }

      df <- phenotype_wide

    } else {

      df <- .preserve_export_id_text(
        read_tbl(
          table
        )
      )
    }

    out_stem <- table

    appendable <- table %in% c(
      "gene_count",
      "protein_count",
      "struct"
    ) ||
      grepl(
        "^protein_.+_count$",
        table
      )

    if (
      identical(
        amr_phenotype_mode,
        "append"
      ) &&
      isTRUE(appendable) &&
      !is.null(phenotype_wide) &&
      "genome_id" %in% names(df)
    ) {

      df <- dplyr::left_join(
        df,
        phenotype_wide,
        by = "genome_id"
      )

      df <- .preserve_export_id_text(
        df
      )

      out_stem <- paste0(
        table,
        "_with_phenotypes"
      )
    }

    if (isTRUE(export_tables)) {

      write_one(
        df,
        out_stem
      )

      exported <- c(
        exported,
        out_stem
      )

      if (isTRUE(verbose)) {
        message(
          "Exported: ",
          out_stem
        )
      }
    }
  }

  manifest_path <- .manifest_find_latest(
    duckdb_path
  )

  invisible(list(
    requested_duckdb_path = requested_duckdb_path,
    duckdb_path = duckdb_path,
    output_path = output_path,
    tables = exported,
    amr_phenotype_mode = amr_phenotype_mode,
    export_formats = export_formats,
    export_sequences = isTRUE(export_sequences),
    hmmer_databases = hmmer_databases,
    manifest_path = manifest_path
  ))
}

#' Remove local amRdata dataset files
#'
#' Removes disposable intermediate files from an amRdata dataset while
#' retaining the final ORB and any human-readable exports.
#'
#' If `dataset_path` is not supplied in an interactive session, registered
#' amRdata dataset manifests are used to identify available datasets. A single
#' available dataset is selected automatically; when multiple datasets are
#' available, an interactive menu is shown.
#'
#' With `complete_remove = TRUE`, the entire dataset directory is permanently
#' removed after strict structural validation and interactive confirmation. If your
#' directory structure does not conform to expectations, you will have to manually
#' delete your data. This is a safety mechanism so you don't accidentally delete
#' the photos of your children on your desktop.
#'
#' @param dataset_path Character scalar or `NULL`. Path to the amRdata dataset
#'   directory, for example `"data/Shigella_flexneri"`. If `NULL`
#'   interactively, registered dataset manifests are used to select a dataset.
#' @param complete_remove Logical. If `FALSE` (default), remove only disposable
#'   build directories and retain `orb/` and `exports/`. If `TRUE`, permanently
#'   remove the entire validated amRdata dataset directory.
#' @param verbose Logical. Print cleanup messages. Default `TRUE`.
#'
#' @return Invisibly returns the removed paths. For a cancelled selection or
#'   complete removal, invisibly returns `FALSE`.
#'
#' @export
removeLocalFiles <- function(
    dataset_path = NULL,
    complete_remove = FALSE,
    verbose = TRUE
) {

  if (is.null(dataset_path)) {
  candidates <- .amr_processing_datasets()

  if (!nrow(candidates)) {
    stop(
      "No prepared amRdata datasets were found.",
      call. = FALSE
    )
  }

  if (nrow(candidates) == 1L) {
    duckdb_path <- candidates$duckdb_path[[1]]
  } else {
    if (!interactive()) {
      stop(
        "`dataset_path` must be supplied when multiple datasets are available ",
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
      title = "Select an amRdata dataset to clean:"
    )

    if (selection == 0L) {
      return(invisible(FALSE))
    }

    duckdb_path <- candidates$duckdb_path[[selection]]
  }

  dataset_path <- dirname(
    dirname(duckdb_path)
  )
}

  if (
    length(complete_remove) != 1L ||
    is.na(complete_remove) ||
    !is.logical(complete_remove)
  ) {
    stop(
      "`complete_remove` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  dataset_path <- normalizePath(
    dataset_path,
    mustWork = TRUE
  )

  if (!dir.exists(dataset_path)) {
    stop(
      "`dataset_path` must be an amRdata dataset directory.",
      call. = FALSE
    )
  }

  data_dir <- dirname(dataset_path)

  if (!identical(basename(data_dir), "data")) {
    stop(
      "Removal refused: dataset directory is not directly ",
      "inside a directory named 'data/'.",
      call. = FALSE
    )
  }

  orb_dir <- file.path(
    dataset_path,
    "orb"
  )

  if (!dir.exists(orb_dir)) {
    stop(
      "Dataset does not contain an 'orb/' directory: ",
      dataset_path,
      call. = FALSE
    )
  }

  # Calculate disk use without relying on platform-specific shell commands
  path_bytes <- function(path) {
    if (
      !file.exists(path) &&
      !dir.exists(path)
    ) {
      return(0)
    }

    if (!dir.exists(path)) {
      size <- file.info(path)$size

      if (is.na(size)) {
        return(0)
      }

      return(
        as.numeric(size)
      )
    }

    files <- list.files(
      path,
      all.files = TRUE,
      full.names = TRUE,
      recursive = TRUE,
      include.dirs = FALSE,
      no.. = TRUE
    )

    if (!length(files)) {
      return(0)
    }

    info <- file.info(files)

    sum(as.numeric(info$size[!info$isdir &
                               !is.na(info$size)]), na.rm = TRUE)
  }

  format_bytes <- function(bytes) {
    bytes <- as.numeric(bytes)

    if (
      length(bytes) != 1L ||
      is.na(bytes) ||
      !is.finite(bytes) ||
      bytes < 0
    ) {
      return("unknown")
    }

    units <- c(
      "B",
      "KiB",
      "MiB",
      "GiB",
      "TiB"
    )

    if (bytes == 0) {
      return("0 B")
    }

    unit_index <- min(
      floor(
        log(
          bytes,
          base = 1024
        )
      ) + 1L,
      length(units)
    )

    value <- bytes / 1024^(unit_index - 1L)

    if (unit_index == 1L) {
      sprintf(
        "%.0f %s",
        value,
        units[[unit_index]]
      )
    } else {
      sprintf(
        "%.2f %s",
        value,
        units[[unit_index]]
      )
    }
  }

  if (isTRUE(complete_remove)) {
    # Reject filesystem roots and especially important directories even if
    # someone has created an unusual path named "data"
    protected_paths <- unique(c(
      normalizePath(path.expand("~"), mustWork = FALSE),
      normalizePath(getwd(), mustWork = FALSE),
      normalizePath(data_dir, mustWork = FALSE)
    ))

    if (identical(dirname(dataset_path), dataset_path) ||
        dataset_path %in% protected_paths) {
      stop("Complete removal refused: protected path.", call. = FALSE)
    }

    # Do not perform destructive recursive deletion through a symlinked
    # dataset root
    root_link <- Sys.readlink(dataset_path)

    if (length(root_link) &&
        !is.na(root_link) &&
        nzchar(root_link)
    ) {
      stop("Complete removal refused: dataset directory is a symbolic link.",
           call. = FALSE)
    }

    # The top level must look exactly like an amRdata dataset. Any foreign
    # material causes a hard stop rather than being silently deleted.
    allowed_entries <- c(
      "genomes",
      "panaroo",
      "cd-hit",
      "hmmer",
      "work",
      "orb",
      "exports",
      ".DS_Store"
    )

    top_entries <- list.files(
      dataset_path,
      all.files = TRUE,
      no.. = TRUE
    )

    unexpected_entries <- setdiff(
      top_entries,
      allowed_entries
    )

    if (length(unexpected_entries)) {
      stop(
        "Complete removal refused: unexpected top-level file(s) or ",
        "directory/directories were found:\n",
        paste(
          unexpected_entries,
          collapse = "\n"
        ),
        call. = FALSE
      )
    }

    # Reject symlinked top-level dataset components as an additional safety
    top_paths <- file.path(
      dataset_path,
      setdiff(
        top_entries,
        ".DS_Store"
      )
    )

    if (length(top_paths)) {
      link_targets <- Sys.readlink(
        top_paths
      )

      linked_entries <- basename(
        top_paths[
          !is.na(link_targets) &
            nzchar(link_targets)
        ]
      )

      if (length(linked_entries)) {
        stop(
          "Complete removal refused: symbolic link(s) found inside ",
          "the dataset root:\n",
          paste(
            linked_entries,
            collapse = "\n"
          ),
          call. = FALSE
        )
      }
    }

    # A real amRdata dataset must contain at least one valid amR dataset
    # manifest. For complete removal we intentionally do not require a
    # successful run, because failed/incomplete builds also need to be
    # removable
    manifest_paths <- list.files(
      orb_dir,
      pattern = "^manifest_.*\\.json$",
      full.names = TRUE
    )

    if (!length(manifest_paths)) {
      stop(
        "Complete removal refused: no amRdata manifest was found in 'orb/'.",
        call. = FALSE
      )
    }

    manifest_paths <- manifest_paths[
      order(
        file.info(manifest_paths)$mtime,
        decreasing = TRUE
      )
    ]

    manifest <- NULL
    manifest_path <- NULL

    for (candidate_path in manifest_paths) {
      candidate <- tryCatch(
        jsonlite::read_json(
          candidate_path,
          simplifyVector = FALSE
        ),
        error = function(e) {
          NULL
        }
      )

      if (is.null(candidate)) {
        next
      }

      valid_manifest <- tryCatch(
        {
          .manifest_validate(
            candidate
          )

          TRUE
        },
        error = function(e) {
          FALSE
        }
      )

      if (
        isTRUE(valid_manifest) &&
        identical(
          candidate$manifest_type %||% "",
          "amR_dataset"
        )
      ) {
        manifest <- candidate
        manifest_path <- candidate_path
        break
      }
    }

    if (is.null(manifest)) {
      stop(
        "Complete removal refused: no valid amRdata dataset manifest ",
        "was found.",
        call. = FALSE
      )
    }

    dataset_id <- as.character(
      manifest$dataset_id %||% ""
    )

    if (!nzchar(dataset_id)) {
      stop(
        "Complete removal refused: manifest does not contain a dataset ID.",
        call. = FALSE
      )
    }

    user_bacs <- unlist(
      manifest$dataset$selection$user_bacs %||% character(),
      use.names = FALSE
    )

    if (!length(user_bacs)) {
      stop(
        "Complete removal refused: manifest does not identify the ",
        "dataset selection.",
        call. = FALSE
      )
    }

    expected_dataset_name <- paste(
      user_bacs,
      collapse = "__"
    ) |>
      stringr::str_replace_all(
        "\\s+",
        "_"
      ) |>
      stringr::str_replace_all(
        "[^A-Za-z0-9._-]",
        ""
      )

    if (!identical(
      basename(dataset_path),
      expected_dataset_name
    )) {
      stop(
        "Complete removal refused: dataset directory name does not match ",
        "the dataset recorded in the manifest.\n",
        "Expected: ",
        expected_dataset_name,
        "\nFound: ",
        basename(dataset_path),
        call. = FALSE
      )
    }

    expected_dataset_id <- .generateDBname(
      user_bacs
    )

    if (!identical(
      dataset_id,
      expected_dataset_id
    )) {
      stop(
        "Complete removal refused: manifest dataset ID does not match ",
        "the dataset selection.",
        call. = FALSE
      )
    }

    # Any DuckDB files that exist must use the canonical names
    work_dir <- file.path(
      dataset_path,
      "work"
    )

    work_duckdbs <- if (dir.exists(work_dir)) {
      list.files(
        work_dir,
        pattern = "\\.duckdb$",
        full.names = FALSE
      )
    } else {
      character()
    }

    orb_duckdbs <- list.files(
      orb_dir,
      pattern = "\\.duckdb$",
      full.names = FALSE
    )

    unexpected_work_duckdbs <- setdiff(
      work_duckdbs,
      paste0(
        dataset_id,
        ".duckdb"
      )
    )

    unexpected_orb_duckdbs <- setdiff(
      orb_duckdbs,
      paste0(
        dataset_id,
        "_parquet.duckdb"
      )
    )

    unexpected_duckdbs <- c(
      unexpected_work_duckdbs,
      unexpected_orb_duckdbs
    )

    if (length(unexpected_duckdbs)) {
      stop(
        "Complete removal refused: unexpected DuckDB file(s) were found:\n",
        paste(
          unexpected_duckdbs,
          collapse = "\n"
        ),
        call. = FALSE
      )
    }

    dataset_bytes <- path_bytes(
      dataset_path
    )

    # Requiring interactive confirmation protects against accidental execution
    # from command history and prevents unattended scripts from recursively
    # deleting complete datasets
    if (!interactive()) {
      stop(
        "Complete removal requires an interactive R session so deletion ",
        "can be explicitly confirmed.",
        call. = FALSE
      )
    }

    cat(
      "\n",
      "PERMANENT DATASET REMOVAL\n",
      "=========================\n",
      "This will permanently delete the entire amRdata dataset,\n",
      "including its ORB, Parquet files, manifests, exports, all\n",
      "genomes, and intermediate build files.\n\n",
      "Dataset:\n",
      dataset_path,
      "\n\n",
      "Disk space to reclaim: ",
      format_bytes(dataset_bytes),
      "\n\n",
      sep = ""
    )

    confirmation <- readline(
      paste0(
        "Are you sure? Type '",
        basename(dataset_path),
        "' to confirm: "
      )
    )

    if (!identical(
      trimws(confirmation),
      basename(dataset_path)
    )) {
      if (isTRUE(verbose)) {
        message(
          "Complete removal cancelled. Nothing was deleted."
        )
      }

      return(
        invisible(FALSE)
      )
    }

    unlink(
      dataset_path,
      recursive = TRUE,
      force = TRUE
    )

    if (dir.exists(dataset_path)) {
      stop(
        "Complete removal was requested, but the dataset directory ",
        "could not be fully removed:\n",
        dataset_path,
        call. = FALSE
      )
    }

    if (isTRUE(verbose)) {
      message(
        "Completely removed amRdata dataset: ",
        dataset_path,
        "\nReclaimed ",
        format_bytes(dataset_bytes),
        " of disk space."
      )
    }

    return(
      invisible(dataset_path)
    )
  }

  # Basic default cleanup, still retain ORB and exports
  parquet_duckdb <- list.files(
    orb_dir,
    pattern = "_parquet\\.duckdb$",
    full.names = TRUE
  )

  if (length(parquet_duckdb) != 1L) {
    stop(
      "Expected exactly one '*_parquet.duckdb' file in: ",
      orb_dir,
      "\nFound: ",
      length(parquet_duckdb),
      call. = FALSE
    )
  }

  parquet_duckdb <- normalizePath(
    parquet_duckdb[[1]],
    mustWork = TRUE
  )

  paths <- .amr_paths_from_dataset_db(
    parquet_duckdb
  )

  manifest_path <- .manifest_find_latest(
    parquet_duckdb,
    require_success = TRUE
  )

  if (is.null(manifest_path)) {
    stop(
      "No successful dataset manifest was found. ",
      "Local build files will not be removed."
    )
  }

  manifest <- jsonlite::read_json(
    manifest_path,
    simplifyVector = FALSE
  )

  .manifest_validate(
    manifest
  )

  artifact <- manifest$artifacts$amRml_input %||% NULL

  if (
    is.null(artifact) ||
    !identical(
      artifact$status,
      "ready"
    )
  ) {
    stop(
      "The dataset is not marked as a ready amRml input. ",
      "Local build files will not be removed."
    )
  }

  producer_run_id <- artifact$producer_run_id %||% NULL

  if (
    is.null(producer_run_id) ||
    !nzchar(producer_run_id)
  ) {
    stop(
      "The ready ORB does not identify its producing run. ",
      "Local build files will not be removed."
    )
  }

  run_ids <- purrr::map_chr(
    manifest$runs %||% list(),
    ~ .x$run_id %||% ""
  )

  run_index <- which(
    run_ids == producer_run_id
  )

  if (length(run_index) != 1L) {
    stop(
      "Could not uniquely identify the run that produced the ORB. ",
      "Local build files will not be removed."
    )
  }

  producer_run <- manifest$runs[[run_index]]

  if (!identical(
    producer_run$status,
    "success"
  )) {
    stop(
      "The run that produced the ORB is not marked successful. ",
      "Local build files will not be removed."
    )
  }

  hmmer_stage <- NULL

  for (run in rev(
    manifest$runs[
      seq_len(run_index)
    ]
  )) {
    matches <- purrr::keep(
      run$stages %||% list(),
      ~ identical(
        .x$name,
        "hmmer"
      ) &&
        identical(
          .x$status,
          "success"
        )
    )

    if (length(matches)) {
      hmmer_stage <- matches[[length(matches)]]
      break
    }
  }

  if (is.null(hmmer_stage)) {
    stop(
      "No successful HMMER stage was found at or before the ",
      "ORB-producing run. Local build files will not be removed."
    )
  }

  hmmer_databases <- unlist(
    hmmer_stage$parameters$databases %||% character(),
    use.names = FALSE
  )

  hmmer_databases <- unique(
    as.character(
      hmmer_databases
    )
  )

  required_orb_files <- c(
    paths$parquet_duckdb,

    file.path(
      paths$orb,
      c(
        "metadata.parquet",
        "metadata_qc.parquet",
        "metadata_qc_rejections.parquet",
        "selected_genomes.parquet",
        "amr_phenotype.parquet",
        "genome_data.parquet",
        "original_metadata.parquet",

        "gene_count.parquet",
        "gene_names.parquet",
        "gene_seqs.parquet",
        "genome_gene_protein.parquet",
        "struct.parquet",

        "protein_count.parquet",
        "protein_names.parquet",
        "protein_seqs.parquet",
        "protein_members.parquet",

        "dyad_feature.parquet"
      )
    ),

    file.path(
      paths$orb,
      paste0(
        "protein_",
        hmmer_databases,
        ".parquet"
      )
    ),

    file.path(
      paths$orb,
      paste0(
        "protein_",
        hmmer_databases,
        "_count.parquet"
      )
    )
  )

  missing_orb_files <- required_orb_files[
    !file.exists(required_orb_files)
  ]

  if (length(missing_orb_files)) {
    stop(
      "ORB validation failed. Local build files will not be removed.\n",
      "Missing expected ORB file(s):\n",
      paste(
        missing_orb_files,
        collapse = "\n"
      )
    )
  }

  remove_paths <- c(
    paths$genomes,
    paths$panaroo,
    paths$cdhit,
    paths$hmmer,
    paths$work
  )

  remove_paths <- remove_paths[
    dir.exists(remove_paths)
  ]

  reclaimed_bytes <- sum(
    vapply(
      remove_paths,
      path_bytes,
      numeric(1)
    )
  )

  retained_bytes <- sum(
  vapply(
    c(
      paths$orb,
      paths$exports
    ),
    path_bytes,
    numeric(1)
  )
)

  if (!length(remove_paths)) {
  if (isTRUE(verbose)) {
    message("No local build files found to remove.")
  }

  return(
    invisible(character())
  )
}

if (!interactive()) {
  stop(
    "Local file cleanup requires an interactive R session so deletion ",
    "can be explicitly confirmed.",
    call. = FALSE
  )
}

cat(
  "\n",
  "LOCAL FILE CLEANUP\n",
  "==================\n",
  "This will remove disposable build files for:\n",
  dataset_path,
  "\n\n",
  "Disk space to reclaim: ",
  format_bytes(reclaimed_bytes),
  "\n",
  "ORB and exports will be retained.\n\n",
  sep = ""
)

confirmation <- readline(
  "Continue? [y/N]: "
)

if (!tolower(trimws(confirmation)) %in% c("y", "yes")) {
  if (isTRUE(verbose)) {
    message(
      "Local file cleanup cancelled. Nothing was deleted."
    )
  }

  return(
    invisible(FALSE)
  )
}

if (length(remove_paths)) {
  purrr::walk(
    remove_paths,
    ~ unlink(
      .x,
      recursive = TRUE,
      force = TRUE
    )
  )
}

  remaining <- remove_paths[
    dir.exists(remove_paths)
  ]

  if (length(remaining)) {
    stop(
      "Some local build directories could not be removed:\n",
      paste(
        remaining,
        collapse = "\n"
      )
    )
  }

  if (isTRUE(verbose)) {
    message(
      "Removed local build files.\n",
      "Reclaimed ",
      format_bytes(reclaimed_bytes),
      " of disk space.\n",
      "ORB and exports retained (",
      format_bytes(retained_bytes),
      ")."
    )
  }

  invisible(
    remove_paths
  )
}
