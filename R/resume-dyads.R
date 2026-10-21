#' Rebuild the dyads in an existing ORB
#'
#' Rebuilds the gene-cluster/protein-cluster pairs and their DuckDB views without
#' rerunning genome processing, annotations, metadata cleaning, or Parquet export.
#' Checks that the other files still match the last successful export, then
#' records the rebuilt ORB in a new run in the existing provenance manifest. This
#' is likely a development-only function.
#'
#' @param orb_dir Directory containing the Parquet-backed ORB and its manifest.
#' @param threads Maximum DuckDB threads for the dyad rebuild (default 4).
#' @return Invisibly, the paths to the rebuilt dyads and manifest.
#' @export
resumeDyads <- function(orb_dir, threads = 4L) {
  orb_dir <- normalizePath(orb_dir, winslash = "/", mustWork = TRUE)
  if (basename(orb_dir) != "orb") {
    stop("`orb_dir` must be an existing dataset 'orb/' directory.", call. = FALSE)
  }

  databases <- list.files(orb_dir, pattern = "_parquet[.]duckdb$", full.names = TRUE)
  if (length(databases) != 1L) {
    stop("Expected exactly one *_parquet.duckdb in the ORB.", call. = FALSE)
  }
  database <- normalizePath(databases[[1L]], winslash = "/", mustWork = TRUE)
  if (file.exists(paste0(database, ".wal"))) {
    stop("Close other ORB DuckDB connections before rebuilding dyads.", call. = FALSE)
  }
  manifest_path <- .manifest_find_latest(database)
  if (is.null(manifest_path)) {
    stop("No successful provenance manifest was found for this ORB.", call. = FALSE)
  }

  previous <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
  .manifest_validate(previous)

  # The final export lists every file amRml expects. Use most recent
  # successful snapshot, not old paths
  exported <- NULL
  for (run in rev(previous$runs %||% list())) {
    if (!identical(run$status, "success")) next
    stage <- .manifest_prior_stage(run, "clean_metadata_and_export")
    if (!is.null(stage) && identical(stage$status, "success")) {
      exported <- stage
      break
    }
  }
  if (is.null(exported) || !length(exported$outputs)) {
    stop("No successful ORB export was recorded in the manifest.", call. = FALSE)
  }

  names_in_manifest <- vapply(exported$outputs, function(x) {
    basename(chartr("\\", "/", x$path))
  }, character(1L))
  if (anyDuplicated(names_in_manifest)) {
    stop("Duplicate ORB output filenames in the manifest.", call. = FALSE)
  }
  required <- c("dyads.parquet", basename(database),
                "genome_gene_protein.parquet", "protein_members.parquet",
                "metadata.parquet", "amr_phenotype.parquet")
  if (!all(required %in% names_in_manifest)) {
    stop("The manifest does not describe a complete ORB export.", call. = FALSE)
  }

  orb_files <- file.path(orb_dir, names_in_manifest)
  if (!all(file.exists(orb_files))) {
    stop("Missing ORB file(s): ", paste(names_in_manifest[!file.exists(orb_files)],
                                        collapse = ", "), call. = FALSE)
  }

  # dyads.parquet and the DuckDB views are the  only files this operation may
  # change
  changed <- !(names_in_manifest %in% c("dyads.parquet", basename(database)))
  for (i in which(changed)) {
    expected <- as.numeric(exported$outputs[[i]]$size_bytes)
    actual <- as.numeric(file.info(orb_files[[i]])$size)
    if (length(expected) != 1L || is.na(expected) || actual != expected) {
      stop("ORB input changed since the previous export: ", names_in_manifest[[i]],
           ". Rebuild the upstream data before updating the manifest.", call. = FALSE)
    }
    if (!is.null(exported$outputs[[i]]$md5) &&
        !identical(as.character(exported$outputs[[i]]$md5),
                   unname(tools::md5sum(orb_files[[i]])))) {
      stop("ORB input checksum changed: ", names_in_manifest[[i]], call. = FALSE)
    }
  }

  # Restore the original dyads and DuckDB if any rebuild step fails. keep the
  # backup if a restore fails, so the original files can be recovered
  backup_dir <- tempfile("dyad_rebuild_backup_", tmpdir = dirname(orb_dir))
  if (!dir.create(backup_dir)) {
    stop("Could not create an ORB backup directory.", call. = FALSE)
  }
  preserved <- c(database, file.path(orb_dir, "dyads.parquet"),
                 file.path(orb_dir, "dyad_feature.parquet"))
  preserved <- preserved[file.exists(preserved)]
  if (!all(file.copy(preserved, backup_dir, overwrite = FALSE))) {
    unlink(backup_dir, recursive = TRUE, force = TRUE)
    stop("Could not back up the existing dyad files and ORB DuckDB.", call. = FALSE)
  }

  manifest <- NULL
  run_id <- NULL
  success <- FALSE
  on.exit({
    if (!success) {
      # Do this before deleting backups or marking the failed manifest run.
      unlink(paste0(database, ".wal"), force = TRUE)
      restored <- vapply(preserved, function(file) {
        file.copy(file.path(backup_dir, basename(file)), file, overwrite = TRUE)
      }, logical(1L))
      if (!is.null(manifest)) {
        try({
          manifest <- .manifest_artifact(manifest, "amRml_input", "failed",
                                         details = list(producer = "amRdata",
                                                        producer_run_id = run_id))
          .manifest_finish(manifest, status = "failed",
                           error = "Dyad rebuild or ORB validation failed; original files restored if possible.")
        }, silent = TRUE)
      }
      if (!all(restored)) {
        warning("Could not restore every ORB file; backups retained at: ",
                backup_dir, call. = FALSE)
      } else {
        unlink(backup_dir, recursive = TRUE, force = TRUE)
      }
    } else {
      unlink(backup_dir, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)

  manifest <- .manifest_resume(
    manifest_path = manifest_path,
    base_dir = dirname(dirname(dirname(orb_dir))),
    hash_files = any(vapply(exported$outputs,
                            function(x) !is.null(x$md5), logical(1L)))
  )
  run_id <- manifest$manifest$runs[[manifest$run_index]]$run_id

  manifest <- .manifest_artifact(manifest, "amRml_input", "building",
                                 details = list(producer = "amRdata",
                                                producer_run_id = run_id))
  manifest <- .manifest_stage(manifest, name = "rebuild_dyads", status = "running",
                               inputs = file.path(orb_dir,
                                                 c("genome_gene_protein.parquet",
                                                   "protein_members.parquet")),
                               parameters = list(threads = threads))

  dyads <- buildDyadFeatureMap(
    duckdb_path = database,
    output_path = orb_dir,
    threads = threads
  )

  # A valid ORB has cluster-pair edges, member and genome provenance, and all
  # the virtual annotation views supported by this dataset's HMMER run
  local({
    con <- .amr_connect_dataset_db(database, read_only = TRUE)
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    relations <- c("dyads", "dyad_feature", "dyad_protein_cluster",
                   "member_dyads", "genome_dyad", "genome_cluster_dyad")
    for (run in rev(previous$runs %||% list())) {
      stage <- .manifest_prior_stage(run, "hmmer")
      if (!is.null(stage) && identical(stage$status, "success")) {
        hmmer_databases <- as.character(unlist(stage$parameters$databases, use.names = FALSE))
        specs <- .dyadHmmerSpecs(hmmer_databases)
        relations <- c(relations, unlist(lapply(specs, function(spec) {
          c(spec$annotation_view, spec$dyad_view)
        }), use.names = FALSE))
        break
      }
    }
    if (!all(relations %in% DBI::dbListTables(con))) {
      stop("Rebuilt ORB is missing dyad relations.", call. = FALSE)
    }
    for (relation in relations) {
      DBI::dbGetQuery(con, paste("SELECT * FROM",
        DBI::dbQuoteIdentifier(con, relation), "LIMIT 0"))
    }
    invalid <- DBI::dbGetQuery(con, "
      SELECT count(*) AS n FROM dyads d
      LEFT JOIN (SELECT DISTINCT cluster FROM protein_members) p
        ON d.protein = p.cluster
      WHERE p.cluster IS NULL OR d.gene IS NULL OR d.protein IS NULL
    ")$n[[1L]]
    n_pairs <- DBI::dbGetQuery(con, "SELECT count(*) AS n FROM dyads")$n[[1L]]
    repeated <- DBI::dbGetQuery(con, "
      SELECT count(*) AS n FROM (
        SELECT protein, gene FROM dyads
        GROUP BY protein, gene HAVING count(*) > 1
      )
    ")$n[[1L]]
    if (invalid > 0L || n_pairs < 1L || repeated > 0L) {
      stop("The dyad graph is empty, duplicated, or contains non-cluster protein IDs.", call. = FALSE)
    }
  })

  manifest <- .manifest_stage(manifest, name = "rebuild_dyads", status = "success",
                               inputs = file.path(orb_dir,
                                                 c("genome_gene_protein.parquet",
                                                   "protein_members.parquet")),
                               outputs = c(dyads, database),
                               parameters = list(threads = threads))

  # Record a full, fresh export snapshot in the stage understood by amRml
  parameters <- exported$parameters
  parameters$resumed_from <- "dyads"
  manifest <- .manifest_stage(manifest,
    name = "clean_metadata_and_export", status = "success",
    parameters = parameters,
    inputs = file.path(orb_dir, c("genome_gene_protein.parquet", "protein_members.parquet")),
    outputs = orb_files,
    message = "Reused unchanged Parquet exports; rebuilt dyads and finalized the ORB."
  )

  manifest <- .manifest_artifact(manifest, "amRml_input", "ready", details = list(
    producer = "amRdata", producer_run_id = run_id,
    path_mode = "relative_to_manifest", directory = ".",
    parquet_duckdb = basename(database),
    metadata_parquet = "metadata.parquet",
    metadata_qc_parquet = "metadata_qc.parquet",
    metadata_qc_rejections_parquet = "metadata_qc_rejections.parquet",
    selected_genomes_parquet = "selected_genomes.parquet",
    dyads_parquet = "dyads.parquet", dyad_representation = "virtual",
    dyad_feature_relation = "dyad_feature",
    dyad_hmmer_mapping = "protein_members.member_to_cluster"
  ))
  manifest <- .manifest_finish(manifest, status = "success")
  success <- TRUE
  message("Rebuilt dyads and updated ORB manifest: ", manifest_path)
  invisible(list(dyads_parquet = dyads, manifest_path = manifest_path))
}
