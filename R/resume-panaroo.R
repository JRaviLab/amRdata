#' Recover a failed Panaroo merge without rerunning its batches
#'
#' Appends a recovery run to the existing dataset manifest, merges explicitly
#' selected completed batches, and imports the merged output into the working
#' DuckDB. Other processing stages are not run. After success, call
#' [runDataProcessing()] with `resume = TRUE` and the original scientific
#' parameters to continue processing. Use the same container images for both.
#'
#' Existing `merge_output` is renamed and preserved before retrying. The original
#' failed run remains unchanged. Batches must come from the same failed Panaroo
#' attempt; this helper validates required files, but cannot reconstruct or prove
#' batch provenance that older manifests did not record. On a different machine,
#' copy the working DuckDB, manifest, and selected batch directories together,
#' preserving the dataset's `work/`, `orb/`, and `panaroo/` layout. Resuming later
#' stages additionally requires the genome files and valid paths in the working
#' database; transferring only the merge inputs is sufficient for this helper.
#'
#' @param duckdb_path Working DuckDB under the dataset's `work/` directory.
#' @param batch_dirs Character vector of at least two completed `panaroo_out_*`
#'   directories directly under the dataset's `panaroo/` directory. Select the
#'   exact batches from the failed attempt, excluding older batches.
#' @param threads Threads for the merge. Default is six.
#' @param manifest_path Dataset manifest to append to. If `NULL`, use the latest
#'   manifest for `duckdb_path`, including failed attempts.
#' @inheritParams prepareGenomes
#' @return Invisibly, the normalized working DuckDB path.
#' @export
resumePanarooMerge <- function(duckdb_path, batch_dirs, threads = 6L,
                               manifest_path = NULL, container_backend = "auto",
                               container_images = NULL, container_cache = NULL) {
  duckdb_path <- normalizePath(duckdb_path, mustWork = TRUE)
  paths <- .amr_paths_from_duckdb(duckdb_path)
  threads <- .resolve_workers(requested = threads)
  if (!is.character(batch_dirs) || length(batch_dirs) < 2L || anyNA(batch_dirs)) {
    stop("Supply at least two completed Panaroo batch directories.", call. = FALSE)
  }
  batch_dirs <- normalizePath(batch_dirs, mustWork = TRUE)
  panaroo <- normalizePath(paths$panaroo, mustWork = TRUE)
  if (anyDuplicated(batch_dirs) || any(dirname(batch_dirs) != panaroo) ||
      any(!grepl("^panaroo_out_", basename(batch_dirs)))) {
    stop("Batches must be unique panaroo_out_* directories directly under the dataset's panaroo/ directory.", call. = FALSE)
  }
  required <- c("final_graph.gml", "gene_data.csv")
  if (!all(vapply(batch_dirs, function(path) {
    files <- file.path(path, required)
    all(file.exists(files)) && all(file.info(files)$size > 0)
  }, logical(1)))) stop("Each batch needs nonempty final_graph.gml and gene_data.csv files.", call. = FALSE)

  if (is.null(manifest_path)) manifest_path <- .manifest_find_latest(duckdb_path, require_success = FALSE)
  if (is.null(manifest_path)) stop("No dataset manifest found.", call. = FALSE)
  manifest_path <- normalizePath(manifest_path, mustWork = TRUE)
  if (dirname(manifest_path) != normalizePath(paths$orb, mustWork = TRUE)) {
    stop("The recovery manifest must belong to this dataset's orb/ directory.", call. = FALSE)
  }
  previous_manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
  .manifest_validate(previous_manifest)
  previous <- tail(previous_manifest$runs, 1L)[[1]]
  prior <- .manifest_prior_stage(previous, "panaroo")
  if (is.null(prior) || !prior$status %in% c("failed", "running") ||
      !isTRUE(prior$parameters$split_jobs)) {
    stop("The latest manifest run must contain a failed/interrupted split Panaroo stage.", call. = FALSE)
  }
  parameter_names <- c("core_threshold", "len_dif_percent", "cluster_threshold", "family_seq_identity")
  if (!all(parameter_names %in% names(prior$parameters))) {
    stop("The failed stage lacks the merge parameters required for safe recovery.", call. = FALSE)
  }
  config <- .prepareContainerImages(
    .createContainerConfig(container_backend, container_images, container_cache),
    "staphb/panaroo:1.7.0"
  )
  # These files are recovery evidence; neither prior provenance nor partial
  # merge outputs are overwritten. tempfile provides collision-free names.
  backup <- tempfile(pattern = paste0(basename(manifest_path), ".before-merge-"),
                     tmpdir = dirname(manifest_path), fileext = ".bak")
  if (!file.copy(manifest_path, backup)) stop("Cannot back up the manifest.", call. = FALSE)
  merged <- file.path(panaroo, "merge_output")
  if (file.exists(merged)) {
    archived <- tempfile(pattern = "merge_output_failed_", tmpdir = panaroo)
    if (!file.rename(merged, archived)) stop("Cannot preserve the previous merge_output.", call. = FALSE)
  }
  manifest <- .manifest_resume(manifest_path)
  config$log_dir <- file.path(paths$orb, "container-logs",
    manifest$manifest$runs[[manifest$run_index]]$run_id)
  log_path <- paste0(config$log_dir, "-panaroo.log")
  parameters <- prior$parameters
  parameters$threads <- threads
  parameters$threads_requested <- parameters$threads_used <- NULL
  success <- FALSE
  failure <- "Panaroo merge recovery exited before completion."
  on.exit({
    .manifest_finish(manifest, status = if (success) "success" else "failed",
                     error = if (success) NULL else failure)
  }, add = TRUE)
  tryCatch({
    manifest <- .manifest_event(manifest,
      message = "Merge-only recovery; completed Panaroo batches reused. Downstream processing remains pending.",
      details = list(source_run_id = previous$run_id, batches = batch_dirs,
                     manifest_backup = backup, log_path = log_path))
    tool <- list(name = "Panaroo", version = "1.7.0", docker_image = "staphb/panaroo:1.7.0")
    inputs <- c(duckdb_path, file.path(batch_dirs, "final_graph.gml"),
                file.path(batch_dirs, "gene_data.csv"))
    manifest <- .manifest_stage(manifest, "panaroo", status = "running",
      parameters = parameters, inputs = inputs, tool = tool, container_config = config)
    do.call(.mergePanaroo, c(list(input_path = panaroo),
      parameters[parameter_names], list(threads = threads, log_path = log_path,
                                       container_config = config, batch_dirs = batch_dirs)))
    .panaroo2duckdb(merged, duckdb_path)
    outputs <- c(merged, duckdb_path, file.path(panaroo, "panaroo_invalid_cds_audit.csv"))
    manifest <- .manifest_stage(manifest, "panaroo", status = "success",
      parameters = parameters, inputs = inputs, outputs = outputs[file.exists(outputs)],
      tool = tool, container_config = config,
      message = "Recovered merge and DuckDB import; batch generation was not repeated.")
    success <- TRUE
  }, error = function(e) {
    failure <<- conditionMessage(e)
    stop(e)
  })
  invisible(duckdb_path)
}
