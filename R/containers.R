# Container execution shared by the data curation and processing workflows.
# Configuration is passed explicitly to future workers, never stored in options.
.findContainerExecutables <- function(x) Sys.which(x)
.executeContainerProcess <- function(...) system2(...)
.createContainerConfig <- function(backend = "auto", images = NULL, cache = NULL) {
  if (inherits(backend, "amr_container_config")) return(backend)
  backend <- match.arg(backend, c("auto", "docker", "apptainer", "singularity"))
  if (!is.null(images) && (!is.character(images) || is.null(names(images)) ||
      anyNA(images) || anyNA(names(images)) || any(!nzchar(images)) || any(!nzchar(names(images))) ||
      anyDuplicated(names(images)))) {
    stop("`container_images` must be a named character vector of image references or SIF paths.", call. = FALSE)
  }
  if (is.null(cache)) cache <- file.path(tools::R_user_dir("amRdata", "cache"), "containers")
  if (!is.character(cache) || length(cache) != 1L || is.na(cache) || !nzchar(cache)) {
    stop("`container_cache` must be a directory path.", call. = FALSE)
  }
  structure(list(backend = backend, images = images, cache = path.expand(cache),
                 executable = NULL, version = NULL, resolved = list(), log_dir = NULL),
            class = "amr_container_config")
}

.selectContainerBackend <- function(requested, available, interactive_session = interactive(),
                                    choose = utils::menu) {
  if (requested != "auto") {
    if (!requested %in% available) stop("Container runtime '", requested, "' is not on PATH.", call. = FALSE)
    return(requested)
  }
  alternatives <- intersect(c("apptainer", "singularity"), available)
  choices <- c(intersect("docker", available), head(alternatives, 1L))
  if (!length(choices)) stop("No container runtime found. Install Docker, Apptainer, or Singularity, or load the appropriate HPC module.", call. = FALSE)
  if (length(choices) == 1L) return(choices)
  if (!interactive_session) stop("Multiple container backends are available. Set `container_backend` explicitly for batch jobs.", call. = FALSE)
  selected <- choose(choices, title = "Choose an amRdata container backend:")
  if (length(selected) != 1L || is.na(selected) || selected < 1L || selected > length(choices)) {
    stop("Container backend selection cancelled.", call. = FALSE)
  }
  choices[[selected]]
}

.resolveContainerRuntime <- function(config) {
  if (!is.null(config$executable)) return(config)
  paths <- .findContainerExecutables(c("docker", "apptainer", "singularity"))
  config$backend <- .selectContainerBackend(config$backend, names(paths)[nzchar(paths)])
  config$executable <- unname(paths[[config$backend]])
  if (config$backend == "docker") {
    status <- suppressWarnings(.executeContainerProcess(config$executable, "info", stdout = FALSE, stderr = FALSE))
    if (!identical(status, 0L)) stop("Docker is on PATH but its daemon cannot be reached. Start Docker or select another `container_backend`.", call. = FALSE)
  }
  version <- suppressWarnings(.executeContainerProcess(config$executable, "--version", stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(version, "status"))) stop("Cannot execute container runtime: ", config$backend, call. = FALSE)
  config$version <- paste(version, collapse = "\n")
  config
}

.resolveContainerImage <- function(config, image) {
  keys <- c("staphb/panaroo:1.7.0" = "panaroo", "weizhongli1987/cdhit:4.8.1" = "cdhit",
            "staphb/hmmer" = "hmmer", "danylmb/bvbrc:5.3" = "bvbrc")
  key <- unname(keys[image])
  if (image %in% names(config$images)) return(unname(config$images[[image]]))
  if (!is.na(key) && key %in% names(config$images)) return(unname(config$images[[key]]))
  image
}

.prepareContainerImages <- function(config, images) {
  config <- .resolveContainerRuntime(config)
  for (image in unique(images)) {
    if (!is.null(config$resolved[[image]])) next
    source <- .resolveContainerImage(config, image)
    local_sif <- grepl("\\.sif$", source, ignore.case = TRUE)
    if (config$backend == "docker") {
      if (local_sif) stop("Local SIF images require Apptainer or Singularity: ", source, call. = FALSE)
      reference <- sub("^docker://", "", source)
      identity <- suppressWarnings(.executeContainerProcess(config$executable,
        c("image", "inspect", "--format", shQuote("{{.Id}}"), shQuote(reference)),
        stdout = TRUE, stderr = FALSE))
      if (!is.null(attr(identity, "status")) || !length(identity)) {
        platform <- if (image %in% c("staphb/panaroo:1.7.0", "weizhongli1987/cdhit:4.8.1")) c("--platform", "linux/amd64") else character()
        status <- .executeContainerProcess(config$executable, c("pull", platform, shQuote(reference)))
        if (!identical(status, 0L)) stop("Failed to pull Docker image: ", reference, call. = FALSE)
        identity <- .executeContainerProcess(config$executable,
          c("image", "inspect", "--format", shQuote("{{.Id}}"), shQuote(reference)),
          stdout = TRUE, stderr = TRUE)
      }
      if (!is.null(attr(identity, "status")) || !length(identity)) stop("Cannot identify Docker image: ", reference, call. = FALSE)
      record <- list(source = source, reference = reference, identity = identity[[1]])
    } else {
      if (local_sif) {
        sif <- normalizePath(path.expand(source), mustWork = TRUE)
      } else {
        dir.create(config$cache, recursive = TRUE, showWarnings = FALSE)
        cache <- normalizePath(config$cache, mustWork = TRUE)
        uri <- if (startsWith(source, "docker://")) source else paste0("docker://", source)
        # Hash the full reference to avoid collisions from sanitized filenames.
        key_file <- tempfile()
        writeLines(uri, key_file)
        key <- unname(tools::md5sum(key_file))
        unlink(key_file)
        sif <- file.path(cache, paste0(substr(gsub("[^A-Za-z0-9_.-]", "_", source), 1L, 80L), "-", key, ".sif"))
        if (!file.exists(sif)) {
          temporary <- tempfile(pattern = ".pull-", tmpdir = cache, fileext = ".sif")
          tryCatch({
            status <- .executeContainerProcess(config$executable, c("pull", shQuote(temporary), shQuote(uri)))
            if (!identical(status, 0L) || !file.exists(temporary)) stop("Failed to pull SIF image: ", uri, call. = FALSE)
            if (!file.rename(temporary, sif)) stop("Cannot publish SIF image: ", sif, call. = FALSE)
          }, finally = unlink(temporary))
        }
      }
      if (file.info(sif)$size <= 0) stop("Empty SIF image: ", sif, call. = FALSE)
      record <- list(source = source, reference = sif, identity = unname(tools::md5sum(sif)))
    }
    config$resolved[[image]] <- record
  }
  config
}

.buildContainerCommand <- function(config, image, command, args = character(),
                              binds = character(), workdir = NULL, platform = NULL) {
  record <- config$resolved[[image]]
  if (is.null(record)) stop("Container image was not prepared: ", image, call. = FALSE)
  if (config$backend == "docker") {
    prefix <- c("run", "--rm")
    if (!is.null(platform)) prefix <- c(prefix, "--platform", platform)
    for (bind in binds) prefix <- c(prefix, "-v", bind)
    if (!is.null(workdir)) prefix <- c(prefix, "-w", workdir)
    # Run the inspected image ID so a changing tag cannot alter a prepared job.
    reference <- record$identity
  } else {
    prefix <- c("exec", "--cleanenv")
    for (bind in binds) prefix <- c(prefix, "--bind", bind)
    if (!is.null(workdir)) prefix <- c(prefix, "--pwd", workdir)
    reference <- record$reference
  }
  list(executable = config$executable, args = c(prefix, reference, command, as.character(args)))
}

.runContainer <- function(config, image, command, args = character(), binds = character(),
                          workdir = NULL, platform = NULL, stdout = TRUE, stderr = TRUE) {
  config <- .prepareContainerImages(config, image)
  invocation <- .buildContainerCommand(config, image, command, args, binds, workdir, platform)
  started <- as.character(Sys.time())
  result <- tryCatch(
    .executeContainerProcess(invocation$executable, shQuote(invocation$args), stdout = stdout, stderr = stderr),
    error = function(e) e
  )
  status <- if (inherits(result, "error")) NA_integer_ else if (is.numeric(result)) {
    as.integer(result)
  } else if (is.null(attr(result, "status"))) 0L else as.integer(attr(result, "status"))
  if (!is.null(config$log_dir)) {
    dir.create(config$log_dir, recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(list(
      started_at = started, finished_at = as.character(Sys.time()),
      backend = config$backend, runtime_version = config$version,
      image = config$resolved[[image]], executable = invocation$executable,
      args = invocation$args, exit_status = status,
      error = if (inherits(result, "error")) conditionMessage(result) else NULL
    ), tempfile(pattern = paste0("command-", Sys.getpid(), "-"), tmpdir = config$log_dir, fileext = ".json"),
    auto_unbox = TRUE, pretty = TRUE, null = "null")
  }
  if (inherits(result, "error")) stop(result)
  result
}

# Adapter for existing Docker-style argument vectors. Tool arguments and mounts
# remain at their existing call sites; unsupported runtime flags fail explicitly.
.runContainerCommand <- function(config, args, stdout = TRUE, stderr = TRUE) {
  stopifnot(length(args) > 1L, args[[1]] == "run")
  i <- 2L
  binds <- character()
  workdir <- platform <- NULL
  while (i <= length(args) && startsWith(args[[i]], "-")) {
    flag <- args[[i]]
    if (flag == "--rm") { i <- i + 1L; next }
    if (!flag %in% c("-v", "--volume", "-w", "--workdir", "--platform") || i == length(args)) {
      stop("Unsupported container option: ", flag, call. = FALSE)
    }
    value <- args[[i + 1L]]
    if (flag %in% c("-v", "--volume")) binds <- c(binds, value)
    if (flag %in% c("-w", "--workdir")) workdir <- value
    if (flag == "--platform") platform <- value
    i <- i + 2L
  }
  if (i + 1L > length(args)) stop("Container image and command are required.", call. = FALSE)
  .runContainer(config, args[[i]], args[[i + 1L]], args[-seq_len(i + 1L)],
                binds, workdir, platform, stdout, stderr)
}

.getContainerProvenance <- function(config, image = NULL) {
  if (is.null(config) || is.null(config$executable)) return(NULL)
  records <- config$resolved
  if (!is.null(image)) records <- records[intersect(image, names(records))]
  list(backend = config$backend, runtime_version = config$version,
       images = records, execution_log_dir = config$log_dir)
}

.isContainerStageReusable <- function(stage, config) {
  if (is.null(config) || is.null(stage$tool$docker_image)) return(TRUE)
  image <- stage$tool$docker_image
  current <- config$resolved[[image]]
  if (is.null(current)) return(FALSE)
  previous <- stage$tool$container
  if (is.null(previous)) {
    # Preserve legacy Docker resumes only when the original reference is used.
    return(config$backend == "docker" && identical(current$source, image))
  }
  old <- previous$images[[image]]
  if (is.null(old)) return(FALSE)
  if (identical(previous$backend, config$backend)) return(identical(old$identity, current$identity))
  # SIF checksums and Docker IDs differ. A shared immutable OCI digest can prove
  # equivalence across runtimes; a mutable tag alone cannot.
  source <- sub("^docker://", "", current$source)
  grepl("@sha256:[a-fA-F0-9]{64}$", source) &&
    identical(source, sub("^docker://", "", old$source))
}
