.createContainerFixture <- function(backend = "apptainer", cache = tempfile()) {
  config <- .createContainerConfig(backend, cache = cache)
  config$executable <- "/fake/runtime"
  config$version <- "test runtime"
  config$resolved[["test/image:1"]] <- list(source = "test/image:1",
    reference = "/images/test image.sif", identity = "sha256:test")
  config
}

test_that("automatic backend selection is deterministic and batch-safe", {
  expect_identical(.selectContainerBackend("auto", "docker", FALSE), "docker")
  expect_identical(.selectContainerBackend("auto", c("singularity", "apptainer"), FALSE), "apptainer")
  expect_identical(.selectContainerBackend("singularity", c("singularity", "docker"), FALSE), "singularity")
  expect_error(.selectContainerBackend("auto", c("docker", "apptainer"), FALSE), "explicitly")
  expect_error(.selectContainerBackend("auto", character(), FALSE), "No container runtime")
  expect_error(.selectContainerBackend("docker", "apptainer", FALSE), "not on PATH")
  expect_identical(.selectContainerBackend("auto", c("docker", "apptainer"), TRUE,
    choose = function(...) 2L), "apptainer")
  expect_error(.selectContainerBackend("auto", c("docker", "apptainer"), TRUE,
    choose = function(...) 0L), "cancelled")
})

test_that("Docker daemon failures do not silently switch runtimes", {
  local_mocked_bindings(
    .findContainerExecutables = function(x) stats::setNames(paste0("/bin/", x), x),
    .executeContainerProcess = function(...) 1L
  )
  expect_error(.resolveContainerRuntime(.createContainerConfig("docker")), "daemon cannot be reached")
})

test_that("runtime commands preserve binds, working directories, and platform rules", {
  config <- .createContainerFixture()
  command <- .buildContainerCommand(config, "test/image:1", "panaroo-merge", c("-t", "6"),
    c("/host dir:/work", "/db:/db:ro"), "/work", "linux/amd64")
  expect_identical(command$args, c("exec", "--cleanenv", "--bind", "/host dir:/work",
    "--bind", "/db:/db:ro", "--pwd", "/work", "/images/test image.sif", "panaroo-merge", "-t", "6"))
  config$backend <- "singularity"
  expect_identical(.buildContainerCommand(config, "test/image:1", "hmmpress")$args[[1]], "exec")
  config$backend <- "docker"
  command <- .buildContainerCommand(config, "test/image:1", "cd-hit", platform = "linux/amd64")
  expect_identical(command$args, c("run", "--rm", "--platform", "linux/amd64", "sha256:test", "cd-hit"))
})

test_that("legacy call-site adapter preserves tool arguments and rejects unknown runtime options", {
  captured <- NULL
  local_mocked_bindings(.runContainer = function(config, image, command, args,
      binds, workdir, platform, stdout, stderr) {
    captured <<- list(image = image, command = command, args = args, binds = binds,
      workdir = workdir, platform = platform, stdout = stdout, stderr = stderr)
    7L
  })
  result <- .runContainerCommand(.createContainerFixture(), c("run", "--platform", "linux/amd64",
    "--rm", "-v", "/host space:/work", "-w", "/work", "test/image:1", "sh", "-lc", "echo '$HOME'"),
    stdout = FALSE, stderr = "/tmp/errors")
  expect_identical(result, 7L)
  expect_identical(captured$args, c("-lc", "echo '$HOME'"))
  expect_identical(captured$binds, "/host space:/work")
  expect_identical(captured$stderr, "/tmp/errors")
  expect_error(.runContainerCommand(.createContainerFixture(), c("run", "--privileged", "image", "cmd")), "Unsupported")
})

test_that("local SIF overrides are reusable and missing images fail before execution", {
  sif <- tempfile(fileext = ".sif")
  writeLines("test image", sif)
  on.exit(unlink(sif), add = TRUE)
  config <- .createContainerFixture()
  config$images <- c(panaroo = sif)
  config <- .prepareContainerImages(config, "staphb/panaroo:1.7.0")
  expect_identical(config$resolved[["staphb/panaroo:1.7.0"]]$reference, normalizePath(sif))
  expect_identical(config$resolved[["staphb/panaroo:1.7.0"]]$identity, unname(tools::md5sum(sif)))
  config$backend <- "docker"
  config$resolved <- list()
  expect_error(.prepareContainerImages(config, "staphb/panaroo:1.7.0"), "Local SIF")
  config$backend <- "apptainer"
  config$images <- c(panaroo = paste0(sif, "-missing.sif"))
  expect_error(.prepareContainerImages(config, "staphb/panaroo:1.7.0"))
})

test_that("SIF pulls happen once and failed pulls do not publish partial images", {
  config <- .createContainerFixture(cache = tempfile())
  config$resolved <- list()
  on.exit(unlink(config$cache, recursive = TRUE), add = TRUE)
  calls <- 0L
  fail <- FALSE
  local_mocked_bindings(.executeContainerProcess = function(command, args, ...) {
    calls <<- calls + 1L
    path <- substring(args[[2]], 2L, nchar(args[[2]]) - 1L)
    writeLines("image bytes", path)
    if (fail) 1L else 0L
  })
  prepared <- .prepareContainerImages(config, rep("test/image:1", 3))
  expect_identical(calls, 1L)
  prepared_again <- .prepareContainerImages(config, "test/image:1")
  expect_identical(calls, 1L)
  expect_identical(prepared$resolved, prepared_again$resolved)
  fail <- TRUE
  expect_error(.prepareContainerImages(config, "test/image:2"), "Failed to pull")
  expect_length(list.files(config$cache, all.files = TRUE, pattern = "pull-"), 0L)
  expect_length(list.files(config$cache, pattern = "\\.sif$"), 1L)
})

test_that("image changes invalidate reuse, including unverified runtime switches", {
  config <- .createContainerFixture()
  stage <- list(tool = list(docker_image = "test/image:1",
    container = .getContainerProvenance(config, "test/image:1")))
  expect_true(.isContainerStageReusable(stage, config))
  other <- config
  other$resolved[["test/image:1"]]$identity <- "changed"
  expect_false(.isContainerStageReusable(stage, other))
  other <- config
  other$backend <- "docker"
  expect_false(.isContainerStageReusable(stage, other))
  stage$tool$container <- NULL
  expect_false(.isContainerStageReusable(stage, config))
  expect_true(.isContainerStageReusable(stage, other))
})

test_that("real process arguments survive shell quoting and failed commands are logged", {
  skip_on_os("windows")
  script <- tempfile(pattern = "runtime with spaces ")
  writeLines(c("#!/bin/sh", "printf '%s\\n' \"$@\"", "exit 7"), script)
  Sys.chmod(script, "0755")
  config <- .createContainerFixture()
  config$executable <- script
  config$log_dir <- tempfile()
  on.exit(unlink(c(script, config$log_dir), recursive = TRUE), add = TRUE)
  arguments <- c("space here", "$(touch SHOULD_NOT_EXIST)", "a;b", "'quoted'", "`literal`")
  output <- suppressWarnings(.runContainer(config, "test/image:1", "tool", arguments))
  expect_identical(tail(as.character(output), length(arguments)), arguments)
  expect_identical(attr(output, "status"), 7L)
  record <- jsonlite::read_json(list.files(config$log_dir, full.names = TRUE)[[1]])
  expect_identical(record$exit_status, 7L)
  expect_identical(record$backend, "apptainer")
})

test_that("resolved configuration runs in a separate R worker without rediscovery", {
  skip_on_os("windows")
  script <- tempfile(pattern = "worker runtime ")
  writeLines(c("#!/bin/sh", "printf '%s\\n' \"$@\""), script)
  Sys.chmod(script, "0755")
  on.exit(unlink(script), add = TRUE)
  config <- .createContainerFixture()
  config$executable <- script
  source_file <- normalizePath(testthat::test_path("..", "..", "R", "containers.R"))
  old <- future::plan()
  on.exit(future::plan(old), add = TRUE)
  future::plan(future::multisession, workers = 2L)
  result <- future::value(future::future({
    env <- new.env(parent = baseenv())
    sys.source(source_file, envir = env)
    env$.runContainer(config, "test/image:1", "worker-tool", "a b")
  }, globals = list(config = config, source_file = source_file)))
  expect_identical(tail(result, 2L), c("worker-tool", "a b"))
})
