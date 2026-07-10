# Regression tests for genome_id_file handling.
#
# PR #29: prepareGenomes() accepted a genome_id_file argument but hardcoded
# genome_id_file = NULL in its internal retrieveMetadata() call, so a
# user-supplied genome list was silently ignored via the wrapper (it worked
# only when calling retrieveMetadata() directly).

test_that("prepareGenomes() forwards genome_id_file to retrieveMetadata()", {
  captured <- new.env()
  gfile <- tempfile(fileext = ".txt")
  writeLines(c("511145.12", "83333.111"), gfile)
  on.exit(unlink(gfile), add = TRUE)

  # Stub everything prepareGenomes() touches so we isolate the pass-through.
  # retrieveMetadata() records the genome_id_file it receives; .filterGenomes()
  # returns NULL so prepareGenomes() exits cleanly right after.
  local_mocked_bindings(
    .ensure_bvbrc_cache = function(...) invisible(NULL),
    retrieveMetadata = function(..., genome_id_file = NULL) {
      captured$gid <- genome_id_file
      invisible(NULL)
    },
    .filterGenomes = function(...) NULL,
    .package = "amRdata"
  )

  res <- prepareGenomes(
    user_bacs = "Shigella flexneri",
    genome_id_file = gfile,
    verbose = FALSE
  )

  expect_null(res)
  # Before the fix this was NULL (the bug); now it is the file path.
  expect_identical(captured$gid, gfile)
})

test_that("retrieveMetadata() errors on a missing genome_id_file", {
  missing <- file.path(tempdir(), "definitely_not_here_9f8a.txt")
  expect_error(
    retrieveMetadata(
      user_bacs = "x",
      genome_id_file = missing,
      verbose = FALSE
    ),
    "does not exist"
  )
})

test_that("genome_id_file parsing trims whitespace and drops blank lines", {
  # Mirrors the readLines/trimws/drop-empty/unique logic in retrieveMetadata().
  gfile <- tempfile(fileext = ".txt")
  writeLines(c("511145.12", "  83333.111  ", "", "   ", "511145.12"), gfile)
  on.exit(unlink(gfile), add = TRUE)

  ids <- readLines(gfile, warn = FALSE)
  ids <- trimws(ids)
  ids <- ids[ids != ""]
  ids <- unique(as.character(ids))

  expect_identical(ids, c("511145.12", "83333.111"))
})
