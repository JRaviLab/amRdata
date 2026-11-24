.createHMMERCacheFixture <- function() {
  root <- tempfile("hmmer-cache-")
  dir.create(root)
  bfc <- BiocFileCache::BiocFileCache(file.path(root, "bfc"), ask = FALSE)
  list(root = root, bfc = bfc, default = file.path(root, "default"))
}

.writeHMMERCacheFixture <- function(root, name = "AMRFinder") {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  hmm <- file.path(root, paste0(name, ".hmm"))
  writeLines(c("HMMER3/f [3.4]", paste("NAME", name), "ACC TEST.1", "GA 20 10;", "//"), hmm)
  Sys.setFileTime(hmm, Sys.time() - 60)
  for (suffix in c(".h3m", ".h3i", ".h3f", ".h3p")) writeLines("pressed", paste0(hmm, suffix))
  normalizePath(hmm)
}

test_that("default HMM preparation discovers a registered custom path without downloads", {
  fixture <- .createHMMERCacheFixture()
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)
  local_mocked_bindings(
    .amr_bfc = function() fixture$bfc,
    .defaultHmmerDbDir = function() fixture$default,
    .prepareContainerImages = function(config, ...) config,
    .amr_download_file = function(...) stop("unexpected download"),
    .runContainerCommand = function(...) stop("unexpected hmmpress")
  )
  hmm <- .writeHMMERCacheFixture(file.path(fixture$root, "project"))
  .amr_bfc_register_hmmer("AMRFinder", hmm)
  log <- file.path(fixture$root, "preparation.log")
  result <- .prepareHmmerDatabases(databases = "AMRFinder", progress = FALSE,
    verbose = FALSE, log_path = log)$AMRFinder
  expect_identical(result$hmm, hmm)
  expect_true(file.exists(result$profiles))
  expect_true(all(c("ga_seq", "tc_dom") %in% names(arrow::read_parquet(result$profiles))))
  expect_match(paste(readLines(log), collapse = "\n"), "Reusing AMRFinder.*BiocFileCache")
  expect_length(list.files(fixture$default, pattern = "\\.hmm$", recursive = TRUE), 0L)
  expect_identical(.amr_bfc_find(fixture$bfc, "amR_hmmer_AMRFinder")$rpath[[1]], hmm)
})

test_that("explicit directories do not silently use a registered database elsewhere", {
  fixture <- .createHMMERCacheFixture()
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)
  local_mocked_bindings(
    .amr_bfc = function() fixture$bfc,
    .defaultHmmerDbDir = function() fixture$default,
    .prepareContainerImages = function(config, ...) config,
    .amr_download_file = function(url, destfile, progress, label, ...) {
      expect_identical(label, "AMRFinder HMM database")
      stop("download into explicit directory")
    }
  )
  registered <- .writeHMMERCacheFixture(file.path(fixture$root, "registered"))
  .amr_bfc_register_hmmer("AMRFinder", registered)
  explicit <- file.path(fixture$root, "explicit")
  expect_error(.prepareHmmerDatabases(explicit, databases = "AMRFinder", progress = FALSE,
    verbose = FALSE), "download into explicit directory")
  local <- .writeHMMERCacheFixture(file.path(explicit, "AMRFinder"))
  result <- .prepareHmmerDatabases(explicit, databases = "AMRFinder", progress = FALSE,
    verbose = FALSE)$AMRFinder
  expect_identical(result$hmm, local)
  expect_true(file.exists(registered))
})

test_that("a database already in the default directory takes precedence over the registry", {
  fixture <- .createHMMERCacheFixture()
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)
  local_mocked_bindings(
    .amr_bfc = function() fixture$bfc,
    .defaultHmmerDbDir = function() fixture$default,
    .prepareContainerImages = function(config, ...) config,
    .amr_download_file = function(...) stop("unexpected download")
  )
  registered <- .writeHMMERCacheFixture(file.path(fixture$root, "registered"))
  local <- .writeHMMERCacheFixture(file.path(fixture$default, "AMRFinder"))
  .amr_bfc_register_hmmer("AMRFinder", registered)
  result <- .prepareHmmerDatabases(databases = "AMRFinder", progress = FALSE, verbose = FALSE)
  expect_identical(result$AMRFinder$hmm, local)
})

test_that("stale, empty, and invalid BFC paths fall back to normal preparation", {
  for (state in c("missing", "empty", "invalid")) {
    fixture <- .createHMMERCacheFixture()
    local_mocked_bindings(
      .amr_bfc = function() fixture$bfc,
      .defaultHmmerDbDir = function() fixture$default,
      .prepareContainerImages = function(config, ...) config,
      .amr_download_file = function(...) stop("normal download fallback")
    )
    hmm <- .writeHMMERCacheFixture(file.path(fixture$root, "registered"))
    .amr_bfc_register_hmmer("AMRFinder", hmm)
    if (state == "missing") unlink(hmm)
    if (state == "empty") writeLines(character(), hmm)
    if (state == "invalid") writeLines("not an HMM database", hmm)
    expect_null(.findRegisteredHMMERDatabase("AMRFinder"))
    expect_error(.prepareHmmerDatabases(databases = "AMRFinder", progress = FALSE,
      verbose = FALSE), "normal download fallback")
    unlink(fixture$root, recursive = TRUE)
  }
})

test_that("missing or stale indexes are rebuilt beside the registered HMM", {
  for (state in c("missing", "stale")) {
    fixture <- .createHMMERCacheFixture()
    presses <- 0L
    local_mocked_bindings(
      .amr_bfc = function() fixture$bfc,
      .defaultHmmerDbDir = function() fixture$default,
      .prepareContainerImages = function(config, ...) config,
      .amr_download_file = function(...) stop("unexpected download"),
      .runContainerCommand = function(config, args, ...) {
        expect_true(all(c("hmmpress", "-f") %in% args))
        presses <<- presses + 1L
        for (suffix in c(".h3m", ".h3i", ".h3f", ".h3p")) writeLines("repaired", paste0(hmm, suffix))
        character()
      }
    )
    hmm <- .writeHMMERCacheFixture(file.path(fixture$root, "registered"))
    .amr_bfc_register_hmmer("AMRFinder", hmm)
    if (state == "missing") unlink(paste0(hmm, ".h3m"))
    if (state == "stale") Sys.setFileTime(paste0(hmm, ".h3m"), Sys.time() - 120)
    result <- .prepareHmmerDatabases(databases = "AMRFinder", progress = FALSE, verbose = FALSE)
    expect_identical(result$AMRFinder$hmm, hmm)
    expect_identical(presses, 1L)
    expect_true(file.exists(hmm))
    unlink(fixture$root, recursive = TRUE)
  }
})

test_that("custom database URLs are not replaced by unrelated registered resources", {
  fixture <- .createHMMERCacheFixture()
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)
  local_mocked_bindings(
    .defaultHmmerDbDir = function() fixture$default,
    .prepareContainerImages = function(config, ...) config,
    .findRegisteredHMMERDatabase = function(...) stop("unexpected registry lookup"),
    .amr_download_file = function(url, ...) stop(paste("requested source", url))
  )
  expect_error(.prepareHmmerDatabases(databases = "Custom", hmmer_db_url = "https://example.invalid/custom.gz",
    progress = FALSE, verbose = FALSE), "requested source https://example.invalid/custom.gz", fixed = TRUE)
})

test_that("download labels identify resources even with verbose output disabled", {
  fixture <- tempfile()
  on.exit(unlink(fixture), add = TRUE)
  local_mocked_bindings(req_perform = function(...) invisible(NULL), .package = "httr2")
  expect_message(.amr_download_file("https://example.invalid/pfam.gz", tempfile(),
    progress = TRUE, label = "Pfam HMM database", log_path = fixture, verbose = FALSE),
    "Downloading Pfam HMM database")
  expect_match(paste(readLines(fixture), collapse = "\n"), "https://example.invalid/pfam.gz", fixed = TRUE)
  expect_silent(.amr_download_file("https://example.invalid/cog.zip", tempfile(),
    progress = FALSE, label = "COG HMM database", log_path = fixture, verbose = FALSE))
  expect_match(paste(readLines(fixture), collapse = "\n"), "Downloading COG HMM database")
})

test_that("DefenseCas reuses both registered components without downloading repositories", {
  fixture <- .createHMMERCacheFixture()
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)
  hits <- file.path(fixture$root, "hits")
  searches <- 0L
  local_mocked_bindings(
    .amr_bfc = function() fixture$bfc,
    .prepareContainerImages = function(config, ...) config,
    .amr_download_file = function(...) stop("unexpected repository download"),
    .runContainerCommand = function(config, args, ...) {
      expect_true("hmmsearch" %in% args)
      searches <<- searches + 1L
      output <- args[[match("--domtblout", args) + 1L]]
      writeLines("# synthetic search output", file.path(hits, basename(output)))
      0L
    },
    .parseHMMEROutput = function(...) tibble::tibble(protein = "protein1"),
    .hmmerHits = function(x) x
  )
  defense <- .writeHMMERCacheFixture(file.path(fixture$root, "shared"), "DefenseFinder")
  cas <- .writeHMMERCacheFixture(file.path(fixture$root, "shared"), "CasFinder")
  .amr_bfc_register_hmmer("DefenseCas", defense, component = "DefenseFinder")
  .amr_bfc_register_hmmer("DefenseCas", cas, component = "CasFinder")
  dir.create(file.path(fixture$root, "work"))
  db <- file.path(fixture$root, "work", "test.duckdb")
  con <- DBI::dbConnect(duckdb::duckdb(), db)
  DBI::dbWriteTable(con, "protein_cluster_seq", data.frame(name = "protein1", sequence = "MKTAYIAK"))
  DBI::dbDisconnect(con, shutdown = TRUE)
  result <- .defenseHMMER(file.path(fixture$default, "DefenseCas"), duckdb_path = db,
    output_path = hits, threads = 1L, reuse_registered = TRUE, verbose = FALSE, progress = FALSE)
  expect_identical(result$databases$DefenseFinder$hmm, defense)
  expect_identical(result$databases$CasFinder$hmm, cas)
  expect_identical(searches, 2L)
  expect_false(dir.exists(file.path(fixture$default, "DefenseCas", "DefenseFinder")))
})
