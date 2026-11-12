.createPanarooRecoveryFixture <- function() {
  root <- tempfile()
  for (part in c("work", "orb", "panaroo")) dir.create(file.path(root, part), recursive = TRUE)
  db <- file.path(root, "work", "test.duckdb")
  file.create(db)
  batches <- file.path(root, "panaroo", paste0("panaroo_out_", 1:3))
  for (batch in batches) {
    dir.create(batch)
    writeLines("graph", file.path(batch, "final_graph.gml"))
    writeLines("genes", file.path(batch, "gene_data.csv"))
  }
  parameters <- list(core_threshold = 0.9, len_dif_percent = 0.95,
    cluster_threshold = 0.95, family_seq_identity = 0.5, split_jobs = TRUE,
    refind_mode = "off", strip_pseudogenes = FALSE, remove_all_x_cds = TRUE)
  original <- list(schema_version = 1L, manifest_type = "amR_dataset", manifest_id = "test",
    dataset_id = "test", dataset = list(duckdb = db), artifacts = list(),
    runs = list(list(run_id = "failed-original", status = "failed", stages = list(
      list(name = "panaroo", status = "failed", parameters = parameters)), events = list())))
  manifest <- file.path(root, "orb", "manifest_test.json")
  jsonlite::write_json(original, manifest, auto_unbox = TRUE)
  config <- .createContainerConfig("apptainer")
  config$executable <- "/fake/apptainer"
  config$version <- "fake"
  config$resolved[["staphb/panaroo:1.7.0"]] <- list(source = "staphb/panaroo:1.7.0",
    reference = "/images/panaroo.sif", identity = "image-checksum")
  list(root = root, db = db, batches = batches, manifest = manifest, config = config,
       parameters = parameters, original = jsonlite::read_json(manifest))
}

test_that("merge recovery selects exact batches and produces a reusable Panaroo stage", {
  fixture <- .createPanarooRecoveryFixture()
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)
  partial <- file.path(fixture$root, "panaroo", "merge_output")
  dir.create(partial)
  writeLines("partial", file.path(partial, "partial.txt"))
  selected <- NULL
  imported <- FALSE
  local_mocked_bindings(
    .mergePanaroo = function(input_path, batch_dirs, container_config, ...) {
      selected <<- batch_dirs
      dir.create(file.path(input_path, "merge_output"))
      expect_identical(container_config$backend, "apptainer")
    },
    .panaroo2duckdb = function(panaroo_output_path, duckdb_path) {
      imported <<- TRUE
      expect_identical(duckdb_path, normalizePath(fixture$db))
    }
  )
  result <- resumePanarooMerge(fixture$db, fixture$batches[1:2],
    manifest_path = fixture$manifest, container_backend = fixture$config)
  expect_identical(result, normalizePath(fixture$db))
  expect_identical(selected, normalizePath(fixture$batches[1:2]))
  expect_true(imported)
  saved <- jsonlite::read_json(fixture$manifest)
  expect_identical(saved$runs[[1]], fixture$original$runs[[1]])
  expect_identical(saved$runs[[2]]$status, "success")
  plan <- .resume_plan(saved$runs[[2]], c("panaroo", "cdhit", "hmmer"),
    list(panaroo = fixture$parameters), container_config = fixture$config)
  expect_identical(plan, c(panaroo = TRUE, cdhit = FALSE, hmmer = FALSE))
  expect_length(list.files(dirname(partial), pattern = "^merge_output_failed_"), 1L)
  expect_length(list.files(dirname(fixture$manifest), pattern = "before-merge-.*bak$"), 1L)
})

test_that("merge and import failures are recorded without altering original provenance", {
  for (phase in c("merge", "import")) {
    fixture <- .createPanarooRecoveryFixture()
    local_mocked_bindings(
      .mergePanaroo = function(input_path, ...) {
        if (phase == "merge") stop("exit status 137")
        dir.create(file.path(input_path, "merge_output"))
      },
      .panaroo2duckdb = function(...) stop("import failed")
    )
    expect_error(resumePanarooMerge(fixture$db, fixture$batches[1:2],
      manifest_path = fixture$manifest, container_backend = fixture$config),
      if (phase == "merge") "137" else "import failed")
    saved <- jsonlite::read_json(fixture$manifest)
    expect_identical(saved$runs[[1]], fixture$original$runs[[1]])
    expect_identical(saved$runs[[2]]$status, "failed")
    expect_identical(saved$runs[[2]]$stages[[1]]$status, "failed")
    expect_true(all(file.exists(file.path(fixture$batches, "final_graph.gml"))))
    unlink(fixture$root, recursive = TRUE)
  }
})

test_that("invalid batch selections are rejected before appending a recovery run", {
  fixture <- .createPanarooRecoveryFixture()
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)
  for (selection in list(fixture$batches[1], rep(fixture$batches[1], 2))) {
    expect_error(resumePanarooMerge(fixture$db, selection, manifest_path = fixture$manifest,
      container_backend = fixture$config))
  }
  unlink(file.path(fixture$batches[1], "gene_data.csv"))
  expect_error(resumePanarooMerge(fixture$db, fixture$batches[1:2],
    manifest_path = fixture$manifest, container_backend = fixture$config), "nonempty")
  expect_length(jsonlite::read_json(fixture$manifest)$runs, 1L)
})

test_that("the merge launcher excludes unselected completed batches", {
  fixture <- .createPanarooRecoveryFixture()
  on.exit(unlink(fixture$root, recursive = TRUE), add = TRUE)
  captured <- NULL
  local_mocked_bindings(.runContainerCommand = function(config, args, ...) {
    captured <<- args
    character()
  })
  .mergePanaroo(file.path(fixture$root, "panaroo"), batch_dirs = fixture$batches[1:2],
    container_config = fixture$config)
  expect_true(all(paste0("/work/", basename(fixture$batches[1:2])) %in% captured))
  expect_false(paste0("/work/", basename(fixture$batches[3])) %in% captured)
  expect_error(.mergePanaroo(file.path(fixture$root, "panaroo"),
    batch_dirs = rep(fixture$batches[1], 2), container_config = fixture$config), "unique")
})
