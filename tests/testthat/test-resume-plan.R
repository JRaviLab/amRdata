test_that(".resume_plan() returns all FALSE when there is no prior run", {
  plan <- amRdata:::.resume_plan(NULL, c("panaroo", "cdhit", "hmmer"))

  expect_equal(plan, c(panaroo = FALSE, cdhit = FALSE, hmmer = FALSE))
})

test_that(".resume_plan() skips a contiguous run of successful stages", {
  prev_run <- list(
    stages = list(
      list(name = "panaroo", status = "success", outputs = list()),
      list(name = "cdhit", status = "success", outputs = list())
    )
  )

  plan <- amRdata:::.resume_plan(prev_run, c("panaroo", "cdhit", "hmmer"))

  expect_equal(plan, c(panaroo = TRUE, cdhit = TRUE, hmmer = FALSE))
})

test_that(".resume_plan() stops at the first stage that did not succeed", {
  prev_run <- list(
    stages = list(
      list(name = "panaroo", status = "success", outputs = list()),
      list(name = "cdhit", status = "failed", outputs = list()),
      list(name = "hmmer", status = "success", outputs = list())
    )
  )

  plan <- amRdata:::.resume_plan(prev_run, c("panaroo", "cdhit", "hmmer"))

  expect_equal(plan, c(panaroo = TRUE, cdhit = FALSE, hmmer = FALSE))
})

test_that(".resume_plan() distrusts a success whose recorded output is missing", {
  missing_path <- file.path(tempdir(), paste0("missing-", as.integer(runif(1, 1, 1e6))))

  prev_run <- list(
    stages = list(
      list(
        name = "panaroo",
        status = "success",
        outputs = list(list(path = missing_path))
      ),
      list(name = "cdhit", status = "success", outputs = list())
    )
  )

  plan <- amRdata:::.resume_plan(prev_run, c("panaroo", "cdhit", "hmmer"))

  expect_equal(plan, c(panaroo = FALSE, cdhit = FALSE, hmmer = FALSE))
})

test_that(".manifest_prior_stage() finds the latest entry for a stage name", {
  prev_run <- list(
    stages = list(
      list(name = "panaroo", status = "running", outputs = list()),
      list(name = "panaroo", status = "success", outputs = list())
    )
  )

  stage <- amRdata:::.manifest_prior_stage(prev_run, "panaroo")

  expect_equal(stage$status, "success")
})

test_that(".manifest_prior_stage() returns NULL when the stage isn't recorded", {
  prev_run <- list(stages = list(list(name = "panaroo", status = "success", outputs = list())))

  expect_null(amRdata:::.manifest_prior_stage(prev_run, "hmmer"))
  expect_null(amRdata:::.manifest_prior_stage(NULL, "panaroo"))
})

test_that(".manifest_migrate_legacy() backfills manifest_type/manifest_id/artifacts", {
  legacy <- list(
    schema_version = 1L,
    dataset_id = "Sfl",
    dataset = list(duckdb = "Sfl.duckdb"),
    runs = list()
  )

  migrated <- amRdata:::.manifest_migrate_legacy(legacy, "/tmp/manifest_Sfl.json")

  expect_identical(migrated$manifest_type, "amR_dataset")
  expect_identical(migrated$manifest_id, "manifest_Sfl")
  expect_identical(migrated$artifacts, list())
})

test_that(".manifest_migrate_legacy() leaves an already-current manifest untouched", {
  current <- list(
    schema_version = 1L,
    manifest_type = "amR_dataset",
    manifest_id = "abc123",
    dataset_id = "Sfl",
    dataset = list(duckdb = "Sfl.duckdb"),
    artifacts = list(foo = "bar"),
    runs = list()
  )

  migrated <- amRdata:::.manifest_migrate_legacy(current, "/tmp/manifest_Sfl.json")

  expect_identical(migrated, current)
})

test_that(".manifest_resume() can resume a manifest written before manifest_type existed", {
  manifest_path <- tempfile(fileext = ".json")
  on.exit(unlink(manifest_path), add = TRUE)

  legacy <- list(
    schema_version = 1L,
    dataset_id = "Sfl",
    dataset = list(duckdb = "Sfl.duckdb", selection = list()),
    runs = list()
  )

  jsonlite::write_json(legacy, manifest_path, auto_unbox = TRUE, pretty = TRUE, null = "null")

  manifest_state <- amRdata:::.manifest_resume(manifest_path, base_dir = tempdir())

  expect_s3_class(manifest_state, "amr_manifest")
  expect_identical(manifest_state$manifest$manifest_type, "amR_dataset")
})
