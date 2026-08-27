test_that("buildDyadFeatureMap() errors when the DuckDB path does not exist", {
  expect_error(
    buildDyadFeatureMap(duckdb_path = tempfile(fileext = ".duckdb"))
  )
})

test_that("buildDyadFeatureMap() errors when no provenance manifest is found", {
  tmp_dir <- file.path(tempdir(), paste0("dyad-test-", as.integer(runif(1, 1, 1e6))))
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  fake_db <- file.path(tmp_dir, "fake.duckdb")
  file.create(fake_db)

  expect_error(
    buildDyadFeatureMap(duckdb_path = fake_db),
    regexp = "provenance manifest"
  )
})

# TODO (Bioconductor): add an end-to-end test that runs buildDyadFeatureMap()
# against a small fixture DuckDB + manifest + annotation Parquets and checks the
# resulting edge list (source/target columns, dyad format, type prefixes).
