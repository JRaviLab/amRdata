.mock_genome_tbl <- function(genome_id, checkm_completeness, checkm_contamination,
                              genome_name = "x", species = "x",
                              genome_length = "4600000", gc_content = "50.5", cds = "4300") {
  tibble::tibble(
    `genome.genome_id` = genome_id,
    `genome.genome_name` = genome_name,
    `genome.species` = species,
    `genome.genome_length` = genome_length,
    `genome.gc_content` = gc_content,
    `genome.cds` = cds,
    `genome.checkm_completeness` = checkm_completeness,
    `genome.checkm_contamination` = checkm_contamination
  )
}

test_that(".apply_metadata_qc() keeps a genome that passes the default CheckM gate", {
  tbl <- .mock_genome_tbl("511145.12", checkm_completeness = "99", checkm_contamination = "1")

  out <- .apply_metadata_qc(tbl)

  expect_identical(out$keep_ids, "511145.12")
  expect_equal(nrow(out$rejections), 0L)
})

test_that(".apply_metadata_qc() drops a genome with contamination above the threshold", {
  tbl <- .mock_genome_tbl("511145.12", checkm_completeness = "99", checkm_contamination = "12")

  out <- .apply_metadata_qc(tbl, checkm_contam = 5, checkm_complete = 95)

  expect_length(out$keep_ids, 0L)
  expect_identical(out$rejections$failed_rule, "checkm_contamination")
  expect_identical(out$rejections$observed, "12")
})

test_that(".apply_metadata_qc() drops a genome with completeness below the threshold", {
  tbl <- .mock_genome_tbl("511145.12", checkm_completeness = "80", checkm_contamination = "1")

  out <- .apply_metadata_qc(tbl, checkm_contam = 5, checkm_complete = 95)

  expect_length(out$keep_ids, 0L)
  expect_identical(out$rejections$failed_rule, "checkm_completeness")
})

test_that(".apply_metadata_qc() drops a genome with missing CheckM fields", {
  tbl <- .mock_genome_tbl("511145.12", checkm_completeness = NA_character_, checkm_contamination = "1")

  out <- .apply_metadata_qc(tbl)

  expect_length(out$keep_ids, 0L)
  expect_identical(out$rejections$failed_rule, "checkm_missing")
})

test_that(".apply_metadata_qc() keeps only the passing genomes in a mixed batch", {
  tbl <- .mock_genome_tbl(
    genome_id = c("good.1", "bad_contam.1", "bad_complete.1"),
    checkm_completeness = c("99", "99", "50"),
    checkm_contamination = c("1", "20", "1")
  )

  out <- .apply_metadata_qc(tbl)

  expect_identical(out$keep_ids, "good.1")
  expect_setequal(out$rejections$genome_id, c("bad_contam.1", "bad_complete.1"))
})
