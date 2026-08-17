# Tests for the BV-BRC Data API download path (R/bvbrc_api.R).
# Pure helpers run offline; the live extractor tests skip when offline/CRAN.

test_that(".bvbrc_chunk splits into size-n groups", {
  expect_length(.bvbrc_chunk(1:10, 3), 4)
  expect_length(.bvbrc_chunk(character(0), 5), 0)
  expect_identical(unname(.bvbrc_chunk(1:3, 10)[[1]]), 1:3)
})

test_that(".bvbrc_prefix_fill fills missing fields, prefixes, coerces to character", {
  df <- data.frame(
    genome_id = "1280.1", antibiotic = "ciprofloxacin",
    stringsAsFactors = FALSE
  )
  out <- .bvbrc_prefix_fill(df, c("genome_id", "antibiotic", "source"), "genome_drug")

  expect_identical(
    names(out),
    c("genome_drug.genome_id", "genome_drug.antibiotic", "genome_drug.source")
  )
  expect_true(is.na(out[["genome_drug.source"]])) # missing field -> NA
  expect_type(out[["genome_drug.genome_id"]], "character")
})

test_that(".extractAMRtable_api returns genome_drug.* columns keyed by genome_id (live)", {
  skip_on_cran()
  skip_if_offline("www.bv-brc.org")

  amr <- .extractAMRtable_api("1280.15865", verbose = FALSE)
  expect_s3_class(amr, "data.frame")
  expect_true("genome_drug.genome_id" %in% names(amr))
  expect_true(all(grepl("^genome_drug[.]", names(amr))))
  expect_gt(nrow(amr), 0)
  expect_true(all(amr[["genome_drug.genome_id"]] == "1280.15865"))
})

test_that(".extractGenomeData_api returns genome.* columns incl. QC fields (live)", {
  skip_on_cran()
  skip_if_offline("www.bv-brc.org")

  gen <- .extractGenomeData_api(
    "1280.15865",
    fields = "species,genome_quality,checkm_completeness,cds",
    verbose = FALSE
  )
  expect_true(all(
    c("genome.genome_id", "genome.species", "genome.checkm_completeness") %in%
      names(gen)
  ))
  expect_equal(nrow(gen), 1L)
})
