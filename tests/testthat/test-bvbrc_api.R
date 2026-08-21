# Tests for the BV-BRC Data API download path (R/bvbrc_api.R).
# Pure helpers run offline; the live extractor tests skip when offline/CRAN.

test_that(".bvbrcChunk splits into size-n groups", {
  expect_length(.bvbrcChunk(1:10, 3), 4)
  expect_length(.bvbrcChunk(character(0), 5), 0)
  expect_identical(unname(.bvbrcChunk(1:3, 10)[[1]]), 1:3)
})

test_that(".bvbrcPrefixFill fills missing fields, prefixes, coerces to character", {
  df <- data.frame(
    genome_id = "1280.1", antibiotic = "ciprofloxacin",
    stringsAsFactors = FALSE
  )
  out <- .bvbrcPrefixFill(df, c("genome_id", "antibiotic", "source"), "genome_drug")

  expect_identical(
    names(out),
    c("genome_drug.genome_id", "genome_drug.antibiotic", "genome_drug.source")
  )
  expect_identical(out[["genome_drug.source"]], "") # missing field -> "" (Docker convention)
  expect_type(out[["genome_drug.genome_id"]], "character")
})

test_that(".extractAMRtableApi returns genome_drug.* columns keyed by genome_id (live)", {
  skip_on_cran()
  skip_if_offline("www.bv-brc.org")

  amr <- .extractAMRtableApi("1280.15865", verbose = FALSE)
  expect_s3_class(amr, "data.frame")
  expect_true("genome_drug.genome_id" %in% names(amr))
  expect_true(all(grepl("^genome_drug[.]", names(amr))))
  expect_gt(nrow(amr), 0)
  expect_true(all(amr[["genome_drug.genome_id"]] == "1280.15865"))
})

test_that(".bvbrcPrefixFill handles a zero-row frame (empty query result)", {
  out <- .bvbrcPrefixFill(data.frame(), c("genome_id", "antibiotic"), "genome_drug")
  expect_identical(names(out), c("genome_drug.genome_id", "genome_drug.antibiotic"))
  expect_equal(nrow(out), 0L)
})

test_that(".resolveGenomeIDsApi resolves a species to valid genome IDs (live)", {
  skip_on_cran()
  skip_if_offline("www.bv-brc.org")

  td <- file.path(tempdir(), paste0("res_", as.integer(runif(1, 1, 1e6))))
  dir.create(td, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  ids <- .resolveGenomeIDsApi(
    base_dir = td, user_bacs = "Morganella morganii",
    overwrite = TRUE, verbose = FALSE
  )
  expect_type(ids, "character")
  expect_gt(length(ids), 0L)
  expect_true(all(grepl("^[0-9]+[.][0-9]+$", ids)))
  expect_false(any(duplicated(ids)))
})

test_that(".resolveGenomeIDsApi warns (not errors) on a zero-match input", {
  skip_on_cran()
  skip_if_offline("www.bv-brc.org")

  td <- file.path(tempdir(), paste0("res0_", as.integer(runif(1, 1, 1e6))))
  dir.create(td, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  expect_warning(
    ids <- .resolveGenomeIDsApi(
      base_dir = td, user_bacs = "Nosuchgenusxyzabc",
      overwrite = TRUE, verbose = FALSE
    ),
    "resolved 0 genomes"
  )
  expect_identical(ids, character(0))
})

test_that(".extractGenomeDataApi returns genome.* columns incl. QC fields (live)", {
  skip_on_cran()
  skip_if_offline("www.bv-brc.org")

  gen <- .extractGenomeDataApi(
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
