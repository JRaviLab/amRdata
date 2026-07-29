test_that(".parse_bvbrc_tsv() returns an empty tibble for empty input", {
  expect_identical(.parse_bvbrc_tsv(character(0)), tibble::tibble())
  expect_identical(.parse_bvbrc_tsv(c("", "   ")), tibble::tibble())
})

test_that(".parse_bvbrc_tsv() parses well-formed TSV lines into a tibble", {
  lines <- c(
    "genome_id\tgenome_name",
    "511145.12\tEscherichia coli",
    "83333.111\tEscherichia coli K-12"
  )

  out <- .parse_bvbrc_tsv(lines)

  expect_s3_class(out, "tbl_df")
  expect_identical(names(out), c("genome_id", "genome_name"))
  expect_identical(out$genome_id, c("511145.12", "83333.111"))
})

test_that(".parse_bvbrc_tsv() drops blank lines and non-TSV noise", {
  lines <- c(
    "genome_id\tgenome_name",
    "",
    "some CLI log line with no tabs",
    "511145.12\tEscherichia coli"
  )

  out <- .parse_bvbrc_tsv(lines)

  expect_equal(nrow(out), 1L)
  expect_identical(out$genome_id, "511145.12")
})
