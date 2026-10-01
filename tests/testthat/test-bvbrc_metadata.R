# Tests for the BV-BRC metadata cleaning helpers (R/bvbrc_metadata.R).

test_that(".cleanBVBRCmetadata maps placeholders to NA in character columns", {
  df <- tibble::tibble(a = c(" x  y ", "unknown", "NA", NA), n = 1:4)
  out <- .cleanBVBRCmetadata(df)

  expect_identical(out$a, c("x y", NA, NA, NA))
  expect_identical(out$n, 1:4)
})

test_that(".cleanBVBRCmetadata flattens list columns without literal \"NA\" strings", {
  df <- tibble::tibble(l = list(c(" p ", "N/A"), NULL, "unknown", c("q", "r")))
  out <- .cleanBVBRCmetadata(df)

  expect_type(out$l, "character")
  expect_identical(out$l, c("p", NA, NA, "q;r"))
})
