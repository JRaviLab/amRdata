make_report_metadata <- function(n_genomes = 30, seed = 1) {
  set.seed(seed)
  drugs <- tibble::tibble(
    genome_drug.antibiotic = c("ampicillin", "ciprofloxacin", "gentamicin", "tetracycline"),
    drug_abbr = c("AMP", "CIP", "GEN", "TET"),
    drug_class = c("penicillins", "fluoroquinolones", "aminoglycosides", "tetracyclines"),
    class_abbr = c("PEN", "FLQ", "AMG", "TET")
  )
  genomes <- tibble::tibble(
    genome.genome_id = sprintf("562.%d", seq_len(n_genomes)),
    genome.species = "Escherichia coli",
    genome.collection_year = sample(c(2015, 2016, 2018, NA), n_genomes, replace = TRUE),
    genome.isolation_country = sample(c("USA", "Kenya", ""), n_genomes, replace = TRUE),
    genome.isolation_source = sample(c("urine", "blood", ""), n_genomes, replace = TRUE),
    genome.host_common_name = sample(c("Human", ""), n_genomes, replace = TRUE),
    genome_drug.pmid = sample(c("111,222", "NULL", ""), n_genomes, replace = TRUE),
    genome_drug.laboratory_typing_method = sample(c("MIC", ""), n_genomes, replace = TRUE)
  )
  out <- merge(genomes, drugs, by = NULL) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      genome_drug.resistant_phenotype = sample(
        c("Resistant", "Susceptible"), dplyr::n(), replace = TRUE, prob = c(0.4, 0.6)
      )
    )
  profile <- out |>
    dplyr::filter(genome_drug.resistant_phenotype == "Resistant") |>
    dplyr::group_by(genome.genome_id) |>
    dplyr::summarise(
      resistant_classes = paste(sort(unique(class_abbr)), collapse = "_"),
      num_resistant_classes = dplyr::n_distinct(class_abbr),
      .groups = "drop"
    )
  out |>
    dplyr::left_join(profile, by = "genome.genome_id") |>
    dplyr::mutate(
      num_resistant_classes = dplyr::coalesce(num_resistant_classes, 0L),
      resistant_classes = dplyr::coalesce(resistant_classes, "Susceptible")
    )
}

new_report_dir <- function() {
  dir <- tempfile("amr-report-")
  dir.create(dir)
  dir
}

# A minimal ORB: a folder holding "<id>_parquet.duckdb" with a `metadata` table.
make_test_orb <- function(df) {
  orb <- file.path(new_report_dir(), "orb")
  dir.create(orb)
  con <- DBI::dbConnect(duckdb::duckdb(), file.path(orb, "test_parquet.duckdb"))
  DBI::dbWriteTable(con, "metadata", as.data.frame(df))
  DBI::dbDisconnect(con, shutdown = TRUE)
  orb
}

test_that(".generateSummary writes a readable Markdown report", {
  md <- make_report_metadata()
  pq <- make_test_orb(md)
  out <- new_report_dir()

  path <- suppressMessages(.generateSummary(pq, out_path = out))
  expect_equal(basename(path), "amr_metadata_summary.md")
  txt <- paste(readLines(path), collapse = "\n")

  expect_match(txt, "# AMR summary report for \\*Escherichia coli\\*")
  expect_match(txt, "## Overview", fixed = TRUE)
  expect_match(txt, "Isolates (unique genomes)", fixed = TRUE)
  expect_match(txt, "30", fixed = TRUE)
  expect_match(txt, "## Resistance by antibiotic", fixed = TRUE)
  expect_match(txt, "## Multidrug resistance", fixed = TRUE)
  # raw column names should not leak into headers
  expect_false(grepl("genome_drug.resistant_phenotype", txt, fixed = TRUE))
  # blank / "NULL" PMIDs are dropped; real ones link to PubMed
  expect_false(grepl("NULL", txt, fixed = TRUE))
  expect_match(txt, "pubmed.ncbi.nlm.nih.gov/111/", fixed = TRUE)
  expect_match(txt, "pubmed.ncbi.nlm.nih.gov/222/", fixed = TRUE)
})

test_that(".generateSummary counts isolates, not genome x drug rows", {
  md <- make_report_metadata(n_genomes = 25)
  pq <- make_test_orb(md)
  out <- new_report_dir()
  path <- suppressMessages(.generateSummary(pq, out_path = out))
  txt <- readLines(path)
  isolates_row <- grep("Isolates (unique genomes)", txt, fixed = TRUE, value = TRUE)
  results_row <- grep("Susceptibility results", txt, fixed = TRUE, value = TRUE)
  expect_match(isolates_row, "\\|\\s*25\\|")
  expect_match(results_row, "\\|\\s*100\\|")
})

test_that(".generateSummary flags percentages based on very few results", {
  md <- make_report_metadata(n_genomes = 3)
  pq <- make_test_orb(md)
  out <- new_report_dir()
  path <- suppressMessages(.generateSummary(pq, out_path = out))
  txt <- paste(readLines(path), collapse = "\n")
  expect_match(txt, "%\\*")
  expect_match(txt, "Fewer than 10 results", fixed = TRUE)
})

test_that("reports work with only the required columns and Intermediate calls", {
  md <- tibble::tibble(
    genome.genome_id = c("1.1", "1.1", "1.2", "1.3"),
    genome_drug.antibiotic = c("ampicillin", "gentamicin", "ampicillin", "gentamicin"),
    genome_drug.resistant_phenotype = c("Resistant", "Susceptible", "Intermediate", "Susceptible")
  )
  pq <- make_test_orb(md)
  out <- new_report_dir()

  expect_no_error(suppressMessages(.generateSummary(pq, out_path = out)))
  txt <- paste(readLines(file.path(out, "amr_metadata_summary.md")), collapse = "\n")
  expect_match(txt, "Intermediate", fixed = TRUE)
  expect_match(txt, "Not recorded", fixed = TRUE)
  # class/MDR sections need cleaned columns and are skipped, not errors
  expect_false(grepl("## Multidrug resistance", txt, fixed = TRUE))

  pdf <- suppressMessages(.generatePlots(pq, out_path = out))
  expect_true(file.exists(pdf))
})

test_that(".generatePlots writes a multi-page PDF", {
  md <- make_report_metadata()
  pq <- make_test_orb(md)
  out <- new_report_dir()

  pdf <- suppressMessages(.generatePlots(pq, out_path = out))
  expect_equal(basename(pdf), "amRdata_exploratory_plots.pdf")
  expect_gt(file.size(pdf), 5000)
  pdf_lines <- readLines(pdf, warn = FALSE, encoding = "latin1")
  pages <- sum(grepl("/Type /Page /Parent", pdf_lines, fixed = TRUE, useBytes = TRUE))
  expect_gt(pages, 3)
})

test_that("isolate-level plots count genomes, not genome x drug rows", {
  md <- .prepReportMetadata(make_report_metadata(n_genomes = 20))
  genomes <- .genomeLevelMetadata(md)
  expect_equal(nrow(genomes), 20)

  mdr <- .plotMdr(genomes)
  # every one of the 20 isolates is counted once (not 20 isolates x 4 drugs)
  expect_equal(sum(mdr$data$n), 20)
})

test_that("plots needing absent data return NULL instead of failing", {
  md <- .prepReportMetadata(tibble::tibble(
    genome.genome_id = c("1.1", "1.2"),
    genome_drug.antibiotic = c("ampicillin", "ampicillin"),
    genome_drug.resistant_phenotype = c("Resistant", "Susceptible")
  ))
  genomes <- .genomeLevelMetadata(md)
  expect_null(.plotDrugYearHeatmap(md))
  expect_null(.plotResultsByYear(md))
  expect_null(.plotMdr(genomes))
  expect_null(.plotProfiles(genomes))
  expect_null(.plotCountryYear(genomes))
  expect_null(.plotIsolateBars(genomes, "country", "x", "y"))
})

test_that("year axis breaks stay inside the observed years", {
  expect_equal(.yearBreaks(c(2014, 2016, 2022)), 2014:2022)
  expect_lte(length(.yearBreaks(c(1950, 2024))), 16)
  expect_true(all(.yearBreaks(c(1950, 2024)) >= 1950))
})

test_that("empty metadata is an error", {
  empty <- tibble::tibble(
    genome.genome_id = character(),
    genome_drug.antibiotic = character(),
    genome_drug.resistant_phenotype = character()
  )
  pq <- make_test_orb(empty)
  out <- new_report_dir()
  expect_error(.generateSummary(pq, out_path = out), "empty")
  expect_error(.generatePlots(pq, out_path = out), "empty")
})

test_that("writeMetadataReports writes both files from the ORB's metadata", {
  orb <- make_test_orb(make_report_metadata())
  logged <- character()
  status <- function(...) logged <<- c(logged, paste0(...))

  res <- writeMetadataReports(orb, status = status)

  # returns the paths written, for the manifest
  expect_equal(basename(res), c("amr_metadata_summary.md", "amRdata_exploratory_plots.pdf"))
  expect_true(all(file.exists(res)))
  expect_length(logged, 2)
})

test_that("writeMetadataReports warns instead of failing the pipeline", {
  orb <- new_report_dir() # no ORB DuckDB here
  logged <- character()
  status <- function(...) logged <<- c(logged, paste0(...))

  warned <- character()
  res <- withCallingHandlers(
    writeMetadataReports(orb, status = status),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warned, 2) # one per report; neither stops the other
  expect_match(warned, "Could not write", all = TRUE)
  expect_length(res, 0) # nothing written, so nothing for the manifest
  expect_true(any(grepl("Could not write", logged)))
})

test_that(".reportFilesMessage lists the reports, and is NULL when there are none", {
  expect_null(.reportFilesMessage(character(0)))
  msg <- .reportFilesMessage(c("/orb/amr_metadata_summary.md", "/orb/amRdata_exploratory_plots.pdf"))
  expect_match(msg, "Metadata reports written to the ORB", fixed = TRUE)
  expect_match(msg, "/orb/amr_metadata_summary.md", fixed = TRUE)
  expect_match(msg, "/orb/amRdata_exploratory_plots.pdf", fixed = TRUE)
})

test_that("the summary counts MDR as three or more classes", {
  md <- make_report_metadata()
  n_mdr <- sum(!duplicated(md$genome.genome_id) & md$num_resistant_classes >= 3)
  path <- suppressMessages(.generateSummary(make_test_orb(md), out_path = new_report_dir()))
  txt <- paste(readLines(path), collapse = "\n")
  expect_match(txt, sprintf("Multidrug resistant isolates \\(3 or more classes\\) \\|\\s+%d \\(", n_mdr))
})

test_that("reports default to the ORB folder", {
  orb <- make_test_orb(make_report_metadata())
  expect_message(path <- .generateSummary(orb), "added to the ORB")
  expect_equal(normalizePath(dirname(path)), normalizePath(orb))
})

test_that("writeMetadataReports reports what it wrote when called directly", {
  orb <- make_test_orb(make_report_metadata())
  expect_message(res <- writeMetadataReports(orb), "Wrote metadata summary report")
  expect_equal(basename(res), c("amr_metadata_summary.md", "amRdata_exploratory_plots.pdf"))
})
