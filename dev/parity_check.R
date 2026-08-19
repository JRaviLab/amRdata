# Parity check: docker vs api on the SAME genomes.
# Prereqs:
#   1. Docker running + `docker pull danylmb/bvbrc:5.3`
#   2. Install the package first (the Docker path uses furrr workers that need
#      it INSTALLED, not just load_all'd):  R CMD INSTALL .
# Then run from the repo root:  Rscript dev/parity_check.R
suppressPackageStartupMessages(library(amRdata))

species <- "Morganella morganii"

# Fix the genome set so BOTH methods pull the same genomes (isolates the
# download comparison from any ID-resolution differences).
ids <- amRdata:::.resolveGenomeIDs_api(
  base_dir = tempfile(), user_bacs = species, overwrite = TRUE, verbose = FALSE
)
ids <- utils::head(ids, 30)
gf <- tempfile(fileext = ".txt")
writeLines(ids, gf)

run <- function(method) {
  td <- file.path(tempdir(), paste0("parity_", method))
  unlink(td, recursive = TRUE)
  dir.create(td, recursive = TRUE)
  invisible(retrieveMetadata(
    user_bacs = species, genome_id_file = gf, method = method,
    base_dir = td, overwrite = TRUE, verbose = FALSE
  ))
  db <- list.files(td, pattern = "[.]duckdb$", recursive = TRUE, full.names = TRUE)[1]
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db, read_only = TRUE)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  list(
    genome = DBI::dbReadTable(con, "genome_data"),
    amr = DBI::dbReadTable(con, "amr_phenotype")
  )
}

d <- run("docker")
a <- run("api")

cat("\n================ PARITY: docker vs api ================\n")
cat(sprintf("genome_data rows   docker=%d  api=%d\n", nrow(d$genome), nrow(a$genome)))
cat(sprintf("amr_phenotype rows docker=%d  api=%d\n", nrow(d$amr), nrow(a$amr)))
cat("same genome set:",
    setequal(d$genome[["genome.genome_id"]], a$genome[["genome.genome_id"]]), "\n")

# AMR: compare the (genome, antibiotic, phenotype) tuples as sets
tup <- function(x) {
  sort(paste(
    x[["genome_drug.genome_id"]],
    x[["genome_drug.antibiotic"]],
    x[["genome_drug.resistant_phenotype"]]
  ))
}
cat("AMR (genome,antibiotic,phenotype) identical:", identical(tup(d$amr), tup(a$amr)), "\n")

# Genome metadata: compare core fields for shared IDs
core <- c("genome.genome_id", "genome.species", "genome.genome_quality",
          "genome.genome_status", "genome.checkm_completeness",
          "genome.checkm_contamination", "genome.cds", "genome.genome_length",
          "genome.gc_content")
ord <- function(x) {
  cols <- intersect(core, names(x))
  x <- x[order(x[["genome.genome_id"]]), cols, drop = FALSE]
  x[] <- lapply(x, as.character)
  as.data.frame(x, stringsAsFactors = FALSE)
}
cat("core genome fields identical:", identical(ord(d$genome), ord(a$genome)), "\n")

cat("\ncolumns only in DOCKER genome_data:",
    paste(setdiff(names(d$genome), names(a$genome)), collapse = ", "), "\n")
cat("columns only in API genome_data:",
    paste(setdiff(names(a$genome), names(d$genome)), collapse = ", "), "\n")
cat("\nExpected known diffs: api fills measurement_unit / computational_method /",
    "source as NA in AMR (absent from genome_amr). Focus on the 'identical' lines.\n")
