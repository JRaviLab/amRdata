# Diagnose the AMR tuple mismatch between cli and api.
# Prereqs: same as dev/parity_check.R (Docker running + package installed).
# Run:  Rscript dev/parity_diff.R
suppressPackageStartupMessages(library(amRdata))

species <- "Morganella morganii"
ids <- amRdata:::.resolveGenomeIDsApi(
  base_dir = tempfile(), user_bacs = species, overwrite = TRUE, verbose = FALSE
)
ids <- utils::head(ids, 30)
gf <- tempfile(fileext = ".txt")
writeLines(ids, gf)

amr <- function(metadata_method) {
  td <- file.path(tempdir(), paste0("pd_", metadata_method))
  unlink(td, recursive = TRUE)
  dir.create(td, recursive = TRUE)
  invisible(retrieveMetadata(
    user_bacs = species, genome_id_file = gf, metadata_method = metadata_method,
    base_dir = td, overwrite = TRUE, verbose = FALSE
  ))
  db <- list.files(td, pattern = "[.]duckdb$", recursive = TRUE, full.names = TRUE)[1]
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db, read_only = TRUE)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  DBI::dbReadTable(con, "amr_phenotype")
}

d <- amr("cli")
a <- amr("api")

tup <- function(x) {
  paste(x[["genome_drug.genome_id"]], x[["genome_drug.antibiotic"]],
        x[["genome_drug.resistant_phenotype"]], sep = " | ")
}
dt <- tup(d)
at <- tup(a)

cat("\n--- tuples in DOCKER but not API (", length(setdiff(dt, at)), ") ---\n", sep = "")
print(utils::head(sort(setdiff(dt, at)), 12))
cat("\n--- tuples in API but not DOCKER (", length(setdiff(at, dt)), ") ---\n", sep = "")
print(utils::head(sort(setdiff(at, dt)), 12))

show <- function(x, col) paste(sort(unique(x[[col]])), collapse = " | ")
cat("\ndistinct ANTIBIOTIC — docker:\n ", show(d, "genome_drug.antibiotic"), "\n")
cat("distinct ANTIBIOTIC — api:\n ", show(a, "genome_drug.antibiotic"), "\n")
cat("\ndistinct PHENOTYPE — docker:", show(d, "genome_drug.resistant_phenotype"), "\n")
cat("distinct PHENOTYPE — api:   ", show(a, "genome_drug.resistant_phenotype"), "\n")

# are duplicates the cause? (same count but different multiplicity)
cat("\ndup tuples docker:", sum(duplicated(dt)), " | dup tuples api:", sum(duplicated(at)), "\n")
