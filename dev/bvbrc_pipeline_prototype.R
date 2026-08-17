#!/usr/bin/env Rscript
# BV-BRC API pipeline prototype (Phase 1)
# ---------------------------------------
# Standalone demo that the direct-API path is a viable, RESILIENT replacement
# for the Docker/p3-* metadata+AMR download in amRdata (issue #30: stochastic
# BV-BRC downloads). NOT wired into the package.
#
# For a species it:
#   1. pulls Good-quality genome metadata (chosen columns) -- key = genome_id
#   2. pulls AMR phenotype rows (paginated past 25k)        -- key = id
#   3. joins on genome_id and writes parquet
# All requests inherit retry-with-backoff on 429/5xx (dev/bvbrc_api_prototype.R),
# so a transient 503 is retried, never parsed as data -> no "more columns" crash.
#
# Requires: dev/bvbrc_api_prototype.R + arrow

source("dev/bvbrc_api_prototype.R")

# minimal genome column set (see docs §9); drop the rest
GENOME_COLS <- paste(
  "genome_id", "assembly_accession", "genbank_accessions", "genome_quality",
  "genome_status", "checkm_completeness", "checkm_contamination", "cds",
  "genome_length", "gc_content", "host_name", "isolation_country",
  "geographic_group", "species", "taxon_id", sep = ","
)
AMR_COLS <- paste(
  "genome_id", "genome_name", "antibiotic", "resistant_phenotype",
  "measurement", "measurement_unit", "laboratory_typing_method", "evidence",
  sep = ",")

# Pull Good-quality genome metadata for a species (sequential; genome_id key).
pull_genome_metadata <- function(species, cols = GENOME_COLS, verbose = TRUE) {
  filt <- sprintf("and(eq(species,%s),eq(genome_quality,Good))", .enc(species))
  if (verbose) message("  genome metadata: ", bvbrc_count("genome", filt), " Good genomes")
  bvbrc_fetch_all("genome", filt, cols, key = "genome_id")
}

# Pull AMR phenotype rows for a species (parallel; UUID id key).
pull_amr <- function(species, cols = AMR_COLS, parallel = TRUE, verbose = TRUE) {
  filt <- sprintf("eq(genome_name,%s)", .enc(species))
  if (verbose) message("  AMR rows: ", bvbrc_count("genome_amr", filt))
  bvbrc_fetch_all("genome_amr", filt, cols, parallel = parallel, key = "id")
}

# End-to-end: pull both, join, write parquet, return a small summary.
run_pipeline <- function(species, out_dir = "dev/out", parallel_amr = TRUE) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  message("Species: ", species)

  t_meta <- system.time(meta <- pull_genome_metadata(species))
  t_amr  <- system.time(amr  <- pull_amr(species, parallel = parallel_amr))

  # keep AMR only for the Good genomes we pulled, then join
  amr_good <- amr[amr$genome_id %in% meta$genome_id, ]
  joined <- merge(meta, amr_good, by = "genome_id", suffixes = c("", ".amr"))

  slug <- gsub("[^A-Za-z0-9]+", "_", species)
  arrow::write_parquet(meta,   file.path(out_dir, paste0(slug, "_genome.parquet")))
  arrow::write_parquet(amr,    file.path(out_dir, paste0(slug, "_amr.parquet")))
  arrow::write_parquet(joined, file.path(out_dir, paste0(slug, "_joined.parquet")))

  summary <- data.frame(
    species              = species,
    good_genomes         = nrow(meta),
    amr_rows             = nrow(amr),
    good_genomes_w_amr   = length(unique(amr_good$genome_id)),
    joined_rows          = nrow(joined),
    meta_secs            = round(t_meta["elapsed"], 1),
    amr_secs             = round(t_amr["elapsed"], 1),
    row.names = NULL
  )
  message("  wrote parquet to ", out_dir, "/")
  summary
}

if (!interactive()) {
  # Enterococcus faecium: ~7.7k Good genomes (1 page) + ~69k AMR rows (keyset,
  # multi-page) -> exercises pagination + parallel without a huge runtime.
  s <- run_pipeline("Enterococcus faecium")
  cat("\n=== pipeline summary ===\n"); print(s, row.names = FALSE)
}
