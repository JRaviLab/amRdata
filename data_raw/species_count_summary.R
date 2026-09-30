# Build data/species_count_summary.rda: BV-BRC bacterial species that qualify
# for exploring amRdata. Run by hand from the package root when the list needs
# refreshing:
#   source("data_raw/species_count_summary.R")
#
# Downloads the full BV-BRC bacterial metadata table (slow) on first run and caches
# it in the BiocFileCache; set refresh <- TRUE to re-download.

devtools::load_all() # internal helpers live in R/bvbrc_metadata.R

refresh <- FALSE

# Selection cutoffs
min_genomes <- 20 # inclusive
max_genomes <- 1000 # inclusive
min_checkm <- 98 # inclusive

bvbrc_clean <- .loadBVBRCmetadata(refresh = refresh)

species_count_summary <- bvbrc_clean |>
  dplyr::distinct(
    species,
    antimicrobial_resistance,
    genome_id,
    antimicrobial_resistance_evidence,
    genome_quality,
    genome_status,
    checkm_completeness
  ) |>
  dplyr::filter(
    antimicrobial_resistance_evidence == "AMR Panel",
    genome_quality == "Good",
    genome_status %in% c("WGS", "Complete"),
    as.numeric(checkm_completeness) >= min_checkm,
    !is.na(species),
    !is.na(antimicrobial_resistance)
  ) |>
  dplyr::distinct(species, genome_id) |>
  dplyr::count(species, name = "genome_count") |>
  dplyr::filter(
    genome_count >= min_genomes,
    genome_count <= max_genomes
  ) |>
  dplyr::arrange(genome_count) |>
  as.data.frame()

usethis::use_data(species_count_summary, overwrite = TRUE)
