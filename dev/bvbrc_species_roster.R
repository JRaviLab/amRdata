#!/usr/bin/env Rscript
# Per-species roster for the "raise the 25k row cap" request to the bvbrc
# Python package. For each species, count (all via header-only requests):
#   total_genomes   : genome collection, eq(species,X)
#   clean_genomes   : + genome_quality = Good   (BV-BRC QC; uses CheckM)
#   clean_with_amr  : + antimicrobial_resistance present
#   amr_rows        : genome_amr rows, eq(genome_name,<species>)  [text match;
#                     aggregates all strains -- genome_amr.taxon_id is unreliable
#                     (mixes species- and strain-level taxa), so name match is
#                     the robust per-species key]
# and flag whether any of those exceed the 25,000 per-request ceiling.
#
# See docs/bvbrc-api-feasibility.md. Requires dev/bvbrc_api_prototype.R.

source("dev/bvbrc_api_prototype.R")

CAP <- 25000L

bvbrc_species_row <- function(species) {
  s <- .enc(species)
  total     <- bvbrc_count("genome", sprintf("eq(species,%s)", s))
  clean     <- bvbrc_count("genome",
                 sprintf("and(eq(species,%s),eq(genome_quality,Good))", s))
  clean_amr <- bvbrc_count("genome",
                 sprintf("and(eq(species,%s),eq(genome_quality,Good),eq(antimicrobial_resistance,*))", s))
  amr_rows  <- bvbrc_count("genome_amr", sprintf("eq(genome_name,%s)", s))
  data.frame(species = species, total_genomes = total,
             clean_genomes = clean, clean_with_amr = clean_amr,
             amr_rows = amr_rows,
             genomes_over_cap = total > CAP, amr_rows_over_cap = amr_rows > CAP)
}

bvbrc_roster <- function(species_vec) {
  rows <- lapply(species_vec, function(sp) {
    message("  ", sp); tryCatch(bvbrc_species_row(sp),
      error = function(e) { message("    ! ", conditionMessage(e)); NULL })
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

# WHO Bacterial Priority Pathogens List (2024) + CDC AR Threats (2019), bacterial
WHO_CDC_PRIORITY <- c(
  "Acinetobacter baumannii", "Pseudomonas aeruginosa", "Escherichia coli",
  "Klebsiella pneumoniae", "Enterobacter cloacae", "Staphylococcus aureus",
  "Enterococcus faecium", "Enterococcus faecalis", "Streptococcus pneumoniae",
  "Salmonella enterica", "Shigella flexneri", "Shigella sonnei",
  "Neisseria gonorrhoeae", "Neisseria meningitidis", "Campylobacter jejuni",
  "Campylobacter coli", "Haemophilus influenzae", "Helicobacter pylori",
  "Mycobacterium tuberculosis", "Clostridioides difficile",
  "Streptococcus pyogenes", "Streptococcus agalactiae", "Serratia marcescens",
  "Proteus mirabilis", "Morganella morganii"
)

if (!interactive()) {
  r <- bvbrc_roster(WHO_CDC_PRIORITY)
  r <- r[order(-r$amr_rows), ]
  out <- "dev/bvbrc_species_roster.csv"
  write.csv(r, out, row.names = FALSE)
  cat("\n=== Per-species roster (WHO/CDC priority) ===\n")
  print(r, row.names = FALSE)
  cat(sprintf("\nSpecies with >25k genomes:     %d / %d\n",
              sum(r$genomes_over_cap), nrow(r)))
  cat(sprintf("Species with >25k genome_amr rows: %d / %d\n",
              sum(r$amr_rows_over_cap), nrow(r)))
  cat("Saved:", out, "\n")
}
