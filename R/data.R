#' Drug class abbreviations
#'
#' A dataset mapping drug classes to their abbreviations.
#'
#' @format A data frame
#' @source Internal reference data
"class_abbr"

#' Clean drug names
#'
#' A dataset mapping original drug names to cleaned/standardized names.
#'
#' @format A data frame
#' @source Internal reference data
"clean_drug"

#' Cleaned BV-BRC country names
#'
#' A dataset mapping raw country entries to cleaned and standardized names.
#'
#' @format A data frame
#' @source Internal reference data
"cleaned_bvbrc_countries"

#' Drug abbreviations
#'
#' A dataset mapping drug names to their abbreviations.
#'
#' @format A data frame
#' @source Internal reference data
"drug_abbr"

#' Drug class mappings
#'
#' A dataset mapping drugs to their drug classes.
#'
#' @format A data frame
#' @source Internal reference data
"drug_class"

#' Species recommended for exploring amRdata
#'
#' BV-BRC bacterial species with enough high-quality genomes and AMR panel
#' phenotypes to explore the package. Built by `data_raw/species_count_summary.R`
#' (AMR Panel evidence, "Good" genome quality, WGS/Complete status, CheckM
#' completeness >= 98, and 20-1,000 qualifying genomes per species).
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{species}{Species name.}
#'   \item{genome_count}{Number of qualifying genomes.}
#' }
#' @source BV-BRC <https://www.bv-brc.org>
"species_count_summary"
