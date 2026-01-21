
#' Generate a summary report for AMR metadata
#'
#' @param metadata_parquet Character string. Path to a Parquet file containing
#'   standardized AMR metadata.
#'
#' @return Prints a structured, human‑readable summary report.
#'
#' @import dplyr
#' @import arrow
#' @import kableExtra
#' 
#' @examples
#' generateSummary(metadata_parquet = "results/metadata.parquet",
#' out_path = "results/")
#' 
#' @export
generateSummary <- function(metadata_parquet, 
                            out_path) {
  
  # Normalize path and load metadata table
  metadata_parquet <- normalizePath(metadata_parquet)
  metadata <- arrow::read_parquet(metadata_parquet)
  md_path  <- file.path(out_path, paste0("amr_metadata_summary.md"))
  # Error if metadata table is empty
  if (nrow(metadata) == 0) {
    stop("The output table is empty. Please check your query or input data.")
  }
  
  # Total and unique genomes
  TotalEntryCount  <- metadata |> dplyr::count()
  CleanEntryCount  <- metadata |> dplyr::distinct(genome_drug.genome_id) |> dplyr::count()
  
  # Antibiotics and classes
  Antibiotics <- metadata |> dplyr::distinct(genome_drug.antibiotic) |> dplyr::pull() |> sort()
  AntibioticClasses <- metadata |> dplyr::distinct(drug_class) |> dplyr::pull() |> sort()
  
  # Lab methods
  LabMethods <- metadata |>
    dplyr::group_by(genome_drug.laboratory_typing_method) |>
    dplyr::count()
  
  # PubMed IDs
  PubMed_ids <- metadata |>
    dplyr::distinct(genome_drug.pmid) |>
    dplyr::filter(!is.na(genome_drug.pmid), genome_drug.pmid != "") |>
    dplyr::pull()
  
  # Phenotype counts
  PhenotypeCount <- metadata |>
    dplyr::group_by(genome_drug.resistant_phenotype) |>
    dplyr::count()
  
  # Phenotype × Drug
  PhenotypebyDrugCount <- metadata |>
    dplyr::group_by(genome_drug.resistant_phenotype, genome_drug.antibiotic) |>
    dplyr::count()
  
  # Resistance proportion per drug
  ResPropbyDrug <- metadata |>
    dplyr::group_by(genome_drug.antibiotic) |>
    dplyr::count(genome_drug.resistant_phenotype) |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::filter(genome_drug.resistant_phenotype == "Resistant") |>
    dplyr::transmute(res_prop = round(prop, 3))
  
  # Phenotype × Antibiotic Class (with consistency filtering)
  PhenotypebyDrugClassCount <- metadata |>
    dplyr::group_by(genome_drug.genome_id, drug_class) |>
    dplyr::filter(!(any(genome_drug.resistant_phenotype == "Resistant") &
                      genome_drug.resistant_phenotype == "Susceptible")) |>
    dplyr::ungroup() |>
    dplyr::group_by(genome_drug.genome_id, drug_class) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::group_by(genome_drug.resistant_phenotype, drug_class) |>
    dplyr::count()
  
  # Resistance proportion × drug class
  ResPropbyDrugClass <- metadata |>
    dplyr::group_by(drug_class) |>
    dplyr::count(genome_drug.resistant_phenotype) |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::filter(genome_drug.resistant_phenotype == "Resistant") |>
    dplyr::transmute(res_prop = round(prop, 3))
  
  # Year summaries
  Year <- metadata |>
    dplyr::distinct(genome.collection_year) |>
    dplyr::filter(!is.na(genome.collection_year)) |>
    dplyr::pull() |> sort()
  
  YearCount <- metadata |>
    dplyr::group_by(genome.collection_year) |>
    dplyr::filter(!is.na(genome.collection_year)) |>
    dplyr::count()
  
  # Country summaries
  Country <- metadata |>
    dplyr::distinct(genome.isolation_country) |>
    dplyr::filter(!is.na(genome.isolation_country), genome.isolation_country != "") |>
    dplyr::pull() |> sort()
  
  CountryCount <- metadata |>
    dplyr::group_by(genome.isolation_country) |>
    dplyr::filter(!is.na(genome.isolation_country), genome.isolation_country != "") |>
    dplyr::count()
  
  # Isolation sources
  Source <- metadata |>
    dplyr::distinct(genome.isolation_source) |>
    dplyr::filter(!is.na(genome.isolation_source), genome.isolation_source != "") |>
    dplyr::pull() |> sort()
  
  SourceCount <- metadata |>
    dplyr::group_by(genome.isolation_source) |>
    dplyr::filter(!is.na(genome.isolation_source), genome.isolation_source != "") |>
    dplyr::count()
  
  # Host names
  Host <- metadata |>
    dplyr::distinct(genome.host_common_name) |>
    dplyr::filter(!is.na(genome.host_common_name), genome.host_common_name != "") |>
    dplyr::pull() |> sort()
  
  # helper to write markdown tables  
  md_tbl <- function(df) {    
    knitr::kable(df, format = "pipe")  
  }
  
  # --------- PRINT SUMMARY ---------
  # ---- Write Markdown (clean tables) ----  
  cat("# AMR Summary Report\n\n", file = md_path)  
  cat(sprintf("- **Entries**: %s\n- **Unique genome IDs**: %s\n\n",              
              TotalEntryCount[[1]], CleanEntryCount[[1]]), 
      file = md_path, append = TRUE)  
  cat(sprintf("- **Publications** (%d): %s\n\n",              
              length(PubMed_ids),              
              if (length(PubMed_ids)) paste(PubMed_ids, collapse = ", ") 
              else "None"),      
      file = md_path, append = TRUE)  
  cat(sprintf("## Antibiotics (%d)\n\n%s\n\n",              
              length(Antibiotics), paste(Antibiotics, collapse = ", ")),      
      file = md_path, append = TRUE)  
  cat(sprintf("## Antibiotic Classes (%d)\n\n%s\n\n",              
              length(AntibioticClasses), 
              paste(AntibioticClasses, collapse = ", ")),      
      file = md_path, append = TRUE)  
  cat("## Phenotype Counts\n\n", 
      md_tbl(PhenotypeCount), "\n\n", 
      file = md_path, append = TRUE)  
  cat("## Phenotype × Antibiotic\n\n", 
      md_tbl(PhenotypebyDrugCount), "\n\n", 
      file = md_path, append = TRUE)  
  cat("## Resistant Proportion per Antibiotic\n\n", 
      md_tbl(ResPropbyDrug), "\n\n", 
      file = md_path, append = TRUE)  
  cat("## Phenotype × Antibiotic Class\n\n", 
      md_tbl(PhenotypebyDrugClassCount), "\n\n", 
      file = md_path, append = TRUE)  
  cat("## Resistant Proportion per Antibiotic Class\n\n", 
      md_tbl(ResPropbyDrugClass), "\n\n", 
      file = md_path, append = TRUE)  
  cat("## Laboratory Methods\n\n", 
      md_tbl(LabMethods), "\n\n", 
      file = md_path, append = TRUE)  
  cat("## Year Counts\n\n", 
      md_tbl(YearCount), "\n\n", 
      file = md_path, append = TRUE)  
  cat("## Country Counts\n\n", 
      md_tbl(CountryCount), "\n\n", 
      file = md_path, append = TRUE)  
  cat("## Isolation Sources\n\n", 
      md_tbl(SourceCount), "\n\n", 
      file = md_path, append = TRUE)  
  if (length(Host)) {    
    cat("## Hosts\n\n", 
        paste("- ", Host, collapse = "\n"), "\n", 
        file = md_path, append = TRUE)  
  }  
  }


#' Write all summary plots to file(s)
#'
#' @param metadata_parquet Character. Path to the Parquet metadata file.
#' @param out_path Character. Output directory for plot files.
#' 
#' @return Invisibly returns a vector of pdf file paths written.
#' @export
generatePlots <- function(metadata_parquet,
                          out_path) {
  
  device <- "pdf"
  if(!dir.exists(out_path)) {dir.create(out_path, showWarnings = FALSE, recursive = TRUE)}
  
  metadata <- arrow::read_parquet(normalizePath(metadata_parquet))
  
  # --------- Build plots (same visuals as your generatePlots) ----------
  # 1) Phenotypes across antibiotics and time
  df_year <- metadata |>
    dplyr::filter(!is.na(genome.collection_year)) |>
    dplyr::select(
      genome_drug.genome_id,
      genome_drug.antibiotic,
      genome_drug.resistant_phenotype,
      genome.isolation_country,
      genome.collection_year
    )
  
  summary_year <- df_year |>
    dplyr::group_by(
      genome_drug.antibiotic,
      genome_drug.resistant_phenotype,
      genome.collection_year
    ) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop")
  
  p1 <- ggplot2::ggplot(
    summary_year,
    ggplot2::aes(x = genome.collection_year,
                 y = count,
                 colour = genome_drug.resistant_phenotype)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ genome_drug.antibiotic, scales = "free_y") +
    ggplot2::labs(
      title = "Resistant phenotypes across antibiotics and time",
      x = "Year", y = "Number of isolates",
      colour = "Phenotype"
    ) +
    ggplot2::scale_color_brewer(palette = "Pastel1") +
    ggplot2::theme_minimal(base_size = 12) +                                 
    ggplot2::theme(text = ggplot2::element_text(colour = "#2D2D2D"),                               
                   legend.position = "bottom",                 
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  # 2) Resistance only over time
  
  
  # 1) Levels for antibiotics (from the same subset used in p2)
  abx_levels <- summary_year |>
    dplyr::filter(genome_drug.resistant_phenotype == "Resistant") |>
    dplyr::distinct(genome_drug.antibiotic) |>
    dplyr::arrange(genome_drug.antibiotic) |>
    dplyr::pull(genome_drug.antibiotic)
  
  # 2) Base Okabe–Ito (CVD-friendly) and pastelizer
  okabe_ito_base <- c(
    "#000000", # black
    "#E69F00", # orange
    "#56B4E9", # sky blue
    "#009E73", # bluish green
    "#F0E442", # yellow
    "#0072B2", # blue
    "#D55E00", # vermillion
    "#CC79A7"  # reddish purple
  )
  
  okabe_ito_pastel <- function(n, lighten = 0.15) {
    # Interpolate if more than 8 needed
    cols <- if (n <= length(okabe_ito_base)) {
      okabe_ito_base[seq_len(n)]
    } else {
      grDevices::colorRampPalette(okabe_ito_base)(n)
    }
    to_rgb <- function(hex) grDevices::col2rgb(hex) / 255
    blend_with_white <- function(hex, a = lighten) {
      rgb <- to_rgb(hex)
      out <- (1 - a) * rgb + a * c(1, 1, 1)
      grDevices::rgb(out[1], out[2], out[3])
    }
    vapply(cols, blend_with_white, character(1), a = lighten)
  }
  
  # 3) Build a NAMED palette aligned to factor levels
  pal_vals <- okabe_ito_pastel(length(abx_levels), lighten = 0.15)
  pal_named <- stats::setNames(pal_vals, abx_levels)
  
  # 4) Plot with factor levels + named palette (no warnings, distinct colors)
  p2 <- ggplot2::ggplot(
    summary_year |>
      dplyr::filter(genome_drug.resistant_phenotype == "Resistant") |>
      dplyr::mutate(
        antibiotic_fac = factor(genome_drug.antibiotic, levels = abx_levels)
      ),
    ggplot2::aes(x = genome.collection_year, y = count,
                 colour = antibiotic_fac)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      title  = "Distribution of resistance data over time",
      x      = "Year",
      y      = "Number of resistant isolates",
      colour = "Antibiotic"
    ) +
    ggplot2::scale_color_manual(values = pal_named, drop = TRUE) +  # <- named palette
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "#2D2D2D"),
      legend.position = "bottom"
    )
  
  # 3) Time × geography × phenotype
  df_country <- metadata |>
    dplyr::filter(genome.isolation_country != "") |>
    dplyr::select(
      genome_drug.genome_id,
      genome_drug.antibiotic,
      genome_drug.resistant_phenotype,
      genome.isolation_country,
      genome.collection_year
    )
  
  summary_country_year <- df_country |>
    dplyr::group_by(
      genome.collection_year,
      genome_drug.resistant_phenotype,
      genome.isolation_country
    ) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop")
  
  p3 <- ggplot2::ggplot(
    summary_country_year,
    ggplot2::aes(x = genome.collection_year,
                 y = genome.isolation_country,
                 size = count,
                 color = genome_drug.resistant_phenotype)
  ) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::scale_size(range = c(3, 15)) +
    ggplot2::scale_color_viridis_d(option = "C", begin = 0.25, end = 0.95) +
    ggplot2::labs(
      title = "AMR isolates across time and geography",
      x = "Year", y = "Country",
      size = "Count", color = "Phenotype"
    ) + 
    ggplot2::theme_minimal(base_size = 12) +                                  
    ggplot2::theme(text = ggplot2::element_text(colour = "#2D2D2D"),                                
                   legend.position = "bottom")
  
  # 4) Phenotype proportion per antibiotic (stacked, normalized)
  p4 <- ggplot2::ggplot(
    metadata,
    ggplot2::aes(x = genome_drug.antibiotic,
                 fill = genome_drug.resistant_phenotype)
  ) +
    ggplot2::geom_bar(position = "fill") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Phenotype proportion per antibiotic",
      x = "Antibiotic", y = "Proportion", fill = "Phenotype"
    ) + ggplot2::scale_fill_brewer(palette = "Pastel1") +                               # <- keep pastel  
    ggplot2::theme_minimal(base_size = 12) +                                  
    ggplot2::theme(text = ggplot2::element_text(colour = "#2D2D2D"),                                
                   legend.position = "bottom",                 
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  # 5) Treemap of isolation sources
  summary_isolation_source <- metadata |>
    dplyr::filter(genome.isolation_source != "") |>
    dplyr::group_by(genome.isolation_source, genome_drug.resistant_phenotype) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop")
  
  p5 <- ggplot2::ggplot(
    summary_isolation_source,
    ggplot2::aes(area = count, fill = genome.isolation_source)
  ) +
    treemapify::geom_treemap() +
    treemapify::geom_treemap_text(
      ggplot2::aes(label = genome_drug.resistant_phenotype),
      color = "grey15", grow = FALSE
    ) +
    ggplot2::labs(
      title = "Distribution of AMR isolates by isolation source",
      fill = "Isolation source"
    ) + ggplot2::scale_fill_brewer(palette = "Pastel2") +                              
    ggplot2::theme_minimal(base_size = 12) +                                  
    ggplot2::theme(text = ggplot2::element_text(colour = "#2D2D2D"),                                
                   legend.position = "bottom",                 
                   plot.title = ggplot2::element_text(face = "bold"))
  
  # 6) Histogram of resistant classes per genome
  p6 <- ggplot2::ggplot(metadata, ggplot2::aes(num_resistant_classes)) +
    ggplot2::geom_histogram(binwidth = 1, fill = "steelblue") +
    ggplot2::labs(
      title = "Distribution of resistant classes per genome",
      x = "# Resistant Classes", y = "Count"
    ) + 
    ggplot2::theme_minimal(base_size = 12) +                                  
    ggplot2::theme(text = ggplot2::element_text(colour = "#2D2D2D"),                                
                   legend.position = "bottom")
  
  plots <- list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6)
  
  # --------- Write to device ----------
  paths <- character(0)
  
  pdf_path <- file.path(out_path, paste0("amRdata_exploratory_plots.pdf"))
  grDevices::pdf(pdf_path, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (nm in names(plots)) {
    print(plots[[nm]])
  }
  paths <- c(paths, pdf_path)
  
  
  invisible(paths)
}

