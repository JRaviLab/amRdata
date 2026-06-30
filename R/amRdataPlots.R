#' Generate a summary report for AMR metadata
#'
#' @param metadata_parquet Character string. Path to a Parquet file containing
#'   standardized AMR metadata.
#' @param out_path Character string. Directory where the Markdown report is written.
#'
#' @return Writes a structured, human‑readable summary report to
#'   "<out_path>/amr_metadata_summary.md".
#'
#' @import dplyr
#' @import arrow
#' @import kableExtra
#'
#' @examples
#' generateSummary(
#'   metadata_parquet = "results/metadata.parquet",
#'   out_path = "results/"
#' )
#'
#' @export
generateSummary <- function(metadata_parquet, out_path) {
  # Little helper to apply distinct + non-empty + sorted vector
  clean_distinct <- function(df, col) {
    df |>
      dplyr::distinct({{ col }}) |>
      dplyr::filter(!is.na({{ col }}), {{ col }} != "") |>
      dplyr::arrange({{ col }}) |>
      dplyr::pull({{ col }})
  }

  # Format for Markdown
  md_tbl <- function(df) {
    knitr::kable(df, format = "pipe")
  }

  # Create a file
  write_new <- function(path, lines) {
    con <- file(path, open = "w", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    writeLines(lines, con = con, sep = "\n", useBytes = TRUE)
  }

  # Add lines to file
  append_lines <- function(path, lines) {
    con <- file(path, open = "a", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    writeLines(lines, con = con, sep = "\n", useBytes = TRUE)
  }

  # Setting paths
  if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
  metadata_parquet <- normalizePath(metadata_parquet)
  metadata <- arrow::read_parquet(metadata_parquet)
  md_path <- file.path(out_path, "amr_metadata_summary.md")

  # Validation (got any data?)
  if (nrow(metadata) == 0) {
    stop("The output table is empty. Please check your query or input data.")
  }

  # Core summaries
  TotalEntryCount <- metadata |> dplyr::count()
  CleanEntryCount <- metadata |>
    dplyr::distinct(genome.genome_id) |>
    dplyr::count()

  Antibiotics <- clean_distinct(metadata, genome_drug.antibiotic)
  AntibioticClasses <- clean_distinct(metadata, drug_class)

  LabMethods <- metadata |>
    dplyr::group_by(genome_drug.laboratory_typing_method) |>
    dplyr::count() |>
    dplyr::ungroup()

  PubMed_ids <- clean_distinct(metadata, genome_drug.pmid)

  PhenotypeCount <- metadata |>
    dplyr::group_by(genome_drug.resistant_phenotype) |>
    dplyr::count() |>
    dplyr::ungroup()

  PhenotypebyDrugCount <- metadata |>
    dplyr::group_by(genome_drug.resistant_phenotype, genome_drug.antibiotic) |>
    dplyr::count() |>
    dplyr::ungroup()

  ResPropbyDrug <- metadata |>
    dplyr::group_by(genome_drug.antibiotic) |>
    dplyr::count(genome_drug.resistant_phenotype) |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::filter(genome_drug.resistant_phenotype == "Resistant") |>
    dplyr::transmute(genome_drug.antibiotic, res_prop = round(prop, 3)) |>
    dplyr::ungroup()

  PhenotypebyDrugClassCount <- metadata |>
    dplyr::group_by(genome.genome_id, drug_class) |>
    dplyr::filter(!(any(genome_drug.resistant_phenotype == "Resistant") &
      genome_drug.resistant_phenotype == "Susceptible")) |>
    dplyr::ungroup() |>
    dplyr::group_by(genome.genome_id, drug_class) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::group_by(genome_drug.resistant_phenotype, drug_class) |>
    dplyr::count() |>
    dplyr::ungroup()

  ResPropbyDrugClass <- metadata |>
    dplyr::group_by(drug_class) |>
    dplyr::count(genome_drug.resistant_phenotype) |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::filter(genome_drug.resistant_phenotype == "Resistant") |>
    dplyr::transmute(drug_class, res_prop = round(prop, 3)) |>
    dplyr::ungroup()

  Year <- metadata |>
    dplyr::distinct(genome.collection_year) |>
    dplyr::filter(!is.na(genome.collection_year)) |>
    dplyr::pull() |>
    sort()

  YearCount <- metadata |>
    dplyr::group_by(genome.collection_year) |>
    dplyr::filter(!is.na(genome.collection_year)) |>
    dplyr::count() |>
    dplyr::ungroup()

  Country <- clean_distinct(metadata, genome.isolation_country)
  CountryCount <- metadata |>
    dplyr::group_by(genome.isolation_country) |>
    dplyr::filter(!is.na(genome.isolation_country), genome.isolation_country != "") |>
    dplyr::count() |>
    dplyr::ungroup()

  Source <- clean_distinct(metadata, genome.isolation_source)
  SourceCount <- metadata |>
    dplyr::group_by(genome.isolation_source) |>
    dplyr::filter(!is.na(genome.isolation_source), genome.isolation_source != "") |>
    dplyr::count() |>
    dplyr::ungroup()

  Host <- clean_distinct(metadata, genome.host_common_name)

  # Header
  write_new(md_path, "# AMR summary report\n")

  # Basic stats
  append_lines(
    md_path,
    c(
      sprintf("- **Entries**: %s", TotalEntryCount[[1]]),
      sprintf("- **Unique genome IDs**: %s", CleanEntryCount[[1]]),
      "",
      sprintf(
        "- **Publications** (%d): %s",
        length(PubMed_ids),
        if (length(PubMed_ids)) paste(PubMed_ids, collapse = ", ") else "None"
      ),
      ""
    )
  )

  # Lists
  append_lines(
    md_path,
    c(
      sprintf("## Antibiotics (%d)", length(Antibiotics)),
      "",
      paste(Antibiotics, collapse = ", "),
      "",
      sprintf("## Antibiotic classes (%d)", length(AntibioticClasses)),
      "",
      paste(AntibioticClasses, collapse = ", "),
      ""
    )
  )

  # Tables!
  append_lines(md_path, c("## Phenotype counts", "", md_tbl(PhenotypeCount), "", ""))
  append_lines(md_path, c("## Phenotype x antibiotic", "", md_tbl(PhenotypebyDrugCount), "", ""))
  append_lines(md_path, c("## Resistant proportion per antibiotic", "", md_tbl(ResPropbyDrug), "", ""))
  append_lines(md_path, c("## Phenotype x antibiotic class", "", md_tbl(PhenotypebyDrugClassCount), "", ""))
  append_lines(md_path, c("## Resistant proportion per antibiotic class", "", md_tbl(ResPropbyDrugClass), "", ""))
  append_lines(md_path, c("## Laboratory methods", "", md_tbl(LabMethods), "", ""))
  append_lines(md_path, c("## Year counts", "", md_tbl(YearCount), "", ""))
  append_lines(md_path, c("## Country counts", "", md_tbl(CountryCount), "", ""))
  append_lines(md_path, c("## Isolation sources", "", md_tbl(SourceCount), "", ""))

  # Hosts as a simple list
  if (length(Host)) {
    append_lines(md_path, c("## Hosts", "", paste0("- ", Host), "", ""))
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
  if (!dir.exists(out_path)) {
    dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  }

  metadata <- arrow::read_parquet(normalizePath(metadata_parquet))

  # --------- Build plots (same visuals as your generatePlots) ----------
  # 1) Phenotypes across antibiotics and time
  df_year <- metadata |>
    dplyr::filter(!is.na(genome.collection_year)) |>
    dplyr::select(
      genome.genome_id,
      drug_abbr,
      genome_drug.resistant_phenotype,
      genome.isolation_country,
      genome.collection_year
    )

  summary_year <- df_year |>
    dplyr::group_by(
      drug_abbr,
      genome_drug.resistant_phenotype,
      genome.collection_year
    ) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop")

  p1 <- ggplot2::ggplot(
    summary_year,
    ggplot2::aes(
      x = genome.collection_year,
      y = count,
      colour = genome_drug.resistant_phenotype
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
<<<<<<< Updated upstream
    ggplot2::facet_wrap(~drug_abbr, scales = "free_y") +
    ggplot2::labs(
      title = "Resistant phenotypes across antibiotics and time",
      x = "Year", y = "Number of isolates",
      colour = "Phenotype"
    ) +
    ggplot2::scale_color_brewer(palette = "Pastel1") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
      axis.text.y = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(colour = "black"),
      panel.grid.minor = element_blank()
    )

  p1
  # 2) Resistance only over time


  # 1) Levels for antibiotics (from the same subset used in p2)
  abx_levels <- summary_year |>
    dplyr::filter(genome_drug.resistant_phenotype == "Resistant") |>
    dplyr::distinct(drug_abbr) |>
    dplyr::arrange(drug_abbr) |>
    dplyr::pull(drug_abbr)
<<<<<<< Updated upstream

  # 2) Base Okabe–Ito (CVD-friendly) and pastelizer
  okabe_ito_base <- c(
    "#000000", # black
    "#E69F00", # orange
    "#56B4E9", # sky blue
    "#009E73", # bluish green
    "#F0E442", # yellow
    "#0072B2", # blue
    "#D55E00", # vermillion
    "#CC79A7" # reddish purple
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
        antibiotic_fac = factor(drug_abbr, levels = abx_levels)
      ),
    ggplot2::aes(
      x = genome.collection_year, y = count,
      colour = antibiotic_fac
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      title  = "Distribution of resistance data over time",
      x      = "Year",
      y      = "Number of resistant isolates",
      colour = "Antibiotic"
    ) +
    ggplot2::scale_color_manual(values = pal_named, drop = TRUE) + # <- named palette
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "#2D2D2D"),
      legend.position = "bottom"
    )

  # 3) Time × geography × phenotype
  df_country <- metadata |>
    dplyr::filter(genome.isolation_country != "") |>
    dplyr::select(
      genome.genome_id,
      drug_abbr,
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
    ggplot2::aes(
      x = genome.collection_year,
      y = genome.isolation_country,
      size = count,
      color = genome_drug.resistant_phenotype
    )
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
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      legend.position = "right"
    )

  # 4) Phenotype proportion per antibiotic (stacked, normalized)
  p4 <- ggplot2::ggplot(
    metadata,
    ggplot2::aes(
      x = drug_abbr,
      fill = genome_drug.resistant_phenotype
    )
  ) +
    ggplot2::geom_bar(position = "fill") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Phenotype proportion per antibiotic",
      x = "Antibiotic", y = "Proportion", fill = "Phenotype"
    ) +
    ggplot2::scale_fill_brewer(palette = "Pastel1") + # <- keep pastel
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  # 5) Treemap of isolation sources
  summary_isolation_source <- metadata |>
    dplyr::filter(genome.isolation_source != "") |>
    dplyr::group_by(genome.isolation_source, genome_drug.resistant_phenotype) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
    dplyr::filter(genome_drug.resistant_phenotype == "Resistant")

  p5 <- ggplot2::ggplot(
    summary_isolation_source,
    ggplot2::aes(area = count, fill = genome.isolation_source)
  ) +
    treemapify::geom_treemap() +
    # treemapify::geom_treemap_text(
    #   ggplot2::aes(label = genome_drug.resistant_phenotype),
    #   color = "grey15", grow = FALSE
    # ) +
    treemapify::geom_treemap_text(
      ggplot2::aes(label = genome.isolation_source),
      color = "grey15", grow = FALSE
    ) +
    ggplot2::labs(
      title = "Distribution of Resistant isolates by isolation source",
      fill = "Isolation source"
    ) +
    ggplot2::scale_fill_manual(
      values = colorRampPalette(RColorBrewer::brewer.pal(8, "Pastel2"))(n_distinct(summary_isolation_source$genome.isolation_source))
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      legend.position = "none"
    )

  # 6) Histogram of resistant classes per genome
  p6 <- ggplot2::ggplot(metadata, ggplot2::aes(num_resistant_classes)) +
    ggplot2::geom_histogram(binwidth = 1, fill = "steelblue") +
    ggplot2::labs(
      title = "Distribution of resistant classes per genome",
      x = "# Resistant Classes", y = "Count"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      legend.position = "bottom",
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(colour = "black"),
      panel.grid.minor = element_blank()
    )

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
