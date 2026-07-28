#' Generate a summary report for AMR metadata
#'
#' @param metadata_parquet Character string. Path to a Parquet file containing
#'   standardized AMR metadata.
#' @param out_path Character string. Directory where the Markdown report is written.
#'
#' @return Writes a structured, human‑readable summary report to
#'   "<out_path>/amr_metadata_summary.md".
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

  # get Species name
  Species_name <- metadata |>
    dplyr::distinct(genome.species) |>
    dplyr::pull()

  # Core summaries
  TotalEntryCount <- metadata |> dplyr::count()
  CleanEntryCount <- metadata |>
    dplyr::distinct(genome.genome_id) |>
    dplyr::count()

  Antibiotics <- clean_distinct(metadata, genome_drug.antibiotic)
  AntibioticClasses <- clean_distinct(metadata, drug_class)

  LabMethods <- metadata |>
    dplyr::mutate(
      genome_drug.laboratory_typing_method = dplyr::case_when(
        is.na(genome_drug.laboratory_typing_method) ~ "Not defined",
        genome_drug.laboratory_typing_method == "" ~ "Not defined",
        TRUE ~ genome_drug.laboratory_typing_method
      )
    ) |>
    dplyr::count(genome_drug.laboratory_typing_method)

  PubMed_ids <- clean_distinct(metadata, genome_drug.pmid)

  PhenotypeCount <- metadata |>
    dplyr::group_by(genome_drug.resistant_phenotype) |>
    dplyr::count() |>
    dplyr::arrange(-n) |>
    dplyr::ungroup()

  ## Stats by drugs
  PhenotypebyDrugCount <- metadata |>
    dplyr::group_by(genome_drug.resistant_phenotype, genome_drug.antibiotic) |>
    dplyr::count() |>
    dplyr::arrange(-n) |>
    dplyr::ungroup()

  ResPropbyDrug <- metadata |>
    dplyr::group_by(genome_drug.antibiotic) |>
    dplyr::count(genome_drug.resistant_phenotype) |>
    dplyr::mutate(prop = n / sum(n)) |>
    dplyr::filter(genome_drug.resistant_phenotype == "Resistant") |>
    dplyr::transmute(genome_drug.antibiotic, res_prop = round(prop, 3)) |>
    dplyr::arrange(-res_prop) |>
    dplyr::ungroup()

 # Collapse to one phenotype per genome x drug_class
drug_class_calls <- metadata |>
  dplyr::group_by(genome.genome_id, drug_class) |>
  dplyr::summarise(
    genome_drug.resistant_phenotype = dplyr::case_when(
      any(genome_drug.resistant_phenotype == "Resistant", na.rm = TRUE) ~ "Resistant",
      any(genome_drug.resistant_phenotype == "Intermediate", na.rm = TRUE) ~ "Intermediate",
      any(genome_drug.resistant_phenotype == "Susceptible", na.rm = TRUE) ~ "Susceptible",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  )

# Stats by drug class
PhenotypebyDrugClassCount <- drug_class_calls |>
  dplyr::count(genome_drug.resistant_phenotype, drug_class) |>
  dplyr::arrange(dplyr::desc(n))

ResPropbyDrugClass <- drug_class_calls |>
  dplyr::group_by(drug_class) |>
  dplyr::count(genome_drug.resistant_phenotype) |>
  dplyr::mutate(prop = n / sum(n)) |>
  dplyr::filter(genome_drug.resistant_phenotype == "Resistant") |>
  dplyr::transmute(drug_class, res_prop = round(prop, 3)) |>
  dplyr::arrange(dplyr::desc(res_prop)) |>
  dplyr::ungroup()
  
  ## Collection years
  Year <- metadata |>
    dplyr::distinct(genome.collection_year) |>
    dplyr::filter(!is.na(genome.collection_year)) |>
    dplyr::pull() |>
    sort()

  YearCount <- metadata |>
    dplyr::group_by(genome.collection_year) |>
    dplyr::filter(!is.na(genome.collection_year)) |>
    dplyr::count() |>
    dplyr::arrange(-n) |>
    dplyr::ungroup()

  ## Geographical regions (countries)
  Country <- clean_distinct(metadata, genome.isolation_country)
  CountryCount <- metadata |>
    dplyr::group_by(genome.isolation_country) |>
    dplyr::filter(!is.na(genome.isolation_country), genome.isolation_country != "") |>
    dplyr::count() |>
    dplyr::arrange(-n) |>
    dplyr::ungroup()

  ## Isolation sources
  Source <- clean_distinct(metadata, genome.isolation_source)
  SourceCount <- metadata |>
    dplyr::group_by(genome.isolation_source) |>
    dplyr::filter(!is.na(genome.isolation_source), genome.isolation_source != "") |>
    dplyr::count() |>
    dplyr::arrange(-n) |>
    dplyr::ungroup()

  Host <- clean_distinct(metadata, genome.host_common_name)
  HostCount <- metadata |>
    dplyr::group_by(genome.host_common_name) |>
    dplyr::filter(!is.na(genome.host_common_name), genome.host_common_name != "") |>
    dplyr::count() |>
    dplyr::arrange(-n) |>
    dplyr::ungroup()

  # Header
  write_new(
    md_path,
    sprintf("# AMR summary report for *%s*", Species_name)
  )

  # Basic stats
  append_lines(
    md_path,
    c(
      sprintf("- **Total no. of observations**: %s", TotalEntryCount[[1]]),
      sprintf("- **Unique genome IDs**: %s", CleanEntryCount[[1]]),
      "",
      sprintf(
        "- **Associated publications** (%d): %s",
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
  append_lines(md_path, c("## Phenotypes x antibiotic(s)", "", md_tbl(PhenotypebyDrugCount), "", ""))
  append_lines(md_path, c("## Resistant proportions per antibiotic", "", md_tbl(ResPropbyDrug), "", ""))
  append_lines(md_path, c("## Phenotypes x antibiotic class(es)", "", md_tbl(PhenotypebyDrugClassCount), "", ""))
  append_lines(md_path, c("## Resistant proportion per antibiotic class", "", md_tbl(ResPropbyDrugClass), "", ""))
  append_lines(md_path, c("## Laboratory methods", "", md_tbl(LabMethods), "", ""))
  append_lines(md_path, c("## Collection years", "", md_tbl(YearCount), "", ""))
  append_lines(md_path, c("## Isolation countries", "", md_tbl(CountryCount), "", ""))
  append_lines(md_path, c("## Isolation sources", "", md_tbl(SourceCount), "", ""))
  append_lines(md_path, c("## Hosts", "", md_tbl(HostCount), "", ""))


  # Hosts as a simple list
  # if (length(Host)) {
  #   append_lines(md_path, c("## Hosts", "", paste0("- ", Host), "", ""))
  # }
}


#' Write all summary plots to file(s)
#'
#' Expects `metadata_parquet` to be the output of `runDataProcessing()`'s
#' `cleanData()` step (or an export of the resulting `metadata` table), since
#' `drug_abbr`, `drug_class`, and `num_resistant_classes` are only populated
#' after that step joins in the reference drug tables.
#'
#' @param metadata_parquet Character. Path to the Parquet metadata file.
#' @param out_path Character. Output directory for plot files.
#'
#' @return Invisibly returns the path to the written PDF (all plots as
#'   separate pages of one multi-page file).
#' @export
generatePlots <- function(metadata_parquet,
                          out_path) {
  if (!dir.exists(out_path)) {
    dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  }

  metadata <- arrow::read_parquet(normalizePath(metadata_parquet))

  if (nrow(metadata) == 0) {
    stop("The input metadata table is empty. Please check your query or input data.")
  }

  ## Generate plots
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
    ggplot2::facet_wrap(~drug_abbr, scales = "free_y") +
    ggplot2::labs(
      title = "Distribution of AMR phenotypes across antibiotics by year",
      x = "Year", y = "Number of isolates",
      colour = "Phenotype"
    ) +
    ggplot2::scale_color_manual(values = PHENOTYPE_COLORS, na.value = "gray70") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
      axis.text.y = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(colour = "black"),
      panel.grid.minor = ggplot2::element_blank()
    )

  # 2) Resistance only over time, colored by antibiotic (named palette aligned
  # to factor levels via meta_palette() so colors stay consistent across runs)
  abx_levels <- summary_year |>
    dplyr::filter(genome_drug.resistant_phenotype == "Resistant") |>
    dplyr::distinct(drug_abbr) |>
    dplyr::arrange(drug_abbr) |>
    dplyr::pull(drug_abbr)

  pal_named <- stats::setNames(meta_palette(length(abx_levels)), abx_levels)

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
      title  = "Distribution of resistant isolates by year",
      x      = "Year",
      y      = "Number of AMR isolates",
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
    ggplot2::scale_color_manual(values = PHENOTYPE_COLORS, na.value = "gray70") +
    ggplot2::labs(
      title = "Distribution of AMR phenotype by year and geography",
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
      title = "Distribution of AMR phenotypes",
      x = "Antibiotic", y = "Proportion", fill = "Phenotype"
    ) +
    ggplot2::scale_fill_manual(values = PHENOTYPE_COLORS, na.value = "gray70") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  # 5) Resistant isolates by isolation source (top 10 sources + Other, same
  # convention as amRviz's makeIsolationSourcesPlot)
  summary_isolation_source <- metadata |>
    dplyr::filter(genome.isolation_source != "") |>
    dplyr::group_by(genome.isolation_source, genome_drug.resistant_phenotype) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
    dplyr::filter(genome_drug.resistant_phenotype == "Resistant")

  top_isolation_sources <- summary_isolation_source |>
    dplyr::slice_max(order_by = count, n = 10) |>
    dplyr::pull(genome.isolation_source)

  summary_isolation_source <- summary_isolation_source |>
    dplyr::mutate(
      isolation_source_grp = ifelse(
        genome.isolation_source %in% top_isolation_sources,
        genome.isolation_source,
        "Other"
      )
    )

  n_sources <- dplyr::n_distinct(summary_isolation_source$isolation_source_grp)

  p5 <- ggplot2::ggplot(
    summary_isolation_source,
    ggplot2::aes(
      x = stats::reorder(isolation_source_grp, count),
      y = count, fill = isolation_source_grp
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Distribution of AMR isolates by isolation source",
      x = "Isolation source", y = "Number of AMR isolates"
    ) +
    ggplot2::scale_fill_manual(values = meta_palette(n_sources)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      legend.position = "none"
    )

  # 6) Histogram of resistant classes per genome
  p6 <- ggplot2::ggplot(metadata, ggplot2::aes(num_resistant_classes)) +
    ggplot2::geom_histogram(binwidth = 1, fill = META_COLORS[[1]]) +
    ggplot2::labs(
      title = "Distribution of AMR classes per isolate",
      x = "# Resistant Classes", y = "Count"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      legend.position = "bottom",
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(colour = "black"),
      panel.grid.minor = ggplot2::element_blank()
    )

  plots <- list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6)

  # Table for drug name - abbreviation mapping
  drug_table <- metadata |>
    dplyr::distinct(genome_drug.antibiotic, drug_abbr) |>
    dplyr::rename(
      "Abbreviation" = "drug_abbr",
      "Antibiotic" = "genome_drug.antibiotic"
    ) |>
    dplyr::arrange(Antibiotic)

  ## Write to device
  pdf_path <- file.path(out_path, "amRdata_exploratory_plots.pdf")
  grDevices::pdf(pdf_path, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (nm in names(plots)) {
    print(plots[[nm]])
  }
  # grid::grid.newpage()

  gridExtra::grid.arrange(
    grid::textGrob(
      "Antibiotic name abbreviations",
      gp = grid::gpar(fontsize = 16, fontface = "bold")
    ),
    gridExtra::tableGrob(drug_table, rows = NULL),
    ncol = 1,
    heights = c(0.08, 0.92)
  )

  invisible(pdf_path)
}
