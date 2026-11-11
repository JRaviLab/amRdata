#########################
#  Shared report prep   #
#########################

# Minimum number of results before a resistance percentage is shown without a
# "low n" flag (summary) or a drug gets its own row in the year heatmap (plots).
.REPORT_MIN_N <- 10L

# Standardise the cleaned `metadata` table for reporting: blank/"NULL" strings
# become NA, the year and class-count columns become numeric, optional columns
# that older exports lack are added as NA so downstream code can test for them,
# and `drug_abbr` falls back to the antibiotic name.
.prepReportMetadata <- function(metadata) {
  required <- c(
    "genome.genome_id", "genome_drug.antibiotic",
    "genome_drug.resistant_phenotype"
  )
  missing_cols <- setdiff(required, names(metadata))
  if (length(missing_cols)) {
    stop(
      "The metadata table is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  metadata <- tibble::as_tibble(metadata)

  chr_cols <- c(
    "drug_abbr", "drug_class", "class_abbr", "resistant_classes",
    "genome.isolation_country", "genome.isolation_source",
    "genome.host_common_name", "genome_drug.laboratory_typing_method",
    "genome_drug.pmid", "genome.species"
  )
  for (col in setdiff(chr_cols, names(metadata))) {
    metadata[[col]] <- NA_character_
  }
  for (col in setdiff(c("genome.collection_year", "num_resistant_classes"), names(metadata))) {
    metadata[[col]] <- NA_real_
  }

  blank_na <- function(x) {
    x <- trimws(as.character(x))
    x[x %in% c("", "NULL", "NA")] <- NA_character_
    x
  }

  metadata |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c(chr_cols, "genome_drug.antibiotic", "genome_drug.resistant_phenotype")),
        blank_na
      ),
      genome.collection_year = suppressWarnings(as.numeric(genome.collection_year)),
      num_resistant_classes = suppressWarnings(as.numeric(num_resistant_classes)),
      drug_abbr = dplyr::coalesce(drug_abbr, genome_drug.antibiotic)
    ) |>
    dplyr::filter(!is.na(genome_drug.resistant_phenotype))
}

# One row per genome, for anything that should count isolates rather than
# genome x drug results (countries, sources, MDR, ...). Genome-level fields are
# constant within a genome, so the first non-missing value is taken.
.genomeLevelMetadata <- function(metadata) {
  first_val <- function(x) x[!is.na(x)][1]
  metadata |>
    dplyr::group_by(genome.genome_id) |>
    dplyr::summarise(
      year = first_val(genome.collection_year),
      country = first_val(genome.isolation_country),
      source = first_val(genome.isolation_source),
      host = first_val(genome.host_common_name),
      n_classes = first_val(num_resistant_classes),
      profile = first_val(resistant_classes),
      any_resistant = any(genome_drug.resistant_phenotype == "Resistant"),
      .groups = "drop"
    )
}

# Collapse to one call per genome x drug class (Resistant > Intermediate >
# Susceptible), so a class counts once per isolate however many of its drugs
# were tested.
.classCalls <- function(metadata) {
  metadata |>
    dplyr::filter(!is.na(drug_class)) |>
    dplyr::group_by(genome.genome_id, drug_class, class_abbr) |>
    dplyr::summarise(
      genome_drug.resistant_phenotype = dplyr::case_when(
        any(genome_drug.resistant_phenotype == "Resistant") ~ "Resistant",
        any(genome_drug.resistant_phenotype == "Intermediate") ~ "Intermediate",
        TRUE ~ "Susceptible"
      ),
      .groups = "drop"
    )
}

# Resistant / Intermediate / Susceptible counts and % resistant per group.
.phenotypeCounts <- function(df, ...) {
  df |>
    dplyr::group_by(...) |>
    dplyr::summarise(
      Resistant = sum(genome_drug.resistant_phenotype == "Resistant"),
      Intermediate = sum(genome_drug.resistant_phenotype == "Intermediate"),
      Susceptible = sum(genome_drug.resistant_phenotype == "Susceptible"),
      Total = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(pct_resistant = Resistant / Total)
}

.speciesLabel <- function(metadata) {
  sp <- sort(unique(stats::na.omit(metadata$genome.species)))
  if (length(sp)) paste(sp, collapse = ", ") else "unknown species"
}

.fmtInt <- function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)
.fmtPct <- function(x, digits = 1) sprintf(paste0("%.", digits, "f%%"), 100 * x)
.pctAxis <- function(x) paste0(round(100 * x), "%")
.tidyLabel <- function(x) gsub("_", " ", x, fixed = TRUE)


#########################
#   .generateSummary    #
#########################

# Internal: called by writeMetadataReports().
#
# Writes a Markdown report that gives an at-a-glance overview of the cleaned
# metadata (isolates, results, antibiotics, multidrug resistance, collection
# years, geography, sources, hosts, laboratory methods, publications), data
# completeness, and per-antibiotic and per-class resistance tables. Results
# with fewer than 10 observations are flagged because their percentages are
# unreliable.
#
# Counts are reported in two units: *isolates* are unique genomes, and
# *results* are genome x antibiotic susceptibility calls. Per-class calls
# collapse a genome's drugs in that class (any resistant call makes the class
# resistant). Multidrug resistant (MDR) means resistant to three or more drug
# classes.
#
# @param duckdb_path Character or NULL. Which ORB to summarize. If NULL
#   (default), completed ORBs are discovered from the dataset manifests and
#   you are prompted to choose when more than one exists (a path is required
#   in non-interactive sessions). May also be an ORB directory or an ORB
#   DuckDB path.
# @param out_path Character or NULL. Directory where the Markdown report is
#   written. Defaults to the ORB folder.
#
# @return Invisibly returns the path to the written Markdown summary report
#   ("<out_path>/amr_metadata_summary.md").
#
.generateSummary <- function(duckdb_path = NULL,
                             out_path = NULL
) {
  resolved <- .readOrbMetadata(duckdb_path)
  to_orb <- is.null(out_path)
  if (to_orb) out_path <- resolved$out_path

  if (nrow(resolved$metadata) == 0) {
    stop("The output table is empty. Please check your query or input data.")
  }
  metadata <- .prepReportMetadata(resolved$metadata)
  if (nrow(metadata) == 0) {
    stop("No susceptibility results with a resistance phenotype were found.")
  }

  if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
  md_path <- file.path(out_path, "amr_metadata_summary.md")

  genomes <- .genomeLevelMetadata(metadata)
  n_genomes <- nrow(genomes)
  n_results <- nrow(metadata)

  md_tbl <- function(df, align = NULL) {
    knitr::kable(df, format = "pipe", align = align, row.names = FALSE)
  }

  # A "name | count | share" table for a genome-level field. Shares are of the
  # isolates that have a value; the remainder beyond `top` is lumped together.
  isolate_tbl <- function(col, header, top = 10L) {
    counts <- genomes |>
      dplyr::filter(!is.na(.data[[col]])) |>
      dplyr::count(.data[[col]], name = "Isolates", sort = TRUE)
    if (nrow(counts) == 0) {
      return(NULL)
    }
    if (nrow(counts) > top) {
      rest <- counts[-seq_len(top), ]
      counts <- dplyr::bind_rows(
        counts[seq_len(top), ],
        stats::setNames(
          tibble::tibble(sprintf("Other (%d values)", nrow(rest)), sum(rest$Isolates)),
          names(counts)
        )
      )
    }
    out <- tibble::tibble(
      counts[[1]],
      .fmtInt(counts$Isolates),
      .fmtPct(counts$Isolates / sum(counts$Isolates))
    )
    names(out) <- c(header, "Isolates", "Share")
    md_tbl(out, align = c("l", "r", "r"))
  }

  # Resistance rate table (per antibiotic / class / year / country ...).
  # `key` is the grouping column in `df`, `header` its displayed name.
  rate_tbl <- function(counts, key_label, extra = NULL) {
    has_i <- any(counts$Intermediate > 0)
    out <- tibble::tibble(key_label)
    if (!is.null(extra)) out <- tibble::tibble(key_label, extra)
    out$Resistant <- .fmtInt(counts$Resistant)
    if (has_i) out$Intermediate <- .fmtInt(counts$Intermediate)
    out$Susceptible <- .fmtInt(counts$Susceptible)
    out$Total <- .fmtInt(counts$Total)
    out$pct <- paste0(
      .fmtPct(counts$pct_resistant),
      ifelse(counts$Total < .REPORT_MIN_N, "*", "")
    )
    out
  }

  has_class <- any(!is.na(metadata$drug_class))
  has_mdr <- any(!is.na(genomes$n_classes))
  has_low_n <- FALSE

  lines <- c(
    sprintf("# AMR summary report for *%s*", .speciesLabel(metadata)),
    "",
    sprintf("_Generated %s by amRdata._", format(Sys.Date(), "%Y-%m-%d")),
    ""
  )

  ## Overview ---------------------------------------------------------------
  n_res <- sum(metadata$genome_drug.resistant_phenotype == "Resistant")
  years <- genomes$year[!is.na(genomes$year)]
  n_countries <- dplyr::n_distinct(genomes$country, na.rm = TRUE)
  pmids <- unique(trimws(unlist(strsplit(
    stats::na.omit(metadata$genome_drug.pmid), "[,;[:space:]]+"
  ))))
  pmids <- sort(pmids[nzchar(pmids)])

  overview <- tibble::tribble(
    ~Metric, ~Value,
    "Isolates (unique genomes)", .fmtInt(n_genomes),
    "Susceptibility results (genome x antibiotic)", .fmtInt(n_results),
    "Antibiotics", .fmtInt(dplyr::n_distinct(metadata$genome_drug.antibiotic)),
    "Resistant results", sprintf("%s (%s)", .fmtInt(n_res), .fmtPct(n_res / n_results)),
    "Isolates resistant to at least one antibiotic",
    sprintf(
      "%s (%s)", .fmtInt(sum(genomes$any_resistant)),
      .fmtPct(mean(genomes$any_resistant))
    ),
    "Collection years",
    if (length(years)) sprintf("%d to %d", min(years), max(years)) else "Not recorded",
    "Countries", if (n_countries) .fmtInt(n_countries) else "Not recorded",
    "Publications", .fmtInt(length(pmids))
  )
  if (has_class) {
    overview <- tibble::add_row(
      overview,
      Metric = "Antibiotic classes",
      Value = .fmtInt(dplyr::n_distinct(metadata$drug_class, na.rm = TRUE)),
      .after = 3
    )
  }
  if (has_mdr) {
    n_mdr <- sum(genomes$n_classes >= 3, na.rm = TRUE)
    n_known <- sum(!is.na(genomes$n_classes))
    overview <- tibble::add_row(
      overview,
      Metric = "Multidrug resistant isolates (3 or more classes)",
      Value = sprintf("%s (%s)", .fmtInt(n_mdr), .fmtPct(n_mdr / n_known)),
      .after = which(overview$Metric == "Isolates resistant to at least one antibiotic")
    )
  }
  lines <- c(lines, "## Overview", "", md_tbl(overview, align = c("l", "r")), "")

  ## Completeness -----------------------------------------------------------
  n_with <- c(
    sum(!is.na(genomes$year)), sum(!is.na(genomes$country)),
    sum(!is.na(genomes$source)), sum(!is.na(genomes$host))
  )
  completeness <- tibble::tibble(
    Field = c("Collection year", "Isolation country", "Isolation source", "Host"),
    `Isolates with a value` = .fmtInt(n_with),
    Coverage = .fmtPct(n_with / n_genomes)
  )
  lines <- c(
    lines, "## Data completeness", "",
    "Share of isolates that have each field recorded. Tables below describe only those isolates.",
    "",
    md_tbl(completeness, align = c("l", "r", "r")), ""
  )

  ## Per antibiotic ----------------------------------------------------------
  drug_counts <- .phenotypeCounts(metadata, genome_drug.antibiotic, drug_abbr, drug_class) |>
    dplyr::arrange(dplyr::desc(Total), genome_drug.antibiotic)
  has_low_n <- has_low_n || any(drug_counts$Total < .REPORT_MIN_N)
  drug_out <- rate_tbl(
    drug_counts,
    key_label = .tidyLabel(drug_counts$genome_drug.antibiotic),
    extra = drug_counts$drug_abbr
  )
  names(drug_out)[1:2] <- c("Antibiotic", "Abbrev")
  if (has_class) drug_out$Class <- .tidyLabel(drug_counts$drug_class)
  names(drug_out)[names(drug_out) == "pct"] <- "% resistant"
  lines <- c(
    lines, "## Resistance by antibiotic", "",
    "Sorted by number of results. Percentages use genome x antibiotic results.", "",
    md_tbl(drug_out, align = c("l", "l", rep("r", ncol(drug_out) - 2 - has_class), if (has_class) "l")),
    ""
  )

  ## Per class ----------------------------------------------------------------
  if (has_class) {
    calls <- .classCalls(metadata)
    class_counts <- .phenotypeCounts(calls, drug_class, class_abbr) |>
      dplyr::arrange(dplyr::desc(Total), drug_class)
    has_low_n <- has_low_n || any(class_counts$Total < .REPORT_MIN_N)
    class_out <- rate_tbl(
      class_counts,
      key_label = .tidyLabel(class_counts$drug_class),
      extra = class_counts$class_abbr
    )
    names(class_out)[1:2] <- c("Drug class", "Abbrev")
    names(class_out)[names(class_out) == "Total"] <- "Isolates"
    names(class_out)[names(class_out) == "pct"] <- "% resistant"
    lines <- c(
      lines, "## Resistance by antibiotic class", "",
      paste(
        "One call per isolate and class: an isolate is resistant to a class if",
        "it is resistant to any tested drug in it. Counts here are isolates."
      ),
      "",
      md_tbl(class_out, align = c("l", "l", rep("r", ncol(class_out) - 2))),
      ""
    )
  }

  ## Multidrug resistance -------------------------------------------------------
  if (has_mdr) {
    mdr <- genomes |>
      dplyr::filter(!is.na(n_classes)) |>
      dplyr::count(n_classes, name = "Isolates") |>
      dplyr::arrange(n_classes)
    mdr_out <- tibble::tibble(
      `Resistant classes` = ifelse(
        mdr$n_classes == 0, "0 (susceptible to all tested)", as.character(mdr$n_classes)
      ),
      Isolates = .fmtInt(mdr$Isolates),
      Share = .fmtPct(mdr$Isolates / sum(mdr$Isolates))
    )
    lines <- c(
      lines, "## Multidrug resistance", "",
      md_tbl(mdr_out, align = c("l", "r", "r")), ""
    )

    profiles <- genomes |>
      dplyr::filter(!is.na(profile), n_classes > 0) |>
      dplyr::count(profile, name = "Isolates", sort = TRUE) |>
      utils::head(10)
    if (nrow(profiles)) {
      prof_out <- tibble::tibble(
        `Resistance profile (class abbreviations)` = gsub("_", " + ", profiles$profile, fixed = TRUE),
        Isolates = .fmtInt(profiles$Isolates),
        Share = .fmtPct(profiles$Isolates / sum(genomes$n_classes > 0, na.rm = TRUE))
      )
      lines <- c(
        lines, "### Most common resistance profiles", "",
        "Share is of isolates resistant to at least one class.", "",
        md_tbl(prof_out, align = c("l", "r", "r")), ""
      )
    }
  }

  ## Time, geography, source, host ---------------------------------------------
  year_counts <- metadata |>
    dplyr::filter(!is.na(genome.collection_year)) |>
    dplyr::group_by(genome.collection_year) |>
    dplyr::summarise(
      Isolates = dplyr::n_distinct(genome.genome_id),
      Entries = dplyr::n(),
      Resistant = sum(genome_drug.resistant_phenotype == "Resistant"),
      .groups = "drop"
    ) |>
    dplyr::arrange(genome.collection_year)
  if (nrow(year_counts)) {
    year_out <- tibble::tibble(
      Year = as.character(year_counts$genome.collection_year),
      Isolates = .fmtInt(year_counts$Isolates),
      Entries = .fmtInt(year_counts$Entries),
      `% resistant` = .fmtPct(year_counts$Resistant / year_counts$Entries)
    )
    lines <- c(
      lines, "## Collection years", "",
      md_tbl(year_out, align = c("l", "r", "r", "r")), ""
    )
  }

  country_counts <- metadata |>
    dplyr::filter(!is.na(genome.isolation_country)) |>
    dplyr::group_by(genome.isolation_country) |>
    dplyr::summarise(
      Isolates = dplyr::n_distinct(genome.genome_id),
      Entries = dplyr::n(),
      Resistant = sum(genome_drug.resistant_phenotype == "Resistant"),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(Isolates), genome.isolation_country)
  if (nrow(country_counts)) {
    shown <- utils::head(country_counts, 15)
    country_out <- tibble::tibble(
      Country = shown$genome.isolation_country,
      Isolates = .fmtInt(shown$Isolates),
      Entries = .fmtInt(shown$Entries),
      `% resistant` = .fmtPct(shown$Resistant / shown$Entries)
    )
    lines <- c(lines, "## Isolation countries", "")
    if (nrow(country_counts) > nrow(shown)) {
      lines <- c(lines, sprintf(
        "Top %d of %d countries; the other %d account for %s isolates.", nrow(shown),
        nrow(country_counts), nrow(country_counts) - nrow(shown),
        .fmtInt(sum(country_counts$Isolates[-seq_len(nrow(shown))]))
      ), "")
    }
    lines <- c(lines, md_tbl(country_out, align = c("l", "r", "r", "r")), "")
  }

  for (spec in list(
    list("source", "Isolation source", "## Isolation sources"),
    list("host", "Host", "## Hosts")
  )) {
    tbl <- isolate_tbl(spec[[1]], spec[[2]])
    if (!is.null(tbl)) lines <- c(lines, spec[[3]], "", tbl, "")
  }

  ## Lab methods -----------------------------------------------------------------
  lab <- metadata |>
    dplyr::mutate(
      method = dplyr::coalesce(genome_drug.laboratory_typing_method, "Not defined")
    ) |>
    dplyr::count(method, name = "Entries", sort = TRUE)
  lab_out <- tibble::tibble(
    `Laboratory method` = lab$method,
    Entries = .fmtInt(lab$Entries),
    Share = .fmtPct(lab$Entries / sum(lab$Entries))
  )
  lines <- c(
    lines, "## Laboratory methods", "",
    md_tbl(lab_out, align = c("l", "r", "r")), ""
  )

  ## Publications ----------------------------------------------------------------
  if (length(pmids)) {
    lines <- c(
      lines, sprintf("## Associated publications (%d)", length(pmids)), "",
      paste(
        sprintf("[%s](https://pubmed.ncbi.nlm.nih.gov/%s/)", pmids, pmids),
        collapse = ", "
      ),
      ""
    )
  }

  if (has_low_n) {
    lines <- c(
      lines, "---", "",
      sprintf("\\* Fewer than %d results, so the percentage is unreliable.", .REPORT_MIN_N),
      ""
    )
  }

  writeLines(enc2utf8(lines), md_path, useBytes = TRUE)

  if (to_orb) {
    message("Summary report generated into a md file and added to the ORB.")
  } else {
    message("Summary report generated into a md file: ", md_path)
  }
  invisible(md_path)
}


#########################
#     .generatePlots    #
#########################

# Shared look for every report plot. ASCII-only text throughout so the PDF
# device's default encoding can render it.
.amr_theme <- function(base_size = 11, hbars = FALSE) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      # horizontal bars carry value labels at the bar ends; gridlines would
      # run through them
      panel.grid.major.y = if (hbars) ggplot2::element_blank(),
      text = ggplot2::element_text(colour = "#2D2D2D"),
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.25)),
      plot.title.position = "plot",
      plot.subtitle = ggplot2::element_text(colour = "#5a5a5a"),
      plot.caption = ggplot2::element_text(colour = "#7a7a7a", hjust = 0),
      plot.caption.position = "plot",
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

# Integer breaks spanning the observed years (at most ~15 labels). Computed from
# the data rather than a break function, because ggplot hands those the
# expanded axis range and would label years beyond the data.
.yearBreaks <- function(years) {
  all_years <- seq(min(years), max(years))
  all_years[seq(1, length(all_years), by = max(1, ceiling(length(all_years) / 15)))]
}

# Whole-number breaks for count axes (avoids ticks like 2.5 isolates).
.countBreaks <- function(limits) {
  b <- pretty(limits)
  b[b == round(b)]
}

# Keep the `top` most frequent values of `x` and pool the rest as "Other".
.lumpTop <- function(x, top) {
  keep <- names(sort(table(x), decreasing = TRUE))
  keep <- utils::head(keep, top)
  ifelse(x %in% keep, x, "Other")
}

# Stacked 100% bars of Resistant/Intermediate/Susceptible per group with the
# number of results printed at the end of each bar. Highest resistance on top.
.plotResistanceBars <- function(counts, label_col, title, unit, top = 30L) {
  if (nrow(counts) == 0) {
    return(NULL)
  }
  counts <- counts |>
    dplyr::arrange(dplyr::desc(Total)) |>
    utils::head(top)
  counts$label <- .tidyLabel(counts[[label_col]])
  counts$label <- factor(
    counts$label,
    levels = counts$label[order(counts$pct_resistant, counts$Total)]
  )
  long <- counts |>
    tidyr::pivot_longer(
      c("Resistant", "Intermediate", "Susceptible"),
      names_to = "phenotype", values_to = "n"
    ) |>
    dplyr::filter(n > 0) |>
    dplyr::mutate(phenotype = factor(phenotype, c("Resistant", "Intermediate", "Susceptible")))

  ggplot2::ggplot(long, ggplot2::aes(x = n, y = label, fill = phenotype)) +
    # Resistant first level + reverse: amber starts at 0%, so % resistant can
    # be read straight off the axis.
    ggplot2::geom_col(position = ggplot2::position_fill(reverse = TRUE), width = 0.75) +
    ggplot2::geom_text(
      data = counts,
      ggplot2::aes(x = 1, y = label, label = paste0("n=", .fmtInt(Total))),
      inherit.aes = FALSE, hjust = -0.15, size = 3, colour = "#5a5a5a"
    ) +
    ggplot2::scale_x_continuous(
      labels = .pctAxis, breaks = seq(0, 1, 0.25),
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::scale_fill_manual(values = PHENOTYPE_COLORS, name = NULL, drop = TRUE) +
    ggplot2::labs(
      title = title, x = paste("Share of", unit), y = NULL,
      caption = if (nrow(counts) == top) paste("Top", top, "by number of", unit) else NULL
    ) +
    .amr_theme(hbars = TRUE)
}

# Antibiotic x collection year heatmap: colour = % resistant, label = n.
.plotDrugYearHeatmap <- function(metadata) {
  dated <- metadata |>
    dplyr::filter(!is.na(genome.collection_year))
  if (dplyr::n_distinct(dated$genome.collection_year) < 2) {
    return(NULL)
  }
  totals <- .phenotypeCounts(dated, drug_abbr)
  keep <- totals$drug_abbr[totals$Total >= .REPORT_MIN_N]
  if (length(keep) == 0) {
    keep <- utils::head(totals$drug_abbr[order(-totals$Total)], 12)
  }
  cells <- dated |>
    dplyr::filter(drug_abbr %in% keep) |>
    .phenotypeCounts(drug_abbr, genome.collection_year)
  drug_order <- totals |>
    dplyr::filter(drug_abbr %in% keep) |>
    dplyr::arrange(pct_resistant, Total)
  cells$drug_abbr <- factor(cells$drug_abbr, levels = drug_order$drug_abbr)
  n_years <- dplyr::n_distinct(cells$genome.collection_year)
  n_hidden <- length(setdiff(totals$drug_abbr, keep))

  ggplot2::ggplot(
    cells,
    ggplot2::aes(genome.collection_year, drug_abbr, fill = pct_resistant)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    {
      if (n_years <= 25) {
        ggplot2::geom_text(ggplot2::aes(label = Total), size = 2.8, colour = "#2D2D2D")
      }
    } +
    ggplot2::scale_x_continuous(
      breaks = .yearBreaks(cells$genome.collection_year),
      expand = ggplot2::expansion(add = 0.5)
    ) +
    ggplot2::scale_fill_gradient2(
      low = PHENOTYPE_COLORS[["Susceptible"]], mid = "#f4f1ea",
      high = PHENOTYPE_COLORS[["Resistant"]], midpoint = 0.5, limits = c(0, 1),
      labels = .pctAxis, name = "% resistant"
    ) +
    ggplot2::labs(
      title = "Resistance by antibiotic and collection year",
      x = "Collection year", y = NULL,
      caption = paste0(
        "Numbers in tiles are entries; blank = no entries that year.",
        if (n_hidden > 0) {
          sprintf(" %d antibiotic(s) with fewer than %d entries not shown.", n_hidden, .REPORT_MIN_N)
        } else {
          ""
        }
      )
    ) +
    .amr_theme() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.key.width = grid::unit(1.6, "cm")
    )
}

# Susceptibility results per collection year, split by phenotype.
.plotResultsByYear <- function(metadata) {
  dated <- metadata |>
    dplyr::filter(!is.na(genome.collection_year))
  if (nrow(dated) == 0) {
    return(NULL)
  }
  counts <- dated |>
    dplyr::count(genome.collection_year, genome_drug.resistant_phenotype, name = "n") |>
    dplyr::mutate(phenotype = factor(
      genome_drug.resistant_phenotype, c("Resistant", "Intermediate", "Susceptible")
    ))
  totals <- counts |>
    dplyr::group_by(genome.collection_year) |>
    dplyr::summarise(total = sum(n), .groups = "drop")
  n_missing <- sum(is.na(metadata$genome.collection_year))

  ggplot2::ggplot(
    counts, ggplot2::aes(genome.collection_year, n, fill = phenotype)
  ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_text(
      data = totals,
      ggplot2::aes(genome.collection_year, total, label = .fmtInt(total)),
      inherit.aes = FALSE, vjust = -0.5, size = 3, colour = "#5a5a5a"
    ) +
    ggplot2::scale_x_continuous(breaks = .yearBreaks(counts$genome.collection_year)) +
    ggplot2::scale_y_continuous(
      labels = function(x) .fmtInt(x), expand = ggplot2::expansion(mult = c(0, 0.1))
    ) +
    ggplot2::scale_fill_manual(values = PHENOTYPE_COLORS, name = NULL, drop = TRUE) +
    ggplot2::labs(
      title = "Susceptibility results by collection year",
      x = "Collection year", y = "Results (genome x antibiotic)",
      caption = if (n_missing > 0) {
        sprintf("%s result(s) without a collection year not shown.", .fmtInt(n_missing))
      }
    ) +
    .amr_theme() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
}

# Isolates per category (country, isolation source, ...), coloured by whether
# the isolate is resistant to at least one antibiotic. Counts isolates, not
# genome x drug rows.
.plotIsolateBars <- function(genomes, col, title, ylab, top = 12L) {
  df <- genomes |>
    dplyr::filter(!is.na(.data[[col]]))
  if (nrow(df) == 0) {
    return(NULL)
  }
  df$key <- .lumpTop(df[[col]], top)
  df$status <- factor(
    ifelse(df$any_resistant, "Resistant to at least one antibiotic", "Susceptible to all tested"),
    c("Resistant to at least one antibiotic", "Susceptible to all tested")
  )
  counts <- df |>
    dplyr::count(key, status, name = "n")
  totals <- counts |>
    dplyr::group_by(key) |>
    dplyr::summarise(total = sum(n), .groups = "drop") |>
    dplyr::mutate(sort_key = ifelse(key == "Other", -1, total))
  counts$key <- factor(counts$key, levels = totals$key[order(totals$sort_key)])
  totals$key <- factor(totals$key, levels = levels(counts$key))
  n_missing <- nrow(genomes) - nrow(df)
  fills <- c(
    "Resistant to at least one antibiotic" = PHENOTYPE_COLORS[["Resistant"]],
    "Susceptible to all tested" = PHENOTYPE_COLORS[["Susceptible"]]
  )

  ggplot2::ggplot(counts, ggplot2::aes(n, key, fill = status)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_text(
      data = totals, ggplot2::aes(total, key, label = .fmtInt(total)),
      inherit.aes = FALSE, hjust = -0.2, size = 3, colour = "#5a5a5a"
    ) +
    ggplot2::scale_x_continuous(
      breaks = .countBreaks, labels = function(x) .fmtInt(x),
      expand = ggplot2::expansion(mult = c(0, 0.1))
    ) +
    ggplot2::scale_fill_manual(values = fills, name = NULL, drop = TRUE) +
    ggplot2::labs(
      title = title, x = "Isolates", y = ylab,
      caption = if (n_missing > 0) {
        sprintf("%s isolate(s) with no value recorded not shown.", .fmtInt(n_missing))
      }
    ) +
    .amr_theme(hbars = TRUE)
}

# Country x collection year tiles showing how many isolates were sampled.
.plotCountryYear <- function(genomes, top = 15L) {
  df <- genomes |>
    dplyr::filter(!is.na(country), !is.na(year))
  if (dplyr::n_distinct(df$country) < 2 || dplyr::n_distinct(df$year) < 2) {
    return(NULL)
  }
  df$country <- .lumpTop(df$country, top)
  cells <- df |>
    dplyr::count(country, year, name = "n")
  order_tbl <- cells |>
    dplyr::group_by(country) |>
    dplyr::summarise(total = sum(n), .groups = "drop") |>
    dplyr::mutate(sort_key = ifelse(country == "Other", -1, total))
  cells$country <- factor(cells$country, levels = order_tbl$country[order(order_tbl$sort_key)])

  ggplot2::ggplot(cells, ggplot2::aes(year, country, fill = n)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = n, colour = n > max(n) * 0.6), size = 2.8, show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "#2D2D2D")) +
    ggplot2::scale_x_continuous(
      breaks = .yearBreaks(cells$year), expand = ggplot2::expansion(add = 0.5)
    ) +
    ggplot2::scale_fill_gradient(
      low = "#e8eef4", high = META_COLORS[[1]], name = "Isolates"
    ) +
    ggplot2::labs(
      title = "Where and when isolates were collected",
      x = "Collection year", y = NULL,
      caption = "Isolates with both a country and a collection year."
    ) +
    .amr_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

# Isolates by number of resistant drug classes; MDR (3+) highlighted.
.plotMdr <- function(genomes) {
  df <- genomes |>
    dplyr::filter(!is.na(n_classes))
  if (nrow(df) == 0) {
    return(NULL)
  }
  levels_n <- seq(0, max(df$n_classes))
  counts <- tibble::tibble(n_classes = levels_n) |>
    dplyr::left_join(dplyr::count(df, n_classes, name = "n"), by = "n_classes") |>
    dplyr::mutate(
      n = dplyr::coalesce(n, 0L),
      group = factor(
        ifelse(n_classes >= 3, "MDR (3 or more classes)", "Not MDR"),
        c("Not MDR", "MDR (3 or more classes)")
      ),
      x = factor(n_classes, levels_n)
    )
  n_mdr <- sum(df$n_classes >= 3)

  ggplot2::ggplot(counts, ggplot2::aes(x, n, fill = group)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(n > 0, .fmtInt(n), "")),
      vjust = -0.5, size = 3, colour = "#5a5a5a"
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) .fmtInt(x), expand = ggplot2::expansion(mult = c(0, 0.1))
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Not MDR" = PHENOTYPE_COLORS[["Susceptible"]],
        "MDR (3 or more classes)" = PHENOTYPE_COLORS[["Resistant"]]
      ),
      name = NULL
    ) +
    ggplot2::labs(
      title = "Multidrug resistance",
      subtitle = sprintf(
        "%s of %s isolates (%s) are resistant to 3 or more drug classes",
        .fmtInt(n_mdr), .fmtInt(nrow(df)), .fmtPct(n_mdr / nrow(df))
      ),
      x = "Number of resistant drug classes", y = "Isolates"
    ) +
    .amr_theme() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
}

# Most common combinations of resistant drug classes.
.plotProfiles <- function(genomes, top = 15L) {
  df <- genomes |>
    dplyr::filter(!is.na(profile), !is.na(n_classes), n_classes > 0)
  if (nrow(df) == 0) {
    return(NULL)
  }
  counts <- df |>
    dplyr::count(profile, name = "n", sort = TRUE) |>
    utils::head(top)
  counts$label <- factor(
    gsub("_", " + ", counts$profile, fixed = TRUE),
    levels = rev(gsub("_", " + ", counts$profile, fixed = TRUE))
  )

  ggplot2::ggplot(counts, ggplot2::aes(n, label)) +
    ggplot2::geom_col(fill = PHENOTYPE_COLORS[["Resistant"]], width = 0.75) +
    ggplot2::geom_text(
      ggplot2::aes(label = .fmtInt(n)), hjust = -0.2, size = 3, colour = "#5a5a5a"
    ) +
    ggplot2::scale_x_continuous(
      breaks = .countBreaks, labels = function(x) .fmtInt(x),
      expand = ggplot2::expansion(mult = c(0, 0.1))
    ) +
    ggplot2::labs(
      title = "Most common resistance profiles",
      subtitle = sprintf(
        "Top %d of %s distinct combinations of resistant drug classes",
        nrow(counts), .fmtInt(dplyr::n_distinct(df$profile))
      ),
      x = "Isolates", y = "Resistant classes (abbreviations on the last page)"
    ) +
    .amr_theme(hbars = TRUE)
}

# Paginated reference table(s) (landscape page, top-aligned).
.drawReferenceTable <- function(tbl, title, rows_per_page = 26L) {
  pages <- split(seq_len(nrow(tbl)), ceiling(seq_len(nrow(tbl)) / rows_per_page))
  for (i in seq_along(pages)) {
    grob <- gridExtra::tableGrob(
      tbl[pages[[i]], , drop = FALSE], rows = NULL,
      theme = gridExtra::ttheme_minimal(base_size = 10)
    )
    grob$vp <- grid::viewport(
      y = grid::unit(0.9, "npc"), just = "top",
      width = sum(grob$widths), height = sum(grob$heights)
    )
    grid::grid.newpage()
    grid::grid.text(
      if (length(pages) > 1) sprintf("%s (%d of %d)", title, i, length(pages)) else title,
      x = 0.04, y = 0.95, just = c("left", "top"),
      gp = grid::gpar(fontsize = 16, fontface = "bold", col = "#2D2D2D")
    )
    grid::grid.draw(grob)
  }
}

# Internal: called by writeMetadataReports().
#
# Writes a landscape multi-page PDF of exploratory plots: resistance by
# antibiotic and by drug class (with sample sizes), resistance by antibiotic
# and collection year, results per year, isolates by country, by country and
# year, and by isolation source, multidrug resistance, the most common
# resistance profiles, and reference tables for antibiotic and class
# abbreviations. Isolate-level plots count unique genomes, not genome x
# antibiotic rows. A plot is skipped (with no error) when the data it needs
# is absent, e.g. no collection years.
#
# Expects the ORB `metadata` table to be the output of
# `cleanMetaData()` step (or an export of the resulting `metadata` table), since
# `drug_abbr`, `drug_class`, and `num_resistant_classes` are only populated
# after that step joins in the reference drug tables.
#
# @param duckdb_path Character or NULL. Which ORB to plot. If NULL (default),
#   completed ORBs are discovered from the dataset manifests and you are
#   prompted to choose when more than one exists (a path is required in
#   non-interactive sessions). May also be an ORB directory or an ORB
#   DuckDB path.
# @param out_path Character or NULL. Output directory for plot files. Defaults
#   to the ORB folder.
#
# @return Invisibly returns the path to the written PDF (all plots as
#   separate pages of one multi-page file).
#
.generatePlots <- function(duckdb_path = NULL,
                           out_path = NULL
) {
  resolved <- .readOrbMetadata(duckdb_path)
  to_orb <- is.null(out_path)
  if (to_orb) out_path <- resolved$out_path

  if (nrow(resolved$metadata) == 0) {
    stop("The input metadata table is empty. Please check your query or input data.")
  }
  metadata <- .prepReportMetadata(resolved$metadata)
  if (nrow(metadata) == 0) {
    stop("No susceptibility results with a resistance phenotype were found.")
  }

  if (!dir.exists(out_path)) {
    dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  }

  genomes <- .genomeLevelMetadata(metadata)
  subtitle <- sprintf(
    "%s | %s isolates, %s susceptibility results",
    .speciesLabel(metadata), .fmtInt(nrow(genomes)), .fmtInt(nrow(metadata))
  )

  has_class <- any(!is.na(metadata$drug_class))
  plots <- list(
    drug = .plotResistanceBars(
      .phenotypeCounts(metadata, drug_abbr),
      "drug_abbr", "Resistance by antibiotic", "results"
    ),
    class = if (has_class) {
      .plotResistanceBars(
        .phenotypeCounts(.classCalls(metadata), drug_class),
        "drug_class", "Resistance by antibiotic class", "isolates"
      )
    },
    drug_year = .plotDrugYearHeatmap(metadata),
    year = .plotResultsByYear(metadata),
    country = .plotIsolateBars(genomes, "country", "Isolates by country", "Country"),
    country_year = .plotCountryYear(genomes),
    source = .plotIsolateBars(
      genomes, "source", "Isolates by isolation source", "Isolation source", top = 10L
    ),
    mdr = .plotMdr(genomes),
    profiles = .plotProfiles(genomes)
  )
  plots <- plots[!vapply(plots, is.null, logical(1))]

  ## Write to device: landscape pages. Plain pdf() rather than cairo_pdf(),
  ## because capabilities("cairo") can report TRUE on macOS without XQuartz
  ## even though the cairo library then fails to load.
  pdf_path <- file.path(out_path, "amRdata_exploratory_plots.pdf")
  grDevices::pdf(pdf_path, width = 11, height = 8.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (nm in names(plots)) {
    tryCatch(
      print(plots[[nm]] + ggplot2::labs(subtitle = plots[[nm]]$labels$subtitle %||% subtitle)),
      error = function(e) {
        warning("Skipped the '", nm, "' plot: ", conditionMessage(e), call. = FALSE)
      }
    )
  }

  # Reference tables: antibiotic (and class) abbreviations
  drug_table <- metadata |>
    dplyr::distinct(genome_drug.antibiotic, drug_abbr, drug_class) |>
    dplyr::arrange(genome_drug.antibiotic) |>
    dplyr::transmute(
      Abbreviation = drug_abbr,
      Antibiotic = .tidyLabel(genome_drug.antibiotic),
      Class = .tidyLabel(dplyr::coalesce(drug_class, ""))
    )
  .drawReferenceTable(drug_table, "Antibiotic name abbreviations")

  if (has_class && any(!is.na(metadata$class_abbr))) {
    class_table <- metadata |>
      dplyr::filter(!is.na(drug_class), !is.na(class_abbr)) |>
      dplyr::distinct(class_abbr, drug_class) |>
      dplyr::arrange(drug_class) |>
      dplyr::transmute(Abbreviation = class_abbr, `Drug class` = .tidyLabel(drug_class))
    .drawReferenceTable(class_table, "Antibiotic class abbreviations")
  }

  if (to_orb) {
    message("Exploratory plots generated into a PDF and added to the ORB.")
  } else {
    message("Exploratory plots generated into a PDF: ", pdf_path)
  }
  invisible(pdf_path)
}


#########################
#   Pipeline hook       #
#########################

#' Write the metadata summary and exploratory plots for an ORB
#'
#' Writes two reports of the cleaned `metadata` table of an ORB into the ORB
#' folder: a Markdown summary (`amr_metadata_summary.md`) and a multi-page PDF of
#' exploratory plots (`amRdata_exploratory_plots.pdf`).
#'
#' The summary gives an at-a-glance overview of the cleaned metadata (isolates,
#' results, antibiotics, multidrug resistance, collection years, geography,
#' sources, hosts, laboratory methods, publications), data completeness, and
#' per-antibiotic and per-class resistance tables. Results with fewer than 10
#' observations are flagged because their percentages are unreliable. The plots
#' cover resistance by antibiotic and class, by collection year, isolates by
#' country and source, multidrug resistance and the most common resistance
#' profiles, plus reference tables of the antibiotic and class abbreviations.
#'
#' Counts are reported in two units: *isolates* are unique genomes, and
#' *results* are genome x antibiotic susceptibility calls. Per-class calls
#' collapse a genome's drugs in that class (any resistant call makes the class
#' resistant). Multidrug resistant (MDR) means resistant to three or more drug
#' classes.
#'
#' The reports are a convenience, so a report that cannot be written is warned
#' about instead of stopping the call (or a long [runDataProcessing()] run). The
#' metadata must be the output of [cleanMetaData()], since `drug_abbr`,
#' `drug_class` and `num_resistant_classes` are only filled in by that step.
#'
#' @param duckdb_path Character or NULL. Which ORB to report on. If NULL
#'   (default), completed ORBs are discovered from the dataset manifests and you
#'   are prompted to choose when more than one exists (a path is required in
#'   non-interactive sessions). May also be an ORB directory or an ORB DuckDB
#'   path.
#' @param out_path Character or NULL. Directory the reports are written to.
#'   Defaults to the ORB folder.
#' @param status A function called with progress messages. Defaults to
#'   [message()]; [runDataProcessing()] passes its own logger.
#'
#' @return Invisibly returns the paths of the files actually written
#'   (`character(0)` if none), so a caller can record just those and not a
#'   stale file left by an earlier run.
#'
#' @examples
#' \dontrun{
#' writeMetadataReports()  # pick from available ORBs
#' writeMetadataReports("data/Shigella_flexneri/orb")
#' }
#'
#' @export
writeMetadataReports <- function(duckdb_path = NULL,
                                 out_path = NULL,
                                 status = function(...) message(...)) {
  report <- function(label, fn) {
    tryCatch(
      {
        path <- suppressMessages(fn(duckdb_path, out_path = out_path))
        status("Wrote ", label, ".")
        normalizePath(path, mustWork = TRUE)
      },
      error = function(e) {
        status("Could not write ", label, ": ", conditionMessage(e))
        warning(
          "Could not write the ", label, ": ", conditionMessage(e),
          call. = FALSE
        )
        NA_character_
      }
    )
  }
  written <- c(
    report("metadata summary report", .generateSummary),
    report("exploratory plots", .generatePlots)
  )
  invisible(written[!is.na(written)])
}

# End-of-run notice pointing at the reports (NULL when none were written).
.reportFilesMessage <- function(report_files) {
  if (!length(report_files)) {
    return(NULL)
  }
  paste0(
    "Metadata reports written to the ORB:\n",
    paste0("  ", report_files, collapse = "\n")
  )
}
