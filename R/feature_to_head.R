#' Describe HMMER-derived dyad feature views
#'
#' @param databases Character vector of HMMER database names.
#'
#' @return A named list describing view names, feature prefixes, SQL transforms,
#'   and export column names for each database.
#' @keywords internal
.dyadHmmerSpecs <- function(databases) {
  databases <- unique(as.character(databases))

  purrr::set_names(databases) |>
    purrr::map(function(database) {
      known <- switch(
        database,
        Pfam = list(
          key = "pfam",
          prefix = "pfam",
          feature_expr = "REPLACE(query_name, '-', '.')",
          output_column = "Pfam"
        ),
        COG = list(
          key = "cog",
          prefix = "cog",
          feature_expr = "query_name",
          output_column = "COG"
        ),
        AMRFinder = list(
          key = "amr",
          prefix = "amr",
          feature_expr = "REPLACE(REPLACE(query_name, '-NCBIFAM', ''), '-', '.')",
          output_column = "ARG"
        ),
        DefenseCas = list(
          key = "defense",
          prefix = "defense",
          feature_expr = "REPLACE(query_name, '-', '.')",
          output_column = "DefenseCas"
        ),
        NULL
      )

      if (is.null(known)) {
        key <- tolower(gsub("[^A-Za-z0-9]+", "_", database))

        key <- gsub("^_+|_+$", "", key)

        if (!nzchar(key)) {
          key <- "feature"
        }

        if (grepl("^[0-9]", key)) {
          key <- paste0("x_", key)
        }

        known <- list(
          key = key,
          prefix = key,
          feature_expr = "query_name",
          output_column = database
        )
      }

      c(list(database = database),
        known,
        list(
          annotation_view = paste0("dyad_ann_", known$key),
          dyad_view = paste0("dyad_", known$key)
        ))
    })
}

#' Build virtual protein-gene dyad feature relations using DuckDB
#'
#' Creates a compact protein-gene dyad Parquet and registers virtual graph
#' relations in the Parquet-backed ORB DuckDB. Expanded dyad-feature edges are
#' not materialized. Structural and HMMER feature edges are reconstructed on
#' demand from the canonical ORB Parquets.
#'
#' HMMER annotations are propagated from each individual protein through its
#' CD-HIT representative cluster before joining to profile annotations.
#'
#' @param duckdb_path Character. Path to the source per-selection DuckDB. The
#'   associated Parquet files and provenance manifest are expected to live
#'   alongside it. This is the same `duckdb_path` produced by
#'   \code{\link{runDataProcessing}}.
#' @param additional_feature_scales Character vector of optional feature types to
#'   include. If `NULL`, all HMMER databases recorded in the manifest are used,
#'   plus `struct` to include pangenome graph triplets. If a character vector is
#'   supplied, individual feature scales can be selected.
#' @param output_path Character or NULL. Directory where `dyads.parquet` is
#'   written. If NULL, the dataset ORB directory is used.
#' @param threads Integer. Maximum DuckDB threads used while building the compact
#'   dyad dimension. Default `4`.
#'
#' @details
#' The canonical graph representation stores only unique `(protein, gene)` dyads
#' in `dyads.parquet`. The ORB DuckDB then exposes views for `genome_dyad`,
#' `struct_gene`, `dyad_struct`, each configured HMMER feature scale, and the
#' compatibility relation `dyad_feature`.
#'
#' `dyad_feature` retains the historical two-column API (`dyad`, `feature`) but
#' creates those strings only when queried. Protein and gene self-edges are also
#' reconstructed virtually rather than persisted.
#'
#' @return Invisibly returns the path to the generated `dyads.parquet` file.
#'
#' @examples
#' \dontrun{
#' buildDyadFeatureMap(
#'   duckdb_path = "data/Staphylococcus_argenteus/work/Sar.duckdb"
#' )
#' }
#'
#' @import DBI duckdb
#' @export
buildDyadFeatureMap <- function(
    duckdb_path,
    additional_feature_scales = NULL,
    output_path = NULL,
    threads = 4L
) {

  duckdb_path <- normalizePath(duckdb_path, mustWork = TRUE)

  paths <- .amr_paths_from_duckdb(duckdb_path)

  manifest_path <- .manifest_find_latest(duckdb_path)

  if (is.null(manifest_path)) {
    stop(
      "No provenance manifest found for: ",
      duckdb_path
    )
  }

  manifest <- jsonlite::read_json(
    manifest_path,
    simplifyVector = FALSE
  )

  hmmer_stage <- NULL

  for (run in rev(manifest$runs %||% list())) {
    matches <- purrr::keep(
      run$stages %||% list(),
      ~ identical(.x$name, "hmmer") &&
        identical(.x$status, "success")
    )

    if (length(matches)) {
      hmmer_stage <- matches[[length(matches)]]
      break
    }
  }

  if (is.null(hmmer_stage)) {
    stop(
      "No successful HMMER stage found in manifest: ",
      manifest_path
    )
  }

  hmmer_databases <- unique(
    as.character(
      unlist(
        hmmer_stage$parameters$databases %||% character(),
        use.names = FALSE
      )
    )
  )

  if (!length(hmmer_databases)) {
    stop(
      "HMMER stage in manifest does not contain any databases."
    )
  }

  allowed_features <- c(
    "struct",
    hmmer_databases
  )

  if (is.null(additional_feature_scales)) {
    additional_feature_scales <- allowed_features
  } else {
    additional_feature_scales <- unique(
      as.character(
        additional_feature_scales
      )
    )

    unknown_features <- setdiff(
      additional_feature_scales,
      allowed_features
    )

    if (length(unknown_features)) {
      stop(
        "Unsupported feature scale(s): ",
        paste(unknown_features, collapse = ", "),
        ". Available features: ",
        paste(allowed_features, collapse = ", ")
      )
    }
  }

  out_dir <- if (is.null(output_path)) {
    paths$orb
  } else {
    normalizePath(
      output_path,
      winslash = "/",
      mustWork = FALSE
    )
  }

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  out_dir <- normalizePath(
    out_dir,
    winslash = "/",
    mustWork = TRUE
  )

  if (!file.exists(paths$parquet_duckdb)) {
    stop(
      "Parquet-backed ORB DuckDB was not found: ",
      paths$parquet_duckdb,
      ". Run cleanData() before buildDyadFeatureMap()."
    )
  }

  dyads_parquet <- file.path(
    out_dir,
    "dyads.parquet"
  )

  legacy_dyad_parquet <- file.path(
    paths$orb,
    "dyad_feature.parquet"
  )

  con <- .amr_connect_dataset_db(
    paths$parquet_duckdb
  )

  on.exit(
    try(
      DBI::dbDisconnect(con),
      silent = TRUE
    ),
    add = TRUE
  )

  threads <- as.integer(threads)

  if (
    length(threads) != 1L ||
    is.na(threads) ||
    threads < 1L
  ) {
    stop(
      "`threads` must be a positive integer.",
      call. = FALSE
    )
  }

  threads <- min(
    threads,
    as.integer(parallelly::availableCores())
  )

  DBI::dbExecute(
    con,
    paste0(
      "SET threads = ",
      threads
    )
  )

  DBI::dbExecute(
    con,
    "SET preserve_insertion_order = false"
  )

  duckdb_temp_dir <- tempfile(
    pattern = "dyad_duckdb_",
    tmpdir = paths$work
  )

  dir.create(
    duckdb_temp_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  on.exit(
    unlink(
      duckdb_temp_dir,
      recursive = TRUE,
      force = TRUE
    ),
    add = TRUE
  )

  DBI::dbExecute(
    con,
    paste(
      "SET temp_directory =",
      DBI::dbQuoteString(
        con,
        normalizePath(
          duckdb_temp_dir,
          winslash = "/",
          mustWork = TRUE
        )
      )
    )
  )

  available_relations <- DBI::dbListTables(
    con
  )

  required_relations <- c(
    "genome_gene_protein",
    "protein_members"
  )

  missing_required <- setdiff(
    required_relations,
    available_relations
  )

  if (length(missing_required)) {
    stop(
      "Required ORB relation(s) not found: ",
      paste(missing_required, collapse = ", ")
    )
  }

  specs <- .dyadHmmerSpecs(
    hmmer_databases
  )

  stale_views <- c(
    "dyad_feature",
    "dyads",
    "genome_dyad",
    "dyad_protein_cluster",
    "struct_gene",
    "dyad_struct",
    unlist(
      purrr::map(
        specs,
        ~ c(
          .x$annotation_view,
          .x$dyad_view
        )
      ),
      use.names = FALSE
    )
  )

  purrr::walk(
    unique(stale_views),
    function(view_name) {
      DBI::dbExecute(
        con,
        paste(
          "DROP VIEW IF EXISTS",
          DBI::dbQuoteIdentifier(
            con,
            view_name
          )
        )
      )
    }
  )

  if (file.exists(legacy_dyad_parquet)) {
    unlink(
      legacy_dyad_parquet,
      force = TRUE
    )
  }

  if (file.exists(dyads_parquet)) {
    unlink(
      dyads_parquet,
      force = TRUE
    )
  }

  dyads_sql <- DBI::dbQuoteString(
    con,
    normalizePath(
      dyads_parquet,
      winslash = "/",
      mustWork = FALSE
    )
  )

  DBI::dbExecute(
    con,
    paste0(
      "COPY (",
      "SELECT DISTINCT ",
      "CAST(protein_ids AS VARCHAR) AS protein, ",
      "CAST(REPLACE(Gene, '~', '.') AS VARCHAR) AS gene ",
      "FROM genome_gene_protein ",
      "WHERE protein_ids IS NOT NULL ",
      "  AND protein_ids <> '' ",
      "  AND Gene IS NOT NULL ",
      "  AND Gene <> ''",
      ") TO ",
      dyads_sql,
      " (FORMAT PARQUET, COMPRESSION ZSTD, COMPRESSION_LEVEL 3, PRESERVE_ORDER false)"
    )
  )

  dyads_view_path <- if (identical(
    normalizePath(
      out_dir,
      mustWork = TRUE
    ),
    normalizePath(
      paths$orb,
      mustWork = TRUE
    )
  )) {
    basename(dyads_parquet)
  } else {
    normalizePath(
      dyads_parquet,
      winslash = "/",
      mustWork = TRUE
    )
  }

  DBI::dbExecute(
    con,
    paste0(
      "CREATE OR REPLACE VIEW dyads AS ",
      "SELECT * FROM read_parquet(",
      DBI::dbQuoteString(
        con,
        dyads_view_path
      ),
      ")"
    )
  )

  DBI::dbExecute(
    con,
    "CREATE OR REPLACE VIEW genome_dyad AS
     SELECT DISTINCT
       CAST(genome_ids AS VARCHAR) AS genome_id,
       CAST(protein_ids AS VARCHAR) AS protein,
       CAST(REPLACE(Gene, '~', '.') AS VARCHAR) AS gene
     FROM genome_gene_protein
     WHERE genome_ids IS NOT NULL
       AND genome_ids <> ''
       AND protein_ids IS NOT NULL
       AND protein_ids <> ''
       AND Gene IS NOT NULL
       AND Gene <> ''"
  )

  DBI::dbExecute(
    con,
    "CREATE OR REPLACE VIEW dyad_protein_cluster AS
     SELECT DISTINCT
       CAST(member AS VARCHAR) AS protein,
       CAST(cluster AS VARCHAR) AS cluster
     FROM protein_members
     WHERE member IS NOT NULL
       AND member <> ''
       AND cluster IS NOT NULL
       AND cluster <> ''"
  )

  ambiguous_members <- DBI::dbGetQuery(
    con,
    "SELECT count(*) AS n
     FROM (
       SELECT protein
       FROM dyad_protein_cluster
       GROUP BY protein
       HAVING count(DISTINCT cluster) > 1
     )"
  )$n[[1L]]

  if (ambiguous_members > 0L) {
    stop(
      "CD-HIT membership is ambiguous for ",
      ambiguous_members,
      " protein(s); each protein must map to at most one representative cluster."
    )
  }

  created_edge_views <- character()
  compatibility_parts <- c(
    "SELECT protein || '|' || gene AS dyad, 'protein:' || protein AS feature FROM dyads",
    "SELECT protein || '|' || gene AS dyad, 'gene:' || gene AS feature FROM dyads"
  )

  if (
    "struct" %in% additional_feature_scales &&
    "struct" %in% available_relations
  ) {
    DBI::dbExecute(
      con,
      "CREATE OR REPLACE VIEW struct_gene AS
       WITH unique_structs AS (
         SELECT DISTINCT CAST(struct AS VARCHAR) AS struct
         FROM struct
         WHERE value = 1
       )
       SELECT DISTINCT
         s.struct,
         CAST(REPLACE(gene_raw, '~', '.') AS VARCHAR) AS gene
       FROM unique_structs s
       CROSS JOIN UNNEST(
         string_split(REPLACE(s.struct, '.', '-'), '-')
       ) AS t(gene_raw)"
    )

    DBI::dbExecute(
      con,
      "CREATE OR REPLACE VIEW dyad_struct AS
       SELECT
         d.protein,
         d.gene,
         sg.struct AS feature
       FROM dyads d
       JOIN struct_gene sg
         ON d.gene = sg.gene"
    )

    created_edge_views <- c(
      created_edge_views,
      "dyad_struct"
    )

    compatibility_parts <- c(
      compatibility_parts,
      "SELECT protein || '|' || gene AS dyad, 'struct:' || feature AS feature FROM dyad_struct"
    )
  }

  for (database in intersect(
    hmmer_databases,
    additional_feature_scales
  )) {
    spec <- specs[[database]]
    source_relation <- paste0(
      "protein_",
      database
    )

    if (!source_relation %in% available_relations) {
      message(
        "Skipping ",
        database,
        ": ",
        source_relation,
        " was not found in the ORB."
      )
      next
    }

    annotation_view_sql <- DBI::dbQuoteIdentifier(
      con,
      spec$annotation_view
    )
    dyad_view_sql <- DBI::dbQuoteIdentifier(
      con,
      spec$dyad_view
    )
    source_relation_sql <- DBI::dbQuoteIdentifier(
      con,
      source_relation
    )

    DBI::dbExecute(
      con,
      paste0(
        "CREATE OR REPLACE VIEW ",
        annotation_view_sql,
        " AS ",
        "SELECT DISTINCT ",
        "CAST(protein AS VARCHAR) AS cluster, ",
        "CAST(",
        spec$feature_expr,
        " AS VARCHAR) AS feature ",
        "FROM ",
        source_relation_sql,
        " WHERE protein IS NOT NULL ",
        "   AND query_name IS NOT NULL"
      )
    )

    DBI::dbExecute(
      con,
      paste0(
        "CREATE OR REPLACE VIEW ",
        dyad_view_sql,
        " AS ",
        "SELECT d.protein, d.gene, a.feature ",
        "FROM dyads d ",
        "JOIN dyad_protein_cluster pc ",
        "  ON d.protein = pc.protein ",
        "JOIN ",
        annotation_view_sql,
        " a ON pc.cluster = a.cluster"
      )
    )

    created_edge_views <- c(
      created_edge_views,
      spec$dyad_view
    )

    compatibility_parts <- c(
      compatibility_parts,
      paste0(
        "SELECT protein || '|' || gene AS dyad, '",
        spec$prefix,
        ":' || feature AS feature FROM ",
        dyad_view_sql
      )
    )
  }

  DBI::dbExecute(
    con,
    paste0(
      "CREATE OR REPLACE VIEW dyad_feature AS ",
      paste(
        compatibility_parts,
        collapse = " UNION ALL "
      )
    )
  )

  if (length(created_edge_views)) {
    available_relations <- unique(c(
      available_relations,
      created_edge_views
    ))
  }

  missing_proteins <- DBI::dbGetQuery(
    con,
    "SELECT count(DISTINCT d.protein) AS n
     FROM dyads d
     LEFT JOIN dyad_protein_cluster pc
       ON d.protein = pc.protein
     WHERE pc.protein IS NULL"
  )$n[[1L]]

  if (missing_proteins > 0L) {
    message(
      "Dyad/CD-HIT mapping: ",
      missing_proteins,
      " protein(s) have no protein_members mapping and therefore cannot inherit HMMER annotations."
    )
  }

  invisible(
    dyads_parquet
  )
}
