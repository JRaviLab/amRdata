#' Build a protein-gene dyad feature network using DuckDB
#'
#' Constructs a two-mode (bipartite) network in which every node is either a
#' protein-gene dyad or a biological feature, and every edge links a dyad to one
#' of its features. Features include the dyad's own protein and gene, structural
#' (pangenome graph) gene arrangements, and HMMER annotations: protein domains
#' (Pfam), COGs, phage-defense systems (including Cas), and antimicrobial
#' resistance genes.
#'
#' @param duckdb_path Character. Path to the source per-selection DuckDB. The
#'   associated Parquet files and provenance manifest are expected to live
#'   alongside it. This is the same `duckdb_path` produced by
#'   \code{\link{runDataProcessing}}.
#' @param additional_feature_scales Character vector of optional feature types to
#'   include. If `NULL`, all HMMER databases recorded in the manifest are used,
#'   plus `struct` to include pangenome graph triplets. If a character vector is
#'   supplied, individual feature scales can be selected.
#' @param output_path Character or NULL. Directory where the output Parquet file
#'   will be written. If NULL, the output is written alongside the source DuckDB.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Creates a protein-gene mapping from
#'     \code{genome_gene_protein.parquet}.
#'   \item Constructs a unique protein-gene dyad identifier of the form
#'     \code{"protein|gene"}.
#'   \item Optionally loads structural and HMMER annotations.
#'   \item Generates a network edge list linking each dyad to its associated
#'     features.
#'   \item Writes the resulting edge list to a compressed Parquet file.
#' }
#'
#' Every edge runs from a `protein|gene` dyad to a type-prefixed feature node:
#'
#' \preformatted{
#' protein|gene --> protein:PROTEIN_ID
#' protein|gene --> gene:GENE_ID
#' protein|gene --> pfam:PFXXXXX
#' protein|gene --> cog:COGXXXX
#' protein|gene --> amr:GENE_NAME
#' protein|gene --> defensecas:DEFENSE_SYSTEM
#' protein|gene --> struct:STRUCTURE
#' }
#'
#' The type prefix on each target keeps feature namespaces separate, so a
#' generic feature name cannot collide across scales downstream.
#'
#' The output edge list contains two columns:
#' \describe{
#'   \item{source}{Protein-gene dyad identifier.}
#'   \item{target}{Feature node identifier prefixed by feature type.}
#' }
#'
#' @return Invisibly returns the path to the generated
#'   \code{dyad_feature.parquet} file.
#'
#' @examples
#' \dontrun{
#' buildDyadFeatureMap(
#'   duckdb_path = "data/Staphylococcus_argenteus/Sar.duckdb"
#' )
#' }
#'
#' @import DBI duckdb
#' @export
buildDyadFeatureMap <- function(
    duckdb_path,
    additional_feature_scales = NULL,
    output_path = NULL
) {

  duckdb_path <- normalizePath(
    duckdb_path,
    mustWork = TRUE
  )

  parquet_dir <- dirname(duckdb_path)

  manifest_path <- .manifest_find_latest(
    duckdb_path
  )

  if (is.null(manifest_path)) {
    stop(
      "No provenance manifest found for: ",
      duckdb_path
    )
  }

  # Read in the `runDataProcessing()` output manifest
  manifest <- jsonlite::read_json(
    manifest_path,
    simplifyVector = FALSE
  )

  hmmer_stage <- NULL

  # Find the most recent run whose HMMER stage completed successfully; its
  # recorded databases determine which feature scales are available.
  for (run in rev(manifest$runs)) {

    stages <- run$stages %||% list()

    matches <- purrr::keep(
      stages,
      ~ identical(.x$name, "hmmer") &&
        identical(.x$status, "success")
    )

    if (length(matches)) {
      hmmer_stage <- matches[[1]]
      break
    }
  }

  # Without a successful HMMER stage the annotation Parquets we join against
  # are not guaranteed to exist, so there is nothing to build from.
  if (is.null(hmmer_stage)) {
    stop(
      "No successful HMMER stage found in manifest: ",
      manifest_path
    )
  }

  hmmer_databases <- unlist(
    hmmer_stage$parameters$databases
  )

  # A successful stage with no databases recorded should not happen; fail loudly.
  if (!length(hmmer_databases)) {
    stop(
      "HMMER stage in manifest does not contain any databases."
    )
  }

  # Feature scales the caller is allowed to request: the structural view plus
  # every HMMER database from the manifest (so new databases are picked up
  # automatically).
  allowed_features <- c(
    "struct",
    hmmer_databases
  )

  # Default to every allowed scale; otherwise reject anything unrecognised.
  if (is.null(additional_feature_scales)) {
    additional_feature_scales <- allowed_features
  } else {
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
    parquet_dir
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

  parquet_path <- file.path(
    out_dir,
    "dyad_feature.parquet"
  )

  con <- DBI::dbConnect(
    duckdb::duckdb(),
    dbdir = ":memory:"
  )

  on.exit(
    DBI::dbDisconnect(
      con,
      shutdown = TRUE
    ),
    add = TRUE
  )

  # Local helpers -------------------------------------------------------------

  # Escape single quotes so a path can be embedded in a SQL string literal.
  .sql_escape <- function(x) {
    gsub(
      "'",
      "''",
      x,
      fixed = TRUE
    )
  }

  # Resolve "<parquet_dir>/<dataset_name>.parquet" to an escaped absolute path,
  # erroring if the expected Parquet file is missing.
  .parquet_dataset_sql <- function(
    parquet_dir,
    dataset_name
  ) {

    path <- file.path(
      parquet_dir,
      paste0(
        dataset_name,
        ".parquet"
      )
    )

    if (!file.exists(path)) {
      stop(
        "Required Parquet file not found: ",
        path
      )
    }

    .sql_escape(
      normalizePath(
        path,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  # =========================
  # Gene -> protein mapping
  # =========================

  # Base protein/gene pairs. Empty strings are excluded alongside NULLs to match
  # the `value != ""` convention in data_processing.R and to avoid emitting
  # bogus "protein|" or "|gene" dyads.
  sql_path <- .parquet_dataset_sql(
    parquet_dir,
    "genome_gene_protein"
  )

  DBI::dbExecute(
    con,
    sprintf(
      "
      CREATE OR REPLACE VIEW protein_gene AS
      SELECT DISTINCT
        protein_ids AS protein,
        REPLACE(Gene, '~', '.') AS gene
      FROM read_parquet('%s')
      WHERE protein_ids IS NOT NULL
        AND protein_ids <> ''
        AND Gene IS NOT NULL
        AND Gene <> ''
      ",
      sql_path
    )
  )

  # =========================
  # Protein-gene -> dyad
  # =========================

  # Collapse each protein/gene pair into a single "protein|gene" dyad id; this
  # is the source node for every edge in the output network.
  DBI::dbExecute(
    con,
    "
    CREATE OR REPLACE VIEW protein_gene_dyad AS
    SELECT DISTINCT
      protein,
      gene,
      CONCAT(protein, '|', gene) AS dyad
    FROM protein_gene
    "
  )

  # Track which optional feature views were successfully created, so the edge
  # queries below only join against views that exist.
  has_struct <- FALSE
  has_pfam <- FALSE
  has_cog <- FALSE
  has_amr <- FALSE
  has_defensecas <- FALSE

  # =========================
  # Structural gene features
  # =========================

  if ("struct" %in% additional_feature_scales) {

    struct_path <- file.path(
      parquet_dir,
      "struct.parquet"
    )

    if (!file.exists(struct_path)) {

      # Same skip-and-continue behaviour as .create_feature_view() uses for the
      # HMMER datasets, kept inline here because the struct view is built with a
      # different (UNNEST) query.
      message(
        "Skipping struct: no parquet found. Generate struct parquet first."
      )

    } else {

      sql_path <- .parquet_dataset_sql(
        parquet_dir,
        "struct"
      )

      DBI::dbExecute(
        con,
        sprintf(
          "
          CREATE OR REPLACE VIEW v_struct_genes AS
          SELECT DISTINCT
            struct,
            gene
          FROM read_parquet('%s') s
          CROSS JOIN UNNEST(
            string_split(s.struct, '.')
          ) AS t(gene)
          WHERE s.value = 1
          ",
          sql_path
        )
      )

      has_struct <- TRUE
    }
  }

  # ==========================
  # Generic HMMER feature view
  # ==========================

  # Local helper: build a `protein -> feature` view from one HMMER annotation
  # Parquet. `feature_expr` is the SQL expression that derives the feature label
  # from `query_name` (identity for most databases, a REPLACE() for AMRFinder).
  # Returns TRUE if the view was created, FALSE if the Parquet was missing.
  .create_feature_view <- function(
    con,
    view_name,
    parquet_dir,
    dataset_name,
    feature_expr
  ) {

    sql_path <- tryCatch(
      .parquet_dataset_sql(
        parquet_dir,
        dataset_name
      ),
      error = function(e) {
        message(
          "Skipping ",
          dataset_name,
          ": ",
          e$message
        )
        NULL
      }
    )

    if (is.null(sql_path)) {
      return(FALSE)
    }

    DBI::dbExecute(
      con,
      sprintf(
        "
        CREATE OR REPLACE VIEW %s AS
        SELECT DISTINCT
          protein,
          %s AS feature
        FROM read_parquet('%s')
        WHERE protein IS NOT NULL
          AND query_name IS NOT NULL
        ",
        view_name,
        feature_expr,
        sql_path
      )
    )

    TRUE
  }

  # =========================
  # HMMER annotation features
  # =========================

  if ("Pfam" %in% additional_feature_scales) {

    has_pfam <- .create_feature_view(
      con = con,
      view_name = "v_pfam",
      parquet_dir = parquet_dir,
      dataset_name = "protein_Pfam",
      feature_expr = "query_name"
    )
  }

  if ("COG" %in% additional_feature_scales) {

    has_cog <- .create_feature_view(
      con = con,
      view_name = "v_cog",
      parquet_dir = parquet_dir,
      dataset_name = "protein_COG",
      feature_expr = "query_name"
    )
  }

  if ("AMRFinder" %in% additional_feature_scales) {

    has_amr <- .create_feature_view(
      con = con,
      view_name = "v_amrfinder",
      parquet_dir = parquet_dir,
      dataset_name = "protein_AMRFinder",
      feature_expr = "REPLACE(REPLACE(query_name, '-NCBIFAM', ''), '-', '.')"
    )
  }

  if ("DefenseCas" %in% additional_feature_scales) {

    has_defensecas <- .create_feature_view(
      con = con,
      view_name = "v_defensecas",
      parquet_dir = parquet_dir,
      dataset_name = "protein_DefenseCas",
      feature_expr = "query_name"
    )
  }

  # =========================
  # Build network edge queries
  # =========================

  # Each entry is a SELECT returning (source, target) rows that are UNIONed into
  # the final edge list. In the joined queries below `pgd` aliases the
  # `protein_gene_dyad` view, so `pgd.dyad` is the "protein|gene" dyad id.
  edge_queries <- c(
    "
    SELECT DISTINCT
      dyad AS source,
      CONCAT('protein:', protein) AS target
    FROM protein_gene_dyad
    ",
    "
    SELECT DISTINCT
      dyad AS source,
      CONCAT('gene:', gene) AS target
    FROM protein_gene_dyad
    "
  )

  # =========================
  # Structure edges
  # =========================

  if (has_struct) {

    edge_queries <- c(
      edge_queries,
      "
      SELECT DISTINCT
        pgd.dyad AS source,
        CONCAT('struct:', sg.struct) AS target
      FROM protein_gene_dyad pgd
      JOIN v_struct_genes sg
        ON pgd.gene = sg.gene
      "
    )
  }

  # =========================
  # Pfam edges
  # =========================

  if (has_pfam) {

    edge_queries <- c(
      edge_queries,
      "
      SELECT DISTINCT
        pgd.dyad AS source,
        CONCAT('pfam:', pf.feature) AS target
      FROM protein_gene_dyad pgd
      JOIN v_pfam pf
        ON pgd.protein = pf.protein
      "
    )
  }

  # =========================
  # COG edges
  # =========================

  if (has_cog) {

    edge_queries <- c(
      edge_queries,
      "
      SELECT DISTINCT
        pgd.dyad AS source,
        CONCAT('cog:', cf.feature) AS target
      FROM protein_gene_dyad pgd
      JOIN v_cog cf
        ON pgd.protein = cf.protein
      "
    )
  }

  # =========================
  # AMRFinder edges
  # =========================

  if (has_amr) {

    edge_queries <- c(
      edge_queries,
      "
      SELECT DISTINCT
        pgd.dyad AS source,
        CONCAT('amr:', af.feature) AS target
      FROM protein_gene_dyad pgd
      JOIN v_amrfinder af
        ON pgd.protein = af.protein
      "
    )
  }

  # =========================
  # DefenseCas edges
  # =========================

  if (has_defensecas) {

    edge_queries <- c(
      edge_queries,
      "
      SELECT DISTINCT
        pgd.dyad AS source,
        CONCAT('defensecas:', df.feature) AS target
      FROM protein_gene_dyad pgd
      JOIN v_defensecas df
        ON pgd.protein = df.protein
      "
    )
  }

  # =========================
  # Final edge list
  # =========================

  DBI::dbExecute(
    con,
    paste0(
      "
      CREATE OR REPLACE VIEW network_edges AS
      ",
      paste(
        edge_queries,
        collapse = "\nUNION\n"
      )
    )
  )

  # =========================
  # Export Parquet
  # =========================

  parquet_sql <- DBI::dbQuoteString(
    con,
    normalizePath(
      parquet_path,
      winslash = "/",
      mustWork = FALSE
    )
  )

  DBI::dbExecute(
    con,
    paste0(
      "
      COPY network_edges
      TO ",
      parquet_sql,
      "
      (FORMAT PARQUET, COMPRESSION ZSTD)
      "
    )
  )

  invisible(parquet_path)
}
