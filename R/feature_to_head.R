#' Build a protein-gene dyad feature network using DuckDB
#'
#' Constructs a bipartite network linking protein-gene dyads to biological
#' features such as proteins, genes, structural gene arrangements, and HMMER
#' annotations for protein domains (Pfam), COGs, Defense including Cas and antimicrobial resistance genes.
#'
#' @param duckdb_path Character. Path to the source dataset DuckDB. The
#'   associated Parquet files and provenance manifest are expected to live
#'   alongside it.
#' @param addtnl_feature_scales Character vector of optional feature types to
#'   include. If `NULL`, all HMMER databases recorded in the manifest are used,
#'   plus `struct` to include pangenome graph triplets. If a character vector is
#'   supplied, individual feature scales can be selected.
#' @param output_path Character or NULL. Directory where the output Parquet file
#'   will be written. If NULL, the output is written alongside the source DuckDB.
#'
#'   #' @details
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
#' The resulting network is bipartite:
#'
#' \preformatted{
#' protein|gene --> protein:PROTEIN_ID
#' protein|gene --> gene:GENE_ID
#' protein|gene --> pfam:PFXXXXX
#' protein|gene --> cog:COGXXXX
#' protein|gene --> amr:GENE_NAME
#' protein|gene --> defensecas:FEATURE
#' protein|gene --> struct:STRUCTURE
#' }
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
#' @export
buildDyadFeatureMap <- function(
    duckdb_path,
    addtnl_feature_scales = NULL,
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

  # Learn what databases were used, and used successfully (pretty important)
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

  # Looks like stuff didn't work? Stop!
  if (is.null(hmmer_stage)) {
    stop(
      "No successful HMMER stage found in manifest: ",
      manifest_path
    )
  }

  hmmer_databases <- unlist(
    hmmer_stage$parameters$databases
  )

  # HMMER didn't run or something? Stop!
  if (!length(hmmer_databases)) {
    stop(
      "HMMER stage in manifest does not contain any databases."
    )
  }

  # Should probably list what's allowed to input as a feature
  allowed_features <- c(
    "struct",
    hmmer_databases # Should be flexible enough to allow for adding HMMER DBs
  )

  # Setdiff to find what was provided vs. what's allowed
  if (is.null(addtnl_feature_scales)) {
    addtnl_feature_scales <- allowed_features
  } else {
    unknown_features <- setdiff(
      addtnl_feature_scales,
      allowed_features
    )

    # If nonsense exists? Shout about it
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

  .sql_escape <- function(x) {
    gsub(
      "'",
      "''",
      x,
      fixed = TRUE
    )
  }

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
  # Gene -> protein features
  # =========================

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
        AND Gene IS NOT NULL
      ",
      sql_path
    )
  )

  # =========================
  # Protein-gene -> dyad
  # =========================

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

  # Feature flags
  has_struct <- FALSE
  has_pfam <- FALSE
  has_cog <- FALSE
  has_amr <- FALSE
  has_defensecas <- FALSE

  # =========================
  # Structural gene features
  # =========================

  if ("struct" %in% addtnl_feature_scales) {

    struct_path <- file.path(
      parquet_dir,
      "struct.parquet"
    )

    if (!file.exists(struct_path)) {

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

  create_feature_view <- function(
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

  if ("Pfam" %in% addtnl_feature_scales) {

    has_pfam <- create_feature_view(
      con = con,
      view_name = "v_pfam",
      parquet_dir = parquet_dir,
      dataset_name = "protein_Pfam",
      feature_expr = "query_name"
    )
  }

  if ("COG" %in% addtnl_feature_scales) {

    has_cog <- create_feature_view(
      con = con,
      view_name = "v_cog",
      parquet_dir = parquet_dir,
      dataset_name = "protein_COG",
      feature_expr = "query_name"
    )
  }

  if ("AMRFinder" %in% addtnl_feature_scales) {

    has_amr <- create_feature_view(
      con = con,
      view_name = "v_amrfinder",
      parquet_dir = parquet_dir,
      dataset_name = "protein_AMRFinder",
      feature_expr = "REPLACE(REPLACE(query_name, '-NCBIFAM', ''), '-', '.')"
    )
  }

  if ("DefenseCas" %in% addtnl_feature_scales) {

    has_defensecas <- create_feature_view(
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
