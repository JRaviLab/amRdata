#' Build a protein cluster-to-feature mapping using DuckDB
#'
#' This function constructs a mapping between protein clusters and functional
#' features (e.g., gene families, domains, COGs, and ARGs) using a DuckDB-backed
#' workflow. The implementation is SQL-first and memory-efficient, leveraging
#' DuckDB views and Parquet output.
#'
#' @param parquet_dir Character. Path to the DuckDB database file with parquet file views
#'   containing the input tables.
#' @param output_path Character or NULL. Directory where the output Parquet file
#'   (\code{cluster_feature.parquet}) will be written. If NULL, the output is written
#'   alongside the DuckDB database.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Creates views for each feature type:
#'     \itemize{
#'       \item Gene → protein features
#'       \item Domain annotations
#'       \item COG annotations
#'       \item Antibiotic resistance gene (ARG) annotations
#'     }
#'   \item Combines all feature mappings into a unified protein–feature view
#'   \item Joins protein–feature mappings to cluster membership
#'   \item Writes the resulting cluster–feature mapping to a compressed Parquet file
#' }
#'
#' All joins and transformations are executed inside DuckDB, ensuring scalability
#' for large datasets without loading data into R memory.
#'
#' @return Invisibly returns the file path to the generated Parquet file.
#'
#'
#' @export
buildFeatureHeadMap <- function(
  parquet_dir,
  feature_scales = c("gene", "struct", "protein", "Pfam", "COG", "AMRFinder"),
  output_path = NULL
) {
  # con <- DBI::dbConnect(
  #   duckdb::duckdb(),
  #   normalizePath(parquet_path)
  # )
  # on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  out_dir <- if (is.null(output_path)) {
    dirname(parquet_dir)
  } else {
    normalizePath(output_path)
  }
  # parquet_path <- file.path(out_dir, "cluster_feature.parquet")

  .sql_escape <- function(x) {
  gsub("'", "''", x, fixed = TRUE)
}

.parquet_dataset_sql <- function(parquet_dir, dataset_name) {
  parquet_dir <- normalizePath(parquet_dir, winslash = "/", mustWork = TRUE)
  path <- file.path(parquet_dir, paste0(dataset_name, "*.parquet"))
  .sql_escape(path)
}
  
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  # =========================
  # Gene → protein features
  # =========================

  if("gene" %in% feature_scales) {

    sql_path <- .parquet_dataset_sql(parquet_dir, "genome_gene_protein")

DBI::dbExecute(
  con,
  sprintf(
    "
    CREATE OR REPLACE VIEW v_gene AS
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
  }

  # helper to create the view 
  create_feature_view <- function(con, view_name, parquet_dir,
                                dataset_name, feature_expr) {

  sql_path <- .parquet_dataset_sql(parquet_dir, dataset_name)

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
      ",
      view_name,
      feature_expr,
      sql_path
    )
  )
}

  # =========================
  # Domain features
  # =========================
  
  if ("Pfam" %in% feature_scales) {
  create_feature_view(
    con = con,
    view_name = "v_pfam",
    parquet_dir = parquet_dir,
    dataset_name = "protein_Pfam",
    feature_expr = "query_name"
  )
}

  # =========================
  # COG features
  # =========================

 if ("COG" %in% feature_scales) {
  create_feature_view(
    con = con,
    view_name = "v_cog",
    parquet_dir = parquet_dir,
    dataset_name = "protein_COG",
    feature_expr = "query_name"
  )
}
  # =========================
  # ARG features
  # =========================

  if ("AMRFinder" %in% feature_scales) {
  create_feature_view(
    con = con,
    view_name = "v_amrfinder",
    parquet_dir = parquet_dir,
    dataset_name = "protein_AMRFinder",
    feature_expr = "REPLACE(REPLACE(query_name, '-NCBIFAM', ''), '-', '.')"
  )
}
  
  # =========================
  # Structural gene features
  # (equivalent to separate_rows + inner_join(gp))
  # =========================
  #  DBI::dbExecute(con, "
  #   CREATE OR REPLACE VIEW v_struct AS
  #    SELECT DISTINCT
  #      gp.protein_ids AS protein_id,
  #     s_gene AS feature
  #   FROM struct s
  #   JOIN genome_gene_protein gp
  #      ON gp.genome_ids = s.genome_id
  #   CROSS JOIN UNNEST(string_split(s.struct, '.')) AS t(s_gene)
  #   WHERE s.value = 1
  #  ")
if("struct" %in% feature_scales) {
  sql_path <- .parquet_dataset_sql(parquet_dir, "struct")

DBI::dbExecute(
  con,
  sprintf(
    "
    CREATE OR REPLACE VIEW v_struct_genes AS
    SELECT DISTINCT
      genome_id,
      s_gene
    FROM read_parquet('%s') s
    CROSS JOIN UNNEST(string_split(s.struct, '.')) AS t(s_gene)
    WHERE s.value = 1
    ",
    sql_path
  )
)
}
  

  # =========================
  # Union: protein → feature
  # =========================
  DBI::dbExecute(con, "
    CREATE OR REPLACE VIEW v_protein_feature AS
    SELECT protein_id, feature FROM v_gene
    UNION
    SELECT protein_id, feature FROM v_domain
    UNION
    SELECT protein_id, feature FROM v_cog
    UNION
    SELECT protein_id, feature FROM v_arg
  ")

  # =========================
  # Cluster → feature mapping
  # =========================
  DBI::dbExecute(con, "
    CREATE OR REPLACE VIEW cluster_feature AS
    SELECT DISTINCT
      cm.cluster,
      pf.feature
    FROM protein_members cm
    JOIN v_protein_feature pf
      ON pf.protein_id = cm.member
    WHERE cm.cluster IS NOT NULL
      AND pf.feature IS NOT NULL
  ")

  # =========================
  # Write Parquet from DuckDB
  # =========================
  DBI::dbExecute(
    con,
    sprintf(
      "COPY cluster_feature TO '%s'
       (FORMAT PARQUET, COMPRESSION ZSTD)",
      parquet_path
    )
  )

  invisible(parquet_path)
}
