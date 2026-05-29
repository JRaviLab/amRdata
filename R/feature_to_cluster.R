#' Build a protein cluster-to-feature mapping using DuckDB
#'
#' This function constructs a mapping between protein clusters and functional
#' features (e.g., gene families, domains, COGs, and ARGs) using a DuckDB-backed
#' workflow. The implementation is SQL-first and memory-efficient, leveraging
#' DuckDB views and Parquet output.
#'
#' @param duckdb_parquet_path Character. Path to the DuckDB database file with parquet file views
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
#' @examples
#' \dontrun{
#' buildClusterFeatureMap(
#'   duckdb_parquet_path = "data/Csp_parquet.duckdb",
#'   output_path = NULL
#' )
#' }
#'
#' @import DBI duckdb
#' @export
buildClusterFeatureMap <- function(
  duckdb_parquet_path,
  output_path = NULL
) {

  con <- DBI::dbConnect(
    duckdb::duckdb(),
    normalizePath(duckdb_parquet_path)
  )
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  out_dir <- if (is.null(output_path)) {
    dirname(duckdb_parquet_path)
  } else {
    normalizePath(output_path)
  }
  parquet_path <- file.path(out_dir, "cluster_feature.parquet")

  # =========================
  # Gene → protein features
  # =========================
  DBI::dbExecute(con, "
    CREATE OR REPLACE VIEW v_gene AS
    SELECT DISTINCT
      protein_ids AS protein_id,
      REPLACE(Gene, '~', '.') AS feature
    FROM genome_gene_protein
    WHERE protein_ids IS NOT NULL
      AND Gene IS NOT NULL
  ")

  # =========================
  # Domain features
  # =========================
  DBI::dbExecute(con, "
    CREATE OR REPLACE VIEW v_domain AS
    SELECT DISTINCT
      AccNum AS protein_id,
      \"DB.ID\" AS feature
    FROM domain_names
    WHERE AccNum IS NOT NULL
      AND \"DB.ID\" IS NOT NULL
  ")

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

  # =========================
  # COG features
  # =========================
  DBI::dbExecute(con, "
    CREATE OR REPLACE VIEW v_cog AS
    SELECT DISTINCT
      query_name AS protein_id,
      name AS feature
    FROM protein_COG
    WHERE query_name IS NOT NULL
      AND name IS NOT NULL
  ")

  # =========================
  # ARG features
  # =========================
  DBI::dbExecute(con, "
    CREATE OR REPLACE VIEW v_arg AS
    SELECT DISTINCT
      query_name AS protein_id,
      REPLACE(REPLACE(name, '-NCBIFAM', ''), '-', '.') AS feature
    FROM protein_ResFinder
    WHERE query_name IS NOT NULL
      AND name IS NOT NULL
  ")

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
