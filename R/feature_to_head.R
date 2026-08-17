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
  addtnl_feature_scales = c("struct", "Pfam", "COG", "AMRFinder"),
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
  path <- file.path(parquet_dir, paste0(dataset_name, ".parquet"))
  .sql_escape(path)
}
  
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  
  # =========================
  # Gene → protein features
  # =========================

sql_path <- .parquet_dataset_sql(parquet_dir, "genome_gene_protein")

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
# protein-gene --> dyad
# =========================

DBI::dbExecute(con, "
CREATE OR REPLACE VIEW protein_gene_dyad AS
SELECT DISTINCT
    protein,
    gene,
    CONCAT(protein, '|', gene) AS dyad
FROM protein_gene
")

# =========================
  # Structural gene features
  # =========================
 
if("struct" %in% addtnl_feature_scales) {
  sql_path <- .parquet_dataset_sql(parquet_dir, "struct")

DBI::dbExecute(
  con,
  sprintf(
    "
    CREATE OR REPLACE VIEW v_struct_genes AS
    SELECT DISTINCT
      struct,
      gene
    FROM read_parquet('%s') s
    CROSS JOIN UNNEST(string_split(s.struct, '.')) AS t(gene)
    WHERE s.value = 1
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
  
  if ("Pfam" %in% addtnl_feature_scales) {
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

 if ("COG" %in% addtnl_feature_scales) {
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

  if ("AMRFinder" %in% addtnl_feature_scales) {
  create_feature_view(
    con = con,
    view_name = "v_amrfinder",
    parquet_dir = parquet_dir,
    dataset_name = "protein_AMRFinder",
    feature_expr = "REPLACE(REPLACE(query_name, '-NCBIFAM', ''), '-', '.')"
  )
}
# ==================================================
# Build network edge queries
# ==================================================

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

# ==================================================
# Structure edges
# ==================================================

if ("struct" %in% addtnl_feature_scales) {

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

# ==================================================
# Pfam edges
# ==================================================

if ("Pfam" %in% addtnl_feature_scales) {

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

# ==================================================
# COG edges
# ==================================================

if ("COG" %in% addtnl_feature_scales) {

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

# ==================================================
# AMRFinder edges
# ==================================================

if ("AMRFinder" %in% addtnl_feature_scales) {

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

# ==================================================
# Final edge list
# ==================================================

DBI::dbExecute(
  con,
  paste0(
    "
    CREATE OR REPLACE VIEW network_edges AS
    ",
    paste(edge_queries, collapse = "\nUNION\n")
  )
)

# ==================================================
# Export parquet
# ==================================================

parquet_path <- file.path(out_dir, "dyad_feature.parquet")

DBI::dbExecute(
  con,
  sprintf(
    "
    COPY network_edges
    TO '%s'
    (FORMAT PARQUET, COMPRESSION ZSTD)
    ",
    parquet_path
  )
)

DBI::dbDisconnect(con, shutdown = TRUE)

invisible(parquet_path)
}
