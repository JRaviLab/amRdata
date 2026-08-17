#' Build a protein-gene dyad feature network using DuckDB
#'
#' Constructs a bipartite network linking protein-gene dyads to biological
#' features such as proteins, genes, structural gene arrangements, protein
#' domains (Pfam), COG annotations, and antimicrobial resistance annotations.
#'
#' The workflow is implemented entirely within DuckDB using SQL views over
#' Parquet datasets, minimizing memory usage and enabling scalable processing
#' of large genomic datasets.
#'
#' @param parquet_dir Character. Directory containing input Parquet datasets.
#'   At minimum, the directory must contain
#'   \code{genome_gene_protein.parquet}.
#' @param addtnl_feature_scales Character vector of optional feature types to
#'   include. Supported values are:
#'   \describe{
#'     \item{\code{"struct"}}{Structural gene arrangements.}
#'     \item{\code{"Pfam"}}{Pfam protein domain annotations.}
#'     \item{\code{"COG"}}{Clusters of Orthologous Groups annotations.}
#'     \item{\code{"AMRFinder"}}{AMRFinder antimicrobial resistance annotations.}
#'   }
#' @param output_path Character or NULL. Directory where the output Parquet file
#'   will be written. If NULL, the output is written to the parent directory of
#'   \code{parquet_dir}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Creates a protein-gene mapping from
#'     \code{genome_gene_protein.parquet}.
#'   \item Constructs a unique protein-gene dyad identifier of the form
#'     \code{"protein|gene"}.
#'   \item Optionally loads structural, Pfam, COG, and AMRFinder annotations.
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
#'   parquet_dir = "data/parquet",
#'   addtnl_feature_scales = c(
#'     "struct",
#'     "Pfam",
#'     "COG",
#'     "AMRFinder"
#'   )
#' )
#' }
#'
#' @export
buildDyadFeatureMap <- function(
  parquet_dir,
  addtnl_feature_scales = c("struct", "Pfam", "COG", "AMRFinder"),
  output_path = NULL
) {
  # con <- DBI::dbConnect(
  #   duckdb::duckdb(),
  #   normalizePath(parquet_path)
  # )
  # on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # out_dir <- if (is.null(output_path)) {
  #   dirname(parquet_dir)
  # } else {
  #   normalizePath(output_path)
  # }
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
 
if ("struct" %in% addtnl_feature_scales) {

  struct_files <- Sys.glob(
    file.path(parquet_dir, "struct*.parquet")
  )

  if (length(struct_files) == 0) {

    message(
      "Skipping struct: no parquet found. Generate struct parquet first."
    )

  } else {

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

  # helper to create the view 
  create_feature_view <- function(
    con,
    view_name,
    parquet_dir,
    dataset_name,
    feature_expr
) {

  parquet_files <- Sys.glob(
    file.path(parquet_dir, paste0(dataset_name, "*.parquet"))
  )

  if (length(parquet_files) == 0) {

    message(
      sprintf(
        "# Skipping %s: no parquet found. Generate %s parquet first.",
        view_name,
        dataset_name
      )
    )

    return(FALSE)
  }

  sql_path <- .parquet_dataset_sql(
    parquet_dir,
    dataset_name
  )

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

  TRUE
}

  # =========================
  # Domain features
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

  # =========================
  # COG features
  # =========================

 if ("COG" %in% addtnl_feature_scales) {
  has_cog <- create_feature_view(
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
  has_amr <- create_feature_view(
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

# ==================================================
# Pfam edges
# ==================================================

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

# ==================================================
# COG edges
# ==================================================

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

# ==================================================
# AMRFinder edges
# ==================================================

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

parquet_path <- file.path(parquet_dir, "dyad_feature.parquet")

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
