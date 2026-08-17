#!/usr/bin/env Rscript
# BV-BRC Data API prototype for amRdata
# -------------------------------------
# Exploratory/benchmark script — NOT wired into the package namespace.
# Demonstrates the three things that fix the current Docker bottleneck:
#   1. bvbrc_count()      -- total matches WITHOUT downloading (Content-Range)
#   2. bvbrc_fetch_all()  -- keyset pagination past the 25k cap, optional
#                            id-space partitioning across a parallel worker pool
#   3. bvbrc_rank_genomes_by_drug_coverage() -- faceted subsampling
#
# Design rationale + limitations: docs/bvbrc-api-feasibility.md
# Requires: httr2, jsonlite, data.table  (+ future, future.apply for parallel)

suppressPackageStartupMessages({
  library(httr2)
  library(jsonlite)
})

BVBRC_BASE <- "https://www.bv-brc.org/api"
PAGE_MAX   <- 25000L # hard per-request AND offset ceiling

# --- helpers ----------------------------------------------------------------

# URL-encode ONE RQL value (BV-BRC wants field names/values encoded individually)
.enc <- function(x) utils::URLencode(as.character(x), reserved = TRUE)

.bvbrc_req <- function(collection, rql, accept = "application/json") {
  url <- sprintf("%s/%s/?%s", BVBRC_BASE, collection, rql)
  request(url) |>
    req_headers(Accept = accept) |>
    req_timeout(120) |>
    req_retry(
      max_tries = 4L,
      is_transient = function(resp) resp_status(resp) %in% c(429, 500, 502, 503, 504)
    )
}

# --- 1. count() -------------------------------------------------------------
# Total matching rows without downloading them. Reads "items 0-1/<N>".
bvbrc_count <- function(collection, rql_filter) {
  resp <- .bvbrc_req(collection, paste0(rql_filter, "&limit(1)")) |> req_perform()
  as.integer(sub(".*/", "", resp_header(resp, "Content-Range")))
}

# --- one keyset page --------------------------------------------------------
# `key` is the unique sort/seek field: "id" for genome_amr/feature/sequence,
# "genome_id" for the genome collection (which has no "id" field).
.bvbrc_page <- function(collection, rql_filter, select, last_id, page_size,
                        key = "id") {
  seek <- if (nzchar(last_id)) {
    sprintf("and(%s,gt(%s,%s))", rql_filter, key, .enc(last_id))
  } else {
    rql_filter
  }
  rql <- sprintf("%s&select(%s,%s)&sort(+%s)&limit(%d)",
                 seek, key, select, key, page_size)
  resp <- .bvbrc_req(collection, rql) |> req_perform()
  fromJSON(resp_body_string(resp), simplifyDataFrame = TRUE)
}

# Walk one result set to exhaustion via keyset (seek) pagination.
# gt(key, X) is INCLUSIVE on this API, so we dedupe the boundary row and judge
# termination on the RAW fetched size, not the post-dedupe size.
.keyset_walk <- function(collection, rql_filter, select, page_size = PAGE_MAX,
                         key = "id") {
  out <- list(); last_id <- ""; k <- 0L
  repeat {
    pg <- .bvbrc_page(collection, rql_filter, select, last_id, page_size, key)
    raw_n <- if (is.data.frame(pg)) nrow(pg) else 0L
    if (raw_n == 0L) break
    if (nzchar(last_id)) pg <- pg[pg[[key]] != last_id, , drop = FALSE]
    if (nrow(pg) > 0L) {
      k <- k + 1L; out[[k]] <- pg
      last_id <- pg[[key]][nrow(pg)]
    }
    if (raw_n < page_size) break
  }
  if (k == 0L) return(NULL)
  data.table::rbindlist(out, fill = TRUE, use.names = TRUE)
}

# --- 2. fetch_all() ---------------------------------------------------------
# Pull an entire result set of ANY size.
#   parallel = FALSE : single keyset walk (simple, good for interactive subsets)
#   parallel = TRUE  : partition the id space into hex buckets and keyset-walk
#                      each in parallel. UUID ids => naturally balanced buckets.
bvbrc_fetch_all <- function(collection, rql_filter, select,
                            page_size = PAGE_MAX, parallel = FALSE, workers = 6L,
                            key = "id") {
  if (!parallel) {
    return(.keyset_walk(collection, rql_filter, select, page_size, key))
  }
  if (key != "id") {
    stop("parallel id-partitioning requires key = 'id' (UUID collections). ",
         "Use parallel = FALSE for the genome collection (key = 'genome_id').")
  }
  # force promises so the VALUES (not unevaluated args) are exported to workers
  force(rql_filter); force(select); force(page_size); force(key)
  hex   <- c(0:9, letters[1:6])      # '0'..'f' (UUID first char)
  edges <- c(hex, "g")               # 'g' > 'f' as an exclusive upper bound
  buckets <- Map(c, edges[-length(edges)], edges[-1])

  worker <- function(b) {
    filt <- sprintf("and(%s,ge(id,%s),lt(id,%s))", rql_filter, b[1], b[2])
    .keyset_walk(collection, filt, select, page_size, key)
  }

  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  parts <- future.apply::future_lapply(
    buckets, worker,
    future.packages = c("httr2", "jsonlite", "data.table"),
    future.globals = c(".keyset_walk", ".bvbrc_page", ".bvbrc_req", ".enc",
                       "BVBRC_BASE", "PAGE_MAX", "rql_filter", "select",
                       "page_size", "key")
  )
  data.table::rbindlist(Filter(Negate(is.null), parts), fill = TRUE, use.names = TRUE)
}

# --- 3. faceted subsampling -------------------------------------------------
# "Give me the N genomes with the most coverage across these drugs."
# Facet fields come back sorted by count desc, so the top-N buckets ARE the
# best-covered genomes.
bvbrc_rank_genomes_by_drug_coverage <- function(species, antibiotics, n = 500L) {
  ab  <- paste(vapply(antibiotics, .enc, character(1)), collapse = ",")
  rql <- sprintf("and(eq(genome_name,%s),in(antibiotic,(%s)))", .enc(species), ab)
  q   <- sprintf(
    "%s&limit(1)&facet((field,genome_id),(mincount,1),(limit,%d))&json(nl,map)",
    rql, n
  )
  resp <- .bvbrc_req("genome_amr", q, accept = "application/solr+json") |> req_perform()
  ff <- fromJSON(resp_body_string(resp))$facet_counts$facet_fields$genome_id
  data.frame(genome_id = names(ff), n_drugs = as.integer(unlist(ff)),
             row.names = NULL)
}

# --- example / benchmark ----------------------------------------------------
# Call bvbrc_demo() explicitly (sourcing this file no longer auto-runs it).
bvbrc_demo <- function() {
  cat("count (S. aureus AMR rows):",
      bvbrc_count("genome_amr", "eq(genome_name,Staphylococcus%20aureus)"), "\n")

  t <- system.time(
    d <- bvbrc_fetch_all("genome_amr", "eq(genome_name,Staphylococcus%20aureus)",
                         "genome_id,antibiotic,resistant_phenotype")
  )
  cat(sprintf("sequential fetch_all: %d rows in %.1fs\n", nrow(d), t["elapsed"]))

  tp <- system.time(
    dp <- bvbrc_fetch_all("genome_amr", "eq(genome_name,Staphylococcus%20aureus)",
                          "genome_id,antibiotic,resistant_phenotype",
                          parallel = TRUE, workers = 6L)
  )
  cat(sprintf("parallel   fetch_all: %d rows in %.1fs\n", nrow(dp), tp["elapsed"]))

  top <- bvbrc_rank_genomes_by_drug_coverage(
    "Staphylococcus aureus", c("ciprofloxacin", "gentamicin", "oxacillin"), n = 5L)
  cat("top genomes by drug coverage:\n"); print(top)
}
