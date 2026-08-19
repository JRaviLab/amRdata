# BV-BRC Data API path (additive; opt-in via retrieveMetadata(method = "api"))
# -----------------------------------------------------------------------------
# Replaces ONLY the Docker/p3-* download of AMR + genome metadata. It returns the
# same two tibbles the Docker path produces -- columns prefixed `genome_drug.*`
# and `genome.*` -- so all downstream QC/join/DuckDB logic is reused unchanged.
#
# Why: the p3-* CLI path is stochastic (issue #30) -- it merges stderr into the
# data stream and never retries, so a transient BV-BRC 503 corrupts a batch. The
# API path uses typed JSON + retry-with-backoff, so a 503 is retried, never
# parsed as data.
#
# @keywords internal

.BVBRC_API_BASE <- "https://www.bv-brc.org/api"
.BVBRC_PAGE_MAX <- 25000L

# One HTTP request with retry on transient failures.
.bvbrc_api_req <- function(collection, rql, accept = "application/json") {
  url <- sprintf("%s/%s/?%s", .BVBRC_API_BASE, collection, rql)
  httr2::request(url) |>
    httr2::req_headers(Accept = accept) |>
    httr2::req_timeout(120) |>
    httr2::req_retry(
      max_tries = 5L,
      is_transient = function(resp) httr2::resp_status(resp) %in%
        c(429, 500, 502, 503, 504)
    )
}

.bvbrc_enc <- function(x) utils::URLencode(as.character(x), reserved = TRUE)

# Keyset-paginate a filtered query past the 25k cap. `key` is the unique
# sort/seek field ("id" for genome_amr, "genome_id" for genome).
.bvbrc_api_fetch <- function(collection, rql_filter, select, key = "id") {
  out <- list()
  last <- ""
  i <- 0L
  repeat {
    seek <- if (nzchar(last)) {
      sprintf("and(%s,gt(%s,%s))", rql_filter, key, .bvbrc_enc(last))
    } else {
      rql_filter
    }
    rql <- sprintf(
      "%s&select(%s,%s)&sort(+%s)&limit(%d)",
      seek, key, select, key, .BVBRC_PAGE_MAX
    )
    resp <- httr2::req_perform(.bvbrc_api_req(collection, rql))
    pg <- jsonlite::fromJSON(httr2::resp_body_string(resp),
      simplifyDataFrame = TRUE
    )
    raw_n <- if (is.data.frame(pg)) nrow(pg) else 0L
    if (raw_n == 0L) break
    if (nzchar(last)) pg <- pg[pg[[key]] != last, , drop = FALSE]
    if (nrow(pg) > 0L) {
      i <- i + 1L
      out[[i]] <- pg
      last <- pg[[key]][nrow(pg)]
    }
    if (raw_n < .BVBRC_PAGE_MAX) break
  }
  if (i == 0L) {
    return(tibble::tibble())
  }
  tibble::as_tibble(data.table::rbindlist(out, fill = TRUE, use.names = TRUE))
}

# Split a vector into chunks of size n (keeps in(...) URLs within length limits).
.bvbrc_chunk <- function(x, n) {
  if (length(x) == 0L) {
    return(list())
  }
  split(x, ceiling(seq_along(x) / n))
}

# Ensure every expected field is present (fill missing with NA), coerce to
# character, order columns, and apply the `prefix.` naming convention that the
# Docker/p3 parser produces.
.bvbrc_prefix_fill <- function(df, expected, prefix) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  n <- nrow(df) # rep() keeps length right when the query returned 0 rows
  for (f in expected) {
    if (!f %in% names(df)) df[[f]] <- rep(NA_character_, n)
  }
  df <- df[, expected, drop = FALSE]
  # coerce to character and use "" for missing, matching the Docker/TSV parser
  # (.parse_bvbrc_tsv yields "" for blank fields, not NA) so the two paths agree.
  df[] <- lapply(df, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  })
  names(df) <- paste0(prefix, ".", expected)
  tibble::as_tibble(df)
}

# --- AMR phenotype (genome_amr) -> genome_drug.* ------------------------------
.extractAMRtable_api <- function(genome_ids, abx = "All",
                                 chunk_size = 150L, verbose = TRUE) {
  # genome_amr collection has these; measurement_unit/computational_method/source
  # are absent and get filled NA to match the Docker column set.
  expected <- c(
    "genome_id", "antibiotic", "computational_method", "evidence",
    "genome_name", "id", "laboratory_typing_method",
    "laboratory_typing_platform", "measurement", "measurement_sign",
    "measurement_unit", "measurement_value", "pmid", "resistant_phenotype",
    "source", "taxon_id", "testing_standard"
  )
  have <- intersect(expected, c(
    "genome_id", "antibiotic", "evidence", "genome_name", "id",
    "laboratory_typing_method", "laboratory_typing_platform", "measurement",
    "measurement_sign", "measurement_value", "pmid", "resistant_phenotype",
    "taxon_id", "testing_standard"
  ))
  # keyset key for genome_amr is "id"; keep genome_id as a data column.
  sel <- paste(setdiff(have, "id"), collapse = ",")

  chunks <- .bvbrc_chunk(genome_ids, chunk_size)
  if (isTRUE(verbose)) {
    message("  [api] AMR: ", length(genome_ids), " genomes in ",
      length(chunks), " chunk(s)")
  }
  parts <- lapply(chunks, function(ids) {
    ab <- if (identical(abx, "All")) {
      ""
    } else {
      sprintf(",in(antibiotic,(%s))",
        paste(vapply(abx, .bvbrc_enc, ""), collapse = ","))
    }
    filt <- sprintf("and(in(genome_id,(%s))%s)",
      paste(vapply(ids, .bvbrc_enc, ""), collapse = ","), ab)
    .bvbrc_api_fetch("genome_amr", filt, sel, key = "id")
  })
  df <- data.table::rbindlist(parts, fill = TRUE, use.names = TRUE)
  .bvbrc_prefix_fill(df, expected, "genome_drug")
}

# --- genome-ID resolution (genome) --------------------------------------------
# API replacement for .retrieveQueryIDs(): resolve species names and/or taxon IDs
# to Good-quality WGS/Complete genome IDs, and write the `bac_data` table (used by
# retrieveMetadata()'s summary). Uses the Data API instead of the Docker-built
# cache, so retrieveMetadata(method = "api") needs no Docker.
.resolveGenomeIDs_api <- function(base_dir = ".", user_bacs,
                                  overwrite = FALSE, verbose = TRUE) {
  sel <- "genome_name,taxon_id,species,strain"
  parts <- lapply(user_bacs, function(ub) {
    ub <- trimws(as.character(ub))
    key_filter <- if (grepl("^[0-9]+$", ub)) {
      sprintf("eq(taxon_lineage_ids,%s)", ub) # taxon ID (any rank)
    } else {
      sprintf("eq(species,%s)", .bvbrc_enc(ub)) # species name
    }
    filt <- sprintf(
      "and(%s,eq(genome_quality,Good),in(genome_status,(WGS,Complete)))",
      key_filter
    )
    if (isTRUE(verbose)) message("  [api] resolving genome IDs for '", ub, "'")
    .bvbrc_api_fetch("genome", filt, sel, key = "genome_id")
  })
  df <- as.data.frame(
    data.table::rbindlist(parts, fill = TRUE, use.names = TRUE),
    stringsAsFactors = FALSE
  )
  if (nrow(df) == 0L) {
    return(character(0))
  }
  df <- df[grepl("^[0-9]+[.][0-9]+$", df$genome_id), , drop = FALSE]
  df <- df[!duplicated(df$genome_id), , drop = FALSE]

  # write bac_data (genome.* columns), mirroring .retrieveQueryIDs()
  paths <- .buildDBpath(base_dir = base_dir, user_bacs = user_bacs, overwrite = overwrite)
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = paths$db_path)
  on.exit(try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)
  bac <- data.frame(
    `genome.genome_id` = df$genome_id,
    `genome.genome_name` = df$genome_name,
    `genome.taxon_id` = df$taxon_id,
    `genome.species` = df$species,
    `genome.strain` = df$strain,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  DBI::dbWriteTable(con, "bac_data", bac, overwrite = TRUE)
  if (isTRUE(verbose)) {
    message("  [api] resolved ", nrow(df), " genome IDs; wrote bac_data")
  }
  unique(df$genome_id)
}

# --- genome metadata (genome) -> genome.* -------------------------------------
.extractGenomeData_api <- function(genome_ids, fields,
                                   chunk_size = 150L, verbose = TRUE) {
  expected <- strsplit(fields, ",", fixed = TRUE)[[1]]
  expected <- unique(c("genome_id", expected))
  sel <- paste(setdiff(expected, "genome_id"), collapse = ",")

  chunks <- .bvbrc_chunk(genome_ids, chunk_size)
  if (isTRUE(verbose)) {
    message("  [api] genome metadata: ", length(genome_ids), " genomes in ",
      length(chunks), " chunk(s)")
  }
  parts <- lapply(chunks, function(ids) {
    filt <- sprintf("in(genome_id,(%s))",
      paste(vapply(ids, .bvbrc_enc, ""), collapse = ","))
    .bvbrc_api_fetch("genome", filt, sel, key = "genome_id")
  })
  df <- data.table::rbindlist(parts, fill = TRUE, use.names = TRUE)
  .bvbrc_prefix_fill(df, expected, "genome")
}
