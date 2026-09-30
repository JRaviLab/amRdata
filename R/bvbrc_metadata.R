# Internal helpers for pulling and cleaning the complete BV-BRC bacterial
# metadata table. Used by data_raw/species_count_summary.R and
# vignettes/BVBRC_stats.Rmd.

#' Fetch the complete BV-BRC bacterial genome metadata table
#'
#' Pages through the BV-BRC `genome` collection (`superkingdom = Bacteria`) with
#' [.bvbrcApiFetch()], selecting every stored field in the collection schema.
#' Unlike [.fetchBVBRCdata()], this is unfiltered and returns all columns.
#'
#' @param verbose Logical. Print progress messages. Default: TRUE.
#'
#' @return A tibble with one row per bacterial genome, character columns coerced
#'   to UTF-8.
#' @keywords internal
#' @noRd
.fetchAllBVBRCmetadataApi <- function(verbose = TRUE) {
  rql_filter <- "eq(superkingdom,Bacteria)"

  # The default response omits fields such as antimicrobial_resistance, so
  # select every stored field listed in the collection schema.
  if (isTRUE(verbose)) message("Discovering BV-BRC genome fields...")
  schema_resp <- httr2::request(paste0(.BVBRC_API_BASE, "/genome/schema")) |>
    httr2::req_timeout(120) |>
    httr2::req_retry(max_tries = 5L) |>
    httr2::req_perform()
  schema_fields <- jsonlite::fromJSON(
    httr2::resp_body_string(schema_resp)
  )$schema$fields
  fields <- schema_fields$name[schema_fields$stored %in% TRUE]
  fields <- setdiff(fields, c("genome_id", "_version_")) # key is always selected
  if (length(fields) == 0L) stop("Unable to determine genome attributes.")

  if (isTRUE(verbose)) {
    message("Fetching BV-BRC bacterial genome metadata (this may take time)...")
  }
  df <- .bvbrcApiFetch(
    "genome", rql_filter,
    select = paste(fields, collapse = ","),
    key = "genome_id"
  )
  if (nrow(df) == 0L) {
    stop("BV-BRC returned no data. The service may be unavailable.")
  }

  df <- dplyr::mutate(
    df,
    dplyr::across(tidyselect::where(is.character), .toUtf8)
  )

  if (isTRUE(verbose)) {
    message(sprintf(
      "Retrieved %s rows x %s columns.",
      format(nrow(df), big.mark = ","), ncol(df)
    ))
  }
  df
}

# Strings that count as missing in text fields
.BVBRC_MISSING_TOKENS <- c(
  "", " ", "NA", "N/A", "n/a", "na", "null", "Null", "NULL",
  "unknown", "Unknown", "UNKNOWN",
  "not provided", "Not provided", "NOT PROVIDED",
  "not reported", "Not reported", "NOT REPORTED",
  "unavailable", "Unspecified", "unspecified", "missing"
)

.isPlaceholderNA <- function(x) {
  if (!is.character(x)) {
    return(rep(FALSE, length(x)))
  }
  stringr::str_trim(x) %in% .BVBRC_MISSING_TOKENS
}

# Clean a list column: squish, placeholder -> NA, all-missing -> NA
.cleanListCol <- function(x) {
  purrr::map(
    x,
    \(z) {
      if (is.null(z) || length(z) == 0) {
        return(NA_character_)
      }
      z <- stringr::str_squish(as.character(unlist(z)))
      z[.isPlaceholderNA(z)] <- NA_character_
      if (all(is.na(z))) {
        return(NA_character_)
      }
      z
    }
  )
}

# Collapse a list column to ';'-separated strings
.flattenListCol <- function(x) {
  purrr::map_chr(
    x,
    \(z) {
      if (is.null(z) || length(z) == 0) {
        return(NA_character_)
      }
      paste(as.character(unlist(z)), collapse = ";")
    }
  )
}

#' Clean the complete BV-BRC metadata table
#'
#' Squishes whitespace in character/factor columns, maps placeholder strings
#' (`"unknown"`, `"not reported"`, ...) to `NA`, and cleans then flattens list
#' columns to `;`-separated strings.
#'
#' @param df Tibble from [.fetchAllBVBRCmetadataApi()].
#' @return A tibble with the same columns, placeholders as `NA`, no list columns.
#' @keywords internal
#' @noRd
.cleanBVBRCmetadata <- function(df) {
  df <- dplyr::mutate(
    df,
    dplyr::across(
      tidyselect::where(~ is.character(.x) || is.factor(.x)),
      ~ {
        x <- stringr::str_squish(as.character(.x))
        x[.isPlaceholderNA(x)] <- NA
        x
      }
    )
  )
  list_cols <- names(df)[vapply(df, is.list, logical(1))]
  dplyr::mutate(
    df,
    dplyr::across(dplyr::all_of(list_cols), .cleanListCol),
    dplyr::across(dplyr::all_of(list_cols), .flattenListCol)
  )
}

# Resource name of the cleaned table in the BiocFileCache
.BVBRC_CACHE_RNAME <- "amRdata_bvbrc_clean"

#' Path of the cached BV-BRC metadata parquet, or NA if not cached yet
#'
#' BiocFileCache stores relative paths, so resolve through `bfcrpath()`. If
#' several entries share the resource name, the newest (highest rid) is used.
#' @noRd
.bvbrcCachedPath <- function(bfc = BiocFileCache::BiocFileCache(ask = FALSE)) {
  hit <- BiocFileCache::bfcquery(
    bfc, .BVBRC_CACHE_RNAME,
    field = "rname", exact = TRUE
  )
  if (nrow(hit) == 0L) {
    return(NA_character_)
  }
  rid <- hit$rid[order(as.integer(sub("^BFC", "", hit$rid)), decreasing = TRUE)][1]
  path <- unname(BiocFileCache::bfcrpath(bfc, rids = rid))
  if (!file.exists(path)) {
    return(NA_character_)
  }
  path
}

#' Load the complete, cleaned BV-BRC bacterial metadata table
#'
#' Reads the table (all columns) from the BiocFileCache if present; otherwise
#' downloads it with [.fetchAllBVBRCmetadataApi()], cleans it with
#' [.cleanBVBRCmetadata()], and caches it as parquet.
#'
#' @param refresh Logical. Re-download even if a cached copy exists.
#' @return A tibble.
#' @keywords internal
#' @noRd
.loadBVBRCmetadata <- function(refresh = FALSE) {
  bfc <- BiocFileCache::BiocFileCache(ask = FALSE)
  cached <- .bvbrcCachedPath(bfc)
  if (!is.na(cached) && !isTRUE(refresh)) {
    return(tibble::as_tibble(arrow::read_parquet(cached)))
  }

  bvbrc_clean <- .cleanBVBRCmetadata(.fetchAllBVBRCmetadataApi())
  path <- if (is.na(cached)) {
    BiocFileCache::bfcnew(bfc, .BVBRC_CACHE_RNAME, ext = ".parquet")
  } else {
    cached
  }
  .write_compressed_parquet(bvbrc_clean, path)
  bvbrc_clean
}
