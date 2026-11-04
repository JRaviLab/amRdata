# Drug name / class / abbreviation tables (generators for data_raw/*.tsv)
# -----------------------------------------------------------------------------
# Two steps (instructions: data_raw/README.md):
#   1. prepareDrugTables()   BV-BRC antibiotic names -> cleaned names, classes
#                            (rules + BV-BRC / WHO AWaRe / AntibioticDB) and
#                            abbreviations, written as four *_review.tsv files
#   ... you edit the review files ...
#   2. finalizeDrugTables()  checks the reviewed files, writes the final TSVs
#
# The final TSVs are only ever written by finalizeDrugTables(), and only when
# every check passes. Everything heuristic (typo fixes, class suggestions,
# generated abbreviations) reaches the final tables only through your review.
#
# @keywords internal

# Name normalisation -----------------------------------------------------------

# Same convention as clean_drug.tsv: lower case, "/" -> "-", spaces -> "_".
# Non-ASCII junk (e.g. the stray "A-circumflex" in some BV-BRC values) is dropped.
.normalizeDrugName <- function(x) {
  x <- gsub("[^\\x01-\\x7F]", "", as.character(x), perl = TRUE)
  x <- tolower(trimws(x))
  x <- sub("^para-", "para_", x) # the one hyphenated single-compound name
  x <- gsub("\\s*/\\s*", "-", x)
  x <- gsub("[[:space:]_]+", "_", x)
  x <- gsub("-{2,}", "-", x)
  gsub("^[-_]+|[-_]+$", "", x)
}

# Reference sources -------------------------------------------------------------

# Distinct `antibiotic` values in genome_amr (what clean_drug is keyed on) with
# their record counts, via a Solr facet.
.fetchAmrDrugNames <- function() {
  resp <- httr2::req_perform(.bvbrcApiReq(
    "genome_amr",
    "keyword(*)&facet((field,antibiotic),(mincount,1),(limit,-1))&limit(0)",
    accept = "application/solr+json"
  ))
  f <- jsonlite::fromJSON(httr2::resp_body_string(resp),
    simplifyVector = FALSE
  )$facet_counts$facet_fields$antibiotic
  tibble::tibble(
    original_drug = unlist(f[seq(1L, length(f), 2L)]),
    n_records = unlist(f[seq(2L, length(f), 2L)])
  )
}

# "Chemical structure" pharmacologic class if present, else the EPC.
.bvbrcPharmClass <- function(x) {
  if (length(x) == 0L) {
    return(NA_character_)
  }
  pick <- function(pattern) {
    hit <- x[grepl(pattern, x)]
    if (length(hit)) sub(".*\\]: ", "", hit[1L]) else NA_character_
  }
  out <- pick("^Chemical/Ingredient")
  if (is.na(out)) out <- pick("^Established Pharmacologic Class")
  out
}

# ATC labels come as concatenated 5-level paths; take the chemical subgroup
# (level 4) of the first systemic anti-infective path.
.bvbrcAtcClass <- function(x) {
  if (length(x) < 5L) {
    return(NA_character_)
  }
  paths <- split(x, ceiling(seq_along(x) / 5L))
  sys <- Filter(function(p) {
    length(p) == 5L && p[1L] == "Antiinfectives for systemic use"
  }, paths)
  if (length(sys)) sys[[1L]][4L] else NA_character_
}

.readWhoEml <- function() {
  who_url <- "https://iris.who.int/server/api/core/bitstreams/abba5c2a-8457-431c-9695-16cb2317dd0e/content"
  destfile <- tempfile(fileext = ".xlsx")
  utils::download.file(who_url, destfile = destfile, mode = "wb", quiet = TRUE)
  # An .xlsx is a zip archive (starts with "PK"); WHO's old bitstream URL
  # returns an HTML page instead, which readxl reports as an unreadable zip.
  if (!identical(readBin(destfile, "raw", 2L), charToRaw("PK"))) {
    stop("WHO EML download did not return an xlsx file; check `who_url`")
  }
  who <- suppressMessages(readxl::read_excel(
    destfile,
    sheet = "AWaRe classification 2023", col_names = FALSE
  )) |>
    tibble::as_tibble() |>
    dplyr::slice(-c(1:3))
  hdr <- unlist(who[1L, ], use.names = FALSE)
  who <- who[-1L, !is.na(hdr)] # drop header row and unnamed (empty) columns
  colnames(who) <- trimws(hdr[!is.na(hdr)])
  if (!all(c("Antibiotic", "Class") %in% colnames(who))) {
    stop("The EML list is missing expected column names and cannot be processed further")
  }
  out <- tibble::tibble(
    name = .normalizeDrugName(who[["Antibiotic"]]),
    who_class = .normalizeDrugName(who[["Class"]])
  )
  out[!duplicated(out$name), ]
}

# AntibioticDB (DOI: 10.1093/jac/dky208). "Drug Class" is free text such as
# "Beta-lactam (carbapenem)"; the parenthetical is the more specific label.
.readAdb <- function(adb_file) {
  adb <- readr::read_csv(adb_file, show_col_types = FALSE, progress = FALSE)
  cls <- adb[["Drug Class"]]
  sub_cls <- ifelse(grepl(" \\(", cls), sub(".* \\(([^)]*)\\).*", "\\1", cls), NA)
  out <- tibble::tibble(
    name = .normalizeDrugName(adb[["Drug Name"]]),
    adb_class = .normalizeDrugName(ifelse(is.na(sub_cls), cls, sub_cls))
  )
  out <- out[!is.na(out$name) & nzchar(out$name), ]
  out[!duplicated(out$name), ]
}

# One fetch of every reference source, reused by all stages.
.loadDrugReferences <- function(adb_file, verbose = TRUE) {
  if (isTRUE(verbose)) message("Loading BV-BRC / WHO / AntibioticDB references")
  bv <- .bvbrcApiFetch(
    "antibiotics", "eq(antibiotic_name,*)",
    "pharmacological_classes,atc_classification,synonyms",
    key = "antibiotic_name"
  )
  bvbrc <- tibble::tibble(
    name = .normalizeDrugName(bv$antibiotic_name),
    bvbrc_class = .normalizeDrugName(
      vapply(bv$pharmacological_classes, .bvbrcPharmClass, "")
    ),
    atc_class = .normalizeDrugName(
      vapply(bv$atc_classification, .bvbrcAtcClass, "")
    )
  )
  # BV-BRC synonyms (brand names, salts, alternate spellings) -> generic name.
  # Only unambiguous aliases: not a drug name in their own right, one target.
  syn <- tibble::tibble(
    alias = .normalizeDrugName(unlist(bv$synonyms)),
    name = rep(.normalizeDrugName(bv$antibiotic_name), lengths(bv$synonyms))
  )
  syn <- syn[nzchar(syn$alias) & syn$alias != syn$name &
    grepl("^[a-z0-9_-]+$", syn$alias) & !syn$alias %in% bvbrc$name, ]
  n_targets <- tapply(syn$name, syn$alias, function(v) length(unique(v)))
  syn <- syn[syn$alias %in% names(n_targets)[n_targets == 1L], ]
  list(
    bvbrc = bvbrc[!duplicated(bvbrc$name), ],
    synonyms = syn[!duplicated(syn$alias), ],
    who = .readWhoEml(),
    adb = .readAdb(adb_file)
  )
}

# Name matching ------------------------------------------------------------------

# Closest known name (or BV-BRC synonym, reported as its generic name) by edit
# distance, only when it looks like a typo: same first letter, distance <= 1,
# or <= 2 for names of 8+ characters. `alias` is a named vector alias -> name.
.nearestDrugName <- function(x, vocab, alias = character(0)) {
  none <- list(name = NA_character_, dist = NA_real_, via = NA_character_)
  keys <- c(vocab, names(alias)) # names first, so they win ties over aliases
  target <- c(vocab, unname(alias))
  keep <- substr(keys, 1L, 1L) == substr(x, 1L, 1L)
  if (!nzchar(x) || !any(keep)) {
    return(none)
  }
  keys <- keys[keep]
  target <- target[keep]
  d <- drop(utils::adist(x, keys))
  j <- which.min(d)
  if (d[j] > if (nchar(x) >= 8L) 2 else 1) {
    return(none)
  }
  list(
    name = target[j], dist = d[j],
    via = if (keys[j] != target[j]) keys[j] else NA_character_
  )
}

# Suggest a cleaned name for `x`, with how it was found: "synonym" (BV-BRC
# alias), "typo", or "combination" -- for "<a>_<b>", "<a> <b>" or a part with a
# typo, the "<a>-<b>" of two known names (known order preferred:
# "piperacillin-tazobactam" over "tazobactam-piperacillin").
.suggestDrugName <- function(x, vocab, alias = character(0)) {
  resolve <- function(part) {
    cand <- unique(c(part, gsub("_", "-", part, fixed = TRUE), gsub("-", "_", part, fixed = TRUE)))
    hit <- cand[cand %in% vocab]
    if (length(hit)) {
      return(list(name = hit[1L], dist = 0, via = NA_character_))
    }
    hit <- cand[cand %in% names(alias)]
    if (length(hit)) {
      return(list(name = alias[[hit[1L]]], dist = 0, via = hit[1L]))
    }
    .nearestDrugName(part, vocab, alias)
  }
  best <- resolve(x)
  best$type <- if (is.na(best$name)) {
    NA_character_
  } else if (!is.na(best$via)) {
    "synonym"
  } else {
    "typo"
  }
  seps <- gregexpr("[_-]", x)[[1L]]
  for (p in seps[seps > 0L]) {
    a <- resolve(substr(x, 1L, p - 1L))
    b <- resolve(substr(x, p + 1L, nchar(x)))
    if (is.na(a$name) || is.na(b$name)) next
    combo <- if (!paste(a$name, b$name, sep = "-") %in% vocab &&
      paste(b$name, a$name, sep = "-") %in% vocab) {
      paste(b$name, a$name, sep = "-")
    } else {
      paste(a$name, b$name, sep = "-")
    }
    d <- a$dist + b$dist
    if (is.na(best$dist) || d < best$dist || (d == best$dist && combo %in% vocab)) {
      best <- list(name = combo, dist = d, via = NA_character_, type = "combination")
    }
  }
  best
}


# Class rules ---------------------------------------------------------------------

# Second component of a combination that makes it "<class>-beta_lactamase_inhibitors".
.BLI_COMPONENTS <- c(
  "clavulanic_acid", "sulbactam", "tazobactam", "avibactam", "relebactam",
  "vaborbactam", "taniborbactam", "enmetazobactam", "durlobactam"
)

# Match a suggested class to the curated class vocabulary (tolerates a
# singular/plural difference); NA when the class is new to the scheme.
.matchScheme <- function(x, scheme) {
  vapply(x, function(v) {
    if (is.na(v) || !nzchar(v)) {
      return(NA_character_)
    }
    # "third-generation-cephalosporins", "cephalosporin,_third_generation" -> cephalosporin(s)
    v <- sub("^(first|second|third|fourth|fifth)-generation-", "", v)
    v <- sub(",_.*$", "", v)
    hit <- scheme[scheme %in% c(v, paste0(v, "s"), sub("s$", "", v))]
    if (length(hit)) hit[1L] else NA_character_
  }, "", USE.NAMES = FALSE)
}


# Abbreviations -------------------------------------------------------------------

# 3-letter code from the first word, extending letter by letter on collision.
# Only letters are used: class codes are pasted with "_" downstream
# (resistant_classes), so they must not contain one.
.makeAbbr <- function(name, used) {
  first <- gsub("[^a-z]", "", strsplit(name, "_", fixed = TRUE)[[1L]][1L])
  if (is.na(first) || !nzchar(first)) first <- "x"
  # cef-/ceph- drugs all start alike; distinguish on the consonants that follow
  # (ceftriaxone CRO, ceftazidime CAZ, cefotaxime CTX ...)
  if (grepl("^ce(f|ph)[a-z]{3,}", first)) {
    rest <- gsub("[aeiou]", "", sub("^ce(f|ph)", "", first))
    if (nchar(rest) >= 2L) first <- paste0("c", rest)
  }
  first <- toupper(first)
  n <- nchar(first)
  # first letter + next two consonants (rifaximin RFX, telavancin TLV), used
  # when the plain first three letters are taken
  cons <- gsub("[AEIOU]", "", substring(first, 2L))
  cands <- unique(c(
    substr(first, 1L, 3L),
    if (nchar(cons) >= 2L) paste0(substr(first, 1L, 1L), substr(cons, 1L, 2L)),
    vapply(seq_len(max(0L, n - 3L)), function(k) {
      paste0(substr(first, 1L, 2L), substr(first, 2L + k, 2L + k))
    }, ""),
    substr(first, 1L, 4L)
  ))
  ok <- cands[nzchar(cands) & !cands %in% used]
  if (length(ok)) ok[1L] else paste0(substr(first, 1L, 3L), length(used))
}

# Abbreviations for `names` not in `existing` (named vector: name -> abbr).
# "<a>-<b>" combinations join the component abbreviations with "-".
.buildAbbrVector <- function(names, existing) {
  used <- unname(existing)
  out <- character(0)
  lookup <- function(nm) {
    if (nm %in% names(existing)) {
      return(existing[[nm]])
    }
    if (nm %in% names(out)) {
      return(out[[nm]])
    }
    a <- .makeAbbr(nm, used)
    used <<- c(used, a)
    out[[nm]] <<- a
    a
  }
  todo <- setdiff(names, names(existing))
  for (nm in todo) {
    parts <- strsplit(nm, "-", fixed = TRUE)[[1L]]
    if (length(parts) > 1L) {
      a <- paste(vapply(parts, lookup, ""), collapse = "-")
      if (a %in% used) a <- paste0(a, length(used))
      used <- c(used, a)
      out[[nm]] <- a
    } else {
      lookup(nm)
    }
  }
  out[todo]
}


# Working tables ----------------------------------------------------------------

.DRUG_TABLES <- list(
  clean_drug = c("original_drug", "cleaned_drug"),
  drug_class = c("drug", "drug_class"),
  drug_abbr = c("drug", "drug_abbr"),
  class_abbr = c("drug_class", "class_abbr")
)

# Cell value for "you still need to decide this". A blank cleaned_drug keeps its
# meaning from the final table: not a drug, exclude.
.TODO <- "TODO"

# Real drugs whose name is also the singular of a class name.
.CLASS_NAMED_DRUGS <- c("penicillin", "tetracycline")

.tsvEol <- function(file) {
  crlf <- file.exists(file) &&
    any(readBin(file, "raw", file.size(file)) == as.raw(13L))
  if (crlf) "\r\n" else "\n"
}

# All-character table, blank -> NA; a missing file is an empty table.
.readDrugTable <- function(file, cols) {
  if (!file.exists(file)) {
    return(tibble::as_tibble(stats::setNames(
      rep(list(character(0)), length(cols)), cols
    )))
  }
  x <- readr::read_tsv(file,
    col_types = readr::cols(.default = "c"), na = "", progress = FALSE
  )
  miss <- setdiff(cols, names(x))
  if (length(miss)) {
    stop(basename(file), " is missing column(s): ", paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  for (cl in names(x)) {
    x[[cl]] <- trimws(x[[cl]])
    x[[cl]][!is.na(x[[cl]]) & !nzchar(x[[cl]])] <- NA_character_
  }
  x
}

.reviewFile <- function(nm, review_dir) file.path(review_dir, paste0(nm, "_review.tsv"))
.finalFile <- function(nm, final_dir) file.path(final_dir, paste0(nm, ".tsv"))

.isDecided <- function(x) !is.na(x) & x != .TODO

# Ensure `x` has all of `cols` (NA where missing) and return them in that order.
.withCols <- function(x, cols) {
  for (cl in setdiff(cols, names(x))) x[[cl]] <- rep(NA_character_, nrow(x))
  x[, cols, drop = FALSE]
}

# Proposals -----------------------------------------------------------------------

# New antibiotic names in genome_amr -> rows for the clean_drug review table.
# An exact match to a known antibiotic is filled in; anything else is TODO with
# a suggestion (BV-BRC synonym, typo fix, or two-drug combination).
.proposeCleanDrug <- function(clean, amr, refs) {
  vocab <- unique(c(
    refs$bvbrc$name, refs$who$name, refs$adb$name,
    clean$cleaned_drug[.isDecided(clean$cleaned_drug)]
  ))
  vocab <- vocab[nzchar(vocab)]
  simple <- vocab[grepl("^[a-z0-9_-]+$", vocab)]
  alias <- stats::setNames(refs$synonyms$name, refs$synonyms$alias)
  alias <- alias[!names(alias) %in% vocab]

  new <- amr[!amr$original_drug %in% clean$original_drug, ]
  new <- new[order(-new$n_records), ]
  norm <- .normalizeDrugName(new$original_drug)
  n <- nrow(new)
  cleaned <- rep(.TODO, n)
  suggestion <- rep(NA_character_, n)
  note <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    if (norm[i] %in% vocab) {
      cleaned[i] <- norm[i]
      note[i] <- "exact match"
    } else if (nzchar(norm[i])) {
      r <- .suggestDrugName(norm[i], simple, alias)
      suggestion[i] <- r$name
      note[i] <- if (is.na(r$type)) {
        "no match"
      } else if (is.na(r$via)) {
        r$type
      } else {
        paste0(r$type, " (", r$via, ")")
      }
    } else {
      note[i] <- "no match"
    }
  }
  cols <- c("original_drug", "cleaned_drug", "n_records", "suggestion", "note")
  rows <- tibble::tibble(
    original_drug = new$original_drug, cleaned_drug = cleaned,
    n_records = as.character(new$n_records), suggestion = suggestion, note = note
  )
  clean <- .withCols(clean, cols)
  clean$n_records <- as.character(amr$n_records[match(clean$original_drug, amr$original_drug)])
  rbind(clean, rows)
}

# Drugs without a class -> rows for the drug_class review table. A class that
# fits the existing scheme (or follows the beta-lactamase-inhibitor rule) is
# filled in; anything else is TODO with the suggestion beside it.
.proposeClasses <- function(todo, cls, refs) {
  cols <- c(
    "drug", "drug_class", "suggestion", "note",
    "who_class", "adb_class", "bvbrc_class", "atc_class"
  )
  cls <- .withCols(cls, cols)
  if (length(todo) == 0L) {
    return(cls)
  }
  scheme <- unique(cls$drug_class[.isDecided(cls$drug_class)])
  srcs <- c("who_class", "adb_class", "bvbrc_class", "atc_class")
  cand <- tibble::tibble(drug = todo)
  cand$who_class <- refs$who$who_class[match(todo, refs$who$name)]
  cand$adb_class <- refs$adb$adb_class[match(todo, refs$adb$name)]
  cand$bvbrc_class <- refs$bvbrc$bvbrc_class[match(todo, refs$bvbrc$name)]
  cand$atc_class <- refs$bvbrc$atc_class[match(todo, refs$bvbrc$name)]
  for (s in srcs) cand[[s]][!is.na(cand[[s]]) & !nzchar(cand[[s]])] <- NA_character_

  # first available label, then upgrade to the first one that fits the scheme
  suggested <- rep(NA_character_, nrow(cand))
  source <- rep(NA_character_, nrow(cand))
  for (s in rev(srcs)) {
    hit <- !is.na(cand[[s]])
    suggested[hit] <- cand[[s]][hit]
    source[hit] <- s
  }
  fits <- rep(NA_character_, nrow(cand))
  for (s in srcs) {
    m <- .matchScheme(cand[[s]], scheme)
    take <- is.na(fits) & !is.na(m)
    fits[take] <- m[take]
    suggested[take] <- m[take]
    source[take] <- s
  }

  # <drug>-<beta-lactamase inhibitor> combinations follow the first component
  base_class <- stats::setNames(
    c(cls$drug_class[.isDecided(cls$drug_class)], suggested),
    c(cls$drug[.isDecided(cls$drug_class)], cand$drug)
  )
  rule <- rep(FALSE, nrow(cand))
  for (i in seq_len(nrow(cand))) {
    parts <- strsplit(cand$drug[i], "-", fixed = TRUE)[[1L]]
    if (length(parts) == 2L && parts[2L] %in% .BLI_COMPONENTS) {
      b <- base_class[parts[1L]]
      if (!is.na(b) && !grepl("beta_lactamase_inhibitors$", b)) {
        suggested[i] <- paste0(b, "-beta_lactamase_inhibitors")
        source[i] <- "rule: <class>-beta_lactamase_inhibitors"
        rule[i] <- TRUE
      }
    }
  }
  fill <- rule | (!is.na(suggested) & suggested %in% scheme)
  rows <- cand
  rows$drug_class <- ifelse(fill, suggested, .TODO)
  rows$suggestion <- ifelse(fill, NA_character_, suggested)
  rows$note <- ifelse(is.na(source), "no match", source)
  rbind(cls, rows[, cols])
}

# Codes for decided drugs / classes that have none (generated, so to be checked).
.proposeAbbr <- function(todo, table, key, value) {
  cols <- c(key, value, "note")
  table <- .withCols(table, cols)
  todo <- setdiff(todo, table[[key]])
  if (length(todo) == 0L) {
    return(table)
  }
  ok <- .isDecided(table[[value]])
  codes <- .buildAbbrVector(todo, stats::setNames(table[[value]][ok], table[[key]][ok]))
  rows <- tibble::tibble(names(codes), unname(codes), "generated")
  names(rows) <- cols
  rbind(table, rows)
}

#' Prepare review tables for drug names, classes and abbreviations (step 1 of 2)
#'
#' Gets every antibiotic name in BV-BRC `genome_amr`, cleans them, and proposes a
#' class and abbreviations for each drug, writing four `*_review.tsv` files. The
#' final tables are not touched. Start from the existing final TSVs (or from an
#' earlier review file, which is continued so your edits are kept; delete the
#' review files to start over). Missing tables are treated as empty.
#'
#' Cells you must decide are `TODO`; a `suggestion` column beside them shows the
#' best guess (BV-BRC synonym, typo fix, combination, or class from WHO / ADB /
#' BV-BRC). Exact name matches, classes that fit the existing scheme,
#' beta-lactamase-inhibitor combinations, and generated abbreviations are filled
#' in but still need a look. A blank `cleaned_drug` means "not a drug, exclude".
#'
#' Class and abbreviation rows appear only for drugs whose cleaned name is
#' decided, so after resolving the `TODO`s in `clean_drug_review.tsv` run this
#' again to get their class and abbreviation rows, then finish those.
#'
#' @param final_dir [chr] directory of the final TSVs
#' @param review_dir [chr] where the review files are written
#' @param adb_file [chr] AntibioticDB `ADB_all_compounds.csv`
#' @param refs reference list from `.loadDrugReferences()` (fetched if `NULL`)
#' @param verbose [lgl] print progress
#' @return invisibly, a list of the four review tables
#' @seealso finalizeDrugTables
#' @keywords internal
prepareDrugTables <- function(final_dir = "data_raw", review_dir = "data_raw",
                              adb_file = file.path(final_dir, "ADB_all_compounds.csv"),
                              refs = NULL, verbose = TRUE) {
  say <- function(...) if (isTRUE(verbose)) message(...)
  if (is.null(refs)) refs <- .loadDrugReferences(adb_file, verbose)

  tbl <- list()
  for (nm in names(.DRUG_TABLES)) {
    rf <- .reviewFile(nm, review_dir)
    src <- if (file.exists(rf)) rf else .finalFile(nm, final_dir)
    say(nm, ": continuing from ", src, if (!file.exists(src)) " (missing, starting empty)")
    tbl[[nm]] <- .readDrugTable(src, .DRUG_TABLES[[nm]])
  }

  tbl$clean_drug <- .proposeCleanDrug(tbl$clean_drug, .fetchAmrDrugNames(), refs)
  decided <- unique(tbl$clean_drug$cleaned_drug[.isDecided(tbl$clean_drug$cleaned_drug)])

  no_class <- setdiff(decided, tbl$drug_class$drug)
  tbl$drug_class <- .proposeClasses(no_class, tbl$drug_class, refs)

  tbl$drug_abbr <- .proposeAbbr(decided, tbl$drug_abbr, "drug", "drug_abbr")
  classes <- unique(tbl$drug_class$drug_class[.isDecided(tbl$drug_class$drug_class)])
  tbl$class_abbr <- .proposeAbbr(classes, tbl$class_abbr, "drug_class", "class_abbr")

  dir.create(review_dir, showWarnings = FALSE, recursive = TRUE)
  for (nm in names(tbl)) {
    readr::write_tsv(tbl[[nm]], .reviewFile(nm, review_dir),
      na = "", eol = .tsvEol(.finalFile(nm, final_dir))
    )
    n_todo <- sum(vapply(tbl[[nm]], function(v) sum(v %in% .TODO), 0L))
    say(nm, "_review.tsv: ", nrow(tbl[[nm]]), " rows, ", n_todo, " TODO")
  }
  say(
    "Next: fix the TODO cells in the review files (see data_raw/README.md), ",
    "run this again if new drugs were decided, then run finalizeDrugTables()."
  )
  invisible(tbl)
}

#' Check the reviewed tables and write the final TSVs (step 2 of 2)
#'
#' Reads the four `*_review.tsv` files, checks them, and writes the final
#' `clean_drug.tsv`, `drug_class.tsv`, `drug_abbr.tsv` and `class_abbr.tsv`
#' (helper columns dropped). Nothing is written unless every check passes; the
#' error lists everything to fix. On success the review files are removed, so
#' the next `prepareDrugTables()` starts from the new final tables. Then run
#' `data_raw/make_data.R` to rebuild `data/*.rda`.
#'
#' Checks: no `TODO` left for a drug that is used; keys and codes unique;
#' cleaned names use only `a-z 0-9 _ -`; codes use only `A-Z 0-9 -` (no `_`, which
#' the resistance summary uses as a separator); every used drug has a class and
#' a code and every used class a code; no drug is named after a class (except
#' `penicillin` and `tetracycline`). Rows for drugs no cleaned name uses are
#' dropped if still `TODO` and kept otherwise.
#'
#' @param review_dir [chr] where the review files are
#' @param final_dir [chr] where the final TSVs are written
#' @param verbose [lgl] print a summary
#' @return invisibly, a list of the four final tables
#' @seealso prepareDrugTables
#' @keywords internal
finalizeDrugTables <- function(review_dir = "data_raw", final_dir = "data_raw",
                               verbose = TRUE) {
  say <- function(...) if (isTRUE(verbose)) message(...)
  files <- vapply(names(.DRUG_TABLES), .reviewFile, "", review_dir = review_dir)
  if (!all(file.exists(files))) {
    stop("Review file(s) not found: ", paste(basename(files[!file.exists(files)]), collapse = ", "),
      ". Run prepareDrugTables() first.",
      call. = FALSE
    )
  }
  tbl <- lapply(names(.DRUG_TABLES), function(nm) {
    .readDrugTable(.reviewFile(nm, review_dir), .DRUG_TABLES[[nm]])
  })
  names(tbl) <- names(.DRUG_TABLES)
  clean <- tbl$clean_drug
  cls <- tbl$drug_class
  dab <- tbl$drug_abbr
  cab <- tbl$class_abbr

  used <- unique(clean$cleaned_drug[.isDecided(clean$cleaned_drug)])
  unused_todo <- function(t, key, value) {
    !(t[[key]] %in% setdiff(t[[key]], used)) | !(t[[value]] %in% .TODO)
  }
  # proposals for names that ended up excluded are dropped, not errors
  cls <- cls[unused_todo(cls, "drug", "drug_class"), ]
  dab <- dab[unused_todo(dab, "drug", "drug_abbr"), ]
  used_classes <- unique(cls$drug_class[cls$drug %in% used])
  cab <- cab[!(!cab$drug_class %in% used_classes & cab$class_abbr %in% .TODO), ]

  errs <- character()
  warns <- character()
  add <- function(msg, x, warn = FALSE) {
    x <- unique(x[!is.na(x)])
    if (!length(x)) {
      return(invisible())
    }
    txt <- paste0(
      msg, ": ", paste(utils::head(x, 12), collapse = ", "),
      if (length(x) > 12) sprintf(" ... (%d in total)", length(x))
    )
    if (warn) warns <<- c(warns, txt) else errs <<- c(errs, txt)
  }
  dup <- function(x) x[duplicated(x) & !is.na(x)]

  add("clean_drug: still TODO", clean$original_drug[clean$cleaned_drug %in% .TODO])
  add("drug_class: blank or TODO", cls$drug[cls$drug_class %in% .TODO | is.na(cls$drug_class)])
  add("drug_abbr: blank or TODO", dab$drug[dab$drug_abbr %in% .TODO | is.na(dab$drug_abbr)])
  add("class_abbr: blank or TODO", cab$drug_class[cab$class_abbr %in% .TODO | is.na(cab$class_abbr)])

  add("clean_drug: duplicated original_drug", dup(clean$original_drug))
  add("drug_class: duplicated drug", dup(cls$drug))
  add("drug_abbr: duplicated drug", dup(dab$drug))
  add("drug_abbr: duplicated code", dup(dab$drug_abbr[.isDecided(dab$drug_abbr)]))
  add("class_abbr: duplicated class", dup(cab$drug_class))
  add("class_abbr: duplicated code", dup(cab$class_abbr[.isDecided(cab$class_abbr)]))

  name_ok <- "^[a-z0-9_-]+$"
  code_ok <- "^[A-Z0-9-]+$"
  add("cleaned_drug: use only a-z 0-9 _ - (\"/\" -> \"-\", space -> \"_\")", used[!grepl(name_ok, used)])
  add("drug_class: bad drug_class name", cls$drug_class[.isDecided(cls$drug_class) & !grepl(name_ok, cls$drug_class)])
  add("drug_abbr: code must be A-Z 0-9 -", dab$drug_abbr[.isDecided(dab$drug_abbr) & !grepl(code_ok, dab$drug_abbr)])
  add(
    "class_abbr: code must be A-Z 0-9 - (no \"_\": the resistance summary splits on it)",
    cab$class_abbr[.isDecided(cab$class_abbr) & !grepl(code_ok, cab$class_abbr)]
  )

  add("used drugs with no class", setdiff(used, cls$drug))
  add("used drugs with no abbreviation", setdiff(used, dab$drug))
  add("used classes with no abbreviation", setdiff(used_classes[.isDecided(used_classes)], cab$drug_class))

  cl_names <- unique(cls$drug_class[.isDecided(cls$drug_class)])
  as_class <- used[used %in% c(cl_names, sub("s$", "", cl_names))]
  add("a drug cannot be named after a class", setdiff(as_class, .CLASS_NAMED_DRUGS))

  hy <- used[grepl("-", used)]
  short <- vapply(strsplit(hy, "-", fixed = TRUE), function(p) any(nchar(p) < 3L), NA)
  add("hyphen is for combinations only; check", hy[short], warn = TRUE)
  add("rows for drugs no cleaned name uses (kept)",
    c(setdiff(cls$drug, used), setdiff(dab$drug, used)),
    warn = TRUE
  )

  for (w in warns) warning(w, call. = FALSE)
  if (length(errs)) {
    stop("Nothing written. Fix these in the review files:\n  ",
      paste(errs, collapse = "\n  "),
      call. = FALSE
    )
  }

  out <- list(clean_drug = clean, drug_class = cls, drug_abbr = dab, class_abbr = cab)
  dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)
  for (nm in names(out)) {
    ff <- .finalFile(nm, final_dir)
    readr::write_tsv(out[[nm]][, .DRUG_TABLES[[nm]], drop = FALSE], ff,
      na = "", eol = .tsvEol(ff)
    )
  }
  file.remove(files)
  say(
    "Wrote ", paste0(names(out), ".tsv", collapse = ", "), " to ", final_dir, ": ",
    length(used), " drugs, ", length(used_classes), " classes. ",
    "Review files removed. Now run data_raw/make_data.R to rebuild data/*.rda."
  )
  invisible(out)
}
