hit_row <- function(protein, query_name, sequence_score = 100, domain_score = 100,
                    sequence_evalue = 1e-30, i_evalue = 1e-30,
                    env_from = 1L, env_to = 100L) {
  tibble::tibble(
    protein = protein,
    query_name = query_name,
    query_accession = NA_character_,
    target_description = NA_character_,
    tlen = 300L,
    qlen = 100L,
    sequence_evalue = sequence_evalue,
    sequence_score = sequence_score,
    domain_num = 1L,
    i_evalue = i_evalue,
    domain_score = domain_score,
    hmm_from = 1L,
    hmm_to = 100L,
    ali_from = env_from,
    ali_to = env_to,
    env_from = env_from,
    env_to = env_to,
    profile_coverage = 1,
    target_coverage = 0.3
  )
}

profile_row <- function(profile_name, ga = c(NA_real_, NA_real_),
                        tc = c(NA_real_, NA_real_)) {
  tibble::tibble(
    profile_name = profile_name,
    profile_accession = NA_character_,
    profile_description = paste("Description of", profile_name),
    ga_seq = ga[[1L]], ga_dom = ga[[2L]],
    tc_seq = tc[[1L]], tc_dom = tc[[2L]],
    nc_seq = NA_real_, nc_dom = NA_real_
  )
}

# Runs .filterHMMER() on hits and profiles given per searched database, and
# returns the annotation table of each requested database
run_filter <- function(hits, profiles, databases = names(hits), ...) {
  dir <- tempfile("hmmer-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  for (name in names(hits)) {
    arrow::write_parquet(
      dplyr::bind_rows(hits[[name]]),
      file.path(dir, paste0("protein_", name, "_hits.parquet"))
    )
  }

  profile_files <- vapply(names(profiles), function(name) {
    path <- file.path(dir, paste0(name, ".profiles.parquet"))
    arrow::write_parquet(dplyr::bind_rows(profiles[[name]]), path)
    path
  }, character(1))

  duckdb_path <- file.path(dir, "test.duckdb")
  con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
  DBI::dbDisconnect(con, shutdown = TRUE)

  parquets <- .filterHMMER(
    duckdb_path, dir, databases, profile_files, verbose = FALSE, ...
  )

  lapply(parquets, arrow::read_parquet)
}

test_that(".parse_hmmer_profiles() reads GA, TC and NC cutoffs", {
  hmm <- tempfile(fileext = ".hmm")
  writeLines(c(
    "HMMER3/f [3.1b2 | February 2015]",
    "NAME  PF00001",
    "ACC   PF00001.1",
    "DESC  First family",
    "GA    27.00 25.50;",
    "TC    27.10 25.60;",
    "NC    26.90 25.40;",
    "HMM          A        C",
    "  COMPO   2.33150  4.35991",
    "//",
    "HMMER3/f [3.1b2 | February 2015]",
    "NAME  COG0001",
    "DESC  No cutoffs",
    "HMM          A        C",
    "//"
  ), hmm)

  profiles <- .parse_hmmer_profiles(hmm)

  expect_equal(profiles$profile_name, c("PF00001", "COG0001"))
  expect_equal(profiles$ga_seq, c(27, NA))
  expect_equal(profiles$ga_dom, c(25.5, NA))
  expect_equal(profiles$tc_seq, c(27.1, NA))
  expect_equal(profiles$nc_dom, c(25.4, NA))
})

test_that(".hmmerHits() adds coverage and keeps the envelope columns", {
  hits <- tibble::tibble(
    protein = "p1", protein_accession = NA_character_, tlen = 400L,
    query_name = "PF1", query_accession = NA_character_, qlen = 200L,
    sequence_evalue = 1e-10, sequence_score = 50, sequence_bias = 0,
    domain_num = 1L, domain_of = 1L, c_evalue = 1e-10, i_evalue = 1e-10,
    domain_score = 50, domain_bias = 0, hmm_from = 11L, hmm_to = 110L,
    ali_from = 101L, ali_to = 300L, env_from = 100L, env_to = 301L,
    acc = 0.9, target_description = NA_character_
  )

  result <- .hmmerHits(hits)

  expect_true(all(c("env_from", "env_to", "tlen", "qlen") %in% names(result)))
  expect_equal(result$profile_coverage, 100 / 200)
  expect_equal(result$target_coverage, 200 / 400)
})

test_that("GA cutoffs keep rounded scores at the cutoff and drop lower ones", {
  hits <- list(Pfam = list(
    hit_row("p_pass", "PF1", sequence_score = 30, domain_score = 30),
    # printed score 23.5 may be a true 23.46, which still rounds to the cutoff
    hit_row("p_tie", "PF1", sequence_score = 30, domain_score = 23.5),
    hit_row("p_low_dom", "PF1", sequence_score = 30, domain_score = 23.4),
    hit_row("p_low_seq", "PF1", sequence_score = 23.4, domain_score = 30)
  ))
  profiles <- list(Pfam = profile_row("PF1", ga = c(23.5, 23.5)))

  result <- run_filter(hits, profiles)

  expect_setequal(result$Pfam$protein, c("p_pass", "p_tie"))
})

test_that("curated cutoffs decide regardless of the E-value", {
  hits <- list(Pfam = list(
    hit_row("p1", "PF1", sequence_score = 40, domain_score = 40,
            sequence_evalue = 0.5, i_evalue = 0.5)
  ))
  profiles <- list(Pfam = profile_row("PF1", ga = c(30, 30)))

  expect_equal(run_filter(hits, profiles)$Pfam$protein, "p1")
})

test_that("AMRFinder uses the TC cutoff", {
  hits <- list(AMRFinder = list(
    hit_row("p_pass", "NF1", sequence_score = 300, domain_score = 300),
    hit_row("p_fail", "NF1", sequence_score = 250, domain_score = 250)
  ))
  profiles <- list(AMRFinder = profile_row("NF1", ga = c(100, 100), tc = c(275, 275)))

  expect_equal(run_filter(hits, profiles)$AMRFinder$protein, "p_pass")
})

test_that("a CasFinder profile without GA uses the fallback", {
  hits <- list(CasFinder = list(
    hit_row("p_strong", "CAS_NOGA"),
    hit_row("p_weak", "CAS_NOGA", sequence_score = 1, domain_score = 1,
            sequence_evalue = 5, i_evalue = 5)
  ))
  profiles <- list(CasFinder = profile_row("CAS_NOGA"))

  expect_setequal(run_filter(hits, profiles)$CasFinder$protein, "p_strong")
})

test_that("overlapping hits to related profiles are all kept for curated databases", {
  hits <- list(Pfam = list(
    hit_row("p1", "PF_A", sequence_score = 90, domain_score = 90),
    hit_row("p1", "PF_B", sequence_score = 40, domain_score = 40)
  ))
  profiles <- list(Pfam = list(
    profile_row("PF_A", ga = c(20, 20)),
    profile_row("PF_B", ga = c(20, 20))
  ))

  expect_setequal(run_filter(hits, profiles)$Pfam$query_name, c("PF_A", "PF_B"))
})

test_that("COG hits are judged on the E-value, which is a parameter", {
  hits <- list(COG = list(
    hit_row("p_weak_domain", "COG1", sequence_evalue = 1e-6, i_evalue = 0.01),
    hit_row("p_weak_seq", "COG1", sequence_evalue = 0.01, i_evalue = 1e-6),
    hit_row("p_ok", "COG1", sequence_evalue = 1e-6, i_evalue = 1e-6),
    # exactly NCBI's cutoff of 1e-4 is accepted, 1e-3 is not
    hit_row("p_ncbi", "COG1", sequence_evalue = 1e-4, i_evalue = 1e-4),
    hit_row("p_loose", "COG1", sequence_evalue = 1e-3, i_evalue = 1e-3)
  ))
  profiles <- list(COG = profile_row("COG1"))

  expect_setequal(run_filter(hits, profiles)$COG$protein, c("p_ok", "p_ncbi"))
  expect_equal(nrow(run_filter(hits, profiles, cog_evalue = 0.1)$COG), 5L)
  expect_error(run_filter(hits, profiles, cog_evalue = 0), "greater than 0")
})

test_that("COG keeps the best hit per region and both COGs of a fused protein", {
  hits <- list(COG = list(
    # one region covered by three COGs: only the lowest E-value is kept
    hit_row("p_same", "COG_A", i_evalue = 1e-40, sequence_evalue = 1e-40, env_from = 1L, env_to = 200L),
    hit_row("p_same", "COG_B", i_evalue = 1e-20, sequence_evalue = 1e-20, env_from = 5L, env_to = 195L),
    hit_row("p_same", "COG_C", i_evalue = 1e-10, sequence_evalue = 1e-10, env_from = 1L, env_to = 200L),
    # two domains sharing 2 residues at the boundary, as in a fused enzyme
    hit_row("p_fused", "COG_D", i_evalue = 1e-50, sequence_evalue = 1e-50, env_from = 1L, env_to = 464L),
    hit_row("p_fused", "COG_E", i_evalue = 1e-30, sequence_evalue = 1e-30, env_from = 463L, env_to = 815L)
  ))
  profiles <- list(COG = lapply(unique(vapply(hits$COG, function(h) h$query_name, character(1))), profile_row))

  result <- run_filter(hits, profiles)$COG

  expect_equal(result$query_name[result$protein == "p_same"], "COG_A")
  expect_setequal(result$query_name[result$protein == "p_fused"], c("COG_D", "COG_E"))
})

test_that("the COG overlap allowance is a parameter", {
  hits <- list(COG = list(
    hit_row("p1", "COG_A", i_evalue = 1e-40, sequence_evalue = 1e-40, env_from = 1L, env_to = 100L),
    # shares 30 of its 100 residues with the better hit
    hit_row("p1", "COG_B", i_evalue = 1e-20, sequence_evalue = 1e-20, env_from = 71L, env_to = 170L)
  ))
  profiles <- list(COG = lapply(unique(vapply(hits$COG, function(h) h$query_name, character(1))), profile_row))

  expect_equal(nrow(run_filter(hits, profiles)$COG), 1L)
  expect_equal(nrow(run_filter(hits, profiles, cog_overlap = 0.3)$COG), 2L)
})

test_that("DefenseCas combines the DefenseFinder and CasFinder searches", {
  hits <- list(
    DefenseFinder = list(hit_row("p1", "DF1", sequence_score = 30, domain_score = 30)),
    CasFinder = list(hit_row("p2", "CAS1", sequence_score = 30, domain_score = 30))
  )
  profiles <- list(
    DefenseFinder = profile_row("DF1", ga = c(20, 20)),
    CasFinder = profile_row("CAS1", ga = c(20, 20))
  )

  result <- run_filter(hits, profiles, databases = "DefenseCas")$DefenseCas

  expect_equal(result$protein, c("p1", "p2"))
  expect_equal(result$database, c("DefenseFinder", "CasFinder"))
  expect_equal(result$description, c("Description of DF1", "Description of CAS1"))
})

test_that(".filterHMMER() writes the annotation columns", {
  hits <- list(Pfam = list(hit_row("p1", "PF1")))
  profiles <- list(Pfam = profile_row("PF1", ga = c(20, 20)))

  result <- run_filter(hits, profiles)$Pfam

  expect_named(
    result,
    c("protein", "query_name", "query_accession", "target_description",
      "i_evalue", "domain_score", "description")
  )
})


test_that("missing and incomplete profile metadata are errors", {
  hits <- list(Pfam = list(hit_row("p1", "PF_MISSING")))
  profiles <- list(Pfam = profile_row("PF_OTHER", ga = c(20, 20)))
  expect_error(run_filter(hits, profiles), "Unknown HMMER profile.*Pfam.*PF_MISSING")
  profiles <- list(Pfam = profile_row("PF_MISSING", ga = c(20, NA_real_)))
  expect_error(run_filter(hits, profiles), "Incomplete curated cutoff pair")
})

test_that("the fallback checks sequence E-value, domain E-value and coverage", {
  good <- hit_row("boundary", "CAS", sequence_evalue = 1e-3, i_evalue = 1e-3)
  good$profile_coverage <- 0.5
  short <- hit_row("short", "CAS")
  short$profile_coverage <- 0.49
  hits <- list(CasFinder = list(good, short,
    hit_row("weak_sequence", "CAS", sequence_evalue = 0.01),
    hit_row("weak_domain", "CAS", i_evalue = 0.01)))
  profiles <- list(CasFinder = profile_row("CAS"))
  expect_equal(run_filter(hits, profiles)$CasFinder$protein, "boundary")
  expect_equal(nrow(run_filter(hits, profiles, fallback_evalue = 0.01,
                              fallback_coverage = 0.49)$CasFinder), 4L)
  good$profile_coverage <- NA_real_
  expect_equal(nrow(run_filter(list(CasFinder = list(good)), profiles)$CasFinder), 0L)
})

test_that("curated cutoffs are not subjected to fallback coverage or E-values", {
  hit <- hit_row("curated", "PF1", sequence_score = 40, domain_score = 40,
                 sequence_evalue = 50, i_evalue = 50)
  hit$profile_coverage <- 0.1
  expect_equal(run_filter(list(Pfam = list(hit)),
    list(Pfam = profile_row("PF1", ga = c(30, 30))))$Pfam$protein, "curated")
})

test_that("fallback settings are validated before workflow side effects", {
  for (value in list(NA_real_, Inf, 0, -1, c(1, 2))) {
    expect_error(runDataProcessing(hmmer_fallback_evalue = value),
                 "hmmer_fallback_evalue")
  }
  for (value in list(NA_real_, Inf, 0, -1, 1.01, c(0.5, 0.6))) {
    expect_error(runDataProcessing(hmmer_fallback_coverage = value),
                 "hmmer_fallback_coverage")
  }
})


test_that("HMMER job requests score-based reporting and saves high-E-value candidates", {
  out <- tempfile("hmmer-command-")
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  captured <- NULL
  run <- .runHmmerJob
  boundary <- new.env(parent = environment(run))
  boundary$system2 <- function(command, args, ...) {
    captured <<- args
    writeLines(paste(
      "p1 - 300 PF1 - 100 50 40 0 1 1 50 50 40 0",
      "1 100 1 100 1 100 0.9 candidate"
    ), file.path(out, "job.tbl"))
    0L
  }
  environment(run) <- boundary
  saved <- run("job", "input.faa", "Pfam", total_proteins = 1000L,
    output_path = out, db_paths = list(Pfam = list(hmm = file.path(out, "db.hmm"))))
  expect_true(all(c("-T", "--domT") %in% captured))
  expect_lt(as.numeric(captured[match("-T", captured) + 1L]), 0)
  expect_lt(as.numeric(captured[match("--domT", captured) + 1L]), 0)
  hits <- arrow::read_parquet(saved)
  expect_equal(hits$sequence_evalue, 50)
  expect_equal(hits$i_evalue, 50)
  profiles <- list(Pfam = profile_row("PF1", ga = c(30, 30)))
  expect_equal(run_filter(list(Pfam = list(hits)), profiles)$Pfam$protein, "p1")
})

test_that("HMMER2duckdb() searches, then filters with the collected profiles", {
  db <- tempfile(fileext = ".duckdb")
  file.create(db)
  out <- tempfile("hmmer-")
  dir.create(out)
  calls <- character()
  filter_args <- NULL
  testthat::local_mocked_bindings(
    .runHMMER = function(..., databases) {
      calls <<- c(calls, paste("search:", paste(databases, collapse = ",")))
      list(databases = list(Pfam = list(profiles = "pfam.parquet")))
    },
    .defenseHMMER = function(...) {
      calls <<- c(calls, "defense")
      list(databases = list(DefenseFinder = list(profiles = "df.parquet")))
    },
    .filterHMMER = function(...) {
      calls <<- c(calls, "filter")
      filter_args <<- list(...)
      c(Pfam = "protein_Pfam.parquet")
    }
  )

  res <- HMMER2duckdb(
    db, output_path = out, databases = c("Pfam", "DefenseCas"),
    hmmer_db_dir = tempdir(), verbose = FALSE
  )

  expect_equal(calls, c("search: Pfam", "defense", "filter"))
  expect_equal(filter_args$databases, c("Pfam", "DefenseCas"))
  expect_equal(
    unname(filter_args$profiles), c("pfam.parquet", "df.parquet")
  )
  expect_equal(names(filter_args$profiles), c("Pfam", "DefenseFinder"))
  expect_equal(res$annotations, c(Pfam = "protein_Pfam.parquet"))
  expect_named(res$databases, c("generic", "DefenseCas"))
})
