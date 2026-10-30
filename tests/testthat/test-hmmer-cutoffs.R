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

test_that("a CasFinder profile without GA accepts any hit", {
  hits <- list(CasFinder = list(
    hit_row("p_strong", "CAS_NOGA"),
    hit_row("p_weak", "CAS_NOGA", sequence_score = 1, domain_score = 1,
            sequence_evalue = 5, i_evalue = 5)
  ))
  profiles <- list(CasFinder = profile_row("CAS_NOGA"))

  expect_setequal(run_filter(hits, profiles)$CasFinder$protein, c("p_strong", "p_weak"))
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
  profiles <- list(COG = profile_row("COG_X"))

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
  profiles <- list(COG = profile_row("COG_X"))

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
