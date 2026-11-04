test_that(".resume_plan() returns all FALSE when there is no prior run", {
  plan <- amRdata:::.resume_plan(NULL, c("panaroo", "cdhit", "hmmer"))

  expect_equal(plan, c(panaroo = FALSE, cdhit = FALSE, hmmer = FALSE))
})

test_that(".resume_plan() skips a contiguous run of successful stages", {
  prev_run <- list(
    stages = list(
      list(name = "panaroo", status = "success", outputs = list()),
      list(name = "cdhit", status = "success", outputs = list())
    )
  )

  plan <- amRdata:::.resume_plan(prev_run, c("panaroo", "cdhit", "hmmer"))

  expect_equal(plan, c(panaroo = TRUE, cdhit = TRUE, hmmer = FALSE))
})

test_that(".resume_plan() stops at the first stage that did not succeed", {
  prev_run <- list(
    stages = list(
      list(name = "panaroo", status = "success", outputs = list()),
      list(name = "cdhit", status = "failed", outputs = list()),
      list(name = "hmmer", status = "success", outputs = list())
    )
  )

  plan <- amRdata:::.resume_plan(prev_run, c("panaroo", "cdhit", "hmmer"))

  expect_equal(plan, c(panaroo = TRUE, cdhit = FALSE, hmmer = FALSE))
})

test_that(".resume_plan() distrusts a success whose recorded output is missing", {
  missing_path <- file.path(tempdir(), paste0("missing-", as.integer(runif(1, 1, 1e6))))

  prev_run <- list(
    stages = list(
      list(
        name = "panaroo",
        status = "success",
        outputs = list(list(path = missing_path))
      ),
      list(name = "cdhit", status = "success", outputs = list())
    )
  )

  plan <- amRdata:::.resume_plan(prev_run, c("panaroo", "cdhit", "hmmer"))

  expect_equal(plan, c(panaroo = FALSE, cdhit = FALSE, hmmer = FALSE))
})

test_that(".manifest_prior_stage() finds the latest entry for a stage name", {
  prev_run <- list(
    stages = list(
      list(name = "panaroo", status = "running", outputs = list()),
      list(name = "panaroo", status = "success", outputs = list())
    )
  )

  stage <- amRdata:::.manifest_prior_stage(prev_run, "panaroo")

  expect_equal(stage$status, "success")
})

test_that(".manifest_prior_stage() returns NULL when the stage isn't recorded", {
  prev_run <- list(stages = list(list(name = "panaroo", status = "success", outputs = list())))

  expect_null(amRdata:::.manifest_prior_stage(prev_run, "hmmer"))
  expect_null(amRdata:::.manifest_prior_stage(NULL, "panaroo"))
})
