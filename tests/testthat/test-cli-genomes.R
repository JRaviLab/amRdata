test_that("CLI downloads keep the first ID and validate zero-exit output files", {
  out <- tempfile()
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  ids <- c("624.102", "624.104")
  emitted <- FALSE
  calls <- 0L
  requested_ids <- ids
  local_mocked_bindings(
    .prepareContainerImages = function(config, ...) config,
    .pick_shell = function(...) "bash",
    .runContainerCommand = function(config, args, ...) {
      calls <<- calls + 1L
      expect_match(tail(args, 1L), "--nohead --col 1", fixed = TRUE)
      expect_identical(readLines(file.path(out, "ids_test.txt")), requested_ids)
      if (emitted) return("service returned no genomes")
      for (id in ids) {
        writeLines(c(paste0(">", id), paste(rep("A", 200), collapse = "")), file.path(out, paste0(id, ".fasta")))
        writeLines(c(paste0(">", id), paste(rep("M", 200), collapse = "")), file.path(out, paste0(id, ".faa")))
        writeLines(paste(rep("x", 200), collapse = ""), file.path(out, paste0(id, ".gto")))
      }
      "downloaded"
    }
  )
  expect_true(.cli_dump_fastas_gto_chunk("image", out, ids, "test", tries = 1L))
  expect_true(all(file.exists(file.path(out, paste0(ids, ".fna")))))
  expect_true(all(file.exists(file.path(out, paste0(ids, ".PATRIC.faa")))))
  expect_true(.cli_dump_fastas_gto_chunk("image", out, ids, "test", tries = 1L))
  expect_identical(calls, 1L)
  unlink(file.path(out, paste0(ids[[1]], ".gto")))
  requested_ids <- ids[[1]]
  emitted <- TRUE
  expect_error(.cli_dump_fastas_gto_chunk("image", out, ids, "test", tries = 1L), "624.102")
  expect_true(any(grepl("service returned no genomes", readLines(file.path(out, "download_chunk_test.log")))))
})

.createGFFFixture <- function() {
  out <- tempfile()
  dir.create(out)
  id <- "624.102"
  jsonlite::write_json(list(contigs = list(list(id = "internal", original_id = "original"))),
    file.path(out, paste0(id, ".gto")), auto_unbox = TRUE)
  exporter <- file.path(out, "fake-exporter")
  writeLines(c(
    "#!/bin/sh",
    'while getopts "i:o:" option; do case "$option" in o) output="$OPTARG";; esac; done',
    "printf '##gff-version 3\\naccn|original\\tRAST\\tCDS\\t1\\t9\\t.\\t+\\t0\\tID=feature1\\n' > \"$output\""
  ), exporter)
  Sys.chmod(exporter, "0755")
  list(out = out, id = id, exporter = exporter)
}

test_that("GFF chunk script exports and remaps a local fixture", {
  skip_on_os("windows")
  fixture <- .createGFFFixture()
  on.exit(unlink(fixture$out, recursive = TRUE), add = TRUE)
  local_mocked_bindings(
    .prepareContainerImages = function(config, ...) config,
    .pick_shell = function(...) "bash",
    .runContainerCommand = function(config, args, ...) {
      script <- tail(args, 1L)
      script <- gsub("/usr/bin/rast-export-genome", fixture$exporter, script, fixed = TRUE)
      script <- gsub("/out/", paste0(fixture$out, "/"), script, fixed = TRUE)
      system2("/bin/bash", c("-lc", shQuote(script)), stdout = TRUE, stderr = TRUE)
    }
  )
  expect_true(.cli_export_gff_chunk("image", fixture$out, fixture$id, "test", tries = 1L))
  gff <- readLines(file.path(fixture$out, paste0(fixture$id, ".PATRIC.gff")))
  expect_match(gff[[2]], "^internal\\t")
})

test_that("GFF chunk failures retain stderr and cannot report success", {
  skip_on_os("windows")
  fixture <- .createGFFFixture()
  on.exit(unlink(fixture$out, recursive = TRUE), add = TRUE)
  writeLines(c("#!/bin/sh", "echo 'annotation service unavailable' >&2", "exit 1"), fixture$exporter)
  local_mocked_bindings(
    .prepareContainerImages = function(config, ...) config,
    .pick_shell = function(...) "bash",
    .runContainerCommand = function(config, args, ...) {
      script <- gsub("/usr/bin/rast-export-genome", fixture$exporter, tail(args, 1L), fixed = TRUE)
      script <- gsub("/out/", paste0(fixture$out, "/"), script, fixed = TRUE)
      suppressWarnings(system2("/bin/bash", c("-lc", shQuote(script)), stdout = TRUE, stderr = TRUE))
    }
  )
  expect_error(.cli_export_gff_chunk("image", fixture$out, fixture$id, "test", tries = 1L), "exit status 1")
  expect_true(any(grepl("annotation service unavailable", readLines(file.path(fixture$out, "gff_chunk_test.log")))))
  expect_match(readLines(file.path(fixture$out, "gff_chunk_test.stderr.txt")), "EXPORT_FAIL_GFF")
  unlink(file.path(fixture$out, paste0(fixture$id, ".gto")))
  expect_error(.cli_export_gff_chunk("image", fixture$out, fixture$id, "missing", tries = 1L), "exit status 1")
  expect_match(readLines(file.path(fixture$out, "gff_chunk_missing.stderr.txt")), "MISSING_GTO")
})
