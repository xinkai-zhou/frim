test_that("split-score workflow loads lowercase frim and dispatches saved standata", {
  script_path <- testthat::test_path("..", "..", "inst", "sample-scores.R")
  if (!file.exists(script_path)) script_path <- system.file("sample-scores.R", package = "frim")
  script <- paste(readLines(script_path,
                            warn = FALSE), collapse = "\n")
  expect_true(grepl("frim:::.sample_scores_part", script, fixed = TRUE))
  expect_false(grepl('package = "FRIM"', script, fixed = TRUE))
  path <- tempfile("frim-split-")
  dir.create(path)
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  writeLines("{}", file.path(path, "standata-part1.json"))

  observed <- new.env(parent = emptyenv())
  mock_fit <- list(draws = function(format) {
    observed$draw_format <- format
    data.frame(`xil1_sc[1,1]` = 0, check.names = FALSE)
  })
  mock_model <- list(sample = function(...) {
    observed$sample_args <- list(...)
    mock_fit
  })
  loader <- function(name, package) {
    observed$model_name <- name
    observed$package <- package
    mock_model
  }
  times <- as.POSIXct(c("2026-11-14 00:00:00", "2026-11-14 00:00:02"), tz = "UTC")
  clock <- local({ index <- 0L; function() { index <<- index + 1L; times[index] } })
  args <- c("sidx=1", paste0("standata_path=", path), "family=gaussian",
            "iter_warmup=10", "iter_sampling=20", "chains=2",
            "num_cores=2", "refresh=0")

  result <- .sample_scores_part(args, model_loader = loader, clock = clock)
  expect_identical(observed$package, "frim")
  expect_identical(observed$model_name, "stanmodel-gaussian")
  expect_identical(observed$sample_args$data,
                   file.path(path, "standata-part1.json"))
  expect_identical(observed$sample_args$iter_warmup, 10L)
  expect_identical(observed$sample_args$iter_sampling, 20L)
  expect_identical(observed$sample_args$chains, 2L)
  expect_identical(observed$sample_args$parallel_chains, 2L)
  expect_identical(observed$draw_format, "df")
  expect_equal(result$runtime$runtime, 2)
  expect_true(all(file.exists(file.path(
    path, c("runtime.csv", "stanfit-part1.rds", "draws-part1.rds")
  ))))
})

test_that("split-score CLI integers use strict canonical lexical form", {
  base <- c(sidx = "1", standata_path = "unused", family = "gaussian",
            iter_warmup = "10", iter_sampling = "20", chains = "2",
            num_cores = "2", refresh = "0")
  for (bad in c("1.5", "1e2", "+1", "01", "-1", "2147483648")) {
    changed <- as.list(base)
    changed$chains <- bad
    expect_error(.sample_scores_integer_value(changed, "chains", 1L),
                 "canonical|integer range")
  }
  expect_identical(.sample_scores_integer_value(as.list(base), "refresh"), 0L)
  duplicate <- c(paste0(names(base), "=", base), "chains=3")
  expect_error(.sample_scores_parse_args(duplicate), "Duplicate argument: chains",
               fixed = TRUE)
})
