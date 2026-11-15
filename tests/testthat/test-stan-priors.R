stan_source <- function(family) {
  path <- testthat::test_path("..", "..", "src", "stan",
                             paste0("stanmodel-", family, ".stan"))
  if (!file.exists(path)) {
    path <- system.file("bin", "stan", paste0("stanmodel-", family, ".stan"),
                        package = "frim")
  }
  if (!nzchar(path) || !file.exists(path)) stop("Installed Stan source is missing")
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

test_that("partial stan.pars are merged with documented defaults", {
  pars <- .normalize_stan_pars(list(chains = 2))

  expect_identical(pars$chains, 2)
  expect_identical(pars$iter_warmup, 1000)
  expect_identical(pars$iter_sampling, 2000)
  expect_identical(pars$refresh, 100)
  expect_identical(pars$prior_sd_scale, 1)

  pars <- .normalize_stan_pars(list(prior_sd_scale = 2.5))
  expect_identical(pars$prior_sd_scale, 2.5)
})

test_that("prior_sd_scale is validated before fitting", {
  invalid <- list(0, -1, Inf, -Inf, NaN, NA_real_, c(1, 2), "1", NULL)
  for (value in invalid) {
    expect_error(
      .normalize_stan_pars(list(prior_sd_scale = value)),
      "finite numeric value > 0",
      fixed = TRUE
    )
  }
  expect_error(.normalize_stan_pars(1), "named list", fixed = TRUE)
  expect_error(.normalize_stan_pars(list(1)), "unique, non-empty names",
               fixed = TRUE)
  expect_error(
    frim(NULL, NULL, NULL, stan.pars = list(prior_sd_scale = 0)),
    "finite numeric value > 0",
    fixed = TRUE
  )
})

test_that("all Stan models use explicit SDs and the shared half-normal scale", {
  for (family in c("gaussian", "binomial", "poisson")) {
    stan <- stan_source(family)
    expect_true(grepl("real<lower=0>           prior_sd_scale;", stan, fixed = TRUE))
    expect_true(grepl("vector<lower=0>[L] sigma_l1", stan, fixed = TRUE))
    expect_true(grepl("vector<lower=0>[M] sigma_l2", stan, fixed = TRUE))
    expect_true(grepl("xil1_sc[,l] = xil1[,l] * sigma_l1[l];", stan, fixed = TRUE))
    expect_true(grepl("xil2_sc[,m] = xil2[,m] * sigma_l2[m];", stan, fixed = TRUE))
    expect_true(grepl("sigma_l1        ~ normal(0, prior_sd_scale);", stan, fixed = TRUE))
    expect_true(grepl("sigma_l2        ~ normal(0, prior_sd_scale);", stan, fixed = TRUE))
    expect_false(grepl("~ inv_gamma", stan, fixed = TRUE))
    expect_false(grepl("~ cauchy", stan, fixed = TRUE))
  }

  gaussian <- stan_source("gaussian")
  expect_true(grepl("real<lower=0>      sigma_e;", gaussian, fixed = TRUE))
  expect_true(grepl("sigma_e         ~ normal(0, prior_sd_scale);", gaussian, fixed = TRUE))
  expect_true(grepl("normal_lpdf(y | eta[start:end], sigma_e)", gaussian,
                    fixed = TRUE))
  expect_true(grepl("reduce_sum(partial_sum, y, grainsize, eta_vec, sigma_e)",
                    gaussian, fixed = TRUE))
  expect_false(grepl("sqrt(sigma_e)", gaussian, fixed = TRUE))

  poisson <- stan_source("poisson")
  expect_true(grepl("poisson_log_lupmf(y | eta[start:end])", poisson,
                    fixed = TRUE))
  expect_false(grepl("poisson_lupmf(y | eta[start:end])", poisson,
                     fixed = TRUE))
  expect_true(grepl("array[ylen]  int        visit;", poisson, fixed = TRUE))
  expect_true(grepl("row_vector[S]           mu_0;", poisson, fixed = TRUE))
  expect_true(grepl("matrix[S, max(nvisits)] mu_j;", poisson, fixed = TRUE))
  expect_true(grepl("eta_mat_i[i, ] = eta_mat_i[i, ] + mu_0;", poisson,
                    fixed = TRUE))
  expect_true(grepl("fixef + (mu_j[, visit[inx]])' +", poisson, fixed = TRUE))
})

test_that("Poisson eta is a log rate in the known-answer likelihood scaffold", {
  y <- c(0L, 1L, 4L)
  eta <- log(c(0.5, 2, 5))
  expected <- sum(stats::dpois(y, lambda = exp(eta), log = TRUE))
  direct <- sum(y * eta - exp(eta) - lgamma(y + 1))
  expect_equal(direct, expected, tolerance = 1e-14)
})

test_that("scaled-score names and compatibility aliases are preserved", {
  gaussian <- stan_source("gaussian")
  expect_true(grepl("matrix[Nsubs, L] xil1_sc;", gaussian, fixed = TRUE))
  expect_true(grepl("matrix[Nfuns, M] xil2_sc;", gaussian, fixed = TRUE))
  expect_true(grepl("tau1[l] = inv_square(sigma_l1[l]);", gaussian, fixed = TRUE))
  expect_true(grepl("tau2[m] = inv_square(sigma_l2[m]);", gaussian, fixed = TRUE))
  expect_true(grepl("taue = inv_square(sigma_e);", gaussian, fixed = TRUE))

  for (family in c("binomial", "poisson")) {
    stan <- stan_source(family)
    expect_true(grepl("sigmal1sq[l] = square(sigma_l1[l]);", stan, fixed = TRUE))
    expect_true(grepl("sigmal2sq[m] = square(sigma_l2[m]);", stan, fixed = TRUE))
  }
})

test_that("frim propagates prior_sd_scale into full and split Stan data", {
  source <- paste(readLines(testthat::test_path("..", "..", "R", "frim.R"),
                            warn = FALSE), collapse = "\n")
  occurrences <- gregexpr("prior_sd_scale = stan.pars$prior_sd_scale", source,
                          fixed = TRUE)[[1]]
  expect_length(occurrences[occurrences > 0], 2)
  expect_true(grepl("stan.pars <- .normalize_stan_pars(stan.pars)", source,
                    fixed = TRUE))
})
