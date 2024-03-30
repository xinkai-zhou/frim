if(!stan_cmdstan_exists()) skip("Tests for fli() skipped because cmdstan is not available.")

test_that("fli runs", {
  # cmdstanr::set_cmdstan_path("~/Documents/cmdstan/")
  set.seed(1)
  dsim <- SimulateLFD(I = 20, Jmean = 2, L = 40)
  dsim.Y <- data.frame(unclass(dsim$data$Y))
  dsim.long <- cbind(dplyr::select(dsim$data, id, visit, X), dsim.Y) %>%
    dplyr::pivot_longer(
      cols = -c(id, visit, X),  # this only allows for one fixed effect variable
      names_to = "time", values_to = "Y"
    ) %>%
    dplyr::mutate(time = as.integer(gsub("X", "", time)))

  # formula = Y ~ X + (1|id) + (1|id:visit)
  # formula.fixef = ~ X
  # data = dsim.long
  # family = "gaussian"
  # argvals = 1:max(dsim.long$time)
  # bin.pars = list(
  #   bin_midpoints = 1:max(dsim.long$time),
  #   bin_len = 5,
  #   cyclic = T
  # )
  # var.names = list(
  #   id_name = "id",
  #   visit_name = "visit",
  #   domain_name = "time"
  # )
  # smoothing.pars = list(
  #   nknots_min = NULL,
  #   smooth_method = "GCV.Cp",
  #   splines = "tp"
  # )
  # mfpca.pars = list(
  #   pve = 0.99,
  #   npc_l1 = 4,
  #   npc_l2 = 4
  # )
  # stan.pars = list(
  #   iter_warmup = 100,
  #   iter_sampling = 200,
  #   chains = 1,
  #   refresh = 100
  # )
  # dir.path = list(
  #   save_intermediate_res = FALSE,
  #   intermediate_res_path = "results/sim/"
  # )
  # ci.level = 0.95
  # parallel = T
  # silent = FALSE
  # num_cores = 8


  fli.res <- fli(
    formula = Y ~ X + (1|id) + (1|id:visit),
    formula.fixef = ~ X,
    data = dsim.long,
    family = "gaussian",
    argvals = 1:max(dsim.long$time),
    bin.pars = list(
      bin_midpoints = 1:max(dsim.long$time),
      bin_len = 5,
      cyclic = T
    ),
    var.names = list(
      id_name = "id",
      visit_name = "visit",
      domain_name = "time"
    ),
    smoothing.pars = list(
      nknots_min = NULL,
      smooth_method = "GCV.Cp",
      splines = "tp"
    ),
    mfpca.pars = list(
      pve = 0.99,
      npc_l1 = 4,
      npc_l2 = 4
    ),
    stan.pars = list(
      iter_warmup = 500,
      iter_sampling = 1000,
      chains = 1,
      refresh = 100
    ),
    dir.path = list(
      save_intermediate_res = FALSE,
      intermediate_res_path = "results/sim/"
    ),
    ci.level = 0.95,
    parallel = T,
    silent = FALSE,
    num_cores = 8
  )

})
