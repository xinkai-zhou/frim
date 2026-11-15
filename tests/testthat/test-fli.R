# Before running tests, compile stan models by running instantiate::stan_package_compile()
# Make sure compilation was successful:
# stanmodel <- instantiate::stan_package_model(
#   name = paste0("stanmodel-gaussian"), package = "fli"
# )


# todo:
# 1. simulate irregularly spaced functional data and make sure the code works

if(!instantiate::stan_cmdstan_exists()) skip("Tests for frim() skipped because cmdstan is not available.")

test_that("frim-gaussian", {
  # cmdstanr::set_cmdstan_path("~/Documents/cmdstan/")
  set.seed(1)

  nsub <- 20 # number of subjects for this analysis
  pmissing <- 0.1 # percent subjects with missing
  mean.duration.missing <- 5

  dsim <- SimulateLFD(I = nsub, Jmean = 7, L = 40)
  dsim.Y <- data.frame(unclass(dsim$data$Y))
  dsim.long <- cbind(dplyr::select(dsim$data, id, visit, X), dsim.Y) %>%
    tidyr::pivot_longer(
      cols = -c(id, visit, X),  # this only allows for one fixed effect variable
      names_to = "time", values_to = "Y"
    ) %>%
    dplyr::mutate(time = as.integer(gsub("X", "", time)))

  ids.missing <- sample(unique(dsim$data$id), size = nsub * pmissing, replace = F)
  visit.missing <- rep(1, length(ids.missing))
  argvals <- unique(dsim.long$time)
  start.missing <- sample(argvals, size = length(ids.missing), replace = T)
  duration.missing <- rpois(length(ids.missing), lambda = mean.duration.missing)
  end.missing <- pmin(start.missing + duration.missing, max(argvals))
  dmissing <- data.frame(
    id = ids.missing,
    visit = visit.missing,
    missing_start = start.missing,
    missing_end   = end.missing
  ) %>% arrange(id, visit)
  datm <- left_join(dsim.long, dmissing, by = c("id", "visit")) %>%
    within({
      Y_new = Y
      Y_new[!is.na(missing_start) & time >= missing_start & time <= missing_end] = NA
    })


  frim.res <- frim(
    formula = Y ~ X + (1|id) + (1|id:visit),
    formula.fixef = ~ X,
    data = datm,
    family = "gaussian",
    argvals = 1:max(datm$time),
    bin.pars = list(
      bin_midpoints = 1:max(datm$time),
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
