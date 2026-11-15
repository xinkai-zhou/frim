.frim_stan_pars_defaults <- function() {
  list(
    iter_warmup = 1000,
    iter_sampling = 2000,
    chains = 4,
    refresh = 100,
    prior_sd_scale = 1
  )
}

.normalize_stan_pars <- function(stan.pars) {
  if (!is.list(stan.pars)) {
    stop("stan.pars must be a named list", call. = FALSE)
  }
  if (length(stan.pars) > 0L &&
      (is.null(names(stan.pars)) || anyNA(names(stan.pars)) ||
       any(names(stan.pars) == "") || anyDuplicated(names(stan.pars)))) {
    stop("stan.pars must have unique, non-empty names", call. = FALSE)
  }

  normalized <- .frim_stan_pars_defaults()
  for (name in names(stan.pars)) {
    normalized[[name]] <- stan.pars[[name]]
  }

  prior_sd_scale <- normalized$prior_sd_scale
  if (!is.numeric(prior_sd_scale) || length(prior_sd_scale) != 1L ||
      !is.finite(prior_sd_scale) || prior_sd_scale <= 0) {
    stop("stan.pars$prior_sd_scale must be one finite numeric value > 0",
         call. = FALSE)
  }
  normalized
}
