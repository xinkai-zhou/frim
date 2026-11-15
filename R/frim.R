#' Functional Random Effect Inference for Generalized Functional Mixed Models
#'
#'
#' The frim approach comprises of the following steps:
#' \enumerate{
#' \item Fit local GLMMs, smooth local fixed effect estimates, and extract local random effects;
#' \item Apply multilevel FPCA on the adjusted local random effects to obtain eigenfunctions;
#' \item Fit a global Bayesian multilevel model to sample scores.
#' }
#'
#' For more information on each step, please refer to the FRIM paper.
#'
#' @param formula Two-sided formula object for the functional mixed model in the lme4 formula syntax.
#' The visit variable should be included as a fixed effect covariate.
#' @param formula.fixef Formula for just the fixed effect.
#' @param data A data frame in the long format.
#' @param family Distribution family of the response. Default to \code{gaussian}.
#' @param argvals A vector containing locations of observations on the
#' functional domain. If not specified, a regular grid across the range of
#' the domain is assumed. Currently only supported for bootstrap (\code{analytic=FALSE}).
#' @param bin.pars A named list specifying the bin midpoints (\code{bin_midpoints}), bin length (\code{bin_len}),
#' and whether bins are cyclic.
#' @param var.names A named list specifying the names of the id, visit, and domain name in the \code{data} object.
#' @param smoothing.pars A named list with components:
#' \code{nknots_min}: the minimal number of knots in the penalized smoothing for the regression coefficients.
#' Defaults to \code{NULL}, which then uses L/2 where L is the dimension of the functional domain.
#' \code{smooth_method}: method for selecting smoothing parameter in step 2. Defaults to \code{"GCV.Cp"}
#' \code{splines}: spline type used for penalized splines smoothing. We use the same syntax as the mgcv package. Defaults to \code{"tp"}
#' @param mfpca.pars A named list with components:
#' \code{pve}: the proportion of variance explained by the eigenfunctions. Defaults to \code{0.99}.
#' \code{npc_l1}: number of eigenfunctions for level 1. Defaults to \code{4}.
#' \code{npc_l2}: number of eigenfunctions for level 2. Defaults to \code{4}.
#' \code{visit.levels}: levels of the visit variable. Defaults to \code{NULL}.
#' \code{twoway}: logical, indicating whether to use two-way MFPCA. Defaults to \code{FALSE}.
#' @param stan.pars A named list with components:
#' \code{iter_warmup}: number of warmup iterations. Defaults to \code{1000}.
#' \code{iter_sampling}: number of sampling iterations. Defaults to \code{2000}.
#' \code{chains}: number of chains. Defaults to \code{4}.
#' \code{refresh}: number of iterations between progress updates. Defaults to \code{100}.
#' \code{prior_sd_scale}: shared scale of the half-normal priors on level-1
#' score SDs, level-2 score SDs, and the Gaussian residual SD. Must be one
#' finite value greater than zero. Defaults to \code{1}.
#' @param dir.path A named list with components:
#' \code{save_intermediate_res}: logical, indicating whether to save intermediate results. Defaults to \code{FALSE}.
#' \code{intermediate_res_path}: path to save intermediate results. Defaults to \code{"results/sim/"}.
#' @param ci.level Confidence level for credible intervals. Defaults to \code{0.95}.
#' @param parallel Logical, indicating whether to fit local GLMMs in parallel. Defaults to \code{TRUE}.
#' @param silent Logical, indicating whether to show descriptions of each step. Defaults to \code{TRUE}.
#' @param num_cores Number of cores for parallel fitting of local GLMMs. Defaults to 1.
#' @param use_complete_visits Logical, indicating whether to use visits with complete observations for estimating eigenfunctions. Defaults to \code{TRUE}.
#' @param skip_STAN Logical, indicating whether to skip STAN sampling. Defaults to \code{FALSE}.
#' @param n.subsets Number of subsets to divide the full data for STAN sampling. Ignored if \code{skip_STAN} is FALSE. Defaults to \code{5}.
#'
#' @return A list containing:
#' \item{argvals}{A vector containing locations of observations on the functional domain.}
#' \item{bin_idx}{A list of bins.}
#' \item{id_visit}{A data frame containing the id and visit variables.}
#' \item{betaTilde}{Matrix of local fixed effect estimates.}
#' \item{betaHat}{Matrix of smoothed fixed effect estimates.}
#' \item{localglmm}{List of local GLMM fits.}
#' \item{evalues}{Eigenvalues of the MFPCA.}
#' \item{efunctions}{List of eigenfunctions.}
#' \item{scores_posterior}{List of posterior samples of scores.}
#' \item{scores_postmean}{List of posterior means of scores.}
#' \item{stanfit}{STAN fit object.}
#' \item{ci.level}{Confidence level for credible intervals.}
#'
#'
#' @import dplyr
#' @import cmdstanr
#'
#' @export frim
#'
#' @examples
#' set.seed(1)
#' fit_nhanes <- frim(value ~ age + gender + (1|id) + (1|id:visit),
#'                  data = nhanes_example)

frim <- function(
    formula,
    formula.fixef,
    data,
    family = "gaussian",
    argvals = NULL,
    bin.pars = list(
      bin_midpoints = NULL,
      bin_len = 20, # given in terms of the number of argvals
      cyclic
    ),
    var.names = list(
      id_name = "id",
      visit_name = "visit",
      domain_name = "time"
    ),
    smoothing.pars = list(
      nknots_min = 35,
      smooth_method = "GCV.Cp",
      splines = "tp"
    ),
    mfpca.pars = list(
      pve = 0.99,
      npc_l1 = 4,
      npc_l2 = 4,
      visit.levels = NULL,
      twoway = F
    ),
    stan.pars = list(
      iter_warmup = 1000,
      iter_sampling = 2000,
      chains = 4,
      refresh = 100,
      prior_sd_scale = 1
    ),
    dir.path = list(
      save_intermediate_res = FALSE,
      intermediate_res_path = "results/sim/"
    ),
    ci.level = 0.95,
    parallel = TRUE,
    silent = FALSE,
    num_cores = 4,
    use_complete_visits = TRUE,
    skip_STAN = FALSE,
    n.subsets = 5
    ){

  # Merge partial user settings with defaults and validate the prior before any
  # local-model fitting or file output.
  stan.pars <- .normalize_stan_pars(stan.pars)

  # If doing parallel computing, set up the number of cores
  if(parallel & !is.integer(num_cores)) num_cores <- parallel::detectCores() - 1

  # Organize the input from the model formula
  model_formula <- as.character(formula)
  stopifnot(model_formula[1] == "~" & length(model_formula) == 3)

  # stop function if there are column names with "." to avoid issues with covariance G() and H() calculations below
  dep_str <- deparse(model_formula[3])
  if(grepl(".", dep_str, fixed = TRUE)){
    # make sure it isn't just a call to all covariates with "Y ~. "
    dep_str_rm <- substr(dep_str, 3, nchar(dep_str)) # remove first character of parsed formula string and check
    if(grepl(".", dep_str_rm, fixed = TRUE)){
      stop('Remove the character "." from all non-functional covariate names and rerun fui() function
           -For example, change the name "X.1" to "X_1"
           -The string "." *should* be kept in the functional outcome names (e.g., "Y.1" *is* proper naming).')
    }
  }
  rm(dep_str)

  if(!is.factor(data[[var.names$visit_name]]) & is.null(mfpca.pars$visit.levels)){
    stop(paste0("Convert ", var.names$visit_name, " to a factor, or specify the visit.levels argument from mfpca.pars."))
  }
  # Convert visit to a factor because MFPCA arranges visit-specific means (eta)
  # as columns in the same order as the levels of the visit variable.
  if(!is.factor(data[[var.names$visit_name]])){
    # data[[var.names$visit_name]] <- as.character(data[[var.names$visit_name]])
    data[[var.names$visit_name]] <-
      factor(data[[var.names$visit_name]], levels = mfpca.pars$visit.levels)
  }

  # create a folder for saving intermediate results
  if(dir.path$save_intermediate_res){
    res.path <- dir.path$intermediate_res_path
    if(!file.exists(res.path)) dir.create(res.path)
  }

  # Filter out NA in the data set
  # data <- filter(data, !is.na(.data[[model_formula[[2]]]]))

  # Convert id to character for merging with random effect df later
  id.unique <- unique(data[[var.names$id_name]])
  data[[var.names$id_name]] <- factor(data[[var.names$id_name]], levels = id.unique)
  # data[[var.names$id_name]] <- as.character(data[[var.names$id_name]])

  # Arrange data in the order of id, visit, time (domain)
  data <- data %>% arrange(
    .data[[var.names$id_name]],
    .data[[var.names$visit_name]],
    .data[[var.names$domain_name]]
  )

  # Get a reference data frame that contains the unique id and visit combo.
  # When local bins do not have data for a specific id-visit, the reference df
  # will be used to mark the corresponding random effect as NA.
  dref <- data %>%
    dplyr::select(all_of(c(var.names$id_name, var.names$visit_name))) %>%
    dplyr::distinct()
  # colnames(dref) <- c("id", "visit")

  ##############################################################################
  ## Step 1 Create data bins
  ##############################################################################
  if(!silent) print("Step 1: Fit Massively Univariate Mixed Models")

  # L: dimension of the functional domain
  L <- max(data[[var.names$domain_name]])
  if(is.null(argvals)){
    argvals <- 1:L
  }else{
    # when only using a subset of the functional domain
    if(max(argvals) > L) stop("Maximum index specified in argvals is greater than the total number of columns for the functional outcome")
    L <- length(argvals)
  }
  if(length(bin.pars$bin_midpoints) > 400){
    message("Your functional data is dense! Consider subsampling along the functional domain (i.e., reduce columns of outcome matrix).")
  }

  # Create bins
  if(is.null(bin.pars$bin_midpoints)){
    bin.pars$bin_midpoints <- c(
      seq(argvals[1], argvals[length(argvals)], by = 10),
      argvals[length(argvals)]
    )
  }
  bin.index <- CreateBins(
    argvals = argvals,
    bin_midpoints = bin.pars$bin_midpoints,
    bin_len = bin.pars$bin_len,
    cyclic = bin.pars$cyclic
  )
  # Count the number of visits per subject
  nvisits <- data %>% dplyr::summarise(
    nvisits = n_distinct(.data[[var.names$visit_name]]),
    .by = all_of(var.names$id_name)
  )

  ##############################################################################
  ## Step 2 Fit mixed models
  ##############################################################################
  ### Create a function that fit a mixed model in bin l
  ### Input: l: index of the list bin.index.
  ###           Each element in the list is a vector containing domain points for that bin.
  ### Output: A list containing point estimates, variance estimates, etc.
  unimm <- function(l){
    dtemp <- data %>% dplyr::filter(.data[[var.names$domain_name]] %in% l)
    if(family == "gaussian"){
      fit_uni <- suppressMessages(
        lme4::lmer(
          formula = formula,
          data = dtemp,
          control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 5000))
        )
      )
    } else {
      fit_uni <- suppressMessages(
        lme4::glmer(
          formula = formula,
          data = dtemp,
          family = family,
          control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 5000))
        )
      )
    }
    betaTilde <- lme4::fixef(fit_uni) ## fixed effects estimates
    re.ls <- lme4::ranef(fit_uni) ## random effects
    id.visit <- strsplit(rownames(re.ls[[paste0(var.names$id_name, ":", var.names$visit_name)]]), ":")
    re.df.id <- data.frame(
      rownames(re.ls[[var.names$id_name]]),
      re.ls[[var.names$id_name]][[1]]
    ) %>% setNames(c(var.names$id_name, "re_id"))
    re.df.id.visit <- data.frame(
      sapply(id.visit, function(x) x[1]),
      sapply(id.visit, function(x) x[2]),
      re.ls[[paste0(var.names$id_name, ":", var.names$visit_name)]][[1]]
    ) %>% setNames(c(var.names$id_name, var.names$visit_name, "re_id_visit"))
    re.df <- dref %>%
      dplyr::left_join(re.df.id, by = var.names$id_name) %>%
      dplyr::left_join(re.df.id.visit, by = c(var.names$id_name, var.names$visit_name))
    # if(!impute_NA){
    #   # Approach 1: if re_id_visit is missing, then re is missing
    #   re.df <- re.df %>% mutate(re = re_id + re_id_visit)
    # } else {
    #   # Approach 2: missing values in re_id_visit are imputed with 0
    #   # If data is completely missing for a visit in a bin, then re_id_visit is NA
    #   # The best prediction for re_id_visit is 0.
    #   # i.e., we use re_id as the subject-visit-specific ranef prediction for that bin.
    re.df <- re.df %>%
        mutate(re = ifelse(is.na(re_id_visit), re_id, re_id + re_id_visit))
    # }


    # se_mat <- summary(fit_uni)$coefficients[, "Std. Error"] ## se of fixed effects
    model_matrix = dref %>%
      left_join(cbind(
        dtemp %>%
          filter(!is.na(.data[[model_formula[[2]]]])) %>%
          dplyr::select(, all_of(c(var.names$id_name, var.names$visit_name))),
        model.matrix(fit_uni)
      ), by = c(var.names$id_name, var.names$visit_name)) %>% dplyr::distinct() %>%
      dplyr::select(- all_of(c(var.names$id_name, var.names$visit_name))) %>%
      as.matrix()
    return(
      list(
        betaTilde = betaTilde,
        # se_mat    = se_mat,
        model_matrix = model_matrix,
        re_names  = paste0(re.df[[var.names$id_name]], ":", re.df[[var.names$visit_name]]),
        #rownames(re.df[[paste0(var.names$id_name, ":", var.names$visit_name)]]),
        re_vals   = re.df$re
        #rep(re.df[[var.names$id_name]][[1]], times = nvisits$nvisits) +
        #re.df[[paste0(var.names$id_name, ":", var.names$visit_name)]][[1]]
      )
    )
  }
  # Fit massively univariate mixed models
  if(parallel == TRUE){
    massmm <- parallel::mclapply(bin.index, unimm, mc.cores = num_cores)
  }else{
    massmm <- lapply(bin.index, unimm)
  }
  if(dir.path$save_intermediate_res) saveRDS(massmm, paste0(res.path, "univariate-mixed-models.rds"))
  # massmm <- read_rds(paste0(res.path, "univariate-mixed-models.rds"))
  # massmm <- read_rds("results/sim/univariate-mixed-models.rds")

  # Obtain betaTilde, fixed effects estimates
  betaTilde <- t(do.call(rbind, lapply(massmm, '[[', "betaTilde")))
  # betaSE <- t(do.call(rbind, lapply(massmm, '[[', "se_mat")))
  colnames(betaTilde) <- names(bin.index)
  p <- nrow(betaTilde)

  # Obtain residuals, AIC, BIC, and random effects estimates (analytic)
  ## AIC/BIC
  # mod_aic <- do.call(c, lapply(massmm, '[[', "aic"))
  # mod_bic <- do.call(c, lapply(massmm, '[[', "bic"))
  # mod_caic <- do.call(c, lapply(massmm, '[[', "caic"))
  # AIC_mat <- cbind(mod_aic, mod_bic, mod_caic)
  # colnames(AIC_mat) <- c("AIC", "BIC", "cAIC")
  # ## residuals
  # resids <- NA
  # if(residuals) resids <- suppressMessages(lapply(massmm, '[[', "residuals") %>% dplyr::bind_cols())
  ## random effects
  randEff <-  t(do.call(rbind, lapply(massmm, '[[', "re_vals")))

  ##############################################################################
  ## Step 3 Smooth fixed effects, adjust random effects
  ##############################################################################
  if(!silent) print("Step 2: Smoothing")

  nknots <- min(round(length(bin.pars$bin_midpoints)/2), smoothing.pars$nknots_min)
  betaHat <- t(apply(
    betaTilde, 1,
    function(x) {
      mgcv::gam(
        x ~ s(bin.pars$bin_midpoints, bs = smoothing.pars$splines, k = (nknots + 1)),
        method = smoothing.pars$smooth_method
      )$fitted.values
    }
  ))
  rownames(betaHat) <- rownames(betaTilde)
  colnames(betaHat) <- bin.pars$bin_midpoints
  if(dir.path$save_intermediate_res) saveRDS(betaHat, paste0(res.path, "smoothed-fixed-effect.rds"))
  # load(paste0(res.path, "smoothed-fixed-effect.rds"))

  # Subtract the smoothed estimates from the fitted values / linear predictors
  # Calculate the difference
  for(i in 1:length(massmm)){
    Xbeta.diff <- massmm[[i]][["model_matrix"]] %*% (betaTilde[, i] - betaHat[, i])
    randEff[, i] <- randEff[, i] + Xbeta.diff
  }

  ##############################################################################
  ## Step 4 MFPCA
  ##############################################################################
  if(!silent) print("Step 3: MFPCA")
  id.visit <- strsplit(massmm[[1]][["re_names"]], ":")
  id.visit <- data.frame(
    id = sapply(id.visit, function(x) as.integer(x[1])),
    visit = sapply(id.visit, function(x) as.integer(x[2]))
  )

  if(use_complete_visits){
    # Approach 2: only use visits with complete observations
    visit.no.missing <- apply(randEff, MARGIN = 1, function(x) all(!is.na(x)))
    mfpca.res <- refund::mfpca.face(
      randEff[visit.no.missing, ],
      id = id.visit$id[visit.no.missing],
      visit = id.visit$visit[visit.no.missing],
      argvals = bin.pars$bin_midpoints,
      pve = mfpca.pars$pve,
      twoway = mfpca.pars$twoway
    )
  } else {
    # Approach 1: use all visits
    mfpca.res <- refund::mfpca.face(
      randEff,
      id = id.visit$id,
      visit = id.visit$visit,
      argvals = bin.pars$bin_midpoints,
      pve = mfpca.pars$pve,
      twoway = mfpca.pars$twoway
    )
  }

  mfpca.res.npc_l1 <- ncol(mfpca.res$efunctions$level1)
  mfpca.res.npc_l2 <- ncol(mfpca.res$efunctions$level2)
  if(dir.path$save_intermediate_res) saveRDS(mfpca.res, paste0(res.path, "mfpca.rds"))
  # mfpca.res <- read_rds(paste0(res.path, "mfpca.rds"))
  # mfpca.res <- read_rds("results/sim/mfpca.rds")

  if(length(bin.pars$bin_midpoints) < length(argvals)){
    # Project eigenfunctions back to the full grid
    deigfun.l1.argvals <- reeval_efunctions(
      knots = nknots,
      argvals_bin = bin.pars$bin_midpoints,
      argvals = argvals,
      efunctions = mfpca.res$efunctions$level1[, 1:min(mfpca.res.npc_l1, mfpca.pars$npc_l1)],
      npc = min(mfpca.res.npc_l1, mfpca.pars$npc_l1)
    )
    deigfun.l2.argvals <- reeval_efunctions(
      knots = nknots,
      argvals_bin = bin.pars$bin_midpoints,
      argvals = argvals,
      efunctions = mfpca.res$efunctions$level2[, 1:min(mfpca.res.npc_l2, mfpca.pars$npc_l2)],
      npc = min(mfpca.res.npc_l2, mfpca.pars$npc_l2)
    )
    # Fix beta in sampling
    betaHat.argvals <- reeval_efunctions(
      knots = nknots,
      argvals_bin = bin.pars$bin_midpoints,
      argvals = argvals,
      efunctions = t(betaHat),
      npc = p # p fixed effects
    )
    betaHat.argvals <- t(betaHat.argvals)
    # mu0.argvals <- reeval_efunctions(
    #   knots = nknots,
    #   argvals_bin = bin.pars$bin_midpoints,
    #   argvals = argvals,
    #   efunctions = data.frame(mfpca.res$mu),
    #   npc = 1
    # )[, 1]
    # muj.argvals <- reeval_efunctions(
    #   knots = nknots,
    #   argvals_bin = bin.pars$bin_midpoints,
    #   argvals = argvals,
    #   efunctions = data.frame(mfpca.res$eta),
    #   npc = ncol(mfpca.res$eta)
    # )
  } else {
    deigfun.l1.argvals <- mfpca.res$efunctions$level1[, 1:min(mfpca.res.npc_l1, mfpca.pars$npc_l1)]
    deigfun.l2.argvals <- mfpca.res$efunctions$level2[, 1:min(mfpca.res.npc_l2, mfpca.pars$npc_l2)]
    betaHat.argvals    <- betaHat
    # mu0.argvals        <- mfpca.res$mu
    # muj.argvals        <- mfpca.res$eta
  }

  if(skip_STAN){
    # For large N (sample size), one can divide the full data into subsets for score sampling.
    # Since all subsets condition on the same set of eigenfunctions from MFPCA,
    # scores are still comparable.
    Nsubjects <- length(id.unique)
    subset.idx <- rep(1:n.subsets, length.out = Nsubjects)

    for(sidx in 1:n.subsets){
      if(!(dir.exists(dir.path$intermediate_res_path))){
        dir.create(dir.path$intermediate_res_path, recursive = T)
      }
      ids.subset <- id.unique[subset.idx == sidx]
      saveRDS(ids.subset, file = paste0(dir.path$intermediate_res_path, "ids-part", sidx, ".rds"))

      data.subset    <- filter(data, .data[[var.names$id_name]] %in% ids.subset)
      nvisits.subset <- filter(nvisits, .data[[var.names$id_name]] %in% ids.subset)
      # This may cause issues if certain factor levels are missing in a subset
      # But it would only happen when sample size is small, when skip_STAN is not necessary.
      fixef.design.mx.subset <- model.matrix(formula.fixef, data.subset)
      nonna.idx.subset <- which(!is.na(data.subset[[model_formula[[2]]]]))

      standata <- list(
        Nsubs     = length(ids.subset),
        Nfuns     = sum(nvisits.subset$nvisits),
        S         = length(argvals),
        L         = min(mfpca.res.npc_l1, mfpca.pars$npc_l1),
        M         = min(mfpca.res.npc_l2, mfpca.pars$npc_l2),
        ylen      = length(nonna.idx.subset), # number of non-NA observations
        p         = p, # dimension of fixed effect
        # visit     = as.integer(data[[var.names$visit_name]][nonna.idx.subset]), # pass factor as integer;
        nvisits   = nvisits.subset$nvisits,
        non_na    = nonna.idx.subset,
        y         = data.subset[[model_formula[[2]]]][nonna.idx.subset], # STAN does not accept NA
        X         = fixef.design.mx.subset[nonna.idx.subset, , drop = F],
        betaHat   = betaHat.argvals,
        # mu_0      = mu0.argvals,
        # mu_j      = muj.argvals,
        efuncs_l1 = t(deigfun.l1.argvals),
        efuncs_l2 = t(deigfun.l2.argvals),
        prior_sd_scale = stan.pars$prior_sd_scale,
        grainsize = 1
      )
      write_stan_json(
        standata,
        file = paste0(dir.path$intermediate_res_path, "standata-part", sidx, ".json"),
        always_decimal = FALSE
      )
    }
    print(paste0("Data for STAN are saved to: ", dir.path$intermediate_res_path))
    print("Run ... to sample scores. ")
    frim.res <- list(
      argvals   = argvals,
      bin_idx   = bin.index,
      id_visit  = id.visit,
      betaTilde = betaTilde,
      betaHat   = betaHat,
      # betaSE    = betaSE,
      localglmm = massmm,
      # mu0_mfpca = mu0.argvals,
      # muj_mfpca = muj.argvals,
      # beta_posterior = beta.posterior,
      # beta_post_mean = beta.post.mean,
      evalues = mfpca.res$evalues,
      efunctions = list(
        level1 = deigfun.l1.argvals,
        level2 = deigfun.l2.argvals
      ),
      scores_posterior = NULL,
      scores_postmean = NULL,
      stanfit = NULL,
      ci.level = ci.level
    )
    return(frim.res)
  }

  ##############################################################################
  ## Step 5 Use STAN to sample scores
  ##############################################################################

  if(!silent) print("Step 4: STAN")
  # deigfun.l1.argvals <- mfpca.res$efunctions$level1[, 1:mfpca.pars$npc_l1]
  # deigfun.l2.argvals <- mfpca.res$efunctions$level2[, 1:mfpca.pars$npc_l2]

  fixef.design.mx <- model.matrix(formula.fixef, data)
  nonna.idx <- which(!is.na(data[[model_formula[[2]]]]))
  standata <- list(
    Nsubs     = length(unique(data[[var.names$id_name]])),
    Nfuns     = sum(nvisits$nvisits),
    S         = length(argvals),
    L         = min(mfpca.res.npc_l1, mfpca.pars$npc_l1),
    M         = min(mfpca.res.npc_l2, mfpca.pars$npc_l2),
    ylen      = length(nonna.idx), # number of non-NA observations
    p         = p, # dimension of fixed effect
    # visit     = as.integer(data[[var.names$visit_name]][nonna.idx]), # pass factor as integer;
    nvisits   = nvisits$nvisits,
    non_na    = nonna.idx,
    y         = data[[model_formula[[2]]]][nonna.idx], # STAN does not accept NA
    X         = fixef.design.mx[nonna.idx, , drop = F],
    betaHat   = betaHat.argvals,
    # mu_0      = mu0.argvals,
    # mu_j      = muj.argvals,
    efuncs_l1 = t(deigfun.l1.argvals),
    efuncs_l2 = t(deigfun.l2.argvals),
    prior_sd_scale = stan.pars$prior_sd_scale,
    grainsize = 1
  )
  if(dir.path$save_intermediate_res){
    write_stan_json(
      standata,
      file = paste0(res.path, "standata.json"),
      always_decimal = FALSE
    )
  }

  # standata <- jsonlite::read_json("results/sim/standata.json")
  # set_cmdstan_path(dir.path$cmdstan_installation_path)
  # stanmodel <- system.file(paste0("stanmodel-", family), package = "frim")
  # model_parallel <- cmdstanr::cmdstan_model(exe_file = stanmodel)
  # if(file.exists(stanmodel)){
  # } else {
  #   model_parallel <-  cmdstan_model(
  #     stan_file = paste0(paste0(dir.path$stanfile_path, "stanmodel-", family), ".stan")
  #   )
  # }

  # new way, using {instantiate}
  stanmodel <- instantiate::stan_package_model(
    name = paste0("stanmodel-", family), package = "frim"
  )
  tic <- Sys.time()
  stanfit <- stanmodel$sample(
    # data = paste0(res.path, "standata.json"),
    data = standata,
    iter_warmup = stan.pars$iter_warmup,
    iter_sampling = stan.pars$iter_sampling,
    chains = stan.pars$chains,
    parallel_chains = num_cores,
    refresh = stan.pars$refresh
  )
  toc <- Sys.time()
  toc - tic
  if(dir.path$save_intermediate_res) saveRDS(stanfit, file = paste0(res.path, "stanfit.rds"))
  # stanfit.summary <- stanfit$summary()
  # rhat.val <- stanfit.summary$rhat
  draws.df <- stanfit$draws(format = "df")

  # if(dir.path$save_intermediate_res) saveRDS(draws.df, file = paste0(res.path, "mcmcdraws.rds"))
  # draws.df <- read_rds(paste0(res.path, "mcmcdraws.rds"))
  # draws.df <- read_rds(paste0("results/sim/mcmcdraws.rds"))
  # sigma.sq <- draws.df[, grep("sigma_epi_sq", colnames(draws.df))]
  # plot(sigma.sq$sigma_epi_sq, type = "l")

  # beta.posterior <- array(NA, c(ncol(fixef.design.mx), stan.pars$iter_sampling, length(argvals)))
  # for(i in 1:ncol(fixef.design.mx)){
  #   # print(i)
  #   for(k in 1:length(argvals)){
  #     temp.varname <- paste0("beta[", i, ",", k, "]")
  #     beta.posterior[i, , k] <- draws.df[[temp.varname]]
  #   }
  # }
  # beta.post.mean <- apply(beta.posterior, c(1,3), mean)

  l1score.posterior <-
    array(NA, c(standata$Nsubs, stan.pars$chains * stan.pars$iter_sampling, min(mfpca.res.npc_l1, mfpca.pars$npc_l1)))
  for(i in 1:standata$Nsubs){
    # print(i)
    for(k in 1:min(mfpca.res.npc_l1, mfpca.pars$npc_l1)){
      temp.varname <- paste0("xil1_sc[", i, ",", k, "]")
      l1score.posterior[i, , k] <- draws.df[[temp.varname]]
    }
  }
  l1score.post.mean <- apply(l1score.posterior, c(1,3), mean)

  l2score.posterior <-
    array(NA, c(standata$Nfuns, stan.pars$chains * stan.pars$iter_sampling, min(mfpca.res.npc_l2, mfpca.pars$npc_l2)))
  for(i in 1:standata$Nfuns){
    # print(i)
    for(k in 1:min(mfpca.res.npc_l2, mfpca.pars$npc_l2)){
      temp.varname <- paste0("xil2_sc[", i, ",", k, "]")
      l2score.posterior[i, , k] <- draws.df[[temp.varname]]
    }
  }
  l2score.post.mean <- apply(l2score.posterior, c(1,3), mean)

  return(
    list(
      argvals   = argvals,
      bin_idx   = bin.index,
      id_visit  = id.visit,
      betaTilde = betaTilde,
      betaHat   = betaHat,
      # betaSE    = betaSE,
      localglmm = massmm,
      # mu0_mfpca = mu0.argvals,
      # muj_mfpca = muj.argvals,
      # beta_posterior = beta.posterior,
      # beta_post_mean = beta.post.mean,
      evalues = mfpca.res$evalues,
      efunctions = list(
        level1 = deigfun.l1.argvals,
        level2 = deigfun.l2.argvals
      ),
      scores_posterior = list(
        level1 = l1score.posterior,
        level2 = l2score.posterior
      ),
      scores_postmean = list(
        level1 = l1score.post.mean,
        level2 = l2score.post.mean
      ),
      stanfit = stanfit,
      ci.level = ci.level
    )
  )
}


# frim.res <- list(
#   argvals  = argvals,
#   id_visit = id.visit,
#   betaHat  = betaHat,
#   betaSE   = betaSE,
#   beta_posterior = beta.posterior,
#   beta_post_mean = beta.post.mean,
#   efunctions = list(
#     level1 = deigfun.l1.argvals,
#     level2 = deigfun.l2.argvals
#   ),
#   scores_posterior = list(
#     level1 = l1score.posterior,
#     level2 = l2score.posterior
#   ),
#   scores_postmean = list(
#     level1 = l1score.post.mean,
#     level2 = l2score.post.mean
#   ),
#   ci.level = ci.level
# )

# format scores
# colQuantile <- function(A, prob){
#   apply(A, 2, function(x) quantile(x, probs = prob))
# }
# l1score.posterior.m1 <- draws.df[, grep("xil1_sc", colnames(draws.df))]
# l1score.post.mean.m1 <- cbind(
#     data.frame(id = unique(id.visit$id)),
#     data.frame(matrix(colMeans(l1score.posterior.m1), ncol = mfpca.pars$npc_l1, byrow = F))
# )
# l2score.posterior.m1 <- draws.df[, grep("xil2_sc", colnames(draws.df))]
# l2score.post.mean.m1 <- cbind(
#   id.visit,
#   data.frame(
#     matrix(colMeans(l2score.posterior.m1), ncol = mfpca.pars$npc_l2, byrow = F)
#   )
# )
# l1score.sd <- cbind(
#   data.frame(id = unique(id.visit$id)),
#   data.frame(
#     matrix(
#       apply(l1score.posterior, 2, sd),
#       ncol = mfpca.pars$npc_l1, byrow = F
#     )
#   )
# )
# l2score.sd <- cbind(
#   id.visit,
#   data.frame(
#     matrix(
#       apply(l2score.posterior, 2, sd),
#       ncol = mfpca.pars$npc_l2, byrow = F
#     )
#   )
# )
# scores_sd = list(
#   ci_level = ci.level,
#   level1 = l1score.sd,
#   level2 = l2score.sd
# )
# scores_quantile = list(
#   ci_level = ci.level,
#   level1_lower = l1score.lower,
#   level1_upper = l1score.upper,
#   level2_lower = l2score.lower,
#   level2_upper = l2score.upper
# )


# frim.obj$id_visit <- id.visit
# check traceplots
# plot(l1score.posterior$`xil1_sc[1,1]`, type = "l")

# frim.res$efunctions <- list(
#   level1 = deigfun.l1.argvals,
#   level2 = deigfun.l2.argvals
# )
# frim.res$scores_postmean <- scores_postmean
# frim.res$scores_quantile <- scores_quantile
# str(frim.res)


# frim.res$scores_sd = list(
#     ci_level = ci.level,
#     level1 = l1score.sd,
#     level2 = l2score.sd
#   )

# frim.res$scores_posterior = list(
#   level1 = l1score.posterior,
#   level2 = l2score.posterior
# )
