#' Simulate longitudinal functional response data
#'
#' @param family distribution of longitudinal functional data, including "gaussian", "binomial", "poisson".
#' @param I number of subjects.
#' @param Jmean mean number of observations per subject.
#' @param L number of grid points on the functional domain.
#' @param beta_true true fixed effects functions.
#'
#' @importFrom funData eFun
#'
#' @export SimulateLFD
#'
#' @return a data frame containing generated predictors and longitudinal functional outcomes.

SimulateLFD <- function(
    family = "gaussian",
    I  = 100, Jmean = 5, L = 100,
    X  = NULL,
    beta_true  = NULL,
    orthonormal = TRUE,
    npc        = list(level1 = 4,    level2 = 4),
    efunctions = list(level1 = NULL, level2 = NULL),
    evalues    = list(level1 = NULL, level2 = NULL),
    scores     = list(level1 = NULL, level2 = NULL),
    epsilon.sd = 1
){
  times <- seq(0, 1, length.out = L)

  ## generate the number of visits for each subject from the Poisson distribution
  Js <- pmax(rpois(I, Jmean), 1)
  nfun <- sum(Js)


  # generate fixed effects
  if(is.null(X)){
    p <- 2
    X <- cbind(1, matrix(rnorm(nfun * (p-1)), ncol = p - 1))
  } else {
    p <- ncol(X)
  }
  if(is.null(beta_true)){
    beta_true <- t(sapply(seq(1, p, by = 1), function(k) sin(k * pi * times)))
  }
  fixef <- X %*% beta_true


  # generate eigenfunctions and scores
  if(is.null(efunctions$level1)){
    efunctions$level1 <- sapply(
      seq(1, npc$level1, by = 1), function(k) sqrt(2) * sin(k * pi * times)
    )
  }
  if(orthonormal){
    efunctions$level2 <- sapply(
      seq(1, npc$level2, by = 1), function(k) sqrt(2) * cos(k * pi * times)
    )
  } else {
    l2temp <- eFun(times, M = L, type = "Poly")
    efunctions$level2 <- t(l2temp@X)
  }

  if(is.null(evalues$level1)){
    evalues$level1 <- 0.5^seq(from = 0, to = npc$level1 - 1, by = 1)
  }
  if(is.null(evalues$level2)){
    evalues$level2 <- 0.5^seq(from = 0, to = npc$level2 - 1, by = 1)
  }
  if(is.null(scores$level1)){
    scores$level1 <- sapply(
      seq(1, npc$level1, by = 1), function(k) rnorm(I, mean = 0, sd = sqrt(evalues$level1[k]))
    )
  }
  if(is.null(scores$level2)){
    scores$level2 <- sapply(
      seq(1, npc$level2, by = 1), function(k) rnorm(nfun, mean = 0, sd = sqrt(evalues$level2[k]))
    )
  }

  ## generate random effects
  subid <- as.factor(rep(1:I, times = Js))
  ranef.design <- model.matrix( ~ 0 + subid)
  ranef.l1 <- scores$level1 %*% t(efunctions$level1)
  ranef <- ranef.design %*% ranef.l1
  ranef <- ranef + scores$level2 %*% t(efunctions$level2)

  ## adjust for relative importance of random effects
  # ranef <- ((sd(fixef) / sd(ranef)) / SNR_B) * ranef

  ## generate linear predictors
  eta <- fixef + ranef

  ## generate longitudinal functional data
  Y <- matrix(NA, nfun, L)
  if(family == "gaussian"){
    # sd_signal <- sd(eta)
    for(i in 1:nfun){
      # Y[i, ] <- eta[i, ] + rnorm(L, mean = 0, sd = sd_signal/SNR_sigma)
      Y[i, ] <- eta[i, ] + rnorm(L, mean = 0, sd = epsilon.sd)
    }
  } else if(family == "binomial"){
    for(i in 1:nfun){
      Y[i, ] <- rbinom(L, size = 1, prob = plogis(eta[i, ]))
    }
  }else if(family == "poisson"){
    for(i in 1:nfun){
      Y[i, ] <- rpois(n = 1, lambda = exp(eta[i, ]))
    }
  }
  visit <- sequence(Js)

  return(
    list(
      data = data.frame(
        id = subid,
        visit = visit,
        X = X[, -1],
        Y = I(Y)
      ),
      beta_true = beta_true,
      ranef = ranef,
      eta = eta,
      npc = npc,
      efunctions = efunctions,
      evalues = evalues,
      scores = scores
    )
  )
}




# SimulateLFD_Easy <- function(
#     family = "gaussian",
#     I,
#     L,
#     Js, # number of visits per subject
#     eta,
#     epsilon.sd
# ){
#   times <- seq(0, 1, length.out = L)
#   nfun <- sum(Js)
#   subid <- as.factor(rep(1:I, times = Js))
#
#   ## generate longitudinal functional data
#   Y <- matrix(NA, nfun, L)
#   if(family == "gaussian"){
#     # sd_signal <- sd(eta)
#     for(i in 1:nfun){
#       # Y[i, ] <- eta[i, ] + rnorm(L, mean = 0, sd = sd_signal/SNR_sigma)
#       Y[i, ] <- eta[i, ] + rnorm(L, mean = 0, sd = epsilon.sd)
#     }
#   } else if(family == "binomial"){
#     for(i in 1:nfun){
#       Y[i, ] <- rbinom(L, size = 1, prob = plogis(eta[i, ]))
#     }
#   }else if(family == "poisson"){
#     for(i in 1:nfun){
#       Y[i, ] <- rpois(n = 1, lambda = exp(eta[i, ]))
#     }
#   }
#
#   visit <- sequence(Js)
#   return(
#     list(
#       data = data.frame(
#         id = subid,
#         visit = visit,
#         X = X[, -1],
#         Y = I(Y)
#       )
#     )
#   )
# }







