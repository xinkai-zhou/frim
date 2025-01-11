#' Plot FRIM results.
#'
#' @param frim.obj An object returned by function \code{frim}.
#' @param id.visit A list specifying the \code{id} and \code{visit} to plot.
#'
#' @export plot.frim

plot.frim <- function(
    frim.obj, id.visit = list(id = NULL, visit = NULL),
    xlab = "Time", ylab = "Random Effect"
){
  # plot the trajectory for a specific subject and visit
  if(is.null(id.visit$id)) id.visit$id <- sample(unique(frim.obj$id_visit$id), size = 1)
  if(is.null(id.visit$visit)) id.visit$visit <- sample(
    filter(frim.obj$id_visit, id == id.visit$id)$visit,
    size = 1
  )

  # calculate subject and subject-visit index
  sub.idx <- which(unique(frim.obj$id_visit$id) == id.visit$id)
  sub.visit.idx <- which(
    frim.obj$id_visit$id == id.visit$id & frim.obj$id_visit$visit == id.visit$visit
  )

  re.temp <- frim.obj$efunctions$level1 %*% t(frim.obj$scores_posterior$level1[sub.idx, , ]) +
    frim.obj$efunctions$level2 %*% t(frim.obj$scores_posterior$level2[sub.visit.idx, , ])

  re.postmean <- apply(re.temp, 1, mean)
  re.lower    <- apply(re.temp, 1, function(x) quantile(x, (1 - frim.obj$ci.level)/2))
  re.upper    <- apply(
    re.temp, 1,
    function(x) quantile(x, frim.obj$ci.level + (1 - frim.obj$ci.level)/2)
  )
  plot(
    x = frim.obj$argvals, y = re.postmean, type = "l",
    ylim = c(min(re.lower), max(re.upper)),
    main = paste0("ID = ", id.visit$id, ", VISIT = ", id.visit$visit),
    xlab = xlab, ylab = ylab
  )
  lines(x = frim.obj$argvals, y = re.lower, lty = 2)
  lines(x = frim.obj$argvals, y = re.upper, lty = 2)
}




