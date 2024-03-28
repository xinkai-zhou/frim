# fli.obj <- fli.res
# id <- 62215
# visit <- 1
# id.visit <- list(id = 62215, visit = 1)

PlotFLI <- function(
    fli.obj, id.visit = list(id = NULL, visit = NULL), 
    xlab = "Time", ylab = "Random Effect"
){
  # plot the trajectory for a specific subject and visit
  if(is.null(id.visit$id)) id.visit$id <- sample(unique(fli.obj$id_visit$id), size = 1)
  if(is.null(id.visit$visit)) id.visit$visit <- sample(
    filter(fli.obj$id_visit, id == id.visit$id)$visit, 
    size = 1
  )
  
  # calculate subject and subject-visit index
  sub.idx <- which(unique(fli.obj$id_visit$id) == id.visit$id)
  sub.visit.idx <- which(
    fli.obj$id_visit$id == id.visit$id & fli.obj$id_visit$visit == id.visit$visit
  )
  
  re.temp <- fli.obj$efunctions$level1 %*% t(fli.obj$scores_posterior$level1[sub.idx, , ]) + 
    fli.obj$efunctions$level2 %*% t(fli.obj$scores_posterior$level2[sub.visit.idx, , ])
  
  re.postmean <- apply(re.temp, 1, mean)
  re.lower    <- apply(re.temp, 1, function(x) quantile(x, (1 - fli.obj$ci.level)/2))
  re.upper    <- apply(
    re.temp, 1, 
    function(x) quantile(x, fli.obj$ci.level + (1 - fli.obj$ci.level)/2)
  )
  plot(
    x = fli.obj$argvals, y = re.postmean, type = "l", 
    ylim = c(min(re.lower), max(re.upper)),
    main = paste0("ID = ", id.visit$id, ", VISIT = ", id.visit$visit),
    xlab = xlab, ylab = ylab
  )
  lines(x = fli.obj$argvals, y = re.lower, lty = 2)
  lines(x = fli.obj$argvals, y = re.upper, lty = 2)
}


# PlotFLI <- function(fli.obj, id.visit = list(id = NULL, visit = NULL)){
#   # plot the trajectory for a specific subject and visit
#   if(is.null(id.visit$id)) id.visit$id <- sample(fli.obj$scores_postmean$level1$id, size = 1)
#   if(is.null(id.visit$visit)) id.visit$visit <- sample(
#     filter(fli.obj$scores_postmean$level2, id == id.visit$id)$visit, 
#     size = 1
#   )
#   # mean curve
#   l1score.mean <- fli.obj$scores_postmean$level1 %>% 
#     filter(id == id.visit$id) %>% 
#     dplyr::select(-id) %>% 
#     as.matrix()
#   l2score.mean <- fli.obj$scores_postmean$level2 %>% 
#     filter(id == id.visit$id & visit == id.visit$visit) %>% 
#     dplyr::select(-c(id, visit)) %>% 
#     as.matrix()
#   re.mean <- fli.obj$efunctions$level1 %*% t(l1score.mean) + 
#     fli.obj$efunctions$level2 %*% t(l2score.mean)
#   
#   # lower 
#   l1score.lower <- fli.obj$scores_quantile$level1_lower %>% 
#     filter(id == id.visit$id) %>% 
#     dplyr::select(-id) %>% 
#     as.matrix()
#   l2score.lower <- fli.obj$scores_quantile$level2_lower %>% 
#     filter(id == id.visit$id & visit == id.visit$visit) %>% 
#     dplyr::select(-c(id, visit)) %>% 
#     as.matrix()
#   re.lower <- fli.obj$efunctions$level1 %*% t(l1score.lower) + 
#     fli.obj$efunctions$level2 %*% t(l2score.lower)
#   
#   # upper
#   l1score.upper <- fli.obj$scores_quantile$level1_upper %>% 
#     filter(id == id.visit$id) %>% 
#     dplyr::select(-id) %>% 
#     as.matrix()
#   l2score.upper <- fli.obj$scores_quantile$level2_upper %>% 
#     filter(id == id.visit$id & visit == id.visit$visit) %>% 
#     dplyr::select(-c(id, visit)) %>% 
#     as.matrix()
#   re.upper <- fli.obj$efunctions$level1 %*% t(l1score.upper) + 
#     fli.obj$efunctions$level2 %*% t(l2score.upper)
#   
#   plot(x = fli.obj$argvals, y = re.mean, type = "l")
#   lines(x = fli.obj$argvals, y = re.lower, lty = 2, col = "blue")
#   lines(x = fli.obj$argvals, y = re.upper, lty = 2, col = "red")
# }






