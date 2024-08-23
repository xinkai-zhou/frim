#' Create bins for regularly spaced functional data.
#'
#'
#' @param argvals A vector containing locations of observations on the
#' functional domain.
#' @param bin_midpoints Centers of the bins.
#' @param bin_len Length of the bin. In the same unit as argvals.
#' @param cyclic
#'
#' @export CreateBins
#'
CreateBins <- function(argvals, bin_midpoints, bin_len, cyclic = F){
  S <- length(argvals)
  argvals.idx <- 1:S

  bin_midpoints <- sort(bin_midpoints)
  nbins <- length(bin_midpoints)
  bmp.idx <- rep(NA, nbins) # bin midpoint index
  for(j in 1:nbins){
    bmp.idx[j] <- which(argvals == bin_midpoints[j])
  }

  bins <- list()
  if(cyclic){
    for(j in 1:nbins){
      left.end  <- bmp.idx[j] - floor(bin_len / 2)
      right.end <- bmp.idx[j] + floor(bin_len / 2)
      if(left.end < argvals.idx[1]){
        if(left.end %% S == 0){
          bins[[j]] <- argvals[c(S, argvals.idx[1]:right.end)]
        } else {
          bins[[j]] <- argvals[c((left.end %% S) : S, argvals.idx[1]:right.end)]
        }
      } else if(right.end > S){
        bins[[j]] <- argvals[c(left.end : S, 1:(right.end %% S))]
      } else {
        bins[[j]] <- argvals[left.end : right.end]
      }
      if(j == 1 & !(argvals[1] %in% bins[[j]])){
        warning("The first bin does not cover the left-most point on the domain.")
      }
      if(j == nbins & !(argvals[S] %in% bins[[j]])){
        warning("The last bin does not cover the right-most point on the domain.")
      }
    }
  } else {
    for(j in 1:nbins){
      left.end  <- bmp.idx[j] - floor(bin_len / 2)
      right.end <- bmp.idx[j] + floor(bin_len / 2)
      if(left.end < argvals.idx[1]){
        bins[[j]] <- argvals[argvals.idx[1] : right.end]
      } else if(left.end + bin_len - 1 > S){
        bins[[j]] <- argvals[left.end : S]
      } else {
        bins[[j]] <- argvals[left.end : (left.end + bin_len - 1)]
      }
      if(j == 1 & !(argvals[1] %in% bins[[j]])){
        warning("The first bin does not cover the left-most point on the domain.")
      }
      if(j == nbins & !(argvals[S] %in% bins[[j]])){
        warning("The last bin does not cover the right-most point on the domain.")
      }
    }
  }
  names(bins) <- bin_midpoints
  return(bins)
}

# # Example
# argvals <- c(1, seq(10, 1440, by = 10))
# bin_len <- 6
# bin_midpoints <- argvals
# CreateBins(argvals = argvals, bin_midpoints = bin_midpoints, bin_len = bin_len, cyclic = T)
# CreateBins(argvals = argvals, bin_midpoints = bin_midpoints, bin_len = bin_len, cyclic = F)
#
#
# argvals <- 1:145
# bin_len <- 5
# bin_midpoints <- argvals
# CreateBins(argvals = argvals, bin_midpoints = bin_midpoints, bin_len = bin_len, cyclic = T)
# CreateBins(argvals = argvals, bin_midpoints = bin_midpoints, bin_len = bin_len, cyclic = F)



