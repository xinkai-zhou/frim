#' Create bins for regularly spaced functional data.
CreateBins <- function(argvals, bin_midpoints, bin_len, cyclic = F){
  S <- length(argvals)
  bin_midpoints <- sort(bin_midpoints)
  nbins <- length(bin_midpoints)
  bin.index <- list()
  if(cyclic){
    for(j in 1:nbins){
      left.end  <- bin_midpoints[j] - floor(bin_len / 2)
      right.end <- bin_midpoints[j] + floor(bin_len / 2)
      if(left.end < argvals[1]){
        if(left.end %% S == 0){
          bin.index[[j]] <- c(S, argvals[1]:right.end)
        } else {
          bin.index[[j]] <- c((left.end %% S) : S, argvals[1]:right.end)
        }
      } else if(right.end > S){
        bin.index[[j]] <- c(left.end : S, 1:(right.end %% S))
      } else {
        bin.index[[j]] <- left.end : right.end
      }
      if(j == 1 & !(argvals[1] %in% bin.index[[j]])){
        warning("The first bin does not cover the left-most point on the domain.")
      }
      if(j == nbins & !(argvals[S] %in% bin.index[[j]])){
        warning("The last bin does not cover the right-most point on the domain.")
      }
    }
  } else {
    for(j in 1:nbins){
      left.end  <- bin_midpoints[j] - floor(bin_len / 2)
      right.end <- bin_midpoints[j] + floor(bin_len / 2)
      if(left.end < argvals[1]){
        bin.index[[j]] <- argvals[1] : right.end
      } else if(left.end + bin_len - 1 > S){
        bin.index[[j]] <- left.end : S
      } else {
        bin.index[[j]] <- left.end : (left.end + bin_len - 1)
      }
      if(j == 1 & !(argvals[1] %in% bin.index[[j]])){
        warning("The first bin does not cover the left-most point on the domain.")
      }
      if(j == nbins & !(argvals[S] %in% bin.index[[j]])){
        warning("The last bin does not cover the right-most point on the domain.")
      }
    }
  }
  names(bin.index) <- bin_midpoints
  return(bin.index)
}

# # Example 
# argvals <- 1:20
# bin_len <- 7
# bin_midpoints <- c(1, 5, 10, 15)
# CreateBins(argvals = argvals, bin_midpoints = bin_midpoints, bin_len = bin_len, cyclic = T)
# CreateBins(argvals = argvals, bin_midpoints = bin_midpoints, bin_len = bin_len, cyclic = F)





