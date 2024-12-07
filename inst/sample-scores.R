

library(tidyverse)
library(magrittr)
library(cmdstanr)

args <- commandArgs(trailingOnly = TRUE)

ParseArgs <- function(args) {
  args_list <- list()
  for (arg in args) {
    key_value <- strsplit(arg, "=")[[1]]
    if (length(key_value) == 2) {
      key <- key_value[1]
      value <- key_value[2]
      args_list[[key]] <- value
    }
  }
  return(args_list)
}
args.parsed <- ParseArgs(args)
print(str(args.parsed))

# subset index
sidx <- as.integer(args.parsed$sidx)
# path to standata
standata.path <- args.parsed$standata_path
family <- args.parsed$family

standatafile <- paste0(standata.path, "standata-part", sidx, ".json")
stanmodel <- instantiate::stan_package_model(
  name = paste0("stanmodel-", family), package = "fli"
)
tic <- Sys.time()
stanfit <- stanmodel$sample(
  # data = paste0(res.path, "standata.json"),
  data = standatafile,
  iter_warmup =  as.integer(args.parsed$iter_warmup),
  iter_sampling =  as.integer(args.parsed$iter_sampling),
  chains =  as.integer(args.parsed$chains),
  parallel_chains =  as.integer(args.parsed$num_cores),
  refresh =  as.integer(args.parsed$refresh)
)
toc <- Sys.time()
toc - tic
saveRDS(stanfit, file = paste0(standata.path, "stanfit-part", sidx, ".rds"))

draws.df <- stanfit$draws(format = "df")
saveRDS(draws.df, file = paste0(standata.path, "draws-part", sidx, ".rds"))

# # stanfit.summary <- stanfit$summary()
# # rhat.val <- stanfit.summary$rhat
# draws.df <- stanfit$draws(format = "df")
#
#
# l1score.posterior <-
#   array(NA, c(dstandata$Nsubs, args.parsed$chains * args.parsed$iter_sampling, dstandata$L))
# for(i in 1:standata$Nsubs){
#   for(k in 1:dstandata$L){
#     temp.varname <- paste0("xil1_sc[", i, ",", k, "]")
#     l1score.posterior[i, , k] <- draws.df[[temp.varname]]
#   }
# }
# l1score.post.mean <- apply(l1score.posterior, c(1,3), mean)
#
# l2score.posterior <-
#   array(NA, c(standata$Nfuns, args.parsed$chains * args.parsed$iter_sampling, dstandata$M))
# for(i in 1:standata$Nfuns){
#   for(k in 1:dstandata$M){
#     temp.varname <- paste0("xil2_sc[", i, ",", k, "]")
#     l2score.posterior[i, , k] <- draws.df[[temp.varname]]
#   }
# }
# l2score.post.mean <- apply(l2score.posterior, c(1,3), mean)


