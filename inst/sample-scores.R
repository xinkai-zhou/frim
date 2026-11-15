# Purpose: execute one split-data FRIM global score sampler.
# Inputs/outputs are validated and implemented by frim:::.sample_scores_part().
# Dependencies: installed frim package, instantiate, and CmdStanR/CmdStan.

frim:::.sample_scores_part(commandArgs(trailingOnly = TRUE))

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


