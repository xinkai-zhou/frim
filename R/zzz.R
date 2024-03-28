# .onLoad <- function(...) {
#   cmdstan_version <- cmdstanr::cmdstan_version(error_on_NA = FALSE)
#   if (is.null(cmdstan_version)) {
#     stop(
#       paste0(
#         "No CmdStan installation found or CmdStan path has not been set. ",
#         "For installation, see https://mc-stan.org/cmdstanr/. ",
#         "For setting CmdStan path, see ?set_cmdstan_path."
#       ),
#       call. = FALSE
#     )
#   }
# }
