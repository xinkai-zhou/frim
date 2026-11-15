.sample_scores_parse_args <- function(args) {
  parsed <- list()
  for (argument in args) {
    key_value <- strsplit(argument, "=", fixed = TRUE)[[1]]
    if (length(key_value) != 2L || key_value[[1]] == "" || key_value[[2]] == "") {
      stop("Each argument must be one nonempty key=value pair", call. = FALSE)
    }
    if (key_value[[1]] %in% names(parsed)) {
      stop("Duplicate argument: ", key_value[[1]], call. = FALSE)
    }
    parsed[[key_value[[1]]]] <- key_value[[2]]
  }
  required <- c("sidx", "standata_path", "family", "iter_warmup",
                "iter_sampling", "chains", "num_cores", "refresh")
  if (!setequal(names(parsed), required)) {
    stop("Required arguments: ", paste(required, collapse = ", "), call. = FALSE)
  }
  parsed
}

.sample_scores_integer_value <- function(parsed, name, minimum = 0L) {
  text <- parsed[[name]]
  if (!is.character(text) || length(text) != 1L ||
      !grepl("^(0|[1-9][0-9]*)$", text)) {
    stop(name, " must use canonical nonnegative-integer notation", call. = FALSE)
  }
  value <- suppressWarnings(as.integer(text))
  if (is.na(value) || value < minimum) {
    stop(name, " must be an integer >= ", minimum, " within R integer range",
         call. = FALSE)
  }
  value
}

.sample_scores_part <- function(
    args,
    model_loader = function(name, package) {
      instantiate::stan_package_model(name = name, package = package)
    },
    clock = Sys.time) {
  parsed <- .sample_scores_parse_args(args)
  integer_value <- function(name, minimum = 0L) {
    .sample_scores_integer_value(parsed, name, minimum)
  }
  sidx <- integer_value("sidx", 1L)
  family <- parsed$family
  if (!family %in% c("gaussian", "binomial", "poisson")) {
    stop("family must be gaussian, binomial or poisson", call. = FALSE)
  }
  standata_path <- parsed$standata_path
  standata_file <- file.path(standata_path, paste0("standata-part", sidx, ".json"))
  if (!file.exists(standata_file)) {
    stop("Missing split Stan data: ", standata_file, call. = FALSE)
  }

  stanmodel <- model_loader(paste0("stanmodel-", family), "frim")
  started <- clock()
  stanfit <- stanmodel$sample(
    data = standata_file,
    iter_warmup = integer_value("iter_warmup"),
    iter_sampling = integer_value("iter_sampling", 1L),
    chains = integer_value("chains", 1L),
    parallel_chains = integer_value("num_cores", 1L),
    refresh = integer_value("refresh")
  )
  finished <- clock()
  runtime <- data.frame(
    runtime = as.numeric(difftime(finished, started, units = "secs")),
    split_idx = sidx,
    step = "step4"
  )
  utils::write.table(
    runtime, file = file.path(standata_path, "runtime.csv"),
    quote = FALSE, sep = ",", append = TRUE, row.names = FALSE,
    col.names = FALSE
  )
  saveRDS(stanfit, file.path(standata_path, paste0("stanfit-part", sidx, ".rds")))
  draws <- stanfit$draws(format = "df")
  saveRDS(draws, file.path(standata_path, paste0("draws-part", sidx, ".rds")))
  invisible(list(fit = stanfit, draws = draws, runtime = runtime,
                 standata_file = standata_file))
}
