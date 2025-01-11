
<!-- README.md is generated from README.Rmd. Please edit that file -->

# frim

<!-- badges: start -->
<!-- badges: end -->

## Installation

- Set environment variables to point to the `cmdstan` installation path.
  Specifically, open the R environment file using

<!-- -->

    usethis::edit_r_environ()

and add the following

    CMDSTAN = "PATH_TO_CMDSTAN/cmdstan/"
    CMDSTAN_INSTALL = "implicit"

Then install the `frim` package using

    remotes::install_github("xinkai-zhou/frim")

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
