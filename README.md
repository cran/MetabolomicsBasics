
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MetabolomicsBasics

<!-- badges: start -->

[![R-CMD-check](https://github.com/janlisec/MetabolomicsBasics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/janlisec/MetabolomicsBasics/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/MetabolomicsBasics)](https://CRAN.R-project.org/package=MetabolomicsBasics)
<!-- badges: end -->

The goal of MetabolomicsBasics is to provide a set of functions to
investigate raw data (a matrix of intensity values) from (metabol)omics
experiments, i.e.  following peak picking and signal deconvolution.
Functions can be used to *i.e.*:

- normalize data
- detect biomarkers
- perform sample classification

A detailed description of best practice usage may be found in the
publication
<https://link.springer.com/protocol/10.1007/978-1-4939-7819-9_20>.

## Installation

You can install the development version of MetabolomicsBasics from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("janlisec/MetabolomicsBasics")
```

## Examples

A typical use case would be to compute a Principal Component Analysis:

``` r
raw <- MetabolomicsBasics::raw
sam <- MetabolomicsBasics::sam
MetabolomicsBasics::RestrictedPCA(dat = raw, sam = sam, group.col = "Group", legend.x = "bottomleft", medsd = TRUE, fmod = "Group")
```

<img src="man/figures/README-example1-1.png" width="100%" /> More
elaborate plots, like the polar coordinate visualization of heterosis
pattern are possible:

``` r
x <- t(raw)
colnames(x) <- sam$GT
MetabolomicsBasics::PolarCoordHeterPlot(x=x, gt=c("B73","B73xMo17","Mo17"), plot_lab="graph", col=1:10, thr=0.5, rev_log=exp(1))
#> Parameter 'col' should be a color vector of length nrow(x)
```

<img src="man/figures/README-example2-1.png" width="100%" />
