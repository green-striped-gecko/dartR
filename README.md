
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dartR <img src='man/figures/dartRlogo.png' align="right" height="180" />

### An accessible genetic analysis platform for conservation, ecology and agriculture

<!-- badges: start -->

Main repository:
[![](https://www.r-pkg.org/badges/version/dartR?color=blue)](https://cran.r-project.org/package=dartR)
[![CRAN
checks](https://cranchecks.info/badges/summary/dartR)](https://cran.r-project.org/web/checks/check_results_dartR.html)
[![R-CMD-check](https://github.com/green-striped-gecko/dartR/workflows/R-CMD-check/badge.svg)](https://github.com/green-striped-gecko/dartR/actions)
[![R-CMD-check-beta](https://github.com/green-striped-gecko/dartR/actions/workflows/R-CMD-check-beta.yaml/badge.svg?branch=beta)](https://github.com/green-striped-gecko/dartR/actions/workflows/R-CMD-check-beta.yaml)
[![](http://cranlogs.r-pkg.org/badges/last-week/dartR?color=orange)](https://cran.r-project.org/package=dartR)
<!-- badges: end -->

Publication:
[![](https://img.shields.io/badge/doi-10.1111/1755--0998.12745-00cccc.svg)](https://doi.org/10.1111/1755-0998.12745)

Zenodo:
[![DOI](https://zenodo.org/badge/86640709.svg)](https://zenodo.org/badge/latestdoi/86640709)

<!-- badges: start -->

Dev repositories:
[![R-CMD-check-dev](https://github.com/green-striped-gecko/dartR/actions/workflows/R-CMD-check-dev.yaml/badge.svg?branch=dev)](https://github.com/green-striped-gecko/dartR/actions/workflows/R-CMD-check-dev.yaml)
[![R-CMD-check-dev_Arthur](https://github.com/green-striped-gecko/dartR/actions/workflows/R-CMD-check-dev_Arthur.yaml/badge.svg?branch=dev_arthur)](https://github.com/green-striped-gecko/dartR/actions/workflows/R-CMD-check-dev_Arthur.yaml)
[![R-CMD-check-dev_Bernd](https://github.com/green-striped-gecko/dartR/actions/workflows/R-CMD-check-dev_Bernd.yaml/badge.svg?branch=dev_bernd)](https://github.com/green-striped-gecko/dartR/actions/workflows/R-CMD-check-dev_Bernd.yaml)
[![R-CMD-check-dev_Luis](https://github.com/green-striped-gecko/dartR/actions/workflows/R-CMD-check-dev_Luis.yaml/badge.svg?branch=dev_luis)](https://github.com/green-striped-gecko/dartR/actions/workflows/R-CMD-check-dev_Luis.yaml)
[![R-CMD-check-dev_Carlo](https://github.com/green-striped-gecko/dartR/actions/workflows/R-CMD-check-dev_Carlo.yaml/badge.svg?branch=dev_carlo)](https://github.com/green-striped-gecko/dartR/actions/workflows/R-CMD-check-dev_Carlo.yaml)
<!-- badges: end -->

## Overview

dartR is a user-friendly R package that provides in the same platform a
wide range of analyses and pipelines along with strong user support
through quality tutorials and documentation.

## Installation

dartR is on CRAN, so to install it simply type:

``` r
install.packages("dartR")
library(dartR)
```

You can use the function gl.install.vanilla.dartR to install all the
packages that are used by dartR.

``` r
gl.install.vanilla.dartR()
```

You can install the development version of dartR from
[GitHub](https://github.com/) with:

``` r
gl.install.vanilla.dartR(flavour = "dev")
```

Please consult [this installation
tutorial](https://github.com/green-striped-gecko/dartR/wiki/Installation-tutorial)
if you run into any problems during setup.

## Usage

<img src='man/figures/Figure_1.png' align="right" width="630"/>

This is a basic example which shows you how to solve a common problem:

``` r
library(dartR)
## basic example code
```

## Getting started

Are you a R rookie? If you want to learn R and RStudio without any fuss,
have a look at our [R-refresher
tutorial](http://georges.biomatix.org/storage/app/media/uploaded-files/Tutorial_1_dartR_RStudio_Refresher_22-Dec-21.pdf).

Let’s get started by reading your genetic data into dartR; if you have
DArT data, follow [this
tutorial](http://georges.biomatix.org/storage/app/media/uploaded-files/tutorial3adartrdatastructuresandinput22-dec-21-2.pdf);
if not, follow [this
one](http://georges.biomatix.org/storage/app/media/uploaded-files/tutorial3bdartrdatastructuresandinputfromsourcesotherthandartlmagv2-2.pdf).

Checking out our [data manipulation
tutorial](http://georges.biomatix.org/storage/app/media/uploaded-files/tutorial4dartrdatamanipulation22-dec-21-3.pdf)
is the easiest way to get your feet wet with dartR.

[This
tutorial](http://georges.biomatix.org/storage/app/media/uploaded-files/tutorial5dartrbasicfiltering22-dec-21-2.pdf)
will provide some pointers on how to filter your data properly, an
important step that depends on making sound threshold assessments.

## Getting help

Q&A forum in support of users can be accessed
[here](https://groups.google.com/g/dartr?pli=1).

There are two main places to get help with R and R-Studio:

1.  The [RStudio community](https://community.rstudio.com/) is a
    friendly place to ask any question.

2.  [Stack Overflow](https://stackoverflow.com/questions/tagged/r) is a
    great source of answers to common questions.

## Contribute

If you want to help shape the future of dartR, [this
tutorial](http://georges.biomatix.org/storage/app/media/uploaded-files/Tutorial_0_dartR_for_the_Developer_2.0_19-Feb-22.pdf)
is for you.

## Citation

Please acknowledge dartR if you use it in your study. Enter the
following in the R console to retrieve the citation information:

``` r
citation("dartR")
```

Check out our articles:

-   [dartR v2: An accessible genetic analysis platform for conservation,
    ecology and agriculture](https://doi.org/10.1111/2041-210X.13918)

-   [dartR: An R package to facilitate analysis of SNP data generated
    from reduced representation genome
    sequencing](https://doi.org/10.1111/1755-0998.12745)

Have fun working with dartR!

Cheers,

Bernd, Arthur, Luis, Carlo & Olly
