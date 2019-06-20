# RVFamSq

Rare Variant-Family-Based Score Test for Quantitative Traits  

## Description

The RVFamSq package provides an efficient approach to examine the association between the region-based rare variants and the quantitative traits in family-based data. The RVFamsSq can be broadly applied to diverse pedigrees with members missing sequence data. In addition, qualitative and quantitative covariates, e.g., age, sex, and body mass index, can be flexibly included. The speed of the package is significantly optimized to analyze large-scale data sets.

## Quick Start

1. Follow the setup instructions below.

2. See the [Quick Start Tutorial]() for an
introduction to `RVFamSq`.

## Setup

To install and load the package,

```R
library("devtools")
install_github('statgenetics/RVFamSq')
library("RVFamSq")
```

## Developer notes

+ To build documentation for the package,

   ```R
   setwd("/path/to/package/root")
   devtools::document()
   ```
  Please **do not** manually edit any `.Rd` files under `man` folder!

+ These are the R commands to build the website (make sure you are
connected to Internet while running these commands):

   ```R
   library(pkgdown)
   build_site(mathjax = FALSE)
   ```
+ To install locally from source code via `devtools`, 

   ```R
   setwd("/path/to/package/root")
   devtools::install()
   ```
