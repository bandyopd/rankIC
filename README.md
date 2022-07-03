# `rankIC`: R package for fitting the partially interval-censored rank regressions

## Overview

`rankIC` is the R package to fit the linear rank regressions when the data are partially interval-censored.

## Installation
```r
devtools::install_github(repo='taehwa015/rankIC')
```

## Status of development

The code provided here is only being provided for research purposes.

## Documentation

Vignette is available at [here](http://htmlpreview.github.io/?https://github.com/taehwa015/rankIC/blob/master/vignettes/rankIC.html).

## Usage

The `rankIC` package provides a rank-based estimating procedures for the partially interval-censored data.
See Choi et al. (2022+) for more detailed description of the method.


Below example is the phase 3 metastatic colorectal cancer clinical trial study.
Since the event times are possibly correlated within same patient, 
we regard patient id as the cluster. 
To adjust informative cluster sizes, we further consider weight function,
motivated by Wang and Zhao (2008).
Larger cluster will be underweighted by letting `alpha = 1`, 
while cluster structure will be ignored `alpha = 0`.
Note that by letting `id = NULL`, we fit the univariate AFT model.
```r
library(PICBayes)
data(mCRC)
dt0 = as.data.frame(mCRC)
d = with(dt0,
         data.frame(U = ifelse(is.na(L), 0, L),
                    V = ifelse(is.na(R), Inf, R),
                    Delta = 1-IC,
                    x1 = TRT_C,
                    x2 = KRAS_C,
                    id = SITE))
U = d$U; V = d$V; X = cbind(d$x1, d$x2); Delta = d$Delta; id = d$id
aft_rank(U = U, V = V, X = X, Delta = Delta, id = id, 
         alpha = 1, type = "gehan", nboot = 10)
aft_rank(U = U, V = V, X = X, Delta = Delta, id = id, 
         alpha = 1, type = "logrank", nboot = 10)
```

## Reference

* Choi, T., Choi, S. and Bandyopadhyay, D. (2022+). 
Rank estimation for the accelerated failure time model with partially interval-censored data, (Under Review)

* Wang, Y. G. and Zhao, Y. (2008). 
Weighted rank regression for clustered data analysis. 
Biometrics, 64(1), 39--45.

* Pan, C. (2021).
PICBayes: Bayesian Models for Partly Interval-Censored Data.
R package.
https://CRAN.R-project.org/package=PICBayes.

