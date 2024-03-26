
# KEPTED

<!-- badges: start -->
<!-- badges: end -->

The goal of KEPTED is to provide an implementation of a kernel-embedding of probability test for elliptical distribution, which has been derived by [Tang and Li (2024)](https://arxiv.org/abs/2306.10594). This is an asymptotic test for elliptical distribution under general alternatives, and the location and shape parameters are assumed to be unknown. Some side-products are posted, including the transformation between rectangular and polar coordinates and two product-type kernel functions.


## Installation

You can install the development version of KEPTED from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tyy20/KEPTED")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(KEPTED)

n=200
d=3

## test under a null distribution
X=matrix(rnorm(d*n),nrow=n,ncol=d)
EllKEPT(X,kerU="Gaussian",kerTheta="Gaussian")
EllKEPT(X,kerU="PIQ",kerTheta="PIQ")

## test under an alternative distribution
X=matrix(rchisq(d*n,2),nrow=n,ncol=d)
EllKEPT(X,kerU="Gaussian",kerTheta="Gaussian")
EllKEPT(X,kerU="PIQ",kerTheta="PIQ")
```

## Reference

Tang, Y. and Li, B. (2024), “A nonparametric test for elliptical distribution based on kernel embedding of probabilities” (https://arxiv.org/abs/2306.10594)
