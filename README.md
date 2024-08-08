
# hdrm

<!-- badges: start -->
<!-- badges: end -->

R package for performing tests on high dimensional repeated measure data
for one group \[1\] or multiple groups \[2\].

## Installation

The current version can be installed with:

``` r
## install devtools package
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
# install package
devtools::install_github("Schnieboli/hdrm", dependencies = TRUE)
```

## One group case

A one group test can be performed by using the function hdrm1. Generic
functions exist for the classes matrix, list and data.frame.

``` r
library(hdrm)
## hdrm1 with matrix
N <- 60
d <- 30
# example data with N = 60 subjects in columns and d = 30 factor levels in rows
M <- matrix(rnorm(d*N), d, N)

# test M for H_0 = mu_1 = mu_2 = mu_3 = 0 (no effect in first three factor levels)
H <- matrix(0,d,d)
H[c(1, d+2, 2*d+3)] <- 1
hdrm1(data = M, hypothesis = H)

## hdrm1 with data.frame in longtable format
head(EEG)
# for 'formula' we are interested in value = value, subject = subject and factor = dimension
formula <- value ~ subject + dimension

# test EEG for equal time profile
hdrm1(data = EEG, formula = formula, hypothesis = "equal")
```

## Multiple group case

A test for multiple groups can be performed by the function hdrm_test.
Here also exist generics for objects of class matrix, list and
data.frame.

``` r
library(hdrm)
## hdrm1 with matrix
N <- 60
d <- 30
# example data with N = 60 subjects in columns and d = 30 factor levels in rows...
M <- matrix(rnorm(d*N), d, N)

# ... divided into a = 6 groups of n_i = 10
group <- rep(1:6, each = 10)

# test M for group differences
#hdrm_test(data = M, hypothesis = "time", group = group)

## hdrm1 with data.frame in longtable format
head(EEG)
# for 'formula' we are interested in value = value, subject = subject and factor = dimension
formula <- value ~ subject + dimension

# test EEG for equal time profile
hdrm1(data = EEG, formula = formula, hypothesis = "equal")
```

<div id="refs" class="references csl-bib-body">

<div id="ref-Pauly2015-jo" class="csl-entry">

<span class="csl-left-margin">\[1\]
</span><span class="csl-right-inline">Pauly M, Ellenberger D and Brunner
E 2015 Analysis of high-dimensional one group repeated measures designs
*Statistics (Ber.)* **49** 1243–61</span>

</div>

<div id="ref-Sattler2018-az" class="csl-entry">

<span class="csl-left-margin">\[2\]
</span><span class="csl-right-inline">Sattler P and Pauly M 2018
Inference for high-dimensional split-plot-designs: A unified approach
for small to large numbers of factor levels *Electron. J. Stat.* **12**
2743–805</span>

</div>

</div>
