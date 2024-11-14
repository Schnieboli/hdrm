This R package is part of my Bachelor's Thesis which is supervised by the authors of the respective Papers (Prof. Dr. Markus Pauly \& Dr. Paavo Sattler). My work has yet to be evaluated!
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

Be aware that the respective
[rtools-version](https://cran.r-project.org/bin/windows/Rtools/) for
your version of R is required for installation.

## One Group Test

A one group test can be performed by using the function
`hdrm_single_longtabe` or `hdrm_single_widetable` depending the format
of your data. Both take a `data.frame` as their first argument. The
package has two data sets included: `birthrates` \[3\] comes in a wide
table format and contains the birthrates of German states from 1990 to
2023. `EEG` \[4\] contains EEG data in 40 dimensions and comes in a long
table format.

``` r
library(hdrm)

### One sample test for data in wide table format with built in data set birthrates
hdrm_single_widetable(data = birthrates,
                      hypothesis = "flat", # test whether time profile is flat
                      )

### One sample test for data in long table format with built in data set EEG
# hypothesis can also be given as a matrix
hypothesis <- matrix(1/40, nrow = 40, ncol = 40) # matrix equivalent to 'flat'
hdrm_single_longtable(data = EEG,
                      hypothesis = hypothesis,
                      value = "value", # can all be given as a character...
                      subject = 4,     # ...or a number
                      dimension = "dimension"
                      )
```

## Multiple Group Test

A test for multiple groups can be performed by using the function
`hdrm_grouped_longtabe` or `hdrm_grouped_widetable` depending the format
of your data. As for the one group test, both take a `data.frame` as
their first argument. The test is performed using bootstraps to estimate
the computationally heaviest estimator. The number of bootstraps can be
controlled via the argument `B`. If your data is large, setting
`bootstrap = TRUE` will call the bootstrap versions for the other
estimators as well.

``` r
library(hdrm)

### Test for multiple groups for data in wide table format with built in data set birthrates

## a vector with the same length as ncol(data) is needed that divides subjects into groups
# divide German states in east and west:
group <- factor(c(1,1,2,2,1,1,1,2,1,1,1,1,2,2,1,2), labels = c("west","east"))
hdrm_grouped_widetable(data = birthrates,
                       hypothesis = "interaction", # test for interaction effect between group and dimension
                       group = group
)

### Test for multiple groups for data in long table format with built in data set EEG
# hypothesis can also be given as a matrix
hypothesis <- list(TW = diag(4) - matrix(1/4, 4, 4),
                   TS =matrix(1/40, 40, 40)
                     ) # list entries equivalent to 'whole'
hdrm_grouped_longtable(data = EEG,
                       hypothesis = hypothesis, # test for time effect in groups
                       group = "group",
                       value = "value", # can all be given as a character...
                       subject = 4,     # ...or a number
                       dimension = "dimension"
)
```

### References

<div id="refs" class="references csl-bib-body">

<div id="ref-Pauly2015" class="csl-entry">

<span class="csl-left-margin">\[1\]
</span><span class="csl-right-inline">Pauly M, Ellenberger D and Brunner
E 2015 [Analysis of high-dimensional one group repeated measures
designs](https://doi.org/10.1080/02331888.2015.1050022) *Statistics*
**49** 1243–61</span>

</div>

<div id="ref-Sattler2018" class="csl-entry">

<span class="csl-left-margin">\[2\]
</span><span class="csl-right-inline">Sattler P and Pauly M 2018
[Inference for high-dimensional split-plot-designs: A unified approach
for small to large numbers of factor
levels](https://doi.org/10.1214/18-ejs1465) *Electronic Journal of
Statistics* **12**</span>

</div>

<div id="ref-birthrates" class="csl-entry">

<span class="csl-left-margin">\[3\]
</span><span class="csl-right-inline">Statistisches Bundesamt (Destatis)
2024 *[Statistischer Bericht - Geburten 2023, Tabelle
12612-09](https://www.destatis.de/DE/Themen/Gesellschaft-Umwelt/Bevoelkerung/Geburten/Publikationen/Downloads-Geburten/statistischer-bericht-geburten-5126104237005.html)*
(Statisitsches Bundesamt)</span>

</div>

<div id="ref-EEG_dataset" class="csl-entry">

<span class="csl-left-margin">\[4\]
</span><span class="csl-right-inline">Höller Y, Bathke A C, Uhl A,
Strobl N, Lang A, Bergmann J, Nardone R, Rossini F, Zauner H, Kirschner
M, Jahanbekam A, Trinka E and Staffen W 2017 [Combining SPECT and
quantitative EEG analysis for the automated differential diagnosis of
disorders with amnestic
symptoms](https://doi.org/10.3389/fnagi.2017.00290) *Front. Aging
Neurosci.* **9** 290</span>

</div>

</div>
