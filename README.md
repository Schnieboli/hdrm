
# hdrm

<!-- badges: start -->
<!-- badges: end -->

R package for performing tests on high dimensional repeated measure data
for one group (Pauly et al., 2015) or multiple groups (Sattler & Pauly,
2018; Sattler & Rosenbaum, 2025).

## Installation

The current version can be installed by:

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
your version of R is required for installation when using Windows.

## Data Sets

The package has two data sets included: `birthrates` (Statistisches
Bundesamt (Destatis), 2024) is a matrix that contains the birthrates of
German states from 1990 to 2023. `EEG` (Höller et al., 2017) is a data
frame that contains EEG data in 40 dimensions.

``` r
data("birthrates")
?birthrates

data("EEG")
?EEG
```

## One Group Test

A one group test can be performed by using the function `hdrm_single`.
The argument `data` that specifies the data can be a matrix or a vector.
If `data` is a matrix, it is expected that subjects are represented by
columns and factor levels are represented by rows. If `data` is a
vector, `subject` must be a vector of the same length that allocates the
measurements in `data` to the subjects.

``` r
library(hdrm)

##### One sample test with built in data set birthrates
hdrm_single(data = birthrates,
            hypothesis = "flat", # test whether time profile is flat
)

##### One sample test with built in data set EEG
# hypothesis can also be given as a matrix
hypothesis <- matrix(1/40, nrow = 40, ncol = 40) 
# matrix equivalent to hypothesis = 'flat'

hdrm_single(data = EEG$value,
            hypothesis = hypothesis,
            subject = EEG$subject,
)
```

## Multiple Group Test

A test for multiple groups can be performed by using the function
`hdrm_grouped`. The argument `data` that specifies the data can be a
matrix or a vector. If `data` is a matrix, it is expected that subjects
are represented by columns and factor levels are represented by rows.
The argument `group` then allocates the subjects to a group. If `data`
is a vector, The arguments `subject` and `group` must be vectors of the
same length that allocates the measurements in `data` to the subjects
and the subjects to a group.

The test is performed using subsampling to estimate the computationally
heaviest estimator. The number of subsamples can be controlled via the
argument `B`. If your data is large, setting `subsampling = TRUE`
results in the subsampling versions for the other estimators being used
as well.

``` r
library(hdrm)

##### Test for multiple groups with built in data set birthrates

## divide German states in east and west:
group <- factor(c(1,1,2,2,1,1,1,2,1,1,1,1,2,2,1,2), labels = c("west", "east"))
hdrm_grouped(data = birthrates,
             ## test for interaction effect between group and dimension
             hypothesis = "interaction",
             group = group
)

##### Test for multiple groups with built in data set EEG
# hypothesis can also be given as a matrix
hypothesis <- list(TW = diag(4) - matrix(1/4, 4, 4),
                   TS = matrix(1/40, 40, 40)
) # list entries equivalent to hypothesis = 'whole'

hdrm_grouped(data = EEG$value,
             hypothesis = hypothesis, # test for time effect in groups
             group = EEG$group,
             subject = EEG$subject,
)
```


The parameter `cov.equal` can be set to `TRUE` to perform the test for
equal covariances between groups if this assumption can be justified
\[5\].
To reduce computational complexity, an alternative hypothesis matrix
can be employed. \[6\] establishes the existence of a companion
matrix with minimal row dimension that yields identical values for the
considered test statistics.The necessary computational tools to determine
these matrices are implemented in the HypoShrink R package \[4\]. 

## References

<div id="refs" class="references csl-bib-body">

<div id="ref-EEG_dataset" class="csl-entry">

<span class="csl-left-margin">[1] </span>
<span class="csl-right-inline">
Höller Y, Bathke A C, Uhl A, Strobl N, Lang A, Bergmann J, Nardone R, Rossini F, Zauner H, Kirschner M, Jahanbekam A, Trinka E and Staffen W (2017)  
<a href="https://doi.org/10.3389/fnagi.2017.00290">Combining SPECT and quantitative EEG analysis for the automated differential diagnosis of disorders with amnestic symptoms</a>.  
<em>Frontiers in Aging Neuroscience</em>, <strong>9</strong>, 290.
</span>

</div>

<br>

<div id="ref-Pauly2015" class="csl-entry">

<span class="csl-left-margin">[2] </span>
<span class="csl-right-inline">
Pauly M, Ellenberger D and Brunner E (2015)  
<a href="https://doi.org/10.1080/02331888.2015.1050022">Analysis of high-dimensional one group repeated measures designs</a>.  
<em>Statistics</em>, <strong>49</strong>, 1243–1261.
</span>

</div>

<br>

<div id="ref-Sattler2018" class="csl-entry">

<span class="csl-left-margin">[3] </span>
<span class="csl-right-inline">
Sattler P and Pauly M (2018)  
<a href="https://doi.org/10.1214/18-ejs1465">Inference for high-dimensional split-plot-designs: A unified approach for small to large numbers of factor levels</a>.  
<em>Electronic Journal of Statistics</em>, <strong>12</strong>.
</span>

</div>

<br>

<div id="ref-Sattler_HypoShrink2025" class="csl-entry">

<span class="csl-left-margin">[4] </span>
<span class="csl-right-inline">
Sattler P and Rosenbaum M (2025)  
<a href="https://doi.org/10.5281/zenodo.17214498">HypoShrink: Optimize your hypotheses. Keep what matters.</a>  
R package version 1.0.0.  
<a href="https://github.com/PSattlerStat/HypoShrink">https://github.com/PSattlerStat/HypoShrink</a>
</span>

</div>

<br>

<div id="ref-Sattler2021" class="csl-entry">

<span class="csl-left-margin">[5] </span>
<span class="csl-right-inline">
Sattler P (2021)  
<a href="https://doi.org/10.1214/21-EJS1865"><span class="nocase">A comprehensive treatment of quadratic-form-based inference in repeated measures designs under diverse asymptotics</span></a>.  
<em>Electronic Journal of Statistics</em>, <strong>15</strong>, 3611–3634.
</span>

</div>

<br>

<div id="ref-Sattler2025" class="csl-entry">

<span class="csl-left-margin">[6] </span>
<span class="csl-right-inline">
Sattler P and Rosenbaum M (2025)  
<a href="https://doi.org/10.1016/j.spl.2025.110356">Choice of the hypothesis matrix for using the anova-type-statistic</a>.  
<em>Statistics & Probability Letters</em>, 110356.
</span>

</div>

<br>

<div id="ref-birthrates" class="csl-entry">

<span class="csl-left-margin">[7] </span>
<span class="csl-right-inline">
Statistisches Bundesamt (Destatis) (2024)  
<a href="https://destatis.de/DE/Themen/Gesellschaft-Umwelt/Bevoelkerung/Geburten/Publikationen/Downloads-Geburten/statistischer-bericht-geburten-5126104237005.html">Statistischer Bericht - Geburten 2023, Tabelle 12612-09</a>.
</span>

</div>

</div>
</div>

<div id="ref-birthrates" class="csl-entry">

Statistisches Bundesamt (Destatis). (2024). *Statistischer Bericht -
Geburten 2023, Tabelle 12612-09*. Statisitsches Bundesamt.
[destatis.de/DE/Themen/Gesellschaft-Umwelt/Bevoelkerung/Geburten/Publikationen/Downloads-Geburten/statistischer-bericht-geburten-5126104237005.html](https://destatis.de/DE/Themen/Gesellschaft-Umwelt/Bevoelkerung/Geburten/Publikationen/Downloads-Geburten/statistischer-bericht-geburten-5126104237005.html)

</div>

</div>
