
# hdrm

<!-- badges: start -->

[![DOI](https://img.shields.io/badge/DOI-10.48550%2FarXiv.2512.17478-blue)](https://doi.org/10.48550/arXiv.2512.17478)
[![R-CMD-check](https://github.com/Schnieboli/hdrm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Schnieboli/hdrm/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `hdrm` package provides inference procedures for high-dimensional
repeated-measures data in one-sample and multiple-group designs (Pauly
et al., 2015; Sattler, 2021; Sattler & Pauly, 2018).

## Installation

The current version can be installed with:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  "Schnieboli/hdrm",
  dependencies = TRUE
)
```

Because `hdrm` contains compiled C++ code, Windows users installing from
source need the `Rtools` version corresponding to their version of `R`.

## Citation

When using `hdrm` in a scientific publication, please cite this article
Sattler & Hichert (2025). The citation can also be obtained directly in
`R`:

``` r
citation("hdrm")
```

## Data formats

Both main functions accept either:

- a numeric matrix with subjects in rows and repeated-measurement
  dimensions in columns, or
- a data.frame with columns `value`, `subject` and `dimension`, giving
  the measurements and the subject and dimension IDs.

In the grouped case, a vector specifying which subject belongs to which
group must be passed to `group`.

Missing values in `data` (and `group`) are not allowed and will result
in an error. At least two dimensions are required. The one-group method
requires at least three subjects, while in the multi-group case, each
group must contain at least six complete subjects.

## One-Group Test

A one-group test can be performed using `hdrm_single()`. The
functionaccepts either a numeric matrix with dimensions in rows and
subjects incolumns, or a data.frame containing columns `value`,
`subject` and `dimension`.

The package includes two example data sets: `birthrates` Statistisches
Bundesamt (Destatis) (2024), a wide-format data set containing birth
rates for the 16 German federal states from 1990 to 2023, and `EEG`
Höller et al. (2017), a long-format data set containing EEG measurements
of 160 individuals in 40 dimensions.

### Wide matrix input

``` r
data("birthrates")

## transform 'birthrates' to matrix and transpose
birthrates_matrix <- t(as.matrix(birthrates))

## test whether the time profile is flat
hdrm_single(
  data = birthrates_matrix,
  hypothesis = "flat"
)

## define hypothesis = "flat" via equivalent projection matrix
d <- ncol(birthrates_matrix)
flat_projection <- diag(d) -
  matrix(1 / d, nrow = d, ncol = d)

## test whether the time profile is flat via custom hypothesis matrix
hdrm_single(
  data = birthrates_matrix,
  hypothesis = flat_projection
)
```

### Long data.frame input

The included `EEG` data contain several diagnostic groups. The following
example selects one group for a one-group analysis.

``` r
data("EEG")
names(EEG) # already contains columns value, dimension and subject

## select only one diagnostic group for one-group analysis
EEG_single <- EEG[EEG$group == "SCC+", ]

## test whether the time profile is flat
hdrm_single(
  data = EEG_single,
  hypothesis = "flat"
)

## define hypothesis = "flat" via equivalent projection matrix
d <- nlevels(EEG_single$dimension)
flat_projection <- diag(d) -
  matrix(1 / d, nrow = d, ncol = d)

## test whether the time profile is flat via custom hypothesis matrix
hdrm_single(
  data = EEG_single,
  hypothesis = flat_projection
)
```

## Multiple-Group inference

`hdrm_grouped()` covers two different covariance settings.

| Setting | Main arguments | Interpretation of `B` |
|----|----|----|
| Heterogeneous covariances with exact available trace estimators | `cov.equal = FALSE`, `subsampling = FALSE` | Base budget; the joint third-trace estimator uses `a * B` joint draws |
| Heterogeneous covariances with additional subsampling | `cov.equal = FALSE`, `subsampling = TRUE` | Base budget; group-specific and pairwise estimators use `B` draws per group or pair, while the joint third-trace estimator uses `a * B` joint draws |
| Equal covariance matrices | `cov.equal = TRUE` | Base budget; the pooled third-trace estimator uses `a * B` total draws distributed across groups |

For the equal-covariance method, `subsampling` has no effect. The
`a * B` third-trace draws are allocated across groups approximately
proportionally to the numbers of available six-subject subsets.

A third-trace quantity is estimated by subsampling in every grouped
analysis. Consequently, `f`, `tau`, and the p-value depend on `B` and
the random seed even when `subsampling = FALSE`. With heterogeneous
covariances and `subsampling = TRUE`, the test statistic itself is seed
dependent as well.

Supplying `seed` makes a call reproducible without permanently changing
the previous `R` random-number state.

### Wide matrix input

``` r
data("birthrates")

birthrates_matrix <- t(as.matrix(birthrates))

# group states into east and west (Berlin as east)
group <- factor(c(1, 1, 2, 2, 1, 1, 1, 2, 1, 1, 1, 1, 2, 2, 1, 2), 
                labels = c("west", "east")) 

## test for interaction effect of group and time
hdrm_grouped(
  data = birthrates_matrix,
  hypothesis = "interaction",
  group = group,
  cov.equal = FALSE,
  subsampling = FALSE,
  B = "100*N",
  seed = 3141
)
```

The equal-covariance method is selected explicitly:

``` r
## test for interaction effect of group and time with equal covariance assumption
hdrm_grouped(
  data = birthrates_matrix,
  hypothesis = "interaction",
  group = group,
  cov.equal = TRUE,
  B = "100*N",
  seed = 3141
)
```

### Long data.frame input

If `data` is a data.frame, it must contain columns `value`, `subject`
and `dimension`.

``` r
data("EEG")
names(EEG) # already contains columns 'value', 'subject' and 'dimension'

## test for differences between the diagnostic groups
hdrm_grouped(
  data = EEG,
  hypothesis = "whole",
  group = EEG$group,
  cov.equal = FALSE,
  subsampling = FALSE,
  B = "100*N",
  seed = 3141
)
```

### Hypotheses

For `hdrm_grouped()`, the predefined hypotheses are:

- `"whole"`: no whole-plot or group main effect,
- `"sub"`: no subplot or dimension main effect,
- `"interaction"`: no group-by-dimension interaction,
- `"identical"`: identical expectation vectors across groups,
- `"flat"`: a flat expectation profile within every group.

A custom grouped hypothesis is supplied as a named list containing
projection matrices `TW` and `TS`:

``` r
## custom hypothesis for testing for effect of time
custom_hypothesis <- list(
  TW = matrix(1 / 2, nrow = 2, ncol = 2),
  TS = diag(34) - matrix(1 / 34, nrow = 34, ncol = 34)
)

## testing for time effect (equivalent to hypothesis = "sub")
hdrm_grouped(
  data = birthrates_matrix,
  hypothesis = custom_hypothesis,
  group = group,
  B = "100*N",
  seed = 3141
)
```

## Reported p-values

Upper-tail probabilities are calculated directly. To avoid reporting an
exact numerical zero, p-values smaller than machine precision are
returned as `.Machine$double.eps`. Such a value should be interpreted as
being no larger than the numerical reporting threshold.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
data-entry-spacing="0" data-line-spacing="2">

<div id="ref-EEG_dataset" class="csl-entry">

Höller, Y., Bathke, A. C., Uhl, A., Strobl, N., Lang, A., Bergmann, J.,
Nardone, R., Rossini, F., Zauner, H., Kirschner, M., Jahanbekam, A.,
Trinka, E., & Staffen, W. (2017). Combining SPECT and quantitative EEG
analysis for the automated differential diagnosis of disorders with
amnestic symptoms. *<span class="nocase">Frontiers in Aging
Neuroscience</span>*, *9*, 290.
<https://doi.org/10.3389/fnagi.2017.00290>

</div>

<div id="ref-Pauly2015" class="csl-entry">

Pauly, M., Ellenberger, D., & Brunner, E. (2015). Analysis of
high-dimensional one group repeated measures designs. *Statistics*,
*49*(6), 1243–1261. <https://doi.org/10.1080/02331888.2015.1050022>

</div>

<div id="ref-Sattler2021" class="csl-entry">

Sattler, P. (2021). A comprehensive treatment of quadratic-form-based
inference in repeated measures designs under diverse asymptotics.
*Electronic Journal of Statistics*, *15*(1), 3611–3634.
<https://doi.org/10.1214/21-EJS1865>

</div>

<div id="ref-SattlerHichert2025hdrm" class="csl-entry">

Sattler, P., & Hichert, N. (2025). *Inference for high dimensional
repeated measure designs with the R package hdrm*. arXiv:2512.17478.
<https://doi.org/10.48550/arXiv.2512.17478>

</div>

<div id="ref-Sattler2018" class="csl-entry">

Sattler, P., & Pauly, M. (2018). Inference for high-dimensional
split-plot-designs: A unified approach for small to large numbers of
factor levels. *Electronic Journal of Statistics*, *12*(2), 2743–2805.
<https://doi.org/10.1214/18-EJS1465>

</div>

<div id="ref-birthrates" class="csl-entry">

Statistisches Bundesamt (Destatis). (2024). *Statistischer Bericht –
Geburten 2023, Tabelle 12612-09*. Zusammengefasste Geburtenziffer nach
Bundesländern (Kinder je Frau).
<https://www.statistischebibliothek.de/mir/receive/DEHeft_mods_00160779>

</div>

</div>
