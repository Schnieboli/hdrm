# hdrm

<!-- badges: start -->

[![DOI](https://img.shields.io/badge/DOI-10.48550%2FarXiv.2512.17478-blue)](https://doi.org/10.48550/arXiv.2512.17478)
<!-- badges: end -->

The `hdrm` package provides inference procedures for expectation vectors in
high-dimensional repeated-measures designs. It implements a one-group method
[1], a multiple-group method allowing heterogeneous covariance matrices [2],
and a multiple-group method under equal covariance matrices [3].

## Installation

The current development version can be installed from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  "Schnieboli/hdrm",
  dependencies = TRUE
)
```

Because `hdrm` contains compiled C++ code, Windows users installing from source
need the Rtools version corresponding to their version of R.

## Citation

When using `hdrm` in a scientific publication, please cite the package paper
[4]. The citation can also be obtained directly in R:

```r
citation("hdrm")
```

## Data formats

Both main functions accept either:

- a numeric matrix with subjects in rows and repeated-measurement dimensions
  in columns, or
- a numeric measurement vector together with subject identifiers.

For vector input, measurements must occur in the same dimensional order for
every subject. In `hdrm_grouped()`, subject labels only need to distinguish
subjects within a group; the same labels may be reused in different groups.

Incomplete subjects are removed as complete blocks and reported through a
warning. At least two dimensions are required. The one-group method requires
at least three complete subjects, while every group in `hdrm_grouped()` must
contain at least six complete subjects.

## One-group inference

`hdrm_single()` implements the one-group procedure. The predefined hypothesis
`"flat"` tests whether the expectation profile is constant over dimensions.
A custom hypothesis can be supplied as a symmetric, idempotent projection
matrix with positive rank.

### Wide matrix input

```r
library(hdrm)

data("birthrates", package = "hdrm")

# The bundled data set stores years in rows and federal states in columns.
# hdrm uses the usual R wide-data convention: subjects in rows.
birthrates_matrix <- t(as.matrix(birthrates))

birthrates_result <- hdrm_single(
  data = birthrates_matrix,
  hypothesis = "flat"
)

birthrates_result
```

### Measurement-vector input

The included `EEG` data contain several diagnostic groups. The following
example selects one group for a one-group analysis.

```r
data("EEG", package = "hdrm")

eeg_single <- droplevels(
  EEG[EEG$group == levels(EEG$group)[1L], ]
)

eeg_single <- eeg_single[
  order(eeg_single$subject, eeg_single$dimension),
]

d <- nlevels(eeg_single$dimension)
flat_projection <- diag(d) -
  matrix(1 / d, nrow = d, ncol = d)

eeg_single_result <- hdrm_single(
  data = eeg_single$value,
  hypothesis = flat_projection,
  subject = eeg_single$subject
)

eeg_single_result
```

## Multiple-group inference

`hdrm_grouped()` covers two different covariance settings.

| Setting | Main arguments | Interpretation of `B` |
|---|---|---|
| Heterogeneous covariances with exact available trace estimators | `cov.equal = FALSE`, `subsampling = FALSE` | Base budget; the joint third-trace estimator uses `a * B` joint draws |
| Heterogeneous covariances with additional subsampling | `cov.equal = FALSE`, `subsampling = TRUE` | Base budget; group-specific and pairwise estimators use `B` draws per group or pair, while the joint third-trace estimator uses `a * B` joint draws |
| Equal covariance matrices | `cov.equal = TRUE` | Base budget; the pooled third-trace estimator uses `a * B` total draws distributed across groups |

For the equal-covariance method, `subsampling` has no effect. The `a * B`
third-trace draws are allocated across groups approximately proportionally to
the numbers of available six-subject subsets.

A third-trace quantity is estimated by subsampling in every grouped analysis.
Consequently, `f`, `tau`, and the p-value depend on `B` and the random seed even
when `subsampling = FALSE`. With heterogeneous covariances and
`subsampling = TRUE`, the test statistic itself is seed dependent as well.

Supplying `seed` makes a call reproducible without permanently changing the
previous R random-number state.

### Wide matrix input

```r
data("birthrates", package = "hdrm")

birthrates_matrix <- t(as.matrix(birthrates))

group <- factor(
  c(1, 1, 2, 2, 1, 1, 1, 2, 1, 1, 1, 1, 2, 2, 1, 2),
  labels = c("west", "east")
)

heterogeneous_result <- hdrm_grouped(
  data = birthrates_matrix,
  hypothesis = "interaction",
  group = group,
  cov.equal = FALSE,
  subsampling = FALSE,
  B = "100*N",
  seed = 3141
)

heterogeneous_result
```

The equal-covariance method is selected explicitly:

```r
equal_covariance_result <- hdrm_grouped(
  data = birthrates_matrix,
  hypothesis = "whole",
  group = group,
  cov.equal = TRUE,
  B = "100*N",
  seed = 3141
)

equal_covariance_result
```

### Measurement-vector input

Subject labels may be reused in different groups. There is no need to construct
globally unique identifiers.

```r
data("EEG", package = "hdrm")

eeg_grouped <- EEG[
  order(EEG$group, EEG$subject, EEG$dimension),
]

eeg_grouped_result <- hdrm_grouped(
  data = eeg_grouped$value,
  hypothesis = "whole",
  group = eeg_grouped$group,
  subject = eeg_grouped$subject,
  cov.equal = FALSE,
  subsampling = FALSE,
  B = "100*N",
  seed = 3141
)

eeg_grouped_result
```

## Hypotheses

For `hdrm_grouped()`, the predefined hypotheses are:

- `"whole"`: no whole-plot or group main effect,
- `"sub"`: no subplot or dimension main effect,
- `"interaction"`: no group-by-dimension interaction,
- `"identical"`: identical expectation vectors across groups,
- `"flat"`: a flat expectation profile within every group.

A custom grouped hypothesis is supplied as a named list containing projection
matrices `TW` and `TS`:

```r
custom_hypothesis <- list(
  TW = matrix(1 / 2, nrow = 2, ncol = 2),
  TS = diag(34) - matrix(1 / 34, nrow = 34, ncol = 34)
)

custom_result <- hdrm_grouped(
  data = birthrates_matrix,
  hypothesis = custom_hypothesis,
  group = group,
  B = "100*N",
  seed = 3141
)
```

## Reported p-values

Upper-tail probabilities are calculated directly. To avoid reporting an exact
numerical zero, p-values smaller than machine precision are returned as
`.Machine$double.eps`. Such a value should be interpreted as being no larger
than the numerical reporting threshold.

## References

1. Pauly M, Ellenberger D and Brunner E (2015). “Analysis of
   high-dimensional one group repeated measures designs.” *Statistics* 49,
   1243–1261. <https://doi.org/10.1080/02331888.2015.1050022>
2. Sattler P and Pauly M (2018). “Inference for high-dimensional
   split-plot-designs: A unified approach for small to large numbers of factor
   levels.” *Electronic Journal of Statistics* 12(2), 2743–2805.
   <https://doi.org/10.1214/18-EJS1465>
3. Sattler P (2021). “A comprehensive treatment of quadratic-form-based
   inference in repeated measures designs under diverse asymptotics.”
   *Electronic Journal of Statistics* 15, 3611–3634.
   <https://doi.org/10.1214/21-EJS1865>
4. Sattler P and Hichert N (2025). “Inference for high dimensional repeated
   measure designs with the R package hdrm.” *arXiv:2512.17478*.
   <https://doi.org/10.48550/arXiv.2512.17478>
5. Statistisches Bundesamt (Destatis) (2024). “Statistischer Bericht –
   Geburten 2023, Tabelle 12612-09.”
6. Höller Y et al. (2017). “Combining SPECT and quantitative EEG analysis for
   the automated differential diagnosis of disorders with amnestic symptoms.”
   *Frontiers in Aging Neuroscience* 9, 290.
   <https://doi.org/10.3389/fnagi.2017.00290>
