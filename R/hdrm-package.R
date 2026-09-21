#' hdrm: Inference for High-Dimensional Repeated Measures
#'
#' The package provides tests for expectation vectors in high-dimensional
#' repeated-measures designs. It implements the one-group procedure described by
#' Pauly et al. (2015), the heterogeneous multiple-group
#' procedure described by Sattler and Pauly (2018), and the
#' equal-covariance multiple-group procedure described by
#' Sattler (2021).
#'
#' The main user-facing functions are:
#' \describe{
#'   \item{[hdrm_single()]}{One-group inference.}
#'   \item{[hdrm_grouped()]}{Multiple-group inference under heterogeneous or
#'   equal covariance matrices.}
#' }
#'
#' Data may be supplied in wide matrix form, with subjects in rows and
#' repeated-measurement dimensions in columns, or as a measurement vector with
#' subject identifiers. See the individual function documentation for the
#' required structure, available hypotheses, and interpretation of the
#' subsampling budget.
#'
#' @references Pauly M, Ellenberger D, Brunner E (2015). “Analysis of high-dimensional one group repeated measures designs.” Statistics, 49(6), 1243–1261. doi: 10.1080/02331888.2015.1050022.
#' @references Sattler P, Pauly M (2018). “Inference for high-dimensional split-plot-designs: A unified approach for small to large numbers of factor levels.” Electronic Journal of Statistics, 12(2), 2743–2805. doi: 10.1214/18-EJS1465.
#' @references Sattler P (2021). “A comprehensive treatment of quadratic-form-based inference in repeated measures designs under diverse asymptotics.” Electronic Journal of Statistics, 15(1), 3611–3634. doi: 10.1214/21-EJS1865.
#' @keywords internal
#' @importFrom Rcpp evalCpp
#' @useDynLib hdrm, .registration = TRUE
"_PACKAGE"
