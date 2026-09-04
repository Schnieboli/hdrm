#' hdrm: Inference for High-Dimensional Repeated Measures
#'
#' The package provides tests for expectation vectors in high-dimensional
#' repeated-measures designs. It implements the one-group procedure of
#' \insertCite{Pauly2015;textual}{hdrm}, the heterogeneous multiple-group
#' procedure of \insertCite{Sattler2018;textual}{hdrm}, and the
#' equal-covariance multiple-group procedure of
#' \insertCite{Sattler2021;textual}{hdrm}.
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
#' @references \insertAllCited
#' @keywords internal
#' @importFrom Rcpp evalCpp
#' @importFrom Rdpack reprompt
#' @useDynLib hdrm, .registration = TRUE
"_PACKAGE"
