#' One-group inference for high-dimensional repeated measures
#'
#' @description
#' Implements the one-group test of
#' \insertCite{Pauly2015;textual}{hdrm}.
#'
#' @param data A numeric vector or matrix. For matrix input, subjects are
#' represented by rows and repeated-measurement dimensions by columns. For
#' vector input, `subject` must identify the subject associated with every
#' measurement.
#' @param hypothesis Either `"flat"` or a finite numeric projection matrix
#' whose dimensions equal the repeated-measurement dimension.
#' @param subject An optional one-dimensional atomic vector identifying
#' subjects. It is required for vector input and ignored with a warning for
#' matrix input.
#' @param AM A single logical value, or alternatively `0` or `1`, specifying
#' whether the compact representation of the hypothesis matrix described by
#' \insertCite{Sattler2025;textual}{hdrm} is used. It may reduce the number of
#' rows used in the calculations without changing the resulting test. The
#' default is `TRUE`.
#'
#' @details
#' For vector input, a missing value in `data` causes the entire affected
#' subject to be removed. The `subject` vector must not contain missing values.
#' For matrix input, every row containing at least one missing value is
#' removed. A warning is issued whenever incomplete subjects are dropped.
#'
#' For vector input, repeated measurements must occur in the same order for
#' every subject. No separate variable identifying the repeated-measurement
#' dimension is supplied, so the within-subject order in `data` determines the
#' component order in the processed data matrix. At least two dimensions and
#' three complete subjects are required.
#'
#' The predefined value `"flat"` tests
#' \deqn{\bm P_d\bm\mu=\bm 0.}
#' Alternatively, `hypothesis` may be a finite numeric
#' \eqn{d\times d} projection matrix. It must be symmetric, idempotent, and
#' have positive rank. Small numerical deviations within the implemented
#' tolerance are accepted.
#'
#' The third trace is estimated by
#' \deqn{
#' B_3=\binom{N}{3}^{-1}
#' \sum_{i<j<k}A_{ij}A_{jk}A_{ki}.
#' }
#' The compiled helper computes the unnormalized sum, and division by
#' \eqn{\binom{N}{3}} is performed in the R wrapper.
#'
#' Upper-tail probabilities are computed directly. Reported p-values are
#' bounded below by `.Machine$double.eps`; a returned value at this boundary
#' should be interpreted as no larger than the numerical reporting threshold.
#'
#' @returns A named list of class `"hdrm_single"` with components:
#' \describe{
#'   \item{data}{The processed data matrix with subjects in rows and
#'   repeated-measurement dimensions in columns.}
#'   \item{f}{The estimated Pearson degrees of freedom.}
#'   \item{statistic}{The standardized test statistic \eqn{W}.}
#'   \item{tau}{The estimated convergence parameter \eqn{\tau=1/f}.}
#'   \item{H}{The hypothesis-matrix representation used in the calculation.}
#'   \item{hypothesis}{`"flat"` or `"custom"`.}
#'   \item{p.value}{The upper-tail p-value, bounded below by
#'   `.Machine$double.eps`.}
#'   \item{dim}{A named numeric vector containing the repeated-measurement
#'   dimension `d` and the number of analyzed subjects `N`.}
#'   \item{removed.cases}{The number of incomplete subjects removed before the
#'   analysis.}
#' }
#'
#' @example man/examples/examples_hdrm_single.R
#'
#' @references \insertAllCited
#'
#' @export
hdrm_single <- function(
    data,
    hypothesis = "flat",
    AM = TRUE
) {
  UseMethod("hdrm_single")
}
