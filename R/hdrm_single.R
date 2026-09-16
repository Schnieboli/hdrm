#' One-group inference for high-dimensional repeated measures
#'
#' @description Implements the one-group test described by Pauly et al. (2015).
#'
#' @param data A matrix or a data.frame. For matrix input, subjects are
#'   represented by rows and repeated-measurement dimensions by columns. For
#'   data.frame input, data must contain columns `value`, `subject` and
#'   `dimension`, giving the observations and the subject and dimension IDs.
#' @param hypothesis Either `"flat"` or a finite numeric projection matrix whose
#'   dimensions equal the repeated-measurement dimension.
#' @param AM A single logical value, specifying whether the compact
#'   representation of the hypothesis matrix described by
#'   Sattler and Rosenbaum (2025) is used. It may reduce the number of
#'   rows used in the calculations without changing the resulting test. The
#'   default is `TRUE`.
#'
#' @details #' The predefined value `"flat"` tests \deqn{\bm P_d\bm\mu=\bm 0.}
#' Alternatively, `hypothesis` may be a finite numeric \eqn{d\times d}
#' projection matrix. It must be symmetric, idempotent, and have positive rank.
#' Small numerical deviations within the implemented tolerance are accepted.
#' 
#' At least two dimensions and three subjects are required. Missing values in 
#' the input data will result in an error.
#' 
#' Upper-tail probabilities are computed directly. Reported p-values are bounded
#' below by `.Machine$double.eps`; a returned value at this boundary should be
#' interpreted as no larger than the numerical reporting threshold.
#' 
#'
#' @returns A named list of class `"hdrm_single"` with components:
#' \describe{
#'   \item{data}{A matrix with the data used. if data was a data.frame, this is
#'   the transformed output.}
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
#' @references Pauly M, Ellenberger D, Brunner E (2015). “Analysis of high-dimensional one group repeated measures designs.” Statistics, 49(6), 1243–1261. doi: 10.1080/02331888.2015.1050022.
#' @references Sattler P, Rosenbaum M (2025). “Choice of the hypothesis matrix for using the Anova-type-statistic.” Statistics & Probability Letters, 219, 110356. doi: 10.1016/j.spl.2025.110356.
#' @export
hdrm_single <- function(data,
                        hypothesis = "flat",
                        AM = TRUE) {
  UseMethod("hdrm_single")
}
