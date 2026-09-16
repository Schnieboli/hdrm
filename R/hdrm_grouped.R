#' Multiple-group inference for high-dimensional repeated measures
#'
#' @description Implements the multiple-group procedure allowing heterogeneous
#'   covariance matrices described by Sattler and Pauly (2018) and the
#'   multiple-group procedure under equal covariance matrices described by
#'   Sattler (2021).
#'
#' @param data A matrix or a data.frame. For matrix input, subjects are
#'   represented by rows and repeated-measurement dimensions by columns. For
#'   data.frame input, data must contain columns value, subject and dimension,
#'   giving the observations and the subject and dimension IDs.
#' @param hypothesis Either one of `"whole"`, `"sub"`, `"interaction"`,
#'   `"identical"`, or `"flat"`, or a named list containing the projection
#'   matrices `TW` and `TS`; see Details.
#' @param group A one-dimensional atomic vector or factor defining the group
#'   allocation. For matrix input it must contain one entry per row. For vector
#'   input it must contain one entry per measurement.
#' @param AM A single logical value, or alternatively `0` or `1`, specifying
#'   whether the compact representation of the hypothesis matrices described by
#'   Sattler and Rosenbaum (2025) is used. It may reduce the number of rows used
#'   in the calculations without changing the resulting test. The default is
#'   `TRUE`.
#' @param cov.equal A single logical value specifying whether the group
#'   covariance matrices are assumed to be equal. The default is `FALSE`.
#' @param subsampling A single logical value specifying whether the subsampling
#'   versions of all available trace estimators are used in the
#'   heterogeneous-covariance procedure. It has no effect when `cov.equal =
#'   TRUE`; see Details.
#' @param B A single numeric value or arithmetic character expression in `N`
#'   defining the subsampling budget. Character expressions may contain only
#'   numeric constants, `N`, parentheses, and the operators `+`, `-`, `*`, `/`,
#'   and `^`. Its interpretation depends on `cov.equal` and `subsampling`; see
#'   Details.
#' @param seed `NULL` or a single integer-valued number used to make stochastic
#'   calculations reproducible. When supplied, the seed is applied locally and
#'   the previous R random-number state is restored after the calculation.
#'
#' @details For vector input, missing values in `data` cause the entire affected
#'   subject to be removed. The vectors `subject` and `group` must not contain
#'   missing values. For matrix input, every row containing at least one missing
#'   value is removed. A warning is issued whenever incomplete subjects are
#'   dropped.
#'
#'   For vector input, repeated measurements must occur in the same order for
#'   every subject. No separate variable identifying the repeated-measurement
#'   dimension is supplied, so the within-subject order in `data` determines the
#'   component order in the processed data matrix. The observations themselves
#'   need not be globally sorted by subject or group.
#'
#'   At least two groups and two repeated-measurement dimensions are required.
#'   Every group must contain at least six complete subjects.
#'
#' The tested hypothesis has the form \deqn{(\bm T_W \otimes \bm T_S)\bm\mu=\bm
#' 0.} The predefined hypotheses are:
#' \itemize{
#'   \item `"whole"`:
#'   \eqn{\bm T_W=\bm P_a} and
#'   \eqn{\bm T_S=\bm J_d/d}; no whole-plot or group main effect.
#'   \item `"sub"`:
#'   \eqn{\bm T_W=\bm J_a/a} and
#'   \eqn{\bm T_S=\bm P_d}; no subplot or dimension main effect.
#'   \item `"interaction"`:
#'   \eqn{\bm T_W=\bm P_a} and
#'   \eqn{\bm T_S=\bm P_d}; no group-by-dimension interaction.
#'   \item `"identical"`:
#'   \eqn{\bm T_W=\bm P_a} and
#'   \eqn{\bm T_S=\bm I_d}; identical expectation vectors across groups.
#'   \item `"flat"`:
#'   \eqn{\bm T_W=\bm I_a} and
#'   \eqn{\bm T_S=\bm P_d}; a flat expectation profile in every group.
#' }
#'
#'   Alternatively, `hypothesis` may be a named list containing `TW` and `TS`.
#'   Both matrices must be finite, symmetric, idempotent projection matrices
#'   with positive rank. `TW` must have one row and column per analyzed group,
#'   and `TS` must have one row and column per repeated-measurement dimension.
#'   Small numerical deviations within the implemented tolerance are accepted.
#'
#'   When `cov.equal = FALSE`, the method of
#'   Sattler and Pauly (2018) is used. The third-trace estimator
#'   entering `f` is always computed by subsampling. Here, `B` is the base
#'   subsampling budget. Group-specific and pairwise subsampling estimators use
#'   `B` draws for each group or group pair, respectively. The joint third-trace
#'   estimator uses \eqn{aB} joint draws, where each draw simultaneously samples
#'   six subjects from every group. When `subsampling = TRUE`, the remaining
#'   available trace estimators are also replaced by their subsampling versions.
#'
#'   When `cov.equal = TRUE`, the method described in Sattler (2021) is used and
#'   `subsampling` has no effect. The pooled third-trace estimator uses an exact
#'   total of \eqn{aB} six-subject draws across all groups. These draws are
#'   allocated approximately proportionally to \eqn{\binom{n_i}{6}} using the
#'   largest-remainder method, with every group receiving at least one draw.
#'
#'   Even when `subsampling = FALSE`, `f`, `tau`, and `p.value` remain seed
#'   dependent because a third-trace quantity is estimated by subsampling. For
#'   heterogeneous covariance matrices with `subsampling = TRUE`, `statistic` is
#'   seed dependent as well.
#'
#'   `B` may be numeric or an arithmetic character expression involving `N`,
#'   such as `"1000*N"` or `"10*(N + 1)"`. The result is rounded up to the next
#'   integer. Functions, assignments, indexing, and additional variable names
#'   are rejected and are never evaluated.
#'
#'   Upper-tail probabilities are computed directly. Reported p-values are
#'   bounded below by `.Machine$double.eps`; a returned value at this boundary
#'   should be interpreted as no larger than the numerical reporting threshold.
#'
#' @returns A named list of class `"hdrm_grouped"` with components:
#' \describe{
#'   \item{data}{The processed data matrix with subjects in rows,
#'   repeated-measurement dimensions in columns, and subjects ordered by group.}
#'   \item{statistic}{The standardized test statistic \eqn{W}.}
#'   \item{f}{The estimated degrees of freedom.}
#'   \item{tau}{The estimated convergence parameter \eqn{\tau=1/f}.}
#'   \item{H}{A named list containing the projection matrices `TW` and `TS`.}
#'   \item{hypothesis}{The selected predefined hypothesis or `"custom"`.}
#'   \item{p.value}{The upper-tail p-value, bounded below by
#'   `.Machine$double.eps`.}
#'   \item{dim}{A named list containing the repeated-measurement dimension `d`
#'   and the number of analyzed subjects `N`.}
#'   \item{groups}{A named list containing the number of analyzed groups `a`
#'   and their subject counts in `table`.}
#'   \item{subsamples}{The evaluated integer base budget `B`. The grouped
#'   third-trace estimators use \eqn{aB} draws; see Details.}
#' }
#'
#' @example man/examples/examples_hdrm_grouped.R
#'
#' @references Sattler P (2021). “A comprehensive treatment of
#'   quadratic-form-based inference in repeated measures designs under diverse
#'   asymptotics.” Electronic Journal of Statistics, 15(1), 3611–3634. doi:
#'   10.1214/21-EJS1865.
#' @references Sattler P, Pauly M (2018). “Inference for high-dimensional
#'   split-plot-designs: A unified approach for small to large numbers of factor
#'   levels.” Electronic Journal of Statistics, 12(2), 2743–2805. doi:
#'   10.1214/18-EJS1465.
#' @references Sattler P, Rosenbaum M (2025). “Choice of the hypothesis matrix
#'   for using the Anova-type-statistic.” Statistics & Probability Letters, 219,
#'   110356. doi: 10.1016/j.spl.2025.110356.
#'
#' @export
hdrm_grouped <- function(data,
                         hypothesis = "whole",
                         group,
                         AM = TRUE,
                         cov.equal = FALSE,
                         subsampling = FALSE,
                         B = "1000*N",
                         seed = NULL) {
  UseMethod("hdrm_grouped")
}
