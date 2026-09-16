#' Multiple-group inference for high-dimensional repeated measures
#'
#' @description Implements the multiple-group procedure allowing heterogeneous
#' covariance matrices of \insertCite{Sattler2018;textual}{hdrm} and the
#' multiple-group procedure under equal covariance matrices of
#' \insertCite{Sattler2021;textual}{hdrm}.
#'
#' @param data A numeric vector or matrix. For matrix input, subjects are
#'   represented by rows and repeated-measurement dimensions by columns. For
#'   data.frame input, colums `subject` and `time` must identify the subject and
#'   time associated with every measurement.
#' @param hypothesis Either one of `"whole"`, `"sub"`, `"interaction"`,
#'   `"identical"`, or `"flat"`, or a named list containing the projection
#'   matrices `TW` and `TS`; see Details.
#' @param AM A single logical value, or alternatively `0` or `1`, specifying
#'   whether the compact representation of the hypothesis matrices described by
#'   \insertCite{Sattler2025;textual}{hdrm} is used. It may reduce the number of
#'   rows used in the calculations without changing the resulting test. The
#'   default is `TRUE`.
#' @param group A one-dimensional atomic vector or factor defining the group
#'   allocation. For matrix input it must contain one entry per row. For vector
#'   input it must contain one entry per measurement.
#' @param subject An optional one-dimensional atomic vector identifying
#'   subjects. It is required for vector input and ignored with a warning for
#'   matrix input. Subject labels need only be unique within groups and may be
#'   reused in different groups.
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
#' subject to be removed. The vectors `subject` and `group` must not contain
#' missing values. For matrix input, every row containing at least one missing
#' value is removed. A warning is issued whenever incomplete subjects are
#' dropped.
#'
#' For vector input, repeated measurements must occur in the same order for
#' every subject. No separate variable identifying the repeated-measurement
#' dimension is supplied, so the within-subject order in `data` determines the
#' component order in the processed data matrix. The observations themselves
#' need not be globally sorted by subject or group.
#'
#' At least two groups and two repeated-measurement dimensions are required.
#' Every group must contain at least six complete subjects.
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
#' Alternatively, `hypothesis` may be a named list containing `TW` and `TS`.
#' Both matrices must be finite, symmetric, idempotent projection matrices with
#' positive rank. `TW` must have one row and column per analyzed group, and `TS`
#' must have one row and column per repeated-measurement dimension. Small
#' numerical deviations within the implemented tolerance are accepted.
#'
#' When `cov.equal = FALSE`, the method of
#' \insertCite{Sattler2018;textual}{hdrm} is used. The third-trace estimator
#' entering `f` is always computed by subsampling. Here, `B` is the base
#' subsampling budget. Group-specific and pairwise subsampling estimators use
#' `B` draws for each group or group pair, respectively. The joint third-trace
#' estimator uses \eqn{aB} joint draws, where each draw simultaneously samples
#' six subjects from every group. When `subsampling = TRUE`, the remaining
#' available trace estimators are also replaced by their subsampling versions.
#'
#' When `cov.equal = TRUE`, the method of \insertCite{Sattler2021;textual}{hdrm}
#' is used and `subsampling` has no effect. The pooled third-trace estimator
#' uses an exact total of \eqn{aB} six-subject draws across all groups. These
#' draws are allocated approximately proportionally to \eqn{\binom{n_i}{6}}
#' using the largest-remainder method, with every group receiving at least one
#' draw.
#'
#' Even when `subsampling = FALSE`, `f`, `tau`, and `p.value` remain seed
#' dependent because a third-trace quantity is estimated by subsampling. For
#' heterogeneous covariance matrices with `subsampling = TRUE`, `statistic` is
#' seed dependent as well.
#'
#' `B` may be numeric or an arithmetic character expression involving `N`, such
#' as `"1000*N"` or `"10*(N + 1)"`. The result is rounded up to the next
#' integer. Functions, assignments, indexing, and additional variable names are
#' rejected and are never evaluated.
#'
#' Upper-tail probabilities are computed directly. Reported p-values are bounded
#' below by `.Machine$double.eps`; a returned value at this boundary should be
#' interpreted as no larger than the numerical reporting threshold.
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
#'   \item{removed.cases}{The number of incomplete subjects removed before the
#'   analysis.}
#'   \item{subsamples}{The evaluated integer base budget `B`. The grouped
#'   third-trace estimators use \eqn{aB} draws; see Details.}
#' }
#'
#' @example man/examples/examples_hdrm_grouped.R
#'
#' @references \insertAllCited
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
  ## checks for AM
  AM <- as.logical(AM)
  if (length(AM) != 1L || is.na(AM)) {
    stop("'AM' must be a single logical value or 0/1.", call. = FALSE)
  }
  ## checks for cov.equal
  cov.equal <- as.logical(cov.equal)
  if (length(cov.equal) != 1L || is.na(cov.equal)) {
    stop("'cov.equal' must be a single non-missing logical value.", 
         call. = FALSE)
  }
  ## checks for subsampling
  subsampling <- as.logical(subsampling)
  if (length(subsampling) != 1L || is.na(subsampling)) {
    stop("'subsampling' must be a single non-missing logical value.",
         call. = FALSE)
  }

  ## checks for seed
  if (!is.null(seed)) {
    if (!is.numeric(seed) ||
        length(seed) != 1L ||
        is.na(seed) ||
        seed != floor(seed) ||
        abs(seed) > .Machine$integer.max) {
      stop("'seed' must be NULL or a single finite integer-valued number.",
           call. = FALSE)
    }

    seed <- as.integer(seed)
  }

  # Identify the two supported data formats
  data_is_df <- is.data.frame(data)
  data_is_matrix <- is.matrix(data) && is.numeric(data)
  
  if (!data_is_df && !data_is_matrix) {
    stop("'data' must be a data frame or matrix.", call. = FALSE)
  }
  
  ## initialize output object
  # Store the processed data in the public N x d orientation while the
  # internal calculations continue to use subjects in columns.
  out <- list(data = data)
  
  ## do all matrix related checks
  if(data_is_matrix){
    data <- t(data)
    if(any(dim(data) == 0)){
      stop("'data' must not be empty.")
    }
    # Check that 'group' is a one-dimensional atomic vector
    if (!is.atomic(group) || length(group) != ncol(data) || !is.null(dim(group))) {
      stop("'group' must be a one-dimensional atomic vector or factor.",
           call. = FALSE)
    }
    d <- nrow(data)
    N <- ncol(data)
    data_list <- lapply(split(seq_len(N), group), function(cols) data[, cols, drop = FALSE])
  }
  
  if(data_is_df){
    if(is.null(data$value) || is.null(data$subject) || is.null(data$time)){
      stop("data must contain columns 'value', 'subject' and 'time'", 
           call. = FALSE)
    }
    
    data <- data.frame(value = data$value, subject = data$subject, time = data$time)
    
    if(nrow(data) < 1){
      stop("'data' must not be empty.", call. = FALSE)
    }
    if(!is.numeric(data$value) || any(!is.finite(data$value))){
      stop("data$value must be numeric and finite", call. = FALSE)
    }
    
    if(!is.atomic(group) || !is.null(dim(group)) || length(group) != nrow(data)){
      stop("'group' must be a vector of length nrow(data)", call. = FALSE)
    }
    
    data$group <- group
    data <- data[order(data$subject, data$time, data$group), ]
    ## reshape data to widetable format
    df_wide <- reshape(
      data,
      idvar = c("subject", "group"),
      timevar = "time",
      direction = "wide"
    )
    ## split df_wide into groups, transform to matrix
    data_list <- lapply(split(df_wide[, -c(1, 2)], df_wide$group), 
                       function(x) unname(as.matrix(t(x))))
    N <- sum(sapply(data_list, ncol))
    d <- nrow(data_list[[1]])
  }
  
  ## check mathematical requirements
  a <- length(data_list)
  n <- sapply(data_list, ncol)
  if(a < 2) stop("there must be at least two groups", call. = FALSE)
  if(d < 2) stop("there must be at least two observations per subject", call. = FALSE)
  if(any(n < 6)) stop("there must be at least six subjects per group", call. = FALSE)
  
  
  
  # Convert B to a positive integer without evaluating arbitrary R code
  reps <- evaluate_subsample_budget(B = B, N = N)
  ## reps <- eval(parse(text = B))
  ## if(!is.finite(reps) && reps < .Machine$integer.max)
  
  
  # The grouped third-trace estimators use a * B draws. Validate the
  # effective budget before any stochastic estimator is evaluated.
  expand_subsample_budget(B = reps, multiplier = a) # TODO kann das weg oder wurde hier vergessen, etwas zuzuweisen?
  

  # Get the hypothesis matrices based on the provided hypothesis
  H <- get_hypothesis_mult(hypothesis, AM, a, d)
  
  
  ### Output
  if (cov.equal) {
    out <- c(
      out,
      hdrm_grouped_eq_cov_internal(
        X_list = data_list,
        H = H,
        B = reps,
        seed = seed
      )
    )
  } else {
    out <- c(
      out,
      hdrm_grouped_internal(
        X_list = data_list,
        H = H,
        subsampling = subsampling,
        B = reps,
        seed = seed
      )
    )
  }

  # Add further output to the result
  out$removed.cases <- N - N
  out$subsamples <- reps
  out$groups = list(a = a, table = table(group))  # Grouping information
  # Description of the hypothesis
  out$hypothesis = ifelse(is.character(hypothesis), hypothesis[1], "custom")
  class(out) <- "hdrm_grouped"
  return(out)

}
