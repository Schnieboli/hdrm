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
  if (
    !(is.logical(AM) || is.numeric(AM)) ||
    length(AM) != 1L ||
    is.na(AM) ||
    !is.finite(AM) ||
    !(AM %in% c(0, 1))
  ) {
    stop(
      "'AM' must be a single logical value or 0/1.",
      call. = FALSE
    )
  }

  AM <- as.logical(AM)
  
  ## checks for cov.equal
  if (
    !is.logical(cov.equal) ||
    length(cov.equal) != 1L ||
    is.na(cov.equal)
  ) {
    stop(
      "'cov.equal' must be a single non-missing logical value.",
      call. = FALSE
    )
  }
  ## checks for subsampling
  if (
    !is.logical(subsampling) ||
    length(subsampling) != 1L ||
    is.na(subsampling)
  ) {
    stop(
      "'subsampling' must be a single non-missing logical value.",
      call. = FALSE
    )
  }

  ## checks for seed
  if (!is.null(seed)) {
    if (
      !is.numeric(seed) ||
      length(seed) != 1L ||
      is.na(seed) ||
      !is.finite(seed) ||
      seed != floor(seed) ||
      abs(seed) > .Machine$integer.max
    ) {
      stop(
        "'seed' must be NULL or a single finite integer-valued number.",
        call. = FALSE
      )
    }

    seed <- as.integer(seed)
  }

  # Identify the two supported data formats
  data_is_vector <- is.numeric(data) && is.null(dim(data))
  data_is_matrix <- is.matrix(data) && is.numeric(data)

  if (!data_is_vector && !data_is_matrix) {
    stop(
      "'data' must be a numeric vector or matrix.",
      call. = FALSE
    )
  }

  # Check that the data contain at least one observation
  if (length(data) == 0L) {
    stop(
      "'data' must not be empty.",
      call. = FALSE
    )
  }

  # Check that 'group' is a one-dimensional atomic vector
  if (!is.atomic(group) || !is.null(dim(group))) {
    stop(
      "'group' must be a one-dimensional atomic vector or factor.",
      call. = FALSE
    )
  }

  # Check that a supplied 'subject' is a one-dimensional atomic vector
  if (
    !is.null(subject) &&
    (!is.atomic(subject) || !is.null(dim(subject)))
  ) {
    stop(
      "'subject' must be a one-dimensional atomic vector or factor.",
      call. = FALSE
    )
  }

  # A subject identifier is required for vector input
  if (data_is_vector && is.null(subject)) {
    stop(
      "'subject' must be provided when 'data' is a vector.",
      call. = FALSE
    )
  }

  # Subject identifiers must be complete
  if (data_is_vector && anyNA(subject)) {
    stop(
      "'subject' must not contain missing values.",
      call. = FALSE
    )
  }

  # Check that the length of 'subject' and 'data' match for vector input
  if (data_is_vector && length(subject) != length(data)) {
    stop(
      "The lengths of 'data' and 'subject' must be equal.",
      call. = FALSE
    )
  }


  # Check that the length of 'group' and 'data' match for vector input
  if (data_is_vector && length(group) != length(data)) {
    stop(
      "The lengths of 'data' and 'group' must be equal.",
      call. = FALSE
    )
  }

  # For matrix input, one group label is required for each subject
  if (data_is_matrix && nrow(data) != length(group)) {
    stop(
      "The length of 'group' must equal the number of rows of 'data' (one group label per subject).",
      call. = FALSE
    )
  }

  # Warn if 'subject' is unnecessarily supplied for matrix input
  if (data_is_matrix && !is.null(subject)) {
    warning(
      "'subject' is ignored when 'data' is a matrix.",
      call. = FALSE
    )
  }
  
  if(any(is.na(data)) || any(is.na(group))){
    stop("'data' and 'group' must not contain missing values", call. = FALSE)
  }



  if (data_is_vector) {  # If 'data' is a vector

    # Create a data frame with 'data', 'subject', and 'group'
    dframe <- data.frame(
      value = data,
      subject = subject,
      whole = group,
      measurement_order = seq_along(data)
    )



    # Convert group labels to a factor
    dframe$whole <- droplevels(as.factor(dframe$whole))

    # Construct subject identifiers that are unique within the full data set.
    # This permits the same subject labels to be reused in different groups.
    dframe$subject <- interaction(
      dframe$whole,
      as.factor(dframe$subject),
      drop = TRUE,
      lex.order = TRUE
    )

    dframe <- dframe[
      order(
        dframe$whole,
        dframe$subject,
        dframe$measurement_order
      ),
      ,
      drop = FALSE
    ]

    ## Store the number of subjects before removing incomplete cases
    N_with_NA <- nlevels(dframe$subject)

    # Identify subjects with at least one missing measurement
    incomplete_subjects <- unique(
      dframe$subject[is.na(dframe$value)]
    )

    # Mark all measurements of incomplete subjects as missing
    if (length(incomplete_subjects) > 0L) {
      dframe$value[
        dframe$subject %in% incomplete_subjects
      ] <- NA_real_
    }

    # Filter out rows with missing values (NA) and remove unnecessary levels
    dframe <- dframe[stats::complete.cases(dframe), ]
    dframe <- droplevels(dframe)

    subject_groups <- unique(
      dframe[c("subject", "whole")]
    )
    group_table <- table(subject_groups$whole)

    a <- nlevels(dframe$whole)

    if (a < 2L) {
      stop(
        "At least two groups must remain after removing incomplete subjects.",
        call. = FALSE
      )
    }

    # Split the data frame by group ('whole')
    L <- split(dframe, dframe$whole)

    # Determine the repeated-measurement dimension in each group
    dimensions <- integer(a)

    for (i in seq_len(a)) {
      subject_dimensions <- unname(table(L[[i]]$subject))
      subject_dimensions <- subject_dimensions[subject_dimensions > 0L]
      unique_dimensions <- unique(subject_dimensions)

      if (length(unique_dimensions) != 1L) {
        stop(
          "All subjects within each group must have the same dimension.",
          call. = FALSE
        )
      }

      dimensions[i] <- unique_dimensions
    }

    # Check that all groups have the same repeated-measurement dimension
    if (length(unique(dimensions)) != 1L) {
      stop(
        "All groups must have the same repeated-measurement dimension.",
        call. = FALSE
      )
    }

    d <- dimensions[[1L]]
    Nv <- integer(a)

    # Get the number of subjects in each group
    for (i in seq_len(a)) {
      Nv[i] <- length(unique(L[[i]]$subject))
    }

    # Calculate the total number of subjects
    N <- sum(Nv)
    # Initialize the data matrix and group vector
    X <- matrix(
      numeric(0),
      nrow = d,
      ncol = 0L
    )
    group <- integer(0)

    # Fill the matrix with data for each group
    for (j in seq_len(a)) {
      temp <- droplevels(L[[j]])
      n_j <- nlevels(temp$subject)

      M <- matrix(
        NA_real_,
        nrow = d,
        ncol = n_j
      )

      group <- c(group, rep.int(j, n_j))

      k <- 1L

      for (i in levels(temp$subject)) {
        M[, k] <- temp$value[temp$subject == i]
        k <- k + 1L
      }

      X <- cbind(X, M)
    }

    # Store the data actually used in the analysis
    out <- list(data = X)
  } else {  # If 'data' is a matrix

    # Remove incomplete subjects from the user-facing N x d matrix and
    # transpose once into the internal d x N representation.
    N_with_NA <- nrow(data)

    complete_subjects <- stats::complete.cases(data)

    group <- group[complete_subjects]
    group <- droplevels(as.factor(group))
    group_table <- table(group)

    X <- t(data[complete_subjects, , drop = FALSE])

    a <- nlevels(group)

    if (a < 2L) {
      stop(
        "At least two groups must remain after removing incomplete subjects.",
        call. = FALSE
      )
    }

    # Get the number of subjects (N) and repeated-measurement dimensions (d)
    # from the internal d x N representation.
    N <- ncol(X)
    d <- nrow(X)



    # Store the processed data in the output
    out <- list(data = t(X))
  }

  # Warning if there were subjects with missing values
  if(N_with_NA > N) warning("Subjects with missing values dropped", call. = FALSE)

  # Convert B to a positive integer without evaluating arbitrary R code
  reps <- evaluate_subsample_budget(
    B = B,
    N = N
  )
  ## reps <- eval(parse(text = B))
  ## if(!is.finite(reps) && reps < .Machine$integer.max)
  
  
  # The grouped third-trace estimators use a * B draws. Validate the
  # effective budget before any stochastic estimator is evaluated.
  expand_subsample_budget(
    B = reps,
    multiplier = a
  ) # TODO kann das weg oder wurde hier vergessen, etwas zuzuweisen?

  # Check the grouping criteria
  check_criteria_grouped(X = X, group = group, hypothesis = hypothesis, reps = reps, subsampling = subsampling)



  ### Output

  # Sort subjects by group while preserving the correspondence
  # between the data columns and the group labels
  group_order <- order(group)
  X_ordered <- X[, group_order, drop = FALSE]
  group_ordered <- as.integer(group[group_order])

  # Store the processed data in the public N x d orientation while the
  # internal calculations continue to use subjects in columns.
  out$data <- t(X_ordered)

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
  out$groups$table <- group_table
  out$removed.cases <- N_with_NA - N
  out$subsamples <- reps
  out$hypothesis = ifelse(is.character(hypothesis), hypothesis[1], "custom")
  class(out) <- "hdrm_grouped"
  return(out)

}
