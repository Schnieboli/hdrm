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
    subject = NULL,
    AM = TRUE
) {
  if ( # test if AM is logical
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

  if (is.character(hypothesis)) { # check if hypothesis as character is correct
    if (
      length(hypothesis) != 1L ||
      is.na(hypothesis)
    ) {
      stop(
        "'hypothesis' must be a single non-missing character value.",
        call. = FALSE
      )
    }

    if (hypothesis != "flat") {
      stop(
        "'hypothesis' must be 'flat' or a numeric projection matrix.",
        call. = FALSE
      )
    }
  } else if (!is.matrix(hypothesis)) {
    stop(
      "'hypothesis' must be 'flat' or a numeric projection matrix.",
      call. = FALSE
    )
  }

  data_is_vector <- is.numeric(data) && is.null(dim(data))
  data_is_matrix <- is.matrix(data) && is.numeric(data)

  if (!data_is_vector && !data_is_matrix) {
    stop(
      "'data' must be a numeric vector or matrix.",
      call. = FALSE
    )
  }

  if (length(data) == 0L) {
    stop(
      "'data' must not be empty.",
      call. = FALSE
    )
  }

  if (
    !is.null(subject) &&
    (!is.atomic(subject) || !is.null(dim(subject)))
  ) {
    stop(
      "'subject' must be a one-dimensional atomic vector or factor.",
      call. = FALSE
    )
  }

  if (data_is_vector && is.null(subject)) {
    stop(
      "'subject' must be provided when 'data' is a vector.",
      call. = FALSE
    )
  }

  if (data_is_vector && anyNA(subject)) {
    stop(
      "'subject' must not contain missing values.",
      call. = FALSE
    )
  }

  if (data_is_vector && length(subject) != length(data)) {
    stop(
      "The lengths of 'data' and 'subject' must be equal.",
      call. = FALSE
    )
  }

  if (data_is_matrix && !is.null(subject)) {
    warning(
      "'subject' is ignored when 'data' is a matrix.",
      call. = FALSE
    )
  }

  if (data_is_vector) {
    subject_values <- as.character(subject)
    subject_factor <- factor(
      subject_values,
      levels = unique(subject_values)
    )

    dframe <- data.frame(
      value = as.numeric(data),
      subject = subject_factor,
      measurement_order = seq_along(data)
    )

    N_with_NA <- nlevels(dframe$subject)

    incomplete_subjects <- unique(
      dframe$subject[is.na(dframe$value)]
    )

    if (length(incomplete_subjects) > 0L) {
      dframe <- dframe[
        !(dframe$subject %in% incomplete_subjects),
        ,
        drop = FALSE
      ]
    }

    dframe <- droplevels(dframe)

    if (nrow(dframe) == 0L) {
      stop(
        "No complete subjects remain after removing missing values.",
        call. = FALSE
      )
    }

    subject_dimensions <- unname(table(dframe$subject))

    if (
      length(subject_dimensions) == 0L ||
      length(unique(subject_dimensions)) != 1L
    ) {
      stop(
        "All complete subjects must have the same number of dimensions.",
        call. = FALSE
      )
    }

    d <- as.integer(subject_dimensions[[1L]])
    subject_levels <- levels(dframe$subject)
    N <- length(subject_levels)

    X <- vapply(
      subject_levels,
      function(current_subject) {
        dframe$value[dframe$subject == current_subject]
      },
      numeric(d)
    )

    X <- matrix(
      X,
      nrow = d,
      ncol = N,
      dimnames = NULL
    )
  } else {
    N_with_NA <- nrow(data)
    complete_subjects <- stats::complete.cases(data)
    X <- t(data[complete_subjects, , drop = FALSE])
    d <- nrow(X)
    N <- ncol(X)
  }

  if (N_with_NA > N) {
    warning(
      "Subjects with missing values dropped",
      call. = FALSE
    )
  }

  check_criteria_single(
    X = X,
    hypothesis = hypothesis
  )

  if (is.character(hypothesis)) {
    T_matrix <- diag(d) - matrix(
      1 / d,
      nrow = d,
      ncol = d
    )
    hypothesis_label <- hypothesis
  } else {
    if (!is.numeric(hypothesis)) {
      stop(
        "The hypothesis matrix must be numeric.",
        call. = FALSE
      )
    }

    if (anyNA(hypothesis) || any(!is.finite(hypothesis))) {
      stop(
        "The hypothesis matrix must contain only finite, non-missing values.",
        call. = FALSE
      )
    }

    if (
      length(dim(hypothesis)) != 2L ||
      any(dim(hypothesis) != c(d, d))
    ) {
      stop(
        paste0(
          "The hypothesis matrix must be a ",
          d,
          " by ",
          d,
          " matrix."
        ),
        call. = FALSE
      )
    }

    T_matrix <- hypothesis
    hypothesis_label <- "custom"
  }

  tol <- sqrt(.Machine$double.eps)

  symmetry_error <- max(
    abs(T_matrix - t(T_matrix))
  )
  idempotence_error <- max(
    abs(T_matrix %*% T_matrix - T_matrix)
  )

  if (symmetry_error > tol) {
    stop(
      paste0(
        "The hypothesis matrix must be symmetric. Maximum deviation: ",
        signif(symmetry_error, 4),
        "."
      ),
      call. = FALSE
    )
  }

  if (idempotence_error > tol) {
    stop(
      paste0(
        "The hypothesis matrix must be idempotent. Maximum deviation: ",
        signif(idempotence_error, 4),
        "."
      ),
      call. = FALSE
    )
  }

  if (qr(T_matrix, tol = tol)$rank == 0L) {
    stop(
      "The hypothesis matrix must have positive rank.",
      call. = FALSE
    )
  }

  computation_matrix <- if (AM) {
    MSrootcompact(T_matrix)
  } else {
    T_matrix
  }

  XT <- computation_matrix %*% X
  transformed_mean <- rowMeans(XT)
  Qn <- N * sum(transformed_mean^2)

  traceNormal <- B0_cpp(XT)
  traceSquare <- B2_cpp(XT)
  traceCubic <- B3_cpp(XT) / choose(N, 3)

  if (
    !is.finite(Qn) ||
    !is.finite(traceNormal) ||
    !is.finite(traceSquare) ||
    traceSquare <= 0 ||
    !is.finite(traceCubic) ||
    traceCubic == 0
  ) {
    stop(
      "The trace estimators are numerically degenerate for these data.",
      call. = FALSE
    )
  }

  W <- (Qn - traceNormal) / sqrt(2 * traceSquare)

  if (!is.finite(W)) {
    stop(
      "The test statistic is numerically undefined for these data.",
      call. = FALSE
    )
  }

  f_raw <- traceSquare^3 / traceCubic^2

  if (!is.finite(f_raw) || f_raw <= 0) {
    stop(
      "The estimated degrees of freedom are numerically undefined.",
      call. = FALSE
    )
  }

  f <- max(1, f_raw)

  p.value <- max(
    stats::pchisq(
      W * sqrt(2 * f) + f,
      df = f,
      lower.tail = FALSE
    ),
    .Machine$double.eps
  )

  out <- list(
    data = t(X),
    f = f,
    statistic = W,
    tau = 1 / f,
    H = computation_matrix,
    hypothesis = hypothesis_label,
    p.value = p.value,
    dim = c(d = d, N = N),
    removed.cases = N_with_NA - N
  )

  class(out) <- "hdrm_single"
  out
}


# Print method ------------------------------------------------------------

#' @method print hdrm_single
#' @export
print.hdrm_single <- function(x, digits = 4, ...) {
  p <- round(x$p.value, digits)

  if (is.na(p)) {
    p_text <- "= NA"
  } else if (p <= 0) {
    p_text <- paste0("< ", 10^(-digits))
  } else {
    p_text <- paste0("= ", p)
  }

  cat(
    "\n",
    "           One Group Repeated Measure\n",
    "Analysis of ", x$dim[["N"]],
    " subjects in ", x$dim[["d"]],
    " dimensions:",
    "\nW = ", round(x$statistic, digits),
    "  f = ", round(x$f, digits),
    "  p.value ", p_text,
    "\nHypothesis type: ", x$hypothesis,
    "\nConvergence parameter \u03c4 = ", round(x$tau, digits),
    "\n",
    sep = ""
  )

  invisible(x)
}
