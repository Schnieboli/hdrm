#' @title Check the prepared inputs for the grouped test
#'
#' @description
#' Checks the assumptions required by the internal grouped procedures after the
#' public wrapper has prepared the data.
#'
#' @param X A numeric data matrix with subjects in columns.
#' @param group An integer vector containing consecutive group codes starting
#' at one.
#' @param hypothesis A predefined hypothesis name or a list containing `TW`
#' and `TS`.
#' @param reps A positive integer giving the requested number of subsamples.
#' @param subsampling A single non-missing logical value.
#'
#' @returns `NULL` invisibly. An error is raised when a condition is violated.
#'
#' @noRd
check_criteria_grouped <- function(
    X,
    group,
    hypothesis,
    reps,
    subsampling
) {
  if (!is.matrix(X) || !is.numeric(X)) {
    stop(
      "'X' must be a numeric matrix.",
      call. = FALSE
    )
  }
  
  if (anyNA(X) || any(!is.finite(X))) {
    stop(
      "'X' must contain only finite, non-missing values.",
      call. = FALSE
    )
  }
  
  d <- nrow(X)
  N <- ncol(X)
  
  if (
    !is.null(dim(group)) ||
    length(group) != N ||
    anyNA(group)
  ) {
    stop(
      "'group' must contain one non-missing group label per column of 'X'.",
      call. = FALSE
    )
  }
  
  if (is.factor(group)) {
    group_codes <- as.integer(group)
  } else {
    if (
      !is.numeric(group) ||
      any(!is.finite(group)) ||
      any(group != floor(group))
    ) {
      stop(
        "'group' must be a factor or an integer-valued vector.",
        call. = FALSE
      )
    }
    
    group_codes <- as.integer(group)
  }
  
  observed_groups <- sort(unique(group_codes))
  a <- length(observed_groups)
  
  if (!identical(observed_groups, seq_len(a))) {
    stop(
      "Internal group codes must be consecutive integers starting at one.",
      call. = FALSE
    )
  }
  
  n <- tabulate(group_codes, nbins = a)
  
  if (a < 2L) {
    stop(
      "At least two groups are required.",
      call. = FALSE
    )
  }
  
  if (d < 2L) {
    stop(
      "At least two repeated-measurement dimensions are required.",
      call. = FALSE
    )
  }
  
  if (any(n < 6L)) {
    stop(
      "All group sizes must be at least 6.",
      call. = FALSE
    )
  }
  
  if (!is.character(hypothesis) && !is.list(hypothesis)) {
    stop(
      "'hypothesis' must be a character value or a list.",
      call. = FALSE
    )
  }
  
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
  
  if (
    !is.numeric(reps) ||
    length(reps) != 1L ||
    is.na(reps) ||
    !is.finite(reps) ||
    reps < 1 ||
    reps != floor(reps) ||
    reps > .Machine$integer.max
  ) {
    stop(
      "The number of subsamples must be a finite positive integer.",
      call. = FALSE
    )
  }
  
  invisible(NULL)
}


################################################################################
# Whole-plot design factor
################################################################################

# Internal helper for the known factor eta_{N,a} in the equal-covariance
# procedure. The calculation is kept separate so that it can be tested
# independently from the stochastic trace estimators.
compute_eta_Na <- function(TW, group_sizes) {
  if (
    !is.matrix(TW) ||
    !is.numeric(TW) ||
    nrow(TW) != ncol(TW) ||
    anyNA(TW) ||
    any(!is.finite(TW))
  ) {
    stop(
      "'TW' must be a finite numeric square matrix.",
      call. = FALSE
    )
  }
  
  a <- nrow(TW)
  
  if (
    !is.numeric(group_sizes) ||
    !is.null(dim(group_sizes)) ||
    length(group_sizes) != a ||
    anyNA(group_sizes) ||
    any(!is.finite(group_sizes)) ||
    any(group_sizes <= 0) ||
    any(group_sizes != floor(group_sizes))
  ) {
    stop(
      paste0(
        "'group_sizes' must contain one positive integer for each ",
        "row of 'TW'."
      ),
      call. = FALSE
    )
  }
  
  symmetry_error <- max(abs(TW - t(TW)))
  
  if (symmetry_error > sqrt(.Machine$double.eps)) {
    stop(
      "'TW' must be symmetric.",
      call. = FALSE
    )
  }
  
  N <- sum(group_sizes)
  D_N <- diag(
    N / group_sizes,
    nrow = a,
    ncol = a
  )
  DNTW <- D_N %*% TW
  
  trace2_whole <- sum(
    diag(DNTW %*% DNTW)
  )
  trace3_whole <- sum(
    diag(DNTW %*% DNTW %*% DNTW)
  )
  
  if (
    !is.finite(trace2_whole) ||
    !is.finite(trace3_whole) ||
    trace2_whole <= 0 ||
    trace3_whole == 0
  ) {
    stop(
      "The whole-plot trace factor is degenerate.",
      call. = FALSE
    )
  }
  
  eta_Na <- trace2_whole^3 / trace3_whole^2
  
  if (
    !is.finite(eta_Na) ||
    eta_Na <= 0
  ) {
    stop(
      "The whole-plot trace factor is not finite and positive.",
      call. = FALSE
    )
  }
  
  as.numeric(eta_Na)
}