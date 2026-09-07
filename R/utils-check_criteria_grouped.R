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
