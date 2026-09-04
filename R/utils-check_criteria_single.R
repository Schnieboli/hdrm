#' @title Check prepared inputs for the one-group test
#'
#' @description
#' Checks the assumptions required by the internal one-group procedure after
#' the public wrapper has prepared the data.
#'
#' @param X A numeric data matrix with subjects in columns.
#' @param hypothesis Either `"flat"` or a numeric matrix.
#'
#' @returns `NULL` invisibly. An error is raised when a condition is violated.
#'
#' @noRd
check_criteria_single <- function(X, hypothesis) {
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

  if (nrow(X) < 2L) {
    stop(
      "At least two repeated-measurement dimensions are required.",
      call. = FALSE
    )
  }

  if (ncol(X) < 3L) {
    stop(
      "At least three complete subjects are required.",
      call. = FALSE
    )
  }

  if (!is.character(hypothesis) && !is.matrix(hypothesis)) {
    stop(
      "'hypothesis' must be 'flat' or a numeric projection matrix.",
      call. = FALSE
    )
  }

  invisible(NULL)
}
