#' @title Build the hypothesis matrices for a multifactorial setting
#'
#' @description
#' Extracts or constructs the whole-plot and subplot hypothesis matrices used
#' by the grouped procedures.
#'
#' @param hypothesis Either one of `"whole"`, `"sub"`, `"interaction"`,
#' `"identical"`, or `"flat"`, or a named list containing `TW` and `TS`.
#' @param a A positive integer giving the number of groups.
#' @param d A positive integer giving the repeated-measurement dimension.
#'
#' @returns A named list with components `TW` and `TS`.
#'
#' @noRd
get_hypothesis_mult <- function(hypothesis, a, d) {
  if (!is.list(hypothesis) && !is.character(hypothesis)) {
    stop(
      "'hypothesis' must be a character value or a list.",
      call. = FALSE
    )
  }
  
  if (is.list(hypothesis)) {
    if (is.null(hypothesis$TW)) {
      stop(
        "No list entry named 'TW' was found in 'hypothesis'.",
        call. = FALSE
      )
    }
    
    if (is.null(hypothesis$TS)) {
      stop(
        "No list entry named 'TS' was found in 'hypothesis'.",
        call. = FALSE
      )
    }
    
    TW <- hypothesis$TW
    TS <- hypothesis$TS
    
    if (!is.matrix(TW) || !is.numeric(TW)) {
      stop(
        "'TW' must be a numeric matrix.",
        call. = FALSE
      )
    }
    
    if (!is.matrix(TS) || !is.numeric(TS)) {
      stop(
        "'TS' must be a numeric matrix.",
        call. = FALSE
      )
    }
    
    if (anyNA(TW) || any(!is.finite(TW))) {
      stop(
        "'TW' must contain only finite, non-missing values.",
        call. = FALSE
      )
    }
    
    if (anyNA(TS) || any(!is.finite(TS))) {
      stop(
        "'TS' must contain only finite, non-missing values.",
        call. = FALSE
      )
    }
    
    if (
      length(dim(TW)) != 2L ||
      any(dim(TW) != c(a, a))
    ) {
      stop(
        paste0(
          "'TW' must be a ",
          a,
          " by ",
          a,
          " matrix."
        ),
        call. = FALSE
      )
    }
    
    if (
      length(dim(TS)) != 2L ||
      any(dim(TS) != c(d, d))
    ) {
      stop(
        paste0(
          "'TS' must be a ",
          d,
          " by ",
          d,
          " matrix."
        ),
        call. = FALSE
      )
    }
    
    tol <- sqrt(.Machine$double.eps)
    
    TW_symmetry_error <- max(abs(TW - t(TW)))
    TW_idempotence_error <- max(abs(TW %*% TW - TW))
    TS_symmetry_error <- max(abs(TS - t(TS)))
    TS_idempotence_error <- max(abs(TS %*% TS - TS))
    
    if (TW_symmetry_error > tol) {
      stop(
        paste0(
          "'TW' must be symmetric. Maximum deviation: ",
          signif(TW_symmetry_error, 4),
          "."
        ),
        call. = FALSE
      )
    }
    
    if (TW_idempotence_error > tol) {
      stop(
        paste0(
          "'TW' must be idempotent. Maximum deviation: ",
          signif(TW_idempotence_error, 4),
          "."
        ),
        call. = FALSE
      )
    }
    
    if (TS_symmetry_error > tol) {
      stop(
        paste0(
          "'TS' must be symmetric. Maximum deviation: ",
          signif(TS_symmetry_error, 4),
          "."
        ),
        call. = FALSE
      )
    }
    
    if (TS_idempotence_error > tol) {
      stop(
        paste0(
          "'TS' must be idempotent. Maximum deviation: ",
          signif(TS_idempotence_error, 4),
          "."
        ),
        call. = FALSE
      )
    }
    
    if (qr(TW, tol = tol)$rank == 0L) {
      stop(
        "'TW' must have positive rank.",
        call. = FALSE
      )
    }
    
    if (qr(TS, tol = tol)$rank == 0L) {
      stop(
        "'TS' must have positive rank.",
        call. = FALSE
      )
    }
  } else {
    if (
      length(hypothesis) != 1L ||
      is.na(hypothesis)
    ) {
      stop(
        "'hypothesis' must be a single non-missing character value.",
        call. = FALSE
      )
    }
    
    if (
      !hypothesis %in%
      c("whole", "sub", "interaction", "identical", "flat")
    ) {
      stop(
        paste0(
          "'hypothesis' must be one of 'whole', 'sub', 'interaction', ",
          "'identical', or 'flat', or a list containing 'TW' and 'TS'."
        ),
        call. = FALSE
      )
    }
    
    P_a <- diag(a) - matrix(1 / a, nrow = a, ncol = a)
    P_d <- diag(d) - matrix(1 / d, nrow = d, ncol = d)
    J_a <- matrix(1 / a, nrow = a, ncol = a)
    J_d <- matrix(1 / d, nrow = d, ncol = d)
    
    if (hypothesis == "whole") {
      TW <- P_a
      TS <- J_d
    } else if (hypothesis == "sub") {
      TW <- J_a
      TS <- P_d
    } else if (hypothesis == "interaction") {
      TW <- P_a
      TS <- P_d
    } else if (hypothesis == "identical") {
      TW <- P_a
      TS <- diag(d)
    } else {
      TW <- diag(a)
      TS <- P_d
    }
  }
  
  list(TW = TW, TS = TS)
}
