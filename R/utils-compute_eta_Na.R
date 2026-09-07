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