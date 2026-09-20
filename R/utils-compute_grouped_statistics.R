#' @keywords internal
compute_grouped_statistic <- function(QN, EW, variance) {
  quantities <- c(QN = QN,
                  EW = EW,
                  variance = variance)
  
  if (any(lengths(list(QN, EW, variance)) != 1L) ||
      !is.numeric(quantities) ||
      anyNA(quantities) ||
      any(!is.finite(quantities))) {
    stop(
      paste0(
        "The quadratic form, its estimated expectation, and its estimated ",
        "variance must be finite numeric scalars."
      ),
      call. = FALSE
    )
  }
  
  if (variance <= 0) {
    stop("The estimated variance of the test statistic must be positive.",
         call. = FALSE)
  }
  
  statistic <- as.numeric((QN - EW) / sqrt(variance))
  
  if (length(statistic) != 1L ||
      is.na(statistic) ||
      !is.finite(statistic)) {
    stop("The standardized test statistic is not finite.", call. = FALSE)
  }
  
  statistic
}