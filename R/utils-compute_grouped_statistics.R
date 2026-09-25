#' @keywords internal
compute_grouped_statistic <- function(QN, EW, variance) {
  quantities <- c(QN = QN,
                  EW = EW,
                  variance = variance)

  statistic <- as.numeric((QN - EW) / sqrt(variance))

  if (any(lengths(list(QN, EW, variance)) != 1L) ||
      !is.numeric(quantities) ||
      anyNA(quantities) ||
      any(!is.finite(quantities)) ||
      variance <= 0 ||
      length(statistic) != 1L ||
      is.na(statistic) ||
      !is.finite(statistic)
      ) {
    stop(
      paste0(
        "Internal error: please contact 'hdrm' package maintainer."
      ),
      call. = FALSE
    )
  }

  statistic
}
