#' @keywords internal
compute_grouped_p_value <- function(statistic, degrees_of_freedom) {
  quantities <- c(statistic = statistic, degrees_of_freedom = degrees_of_freedom)

  if (any(lengths(list(statistic, degrees_of_freedom)) != 1L) ||
      !is.numeric(quantities) ||
      anyNA(quantities) ||
      any(!is.finite(quantities)) ||
      degrees_of_freedom <= 0) {
    stop(
      paste0(
        "Internal error: please contact 'hdrm' package maintainer."
      ),
      call. = FALSE
    )
  }

  p_value <- max(
    stats::pchisq(
      statistic * sqrt(2 * degrees_of_freedom) +
        degrees_of_freedom,
      df = degrees_of_freedom,
      lower.tail = FALSE
    ),
    .Machine$double.eps
  )

  if (length(p_value) != 1L ||
      is.na(p_value) ||
      !is.finite(p_value) ||
      p_value < 0 ||
      p_value > 1) {
    stop("Internal error: please contact 'hdrm' package maintainer.",
         call. = FALSE)
  }

  as.numeric(p_value)
}
