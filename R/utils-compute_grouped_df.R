#' @keywords internal
compute_grouped_df <- function(second_order, third_order, design_factor) {
  quantities <- c(
    second_order = second_order,
    third_order = third_order,
    design_factor = design_factor
  )

  if (any(lengths(list(second_order, third_order, design_factor)) != 1L) ||
  !is.numeric(quantities) ||
  anyNA(quantities) ||
  any(!is.finite(quantities)) ||
  second_order <= 0 ||
  third_order == 0 ||
  design_factor <= 0) {
    stop("Internal error: please contact 'hdrm' package maintainer.",
         call. = FALSE)
  }

  degrees_of_freedom <- as.numeric((second_order^3 / third_order^2) * design_factor)

  if (length(degrees_of_freedom) != 1L ||
      is.na(degrees_of_freedom) ||
      !is.finite(degrees_of_freedom) ||
      degrees_of_freedom <= 0) {
    stop("Internal error: please contact 'hdrm' package maintainer.",
         call. = FALSE)
  }

  max(1, degrees_of_freedom)
}
