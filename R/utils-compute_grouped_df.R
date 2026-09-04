#' @keywords internal
compute_grouped_df <- function(
    second_order,
    third_order,
    design_factor = 1
) {
  quantities <- c(
    second_order = second_order,
    third_order = third_order,
    design_factor = design_factor
  )
  
  if (
    any(lengths(list(second_order, third_order, design_factor)) != 1L) ||
    !is.numeric(quantities) ||
    anyNA(quantities) ||
    any(!is.finite(quantities))
  ) {
    stop(
      paste0(
        "The second-order trace estimate, third-order trace estimate, and ",
        "design factor must be finite numeric scalars."
      ),
      call. = FALSE
    )
  }
  
  if (second_order <= 0) {
    stop(
      "The second-order trace estimate must be positive.",
      call. = FALSE
    )
  }
  
  if (third_order == 0) {
    stop(
      paste0(
        "The third-order trace estimate is zero. Increase 'B' or check ",
        "whether the data are degenerate."
      ),
      call. = FALSE
    )
  }
  
  if (design_factor <= 0) {
    stop(
      "The design factor must be positive.",
      call. = FALSE
    )
  }
  
  degrees_of_freedom <- as.numeric(
    (second_order^3 / third_order^2) * design_factor
  )
  
  if (
    length(degrees_of_freedom) != 1L ||
    is.na(degrees_of_freedom) ||
    !is.finite(degrees_of_freedom) ||
    degrees_of_freedom <= 0
  ) {
    stop(
      paste0(
        "The estimated degrees-of-freedom parameter is not finite and ",
        "positive."
      ),
      call. = FALSE
    )
  }
  
  max(
    1,
    degrees_of_freedom
  )
}