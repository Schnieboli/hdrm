#' Expand a base subsampling budget by an integer multiplier
#'
#' Internal helper used when a grouped third-trace estimator deliberately uses
#' more than `B` draws. It also protects the C++ interface from integer
#' overflow.
#'
#' @param B A positive integer base budget.
#' @param multiplier A positive integer multiplier.
#'
#' @returns The integer product `B * multiplier`.
#'
#' @noRd
expand_subsample_budget <- function(B, multiplier) {
  if (
    !is.numeric(B) ||
    length(B) != 1L ||
    is.na(B) ||
    !is.finite(B) ||
    B != floor(B) ||
    B < 1
  ) {
    stop(
      "'B' must be a finite positive integer.",
      call. = FALSE
    )
  }
  
  if (
    !is.numeric(multiplier) ||
    length(multiplier) != 1L ||
    is.na(multiplier) ||
    !is.finite(multiplier) ||
    multiplier != floor(multiplier) ||
    multiplier < 1
  ) {
    stop(
      "'multiplier' must be a finite positive integer.",
      call. = FALSE
    )
  }
  
  effective_budget <- as.double(B) * as.double(multiplier)
  
  if (
    !is.finite(effective_budget) ||
    effective_budget > .Machine$integer.max
  ) {
    stop(
      paste0(
        "The effective subsampling budget must not exceed ",
        .Machine$integer.max,
        "."
      ),
      call. = FALSE
    )
  }
  as.integer(effective_budget)
}