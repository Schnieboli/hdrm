#' @title Allocate the equal-covariance third-trace budget across groups
#'
#' @description
#' Treats `B` as a base budget per group and allocates the resulting exact total
#' budget \eqn{aB} approximately proportionally to the numbers of available
#' six-subject subsets, \eqn{\binom{n_i}{6}}. Every group receives at least one
#' subsample. The remaining budget is distributed by the largest-remainder
#' method, with ties resolved by group order.
#'
#' @param group_sizes Integer group sizes, all at least six.
#' @param B A positive integer giving the base subsampling budget.
#'
#' @returns An integer vector of group-specific subsample counts whose sum is
#' exactly \eqn{aB}.
#'
#' @noRd
allocate_C1_subsamples <- function(group_sizes, B) {
  if (
    !is.numeric(group_sizes) ||
    length(group_sizes) < 1L ||
    anyNA(group_sizes) ||
    any(!is.finite(group_sizes)) ||
    any(group_sizes != floor(group_sizes)) ||
    any(group_sizes < 6L)
  ) {
    stop(
      "'group_sizes' must contain finite integers of at least 6.",
      call. = FALSE
    )
  }
  
  group_sizes <- as.integer(group_sizes)
  a <- length(group_sizes)
  
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
  
  total_budget <- expand_subsample_budget(
    B = B,
    multiplier = a
  )
  
  combination_counts <- choose(group_sizes, 6)
  
  if (
    any(!is.finite(combination_counts)) ||
    sum(combination_counts) <= 0
  ) {
    stop(
      "The six-subject combination weights could not be computed.",
      call. = FALSE
    )
  }
  
  remaining_budget <- total_budget - a
  
  if (remaining_budget == 0L) {
    return(rep.int(1L, a))
  }
  
  raw_allocation <- remaining_budget *
    combination_counts /
    sum(combination_counts)
  
  integer_allocation <- floor(raw_allocation)
  allocation <- 1L + as.integer(integer_allocation)
  
  remainder <- total_budget - sum(allocation)
  
  if (remainder > 0L) {
    fractional_parts <- raw_allocation - integer_allocation
    
    priority <- order(
      -fractional_parts,
      seq_along(fractional_parts)
    )
    
    selected_groups <- priority[seq_len(remainder)]
    allocation[selected_groups] <-
      allocation[selected_groups] + 1L
  }
  
  if (
    any(allocation < 1L) ||
    sum(allocation) != total_budget
  ) {
    stop(
      "Internal error while allocating the subsampling budget.",
      call. = FALSE
    )
  }
  
  allocation
}
