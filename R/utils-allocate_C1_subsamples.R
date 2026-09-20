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
  group_sizes <- as.integer(group_sizes)
  a <- length(group_sizes)
  
  B_tot <- a * B
  if (B_tot - a <= 0L) {
    return(rep.int(1L, a))
  }
  
  w <- choose(group_sizes, 6)
  B_rem <- B_tot - a
  
  alloc_raw <- B_rem * w / sum(w)
  alloc <- floor(alloc_raw) + 1L
  r <- B_tot - sum(alloc)
  
  if (r > 0L) {
    priority <- head(order(-(alloc_raw %% 1)), r)
    alloc[priority] <- alloc[priority] + 1L
  }
  as.integer(alloc)
}
