# library(testthat)
# load_all()
#
# ## This file tests, whether the new and faster trace estimators produce the same
# ##  result for a randomly generated Matrix X
#
#
# B0_intuitiv <- function(X){
#   S <- 0
#   N <- dim(X)[2]
#   for (i in 1:N) {
#     S <- S + sum(X[,i]^2)
#   }
#   return(S)
# }
#
B2_intuitiv <- function(X){
  S <- 0
  N <- dim(X)[2]
  for (i in 1:N) {
    for(j in 1:N){
      if(i != j) S <- S + sum(X[,i]*X[,j])^2
    }
  }
  return(S)
}
#
# B3_intuitiv <- function(X){
#   N <- ncol(X)
#   S <- 0
#   #browser()
#   for (k in 1:(N-2)) {
#     for (l in (k+1):(N-1)) {
#       for (r in (l+1):N) {
#         S <- S + (sum(X[, k] * X[, l]) * sum(X[, l] * X[, r]) * sum(X[, r] * X[, k]))
#       }
#     }
#   }
#   return(sum(S))
# }
#
#
# # Tests -------------------------------------------------------------------
#
# set.seed(1)
# N <- 30
# d <- 40
# X <- matrix(rnorm(30*40), d, N)

#
# test_that("B0",{
#   expect_equal(B0_intuitiv(X), B0(X))
# })
#
# test_that("B2",{
#   expect_equal(B2_intuitiv(X), B2(X))
# })
#
# test_that("B3",{
#   expect_equal(B3_intuitiv(X), B3(X))
# })
