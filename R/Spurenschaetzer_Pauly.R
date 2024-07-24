#' trace estimator B_0
#'
#' @param X data matrix with individuals in rows and levels in cols multiplied with the hypothesis matrix
#' @return returns the trace estimator for T*V
#' @keywords internal

B0 <- function(X){
  return(sum(X^2)) # wenn ich herausfinde, warum das hier geht, dann kann ich auch herasufinden, warum sich B2 umschreiben lässt...
}



#' trace estimator B_2
#'
#' @param X data matrix with individuals in rows and levels in cols multiplied with the hypothesis matrix
#' @return returns the trace estimator for (T*V)^2
#' @keywords internal
B2 <- function(X){
  #browser()
  N <- ncol(X)
  S <- (colSums(X^2))^2
  M <- sum(X[,1]^2)^2
  for (i in 2:N) {
    Y <-  X[, c(i:N, 1:(i-1))] * X
    S <- S + (colSums(Y))^2
    M <- M + sum(X[, i]^2)^2
  }
  return(sum(S) - M)
}


#' trace estimator B_3
#'
#' @param X data matrix with individuals in rows and levels in cols multiplied with the hypothesis matrix
#' @return returns the trace estimator for (T*V)^3
#' @keywords internal
B3 <- function(X){
  N <- ncol(X)
  S <- 0
  #browser()
  for (k in 1:(N-2)) {
    for (l in (k+1):(N-1)) {
      for (r in (l+1):N) {
        S <- S + (sum(X[, k] * X[, l]) * sum(X[, l] * X[, r]) * sum(X[, r] * X[, k]))
      }
    }
  }
  return(sum(S))
}

# B3_experimentell <- function(X){
#   N <- ncol(X)
#   S <- 0
#   for (k in 1:N) {
#     for (l in 1:N) {
#       for (r in 1:N) {
#         S <- S + (sum(X[, k] * X[, l]) * sum(X[, l] * X[, r]) * sum(X[, r] * X[, k]))
#       }
#     }
#   }
#   return(S)
# }
# ### Stand jetzt produziert B3_ON2 das gleiche Ergebnis wie B3_experimentell.
# ### Beide sind versch Varianten, für alle j,k,r B3 zu berechnen jetzt fehlt
# ### noch, das was zu viel ist abziehen ...
# ### am besten als erstes die elemente identifizieren, die zu viel multipliziert werden...
# ### Vielleicht klappt das auch gar nicht???
# B3_ON2 <- function(X){
# N <- ncol(X)
# S1 <- S2 <- numeric(1)
# E <- numeric(1)
# for (i in 1:N) {
#   for(j in 1:N){
#     S1 <- S1 + sum(colSums(X * X[, c(i:N, 1:(i-1))[1:N]]) * colSums(X[, c(i:N, 1:(i-1))[1:N]] * X[, c(j:N, 1:(j-1))[1:N]]) * colSums(X[, c(j:N, 1:(j-1))[1:N]] * X))
#
#   }
# }
# # for (i in 1:N) {
# #   for(j in 1:N){
# #   M_i <- X[, c(i:N, 1:(i-1))[1:N]]
# #   M_j <- X[, c(j:N, 1:(j-1))[1:N]]
# #     S2 <- S2 + sum(colSums(X * M_i) * colSums(M_i * M_j) * colSums(M_j * X))
# #   }
# # }
# return(sum(S1))
# }
