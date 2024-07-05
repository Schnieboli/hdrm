#' trace estimator B_0
#'
#' @param X data matrix with individuals in rows and levels in cols multiplied with the hypothesis matrix
#' @param N number of individuals
#' @return returns the trace estimator for T*V
#' @keywords internal

B0 <- function(X, N){
  return(sum(rowSums(X^2))/N)
}



#' trace estimator B_2
#'
#' @param X data matrix with individuals in rows and levels in cols multiplied with the hypothesis matrix
#' @param N number of individuals
#' @return returns the trace estimator for (T*V)^2
#' @keywords internal
B2 <- function(X, N){
  S <- (rowSums(X^2))^2
  M <- sum(X[1,]^2)^2
  for (i in 2:N) {
    Y <-  X * X[c(i:N, 1:(i-1)),]
    S <- S + (rowSums(Y))^2
    M <- M + sum(X[i,]^2)^2
  }
  return((sum(S) - M)/(N*(N-1)))
}

#' trace estimator B_3
#'
#' @param X data matrix with individuals in rows and levels in cols multiplied with the hypothesis matrix
#' @param N number of individuals
#' @return returns the trace estimator for (T*V)^3
#' @keywords internal
B3 <- function(X, N){
  S <- 0
  for (k in 1:(N-2)) {
    for (l in (k+1):(N-1)) {
      for (r in (l+1):N) {
        S <- S + (sum(X[k, ] * X[l, ]) * sum(X[l, ] * X[r, ]) * sum(X[r, ] * X[k, ]))
      }
    }
  }
  return(S/choose(N,3))
}
