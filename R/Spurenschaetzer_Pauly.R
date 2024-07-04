#' trace estimator B_0
#'
#' @param X data matrix with individuals in rows and levels in cols multiplied with the hypothesis matrix
#' @param N number of individuals
#' @return returns the trace estimator for T*V
#' @keywords internal

B0 <- function(X, N){
  S <- 0
  for (i in 1:N) {
    S <- S + sum(X[i, ]^2)
  }
  return(S/N)
}

#' trace estimator B_2
#'
#' @param X data matrix with individuals in rows and levels in cols multiplied with the hypothesis matrix
#' @param N number of individuals
#' @return returns the trace estimator for (T*V)^2
#' @keywords internal
B2 <- function(X, N){
  S <- 0
  for (k in 1:N) {
    for(l in 1:N){
      if(k != l) S <- S + (sum(X[k, ] * X[l, ])^2)
    }
  }
  return(S / (N*(N-1)))
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
