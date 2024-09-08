#' trace estimator B_0
#'
#' @param X data matrix with subjects in cols and factor levels in rows multiplied with the hypothesis matrix
#' @return returns the trace estimator for T*V
#' @keywords internal

B0 <- function(X){
  return(sum(X^2) / ncol(X)) # wenn ich herausfinde, warum das hier geht, dann kann ich auch herasufinden, warum sich B2 umschreiben lässt...
}



#' trace estimator B_2
#'
#' @param X data matrix with subjects in cols and factor levels in rows multiplied with the hypothesis matrix
#' @return returns the trace estimator for (T*V)^2
#' @keywords internal
B2 <- function(X){
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
#' @param X data matrix with subjects in colss and factor levels in rows multiplied with the hypothesis matrix T
#' @return returns the trace estimator for (T*V)^3
#' @keywords internal
B3 <- function(X){
  N <- ncol(X)
  S <- 0
  for (k in 1:(N-2)) {
    for (l in (k+1):(N-1)) {
      for (r in (l+1):N) {
        S <- S + (sum(X[, k] * X[, l]) * sum(X[, l] * X[, r]) * sum(X[, r] * X[, k]))
      }
    }
  }
  return(sum(S))
}
