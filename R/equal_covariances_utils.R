
## A1_eq - function that calculates the estimator A1 form the equal covariances paper.
##
## Input: a numeric matrix, where cols represent the individuals of one group
## and rows represent dimensions
##
## Output: a numeric vector that gives the group estimation
A1_eq <- function(X){
  n <- ncol(X)
  S <- 0.0
  for(j in 1:(n-1)){
    for(k in j:n){
      S <- S + crossprod(X[, j], X[, k])
    }
  }
  return(as.numeric(S))
}

## A2_eq - function that calculates the estimator A2 form the equal covariances paper.
##
## Input: a numeric matrix, where cols represent the individuals of one group
## and rows represent dimensions
##
## Output: a numeric vector that gives the group estimation
A2_eq <- function(X){
  n <- ncol(X)
  S <- 0.0
  for (l2 in 1:(n-1)) {
    for(l1 in l2:n){
      for(k2 in 1:(n-1)){
        # only call the loop if this condition is true
        if(k2 != l1 && k2!= l2){
          for(k1 in k2:n){
        # only add to total if this condition is true
            if(k1 != l1 && k1 != l2){
              S <- S + (crossprod(X[, l1] - X[, l2], X[, k1] - X[, k2]))^2
            }
          }
        }
      }
    }
  }

  return(as.numeric(S))
}
