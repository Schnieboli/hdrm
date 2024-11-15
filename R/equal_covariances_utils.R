
## A1_eq - function that calculates the estimator A1 form the equal covariances paper.
##
## Input: a numeric matrix, where cols represent the individuals of one group
## and rows represent dimensions
##
## Output: a numeric vector that gives the group estimation
A1_eq <- function(X){
  n <- ncol(X)
  S <- 0.0
  for(l2 in 1:(n-1)){
    for(l1 in (l2+1):n){
      S <- S + sum((X[, l1] -  X[, l2])^2)
    }
  }
  return(S)
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
    for(l1 in (l2+1):n){
      for(k2 in 1:(n-1)){
        # only call the loop if this condition is true
        if(k2 != l1 && k2!= l2){
          for(k1 in (k2+1):n){
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


## C1_eq - function that calculates the estimator C1 form the equal covariances paper.
##
## Input: a numeric matrix, where cols represent the individuals of one group
## and rows represent dimensions
##
## Output: a numeric vector that gives the group estimation
C1_eq <- function(X){
  n <- ncol(X)
  S <- 0.0
  for(l1 in 1:(n-5)){
    for(l2 in (l1+1):(n-4)){
      for(l3 in (l2+1):(n-3)){
        for(l4 in (l3+1):(n-2)){
          for(l5 in (l4+1):(n-1)){
            for(l6 in (l5+1):(n)){
              S <- S + crossprod(X[, l1] - X[, l2], X[, l3] - X[, l4]) *
                crossprod(X[, l3] - X[, l4], X[, l5] - X[, l6]) *
                crossprod(X[, l5] - X[, l6], X[, l1] - X[, l2])
            }
          }
        }
      }
    }
  }
  return(as.numeric(S))
}


## C1_star_eq - function that calculates the estimator C1* form the equal covariances paper.
##
## Input:
##      X - a numeric matrix, where cols represent the individuals of one group
##          and rows represent dimensions
##      B - an integer, determines the number of subsamplings
##
## Output: a numeric vector that gives the group estimation
C1_star_eq <- function(X, B){
  n <- ncol(X)
  S <- 0.0
  for (i in 1:B) {
    b <- sample.int(n = n, size = 6, replace = FALSE)
    S <- S + crossprod(X[, b[1]] - X[, b[2]], X[, b[3]] - X[, b[4]]) *
      crossprod(X[, b[3]] - X[, b[4]], X[, b[5]] - X[, b[6]]) *
      crossprod(X[, b[5]] - X[, b[6]], X[, b[1]] - X[, b[2]])
  }
  return(as.numeric(S))
}
