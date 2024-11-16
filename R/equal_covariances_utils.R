
## A1_eq - function that calculates the estimator A1 form the equal covariances paper.
##
## Input: a numeric matrix, where cols represent the individuals of one group
## and rows represent dimensions
##
## Output: a numeric vector that gives the group estimation for A1
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

## make_A1_eq - function that combines the calculators A1_eq for all groups.
##
## Input:
##      X - a numeric matrix, where cols represent the individuals of one group
##          and rows represent dimensions
##  group - an integer vector, specifies the group membership
##
## Output: numeric; gives the estimator A1 for equal covariances
make_A1_eq <- function(X, group, version = "cpp"){
  ns <- unname(table(group))
  a <- length(ns)
  out <- 0.0
  prefactor <- 0
    for(i in 1:a){
      n <- sum
      out <- out + A1_eq_cpp(X[, group == i])
      prefactor <- prefactor + (ns[i]* (ns[i] - 1))
    }
  return(out/(prefactor))
}


## A2_eq - function that calculates the estimator A2 form the equal covariances paper.
##
## Input: a numeric matrix, where cols represent the individuals of one group
## and rows represent dimensions
##
## Output: a numeric vector that gives the group estimation for A2
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

## make_A2_eq - function that combines the calculators A2_eq for all groups.
##
## Input:
##      X - a numeric matrix, where cols represent the individuals of one group
##          and rows represent dimensions
##  group - an integer vector, specifies the group membership
##
## Output: numeric; gives the estimator A2 for equal covariances
make_A2_eq <- function(X, group){
  ns <- unname(table(group))
  a <- length(ns)
  out <- 0.0
  prefactor <- 0
    for(i in 1:a){
      n <- sum
      out <- out + A2_eq_cpp(X[, group == i])
      prefactor <- prefactor + (24* choose(ns[i], 4))
    }
  return(out/(prefactor))
}


## C1_eq - function that calculates the estimator C1 form the equal covariances paper.
##
## Input:
##      X - a numeric matrix, where cols represent the individuals of one group
##          and rows represent dimensions
##
## Output: a numeric vector that gives the group estimation
C1_eq <- function(X){
  n <- ncol(X)
  S <- 0.0
  for(l1 in 1:n){
    for(l2 in 1:n){
      if(l1 != l2){
        for(l3 in 1:n){
          if(l2 != l3){
            for(l4 in 1:n){
              if(l3 != l4){
                for(l5 in 1:n){
                  if(l5 != l4){
                    for(l6 in 1:n){
                      if(length(unique(c(l1,l2,l3,l4,l5,l6))) == 6){
                        S <- S + crossprod(X[, l1] - X[, l2], X[, l3] - X[, l4]) *
                          crossprod(X[, l3] - X[, l4], X[, l5] - X[, l6]) *
                          crossprod(X[, l5] - X[, l6], X[, l1] - X[, l2])
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(as.numeric(S))
}

## make_C1_eq - function that combines the calculations C1_eq for all groups.
##
## Input:
##      X - a numeric matrix, where cols represent the individuals of one group
##          and rows represent dimensions
##  group - an integer vector, specifies the group membership
##
## Output: numeric; gives the estimator A2 for equal covariances
make_C1_eq <- function(X, group){
  ns <- unname(table(group))
  a <- length(ns)
  out <- 0.0
  prefactor <- 0
  for(i in 1:a){
    n <- sum
    out <- out + C1_eq_cpp(X[, group == i])
    prefactor <- prefactor + (choose(ns[i], 6))
  }
  return(out/(prefactor * 5760)) ## 5760 = 6! * 8
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

## make_C1_star_eq - function that combines the calculations C1_star_eq for all groups.
##
## Input:
##      X - a numeric matrix, where cols represent the individuals of one group
##          and rows represent dimensions
##      B - an integer, determines the number of subsamplings
##  group - an integer vector, specifies the group membership
##
## Output: numeric; gives the estimator A2 for equal covariances
make_C1_star_eq <- function(X, B, group){
  ns <- unname(table(group))
  a <- length(ns)
  out <- 0.0
  prefactor <- 0
  for(i in 1:a){
    n <- sum
    out <- out + C1_star_eq(X[, group == i], B)
  }
  return(out/(8 * a * B)) ## 8 * a*B considering B is equal for all groups
}
