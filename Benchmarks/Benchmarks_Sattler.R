library(microbenchmark)
library(Rcpp)

source("R/Spurenschaetzer_Sattler.R")
sourceCpp("src/trace_estimators_Sattler.cpp")


# data --------------------------------------------------------------------

N <- 50
d <- c(30, 50, 70, 100, 130)

M1 <- matrix(rnorm(N*d[1]), d[1], N)
M2 <- matrix(rnorm(N*d[2]), d[2], N)
M3 <- matrix(rnorm(N*d[3]), d[3], N)
M4 <- matrix(rnorm(N*d[4]), d[4], N)
M5 <- matrix(rnorm(N*d[5]), d[5], N)



# A2 ----------------------------------------------------------------------





# A3 ----------------------------------------------------------------------

A3_intuitiv <- function(X){
  nX <- ncol(X)
  S <- 0
  comb <- combn(nX, 2)

  for (l1 in 1:(nX-1)) {
    for(l2 in (l1+1):nX){
      for(k2 in 1:(nX-1)){
        for(k1 in (k2+1):nX){
          if(k2 != l1 & k2 != l2 & k1 != l1 & k1 != l2) S <- S + sum((X[, l1] - X[, l2]) * (X[, k1] - X[, k2]))^2
        }
      }
    }
  }
  return(S / (24 * choose(nX,4)))
}


