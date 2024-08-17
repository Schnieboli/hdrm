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
M6 <- matrix(rnorm(100*d[5]), d[5], 100)


# A1 ----------------------------------------------------------------------

### wenn man sehr große Matrizen betrachtet, dann kommt R bis auf 14x an cpp ran -> das muss auch auf jeden Fall in die Arbeit

# d < N
microbenchmark(R = A1(M1),
               cpp = A1_cpp(M1),
               unit = "relative",
               check = "equal")

# d = N
microbenchmark(R = A1(M2),
               cpp = A1_cpp(M2),
               unit = "relative",
               check = "equal")

# d > N
microbenchmark(R = A1(M3),
               cpp = A1_cpp(M3),
               unit = "relative",
               check = "equal")


# d >> N
microbenchmark(R = A1(M4),
               cpp = A1_cpp(M4),
               unit = "relative",
               check = "equal")


# d >>> N
microbenchmark(R = A1(M5),
               cpp = A1_cpp(M5),
               unit = "relative",
               check = "equal")


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


microbenchmark(R = A3(M6),
               cpp = A3_cpp(M6, sum(rowMeans(M6)^2)),
               unit = "relative",
               check = "equal",
               times = 5)


# A2 ----------------------------------------------------------------------

## quasi gleich schnell

n1 <- 580
n2 <- 200
d <- 365
M1 <- matrix(rnorm(n1*d),d,n1)
M2 <- matrix(rnorm(n2*d),d,n2)
# P1 <- diag(n1) - matrix(1/n1,n1,n1)
# P2 <- diag(n2) - matrix(1/n2,n2,n2)
# PM1 <- P1%*%t(M1)
# PM2 <- M2 %*% P2

microbenchmark(R = A2(M1,M2),
               cpp = {
                 P1 <- diag(n1) - matrix(1/n1,n1,n1)
                 P2 <- diag(n2) - matrix(1/n2,n2,n2)
                 sum(A2_cpp(M1, M2, P1, P2))/((n1-1)*(n2-1))
                 },
               unit = "relative",
               check = "equal"
               )

