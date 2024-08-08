library(microbenchmark)

source("R/trace_estimators_Pauly.R")
Rcpp::sourceCpp("src/trace_estimators_Pauly.cpp")



# data --------------------------------------------------------------------

N <- 50
d <- c(30, 50, 70, 100, 130)

M1 <- matrix(rnorm(N*d[1]), d[1], N)
M2 <- matrix(rnorm(N*d[2]), d[2], N)
M3 <- matrix(rnorm(N*d[3]), d[3], N)
M4 <- matrix(rnorm(N*d[4]), d[4], N)
M5 <- matrix(rnorm(N*d[5]), d[5], N)


# B2 ----------------------------------------------------------------------

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

## d < N
microbenchmark(R_int = B2_intuitiv(M1),
               R_opt = B2(M1),
               cpp = B2_cpp(M1),
               unit = "relative",
               check = "equal")
## d = N
microbenchmark(R_int = B2_intuitiv(M2),
               R_opt = B2(M2),
               cpp = B2_cpp(M2),
               unit = "relative",
               check = "equal")

## d > N
microbenchmark(R_int = B2_intuitiv(M3),
               R_opt = B2(M3),
               cpp = B2_cpp(M3),
               unit = "relative",
               check = "equal")



# B3 ----------------------------------------------------------------------

B3_intuitiv <- function(X){
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

## d < N
microbenchmark(
  #R_int = B3_intuitiv(M1),
  R_opt = B3(M1),
  cpp = B3_cpp(M1),
  unit = "relative",
  check = "equal",
  times = 50
)

## d = N
microbenchmark(
  #R_int = B3_intuitiv(M2),
  R_opt = B3(M2),
  cpp = B3_cpp(M2),
  unit = "relative",
  check = "equal",
  times = 50
)

## d > N
microbenchmark(
  #R_int = B3_intuitiv(M3),
  R_opt = B3(M3),
  cpp = B3_cpp(M3),
  unit = "relative",
  check = "equal",
  times = 50
)

## d >> N
microbenchmark(
  #R_int = B3_intuitiv(M4),
  R_opt = B3(M4),
  cpp = B3_cpp(M4),
  unit = "relative",
  check = "equal",
  times = 50
)

## d >>> N
microbenchmark(
  #R_int = B3_intuitiv(M4),
  R_opt = B3(M5),
  cpp = B3_cpp(M5),
  unit = "relative",
  check = "equal",
  times = 50
)

## hier scheint es so, als würde die R Funktion deutlich aufholen

# ### Idee für Simulationsstudie:
#
# N <- c(50, 100, 150, 200, 300, 400, 500)
# M <- matrix(rnorm(max(N)*1.5*max(N)), 1.5*max(N), max(N))
# Erg <- matrix(0, 7, 4)
#
# for(i in 1:length(N)){
#   print(i)
#   d <- floor(N[i]*c(1,1.2,1.4,1.5))
#   for (j in 1:length(d)) {
#     mat <- M[1:d[j], 1:N[i]]
#     temp <- microbenchmark(
#       R_opt = B3(mat),
#       cpp = B3_cpp(mat),
#       times = 1
#     )
#     Erg[i,j] <- median(temp$time[temp$expr == "cpp"]) > median(temp$time[temp$expr == "R_opt"])
#
#   }
#
# }

