library(microbenchmark)
library(Rcpp)

source("R/Spurenschaetzer_Sattler.R")
sourceCpp("src/trace_estimators_Sattler.cpp")
sourceCpp("src/A2_Sattler.cpp")

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
               unit = "relative"
               )

# d = N
microbenchmark(R = A1(M2),
               cpp = A1_cpp(M2),
               unit = "relative"
               )

# d > N
microbenchmark(R = A1(M3),
               cpp = A1_cpp(M3),
               unit = "relative"
               )


# d >> N
microbenchmark(R = A1(M4),
               cpp = A1_cpp(M4),
               unit = "relative"
               )


# d >>> N
microbenchmark(R = A1(M5),
               cpp = A1_cpp(M5),
               unit = "relative"
               )
# Es ergibt Sinn, dass für größere d R cpp aufholt, da


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


microbenchmark(R = A3(M1),
               cpp = A3_cpp(M1, sum(rowMeans(M1)^2)),
               unit = "relative",
               times = 5
               )

microbenchmark(R = A3(M2),
               cpp = A3_cpp(M2, sum(rowMeans(M2)^2)),
               unit = "relative",
               times = 5
)

microbenchmark(R = A3(M3),
               cpp = A3_cpp(M3, sum(rowMeans(M3)^2)),
               unit = "relative",
               times = 5
)

microbenchmark(R = A3(M4),
               cpp = A3_cpp(M4, sum(rowMeans(M4)^2)),
               unit = "relative",
               times = 5
)

microbenchmark(R = A3(M5),
               cpp = A3_cpp(M5, sum(rowMeans(M5)^2)),
               unit = "relative",
               times = 5
)

# A2 ----------------------------------------------------------------------

## quasi gleich schnell

n1 <- c(70, 100, 130, 200)
n2 <- c(50, 90, 150, 250)
d <- 365
M11 <- matrix(rnorm(n1[1]*d),d,n1[1])
M12 <- matrix(rnorm(n2[1]*d),d,n2[1])
M21 <- matrix(rnorm(n1[2]*d),d,n1[2])
M22 <- matrix(rnorm(n2[2]*d),d,n2[2])
M31 <- matrix(rnorm(n1[3]*d),d,n1[3])
M32 <- matrix(rnorm(n2[3]*d),d,n2[3])
M41 <- matrix(rnorm(n1[4]*d),d,n1[4])
M42 <- matrix(rnorm(n2[4]*d),d,n2[4])

microbenchmark(R = A2(M11,M12),
               cpp = {
                 P1 <- diag(n1[1]) - matrix(1/n1[1],n1[1],n1[1])
                 P2 <- diag(n2[1]) - matrix(1/n2[1],n2[1],n2[1])
                 sum(A2_cpp(M11, M12, P1, P2))/((n1[1]-1)*(n2[1]-1))
                 },
               unit = "relative",
               times = 10
               )

microbenchmark(R = A2(M21,M22),
               cpp = {
                 P1 <- diag(n1[2]) - matrix(1/n1[2],n1[2],n1[2])
                 P2 <- diag(n2[2]) - matrix(1/n2[2],n2[2],n2[2])
                 sum(A2_cpp(M21, M22, P1, P2))/((n1[2]-1)*(n2[2]-1))
               },
               unit = "relative",
               times = 10
)

microbenchmark(R = A2(M31,M32),
               cpp = {
                 P1 <- diag(n1[3]) - matrix(1/n1[3],n1[3],n1[3])
                 P2 <- diag(n2[3]) - matrix(1/n2[3],n2[3],n2[3])
                 sum(A2_cpp(M31, M32, P1, P2))/((n1[3]-1)*(n2[3]-1))
               },
               unit = "relative",
               times = 10
)

microbenchmark(R = A2(M41,M42),
               cpp = {
                 P1 <- diag(n1[4]) - matrix(1/n1[4],n1[4],n1[4])
                 P2 <- diag(n2[4]) - matrix(1/n2[4],n2[4],n2[4])
                 sum(A2_cpp(M41, M42, P1, P2))/((n1[4]-1)*(n2[4]-1))
               },
               unit = "relative",
               times = 10
)

# C5 ----------------------------------------------------------------------

##### Versuch Simulation -> es fehlt noch eine Funktion, die nur den Median jeweils ausgibt...

# library(microbenchmark)
# N <- c(144, 288, 576, 1152)
# d <- c(100, 200, 300, 600)
# mat <- matrix(rnorm(max(N)*max(d)), max(d), max(N))
# B <- 10
#
# for(i in N){
#   a <- c(4, 6, 8, 12, 24, 36)
#   group1 <- rep(1:a[1], i/a[1])
#   group2 <- rep(1:a[2], i/a[2])
#   group3 <- rep(1:a[3], i/a[3])
#   group4 <- rep(1:a[4], i/a[4])
#   group5 <- rep(1:a[5], a[5])
#   for (j in d) i/{
#     mat1 <- mat[1:j, 1:i]
#     microbenchmark(R = C5star(mat1, group = group1, TM = diag(a[1]*j), B = B),
#                    cpp = C5star_cpp(mat1, group = group1, B = B, TM = diag(a[1]*j), n = as.integer(table(group1))),
#                    unit = "relative",
#                    times = 10)
#
#     microbenchmark(R = C5star(mat1, group = group2, TM = diag(a[2]*j), B = B),
#                    cpp = C5star_cpp(mat1, group = group2, B = B, TM = diag(a[2]*j), n = as.integer(table(group2))),
#                    unit = "relative",
#                    times = 10)
#
#     microbenchmark(R = C5star(mat1, group = group3, TM = diag(a[3]*j), B = B),
#                    cpp = C5star_cpp(mat1, group = group3, B = B, TM = diag(a[3]*j), n = as.integer(table(group3))),
#                    unit = "relative",
#                    times = 10)
#
#     microbenchmark(R = C5star(mat1, group = group4, TM = diag(a[4]*j), B = B),
#                    cpp = C5star_cpp(mat1, group = group4, B = B, TM = diag(a[4]*j), n = as.integer(table(group4))),
#                    unit = "relative",
#                    times = 10)
#
#     microbenchmark(R = C5star(mat1, group = group5, TM = diag(a[5]*j), B = B),
#                    cpp = C5star_cpp(mat1, group = group5, B = B, TM = diag(a[5]*j), n = as.integer(table(group5))),
#                    unit = "relative",
#                    times = 10)
#
#   }
# }


N <- 144
d <- 50
mat <- matrix(rnorm(N*d), d, N)

group1 <- rep(1:24, 6)
group2 <- rep(1:12, 12)
group3 <- rep(1:8, 18)
group4 <- rep(1:4, 36)
group5 <- rep(1:3, 48)
group6 <- rep(1:2, 72)
a1 <- 24
a2 <- 12
a3 <- 8
a4 <- 4
a5 <- 3
a6 <- 2

B <- 100000

##### für a -> unendlich wird R besser unabhängig von B
microbenchmark(R = C5star(mat, group = group1, TM = diag(a1*d), B = B),
               cpp = C5star_cpp(mat, group = group1, B = B, TM = diag(a1*d)),
               unit = "relative",
               times = 10)

microbenchmark(R = C5star(mat, group = group2, TM = diag(a2*d), B = B),
               cpp = C5star_cpp(mat, group = group2, B = B, TM = diag(a2*d)),
               unit = "relative",
               times = 10)

microbenchmark(R = C5star(mat, group = group3, TM = diag(a3*d), B = B),
               cpp = C5star_cpp(mat, group = group3, B = B, TM = diag(a3*d)),
               unit = "relative",
               times = 10)

microbenchmark(R = C5star(mat, group = group4, TM = diag(a4*d), B = B),
               cpp = C5star_cpp(mat, group = group4, B = B, TM = diag(a4*d)),
               unit = "relative",
               times = 10)

microbenchmark(R = C5star(mat, group = group5, TM = diag(a5*d), B = B),
               cpp = C5star_cpp(mat, group = group5, B = B, TM = diag(a5*d)),
               unit = "relative",
               times = 10)

microbenchmark(R = C5star(mat, group = group6, TM = diag(a6*d), B = B),
               cpp = C5star_cpp(mat, group = group6, B = B, TM = diag(a6*d)),
               unit = "relative",
               times = 5)
