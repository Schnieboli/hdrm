n1 <- 31
n2 <- 35
d <- 40
M1 <- matrix(rnorm(n1*d),d,n1)
M2 <- matrix(rnorm(n2*d),d,n2)
P1 <- diag(n1) - matrix(1/n1,n1,n1)
P2 <- diag(n2) - matrix(1/n2,n2,n2)
PM1 <- P1%*%t(M1)
PM2 <- M2 %*% P2

# A1 ----------------------------------------------------------------------

Erg1_R <- A1(M1)
Erg1_cpp <- A1_cpp(t(M1))
test_that("A1 == A1_cpp",{
  expect_equal(Erg1_cpp, Erg1_R)
})

# A2 ----------------------------------------------------------------------



Erg1 <- sum((t(M1) %*% M2)^2)



S <- 0
for (i in 1:n1) {
  for (j in 1:n2) {
    temp <- 0
    for (k in 1:d) {
      temp <- temp + (M1[k,i] * M2[k,j])
    }
    S <- S + temp^2
  }
}
S
Erg1

Erg2_R <- A2(t(M1),t(M2))
Erg2_cpp <- A2_cpp(PM1, t(PM2))

test_that("A2==A2_cpp",{
  expect_equal(Erg2_R, Erg2_cpp)
})
