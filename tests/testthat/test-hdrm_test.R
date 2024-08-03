library(testthat)
devtools::load_all()

# Dieses Skript erzeugt eine Liste, eine Matrix und einen Data Frame, die das identische Ergebnis ausgeben sollten.
# Anschließend wird getestet, ob die gleichen Ergebnisse rauskommen

## Listenobjekt erstellen
set.seed(1)
a <- 4
N <- 24
d <- 30
L <- list()
for(i in 1:a) L[[i]] = matrix(rnorm(N*d), d, N)


### Matrixobjekt erstellen
set.seed(1)
M <- matrix(NA,d,0)
for (i in 1:a) {
  M <- cbind(M, matrix(rnorm(d*N), d, N))
}

# data.frame erstellen
set.seed(1)
df <- data.frame(value = rnorm(d*N*a))
df$subject <- rep(rep(1:N, each = d), a)
df$whole <- rep(1:a, each = N*d)
df$sub <- rep(1:d, N*a)



Erg_list <- hdrm_test(data = L, hypothesis = "sub")
Erg_matrix <- hdrm_test(M, group = rep(1:a, each = N), hypothesis = "sub")
Erg_data.frame <- hdrm_test(df, hypothesis = "sub")

Erg_list$p.value
Erg_matrix$p.value
Erg_data.frame$p.value

Erg_list$statisitc
Erg_matrix$statisitc
Erg_data.frame$statisitc


# Tests -------------------------------------------------------------------

test_that("equal W",{
  expect_equal(c(Erg_list$statisitc, Erg_data.frame$statisitc), c(Erg_matrix$statisitc, Erg_matrix$statisitc))
})



test_that("equal p.value",{
  expect_equal(c(Erg_list$p.value, Erg_data.frame$p.value), c(Erg_matrix$p.value, Erg_matrix$p.value),
               tolerance = 1e-7)
})


Erg_list <- hdrm_test(data = L, hypothesis = "whole")
Erg_matrix <- hdrm_test(M, group = rep(1:a, each = N), hypothesis = "whole")
Erg_data.frame <- hdrm_test(df, hypothesis = "whole")

test_that("equal W",{
  expect_equal(c(Erg_list$statisitc, Erg_data.frame$statisitc), c(Erg_matrix$statisitc, Erg_matrix$statisitc))
})


Erg_list <- hdrm_test(data = L, hypothesis = "interaction")
Erg_matrix <- hdrm_test(M, group = rep(1:a, each = N), hypothesis = "interaction")
Erg_data.frame <- hdrm_test(df, hypothesis = "interaction")

test_that("equal W",{
  expect_equal(c(Erg_list$statisitc, Erg_data.frame$statisitc), c(Erg_matrix$statisitc, Erg_matrix$statisitc))
})
