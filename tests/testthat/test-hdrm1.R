library(testthat)
devtools::load_all()

# Dieses Skript erzeugt eine Liste, eine Matrix und einen Data Frame, die das identische Ergebnis ausgeben sollten.
# Anschließend wird getestet, ob die gleichen Ergebnisse rauskommen

N <- 30
d <- 40


### Matrixobjekt erstellen
set.seed(1)
M <- matrix(rnorm(N*d), d, N)


# data.frame erstellen
set.seed(1)
df <- data.frame(value = rnorm(d*N))
df$subject <- rep(1:N, each = d)
df$factor <- rep(1:d, N)




#Erg_list <- hdrm1(data = L, hypothesis = "flat")
Erg_matrix <- hdrm1(M, hypothesis = "flat")
Erg_data.frame <- hdrm1(df, formula = value ~ subject + factor, hypothesis = "flat")

#Erg_list$p.value
Erg_matrix$p.value
Erg_data.frame$p.value

#Erg_list$statisitc
Erg_matrix$statisitc
Erg_data.frame$statisitc


# Tests -------------------------------------------------------------------

test_that("equal data",{
  expect_equal(Erg_data.frame$data, Erg_matrix$data)
})


test_that("equal W",{
  expect_equal(Erg_data.frame$statisitc, Erg_matrix$statisitc)
})

# wird halt optimiert, desegen kann es sein, dass das hier mal failed
test_that("equal p.value",{
  expect_equal(Erg_data.frame$p.value, Erg_matrix$p.value



Erg_matrix <- hdrm1(M, hypothesis = "equal")
Erg_data.frame <- hdrm1(df, formula = value ~ subject + factor, hypothesis = "equal")

test_that("equal W",{
  expect_equal(Erg_data.frame$statisitc, Erg_matrix$statisitc)
})


