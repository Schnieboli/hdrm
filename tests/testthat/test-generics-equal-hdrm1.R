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
Erg_matrix_flat <- hdrm1(M, hypothesis = "flat")
Erg_data.frame_flat <- hdrm1(df, formula = value ~ subject + factor, hypothesis = "flat")
Erg_matrix_equal <- hdrm1(M, hypothesis = "flat")
Erg_data.frame_equal <- hdrm1(df, formula = value ~ subject + factor, hypothesis = "flat")



# Tests -------------------------------------------------------------------

test_that("equal data",{
  expect_equal(Erg_data.frame_flat$data, Erg_matrix_flat$data)
  expect_equal(Erg_data.frame_equal$data, Erg_matrix_equal$data)
})

test_that("X in X out data.frame",{
  expect_true(all.equal(Erg_data.frame_flat$data, M))
  expect_true(all.equal(Erg_data.frame_equal$data, M))

})

test_that("X in X out matrix",{
  expect_true(all.equal(Erg_matrix_flat$data, M))
  expect_true(all.equal(Erg_matrix_equal$data, M))
})

test_that("equal W",{
  expect_equal(Erg_data.frame_flat$statisitc, Erg_matrix_flat$statisitc)
  expect_equal(Erg_data.frame_equal$statisitc, Erg_matrix_equal$statisitc)
})

# wird halt optimiert, desegen kann es sein, dass das hier mal failed
test_that("equal p.value",{
  expect_equal(Erg_data.frame_flat$p.value, Erg_matrix_flat$p.value)
  expect_equal(Erg_data.frame_equal$p.value, Erg_matrix_equal$p.value)
})


