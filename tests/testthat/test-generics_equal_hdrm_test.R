# Dieses Skript erzeugt eine Liste, eine Matrix und einen Data Frame, die das identische Ergebnis ausgeben sollten.
# Anschließend wird getestet, ob die gleichen Ergebnisse rauskommen

N <- 30
d <- 40
a <- 3


### Matrixobjekt erstellen
set.seed(1)
M <- matrix(rnorm(N*d, mean = c(rep(0,800), rep(1:d, 10))), d, N)
group <- rep(1:3, each = 10)


# data.frame erstellen
set.seed(1)
df <- data.frame(value = rnorm(N*d, mean = c(rep(0,800), rep(1:d, 10))))
df$group <- rep(1:3, each = d*N/3)
df$subject <- rep(1:N, each = d)
df$factor <- rep(1:d, N)


## globale Hypothesen
set.seed(1)
Erg_matrix_whole <- hdrm_test(M, hypothesis = "whole", group = group)
set.seed(1)
Erg_data.frame_whole <- hdrm_test(df, hypothesis = "whole", "group", "value", "subject","factor")
set.seed(1)
Erg_matrix_sub <- hdrm_test(M, hypothesis = "sub", group = group)
set.seed(1)
Erg_data.frame_sub <- hdrm_test(df, hypothesis = "sub", "group", "value", "subject","factor")

test_that("globale hypothesen gleich",{
  expect_equal(Erg_data.frame_whole$statisitc, Erg_matrix_whole$statisitc)
  expect_equal(Erg_data.frame_sub$statisitc, Erg_matrix_sub$statisitc)
  expect_equal(Erg_data.frame_whole$f, Erg_matrix_whole$f)
  expect_equal(Erg_data.frame_whole$f, Erg_matrix_whole$f)
})

### lokale Hypothesen: gruppen unterschiedlich
test_that("Test auf Gruppenunterschied gleich",{

  e1 <- c(1,0,0)
  e2 <- c(0,1,0)
  e3 <- c(0,0,1)
  P_d <- diag(d) - matrix(1/d,d,d)

  hypothesis <- list(TW = kronecker(e1 %*% t(e1), 1), TS = P_d)
  expect_equal({set.seed(1); hdrm_test(M, hypothesis = hypothesis, group)$p.value},
               {set.seed(1); hdrm_test(df, hypothesis = hypothesis, "group","value","subject","factor")$p.value}
               )

  hypothesis <- list(TW = kronecker(e2 %*% t(e2), 1), TS = P_d)
  expect_equal({set.seed(1); hdrm_test(M, hypothesis = hypothesis, group)$p.value},
               {set.seed(1); hdrm_test(df, hypothesis = hypothesis, "group","value","subject","factor")$p.value}
  )

  hypothesis <- list(TW = kronecker(e3 %*% t(e3), 1), TS = P_d)
  expect_equal({set.seed(1); hdrm_test(M, hypothesis = hypothesis, group)$p.value},
               {set.seed(1); hdrm_test(df, hypothesis = hypothesis, "group","value","subject","factor")$p.value}
  )

})


### lokale Hypothesen: zeitprofil unterschiedlich
# annahme: d ist kreuzung aus zwei faktoren mit 4x25
test_that("lokales Zeitprofil unterschiedlich",{
  a <- 3
  P_a <- diag(a) - matrix(1/a,a,a)
  P_10 <- diag(10) - matrix(1/10,10,10)
  e1 <- c(1,0,0,0)
  e2 <- c(0,1,0,0)
  e3 <- c(0,0,1,0)
  e4 <- c(0,0,0,1)

  hypothesis = list(TW = P_a, TS = kronecker(e1 %*% t(e1), P_10))
  expect_equal({set.seed(1); hdrm_test(M, hypothesis = hypothesis, group)$p.value},
               {set.seed(1); hdrm_test(df, hypothesis = hypothesis, "group","value","subject","factor")$p.value}
  )

  hypothesis = list(TW = P_a, TS = kronecker(e2 %*% t(e2), P_10))
  expect_equal({set.seed(1); hdrm_test(M, hypothesis = hypothesis, group)$p.value},
               {set.seed(1); hdrm_test(df, hypothesis = hypothesis, "group","value","subject","factor")$p.value}
  )

  hypothesis = list(TW = P_a, TS = kronecker(e2 %*% t(e2), P_10))
  expect_equal({set.seed(1); hdrm_test(M, hypothesis = hypothesis, group)$p.value},
               {set.seed(1); hdrm_test(df, hypothesis = hypothesis, "group","value","subject","factor")$p.value}
  )

  hypothesis = list(TW = P_a, TS = kronecker(e2 %*% t(e2), P_10))
  expect_equal({set.seed(1); hdrm_test(M, hypothesis = hypothesis, group)$p.value},
               {set.seed(1); hdrm_test(df, hypothesis = hypothesis, "group","value","subject","factor")$p.value}
  )
})
