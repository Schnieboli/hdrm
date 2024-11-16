
#' @keywords internal
hdrm_grouped_eq_cov_internal <- function(data, group, hypothesis = c("whole","sub","interaction"), B, subsampling){

  # N, n, d, a bestimmen
  N <- ncol(data)
  d <- nrow(data)
  a <- length(table(group))
  n <- as.integer(table(group))

  # Hypothese bestimmen
  H <- get_hypothesis_mult(hypothesis, a, d)
  TW <- H$TW
  TS <- H$TS
  TM <- kronecker(TW, TS)

  # X_TS vorbereiten
  X_TS <- TS %*% data

  # Schaetzer
  A1 <- make_A1_eq(X = X_TS, group = group)
  A2 <- make_A2_eq(X = X_TS, group = group)
  C1 <- make_C1_star_eq(X = X_TS, group = group)

  ### Teststatistik
  X_bar <- numeric(a*d)
  for (i in 1:a) {
  ## data is transposed => rowMeans; colMeans would be more efficient, but its ok because its just a times
    X_bar[1:d + ((i-1)*d)] <- rowMeans(data[, group == i])
  }
  # Erwartungswert, Varianz, Qn und W berechnen
  EW <- sum((N/n) * diag(TW)) * A1
  tmp = 0
  for (i in 1:a) {
    for(r in 1:a){
      tmp = tmp + ( TW[i,r]^2 * (N^2 / (n[i] * n[r])) )
    }
  }
  Var <- 2 * A2 * tmp
  QN <- N * sum(X_bar *(TM %*% X_bar))
  W <- as.numeric((QN - EW) / sqrt(Var))

  # f und p-Wert
  f <- max(1, as.numeric(A4^3 / C5^2))
  p.value <- max(1 - stats::pchisq(W * sqrt(2 * f) + f, df = f), .Machine$double.eps)


  ## Ausgabe
  return(
    list(
      f = f,
      statisitc = W,
      tau = 1 / f,
      H = H,
      # gibt bezeichnung der hypothse aus
      hypothesis = ifelse(is.character(hypothesis), hypothesis[1], "custom"),
      p.value = p.value,
      dim = list(d = d, N = N),
      groups = list(a = a, table = table(group))
    )
  )
}
