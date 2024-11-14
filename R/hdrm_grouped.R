#### Was soll eingegeben werden??
## X = matrix mit dimensionen c(d,N) -> Individuen in Spalten, Dimensionen in Zeilen
## group = vector der Laenge N mit der Gruppenzuweisung
## hypothesis = character %in% c("time","group","interaction") oder liste mit Eintraegen TW und TS mit dim(TW) = c(a,a) und dim(TS) = c(d,d)
#' @keywords internal
hdrm_grouped_internal <- function(data, group, hypothesis = c("whole","sub","interaction"), B, subsampling){

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

  A1 <- A3 <- numeric(a)
  A2 <- matrix(0, a, a)
  C5 <- numeric(1)

  # A1 und A3
  for (i in 1:a) {
    if(subsampling){
      A1[i] <- A1star_cpp(X = X_TS[, group == i], B)
      A3[i] <- A3star_cpp(X = X_TS[, group == i], B)
    }else{
      A1[i] <- A1_cpp(mat = X_TS[, group == i])
      A3[i] <- A3_cpp(mat = X_TS[, group == i], Part6 = sum(rowMeans(X_TS[, group == i])^2))
    }

  }

  # A2
  for (i in 1:(a-1)) {
    for(r in (i+1):a){
      if(subsampling){
        A2[i,r] <- A2star_cpp(X = X_TS[, group == i], Y = X_TS[, group == r], B)
      }else{
        A2[i,r] <- A2(X = X_TS[, group == i], Y = X_TS[, group == r])
      }
    }
  }

  # C5
  C5 <- C5star_cpp(X = data, group = group, TW = TW, TS = TS, B = B)


  # A4
  temp1 <- temp2 <- 0
  for (i in 1:a) {
    temp1 <- temp1 + ((N/n[i])^2 * TW[i,i]^2 * A3[i])
  }
  for (i in 1:(a-1)) {
    for(r in (i+1):a){
      temp2 = temp2 + ( (N^2 / (n[i]*n[r])) * TW[i,r]^2 * A2[i,r])
    }
  }
  A4 <- temp1 + 2*temp2

  ### Teststatistik
  X_bar <- numeric(a*d)
  for (i in 1:a) {
    # hier rowMeans, da data nicht transponiert ist und wir den mean ueber die
    # Zeit haben wollen (-> wenn sich heraus stellt, dass wir mean ueber das
    # individuum haben wollen, dann colMeans)
    X_bar[1:d + ((i-1)*d)] <- rowMeans(data[, group == i])
  }
  # Erwartungswert, Varianz, Qn und W berechnen
  EW <- sum((N/n) * diag(TW) * A1)
  Var <- 2*A4
  QN <- N * sum(X_bar *(TM %*% X_bar))
  W <- as.numeric((QN - EW) / sqrt(Var))

  # f und p-Wert
  f <- max(1, as.numeric(A4^3 / C5^2))
  p.value <- max(1 - stats::pchisq(W * sqrt(2 * f) + f, df = f), .Machine$double.eps)


  ## Ausgabe
  L <- list(f = f,
            statisitc = W,
            tau = 1/f,
            H = H,
            # gibt bezeichnung der hypothse aus
            hypothesis = ifelse(is.character(hypothesis), hypothesis[1], "custom"),
            p.value = p.value,
            dim = list(d = d, N = N),
            groups = list(a = a,
                          table = table(group))
  )
  class(L) <- c("hdrm")
  return(L)
}
