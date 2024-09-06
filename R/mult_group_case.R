#### Was soll eingegeben werden??
## data = matrix mit dimensionen c(d,N) -> Individuen in Spalten, Dimensionen in Zeilen
## group = vector der Länge N mit der Gruppenzuweisung
## hypothesis = character %in% c("time","group","interaction") oder liste mit Einträgen TW und TS mit dim(TW) = c(a,a) und dim(TS) = c(d,d)
#' @keywords internal
hdrm_test_internal_bootstrap <- function(data, group, hypothesis = c("whole","sub","interaction"), B){
  # N,d,a bestimmen
  N_with_NA <- ncol(data)
  group <- group[stats::complete.cases(t(data))] # hier mit factor arbeiten

  data <- t(stats::na.omit(t(data))) # data ist nicht transponiert!!!
  N <- ncol(data)
  d <- nrow(data)
  a <- length(table(group))
  n <- as.integer(table(group))
  stopifnot(length(group) == N)
  ### Stopkriterien
  # mind 2 Gruppen
  # mind 6 Individuen
  stopifnot(a >= 2, all(n >= 6))


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

  # spur1 und spur3
  for (i in 1:a) {
    A1[i] <- A1star_cpp(X = X_TS[, group == i], B = B)
    A3[i] <- A3star_cpp(X = X_TS[, group == i], B = B)
  }

  # spur2
  for (i in 1:(a-1)) {
    for(r in (i+1):a){
      A2[i,r] <- A2star_cpp(X = X_TS[, group == i], Y = X_TS[, group == r], B = B)
    }
  }

  # C5
  C5 <- C5star_cpp(X = data, group = group, TW = TW, TS = TS, B = B)

  # E(Q_N)
  EW <- sum((N/n) * diag(TW) * A1)

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
  rm(temp1, temp2)

  # Var(Q_N)
  Var <- 2*A4


  ### Teststatistik
  X_bar <- NULL
  for (i in 1:a) {
    X_bar <- c(X_bar, rowMeans(data[, group == i])) # hier rowMeans, da data nicht transponiert ist und wir den mean über die Zeit haben wollen (-> wenn sich heraus stellt, dass wir mean über das individuum haben wollen, dann colMeans)
  }
  QN <- N * sum(X_bar *(TM %*% X_bar))
  W <- as.numeric((QN - EW) / sqrt(Var))

  # Test
  f <- as.numeric(A4^3 / C5^2)
  p.value <- 1 - stats::pchisq(W * sqrt(2 * f) + f, df = f)


  ## Ausgabe
  L <- list(f = f,
            statisitc = W,
            tau = 1/f,
            hypothesis = ifelse(is.list(hypothesis), list(TW = hypothesis$TW, TS = hypothesis$TS), hypothesis[1]), # gibt entweder eine liste mit (nur) den matrizen oder hypothesis[1]
            p.value = p.value,
            dim = list(d = d, N = N),
            groups = list(a = a, # hier könnte eine Tabelle stehen, bei der die ursprünglichen namen der gruppen stehen
                          table = table(group)),
            removed.cases = N_with_NA - N
            )
  class(L) <- c("hdrm")
  return(L)
}




#### Was soll eingegeben werden??
## X = matrix mit dimensionen c(d,N) -> Individuen in Spalten, Dimensionen in Zeilen
## group = vector der Länge N mit der Gruppenzuweisung
## hypothesis = character %in% c("time","group","interaction") oder liste mit Einträgen TW und TS mit dim(TW) = c(a,a) und dim(TS) = c(d,d)
#' @keywords internal
hdrm_test_internal <- function(data, group, hypothesis = c("whole","sub","interaction"), B){

  # N,d,a bestimmen
  N_with_NA <- dim(data)[2]
  group <- group[stats::complete.cases(t(data))] # hier mit factor arbeiten

  data <- t(stats::na.omit(t(data))) # data ist nicht transponiert!!!
  N <- dim(data)[2]
  d <- dim(data)[1]
  a <- length(table(group))
  n <- as.integer(table(group))
  stopifnot(length(group) == N)
  ### Stopkriterien
  # mind 2 Gruppen
  # mind 6 Individuen
  stopifnot(a >= 2, all(n >= 6)) # die checks sind ja tatsächlich auch schon in den datenverarbeitungsfunktionen...


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

  # spur1 und spur3
  for (i in 1:a) {
    A1[i] <- A1_cpp(mat = X_TS[, group == i])/(2*choose(n[i],2))
    A3[i] <- A3_cpp(mat = X_TS[, group == i], Part6 = sum(rowMeans(X_TS[, group == i])^2))
  }

  # spur2
  for (i in 1:(a-1)) {
    for(r in (i+1):a){
      A2[i,r] <- A2(X = X_TS[, group == i], Y = X_TS[, group == r])
    }
  }

  # C5
  C5 <- C5star_cpp(X = data, group = group, TW = TW, TS = TS, B = B)

  # E(Q_N)
  EW <- sum((N/n) * diag(TW) * A1)

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
  rm(temp1, temp2)

  # Var(Q_N)
  Var <- 2*A4


  ### Teststatistik
  X_bar <- NULL
  for (i in 1:a) {
    X_bar <- c(X_bar, rowMeans(data[, group == i])) # hier rowMeans, da data nicht transponiert ist und wir den mean über die Zeit haben wollen (-> wenn sich heraus stellt, dass wir mean über das individuum haben wollen, dann colMeans)
  }
  QN <- N * sum(X_bar *(TM %*% X_bar))
  W <- as.numeric((QN - EW) / sqrt(Var))

  # Test
  f <- as.numeric(A4^3 / C5^2)
  p.value <- 1 - stats::pchisq(W * sqrt(2 * f) + f, df = f)


  ## Ausgabe
  L <- list(f = f,
            statisitc = W,
            tau = 1/f,
            hypothesis = ifelse(is.list(hypothesis), list(TW = hypothesis$TW, TS = hypothesis$TS), hypothesis[1]), # gibt entweder eine liste mit (nur) den matrizen oder hypothesis[1]
            p.value = p.value,
            dim = list(d = d, N = N),
            groups = list(a = a, # hier könnte eine Tabelle stehen, bei der die ursprünglichen namen der gruppen stehen
                          table = table(group)),
            removed.cases = N_with_NA - N
  )
  class(L) <- c("hdrm")
  return(L)
}
