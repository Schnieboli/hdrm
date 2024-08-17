#### Was soll eingegeben werden??
## X = matrix mit dimensionen c(d,N) -> Individuen in Zeilen, Dimensionen in Spalten
## group = vector der Länge N mit der Gruppenzuweisung
## hypothesis = character %in% c("time","group","interaction") oder liste mit Einträgen TW und TS mit dim(TW) = c(a,a) und dim(TS) = c(d,d)
#' @export
hdrm_test_internal <- function(data, group, hypothesis = c("whole","sub","interaction"), alpha = 0.05, B = "500*N", subsample = subsample){

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


  # Hypothese bestimmen -> Fälle:
    # character -> hypothesis[1]
    # list -> TW und TS aus der Liste nehmen
  stopifnot(is.list(hypothesis) | is.character(hypothesis))

  if(is.list(hypothesis)){
    TW <- hypothesis$TW
    TS <- hypothesis$TS
    stopifnot(is.matrix(TW), is.matrix(TS), is.numeric(TW), is.numeric(TS))
    stopifnot(dim(TW) == c(a,a), dim(TS) == c(d,d))
    TM <- kronecker(TW, TS)
  }

  if(is.character(hypothesis)){
    if(hypothesis[1] == "whole"){
      TW <- diag(a) - matrix(1/a,a,a)
      TS <- matrix(1,d,d)/d
    }
    if(hypothesis[1] == "sub"){
      TW <- matrix(1,a,a)/a
      TS <- diag(d) - matrix(1/d,d,d)
    }
    if(hypothesis[1] == "interaction"){
      TW <- diag(a) - matrix(1,a,a)/a
      TS <- diag(d) - matrix(1,d,d)/d
    }
    TM <- kronecker(TW, TS) # TM = T
  }

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

  C5 <- C5star(X = data, group = group, TM = TM, B = eval(expr = parse(text = B)))

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
  crit.value <- (stats::qchisq(1-alpha, df = f) - f) / sqrt(2*f)
  pWert_Kf <- function(p, param, statistic) abs( ((stats::qchisq(1-p, param) - param)/sqrt(2*param)) - statistic)
  # hier wird auf[0, 1] das minimum der funktion pWert_Kf gesucht
  # tol = .Machine$double.eps für maximale Accuracy -> sonst kommt ab einer bestimmten Extremität immer der gleiche Wert raus
  p.value <- stats::optimise(pWert_Kf, interval = c(0,1), param = f, statistic = W, tol = .Machine$double.eps)$minimum



  ## Ausgabe
  L <- list(data = data,
            f = f,
            statisitc = W,
            tau = 1/f,
            hypothesis = ifelse(is.list(hypothesis), list(TW = hypothesis$TW, TS = hypothesis$TS), hypothesis[1]), # gibt entweder eine liste mit (nur) den matrizen oder hypothesis[1]
            p.value = p.value,
            critical.value = crit.value,
            dim = list(d = d, N = N),
            groups = list(a = a, # hier könnte eine Tabelle stehen, bei der die ursprünglichen namen der gruppen stehen
                          table = table(group)),
            removed.cases = N_with_NA - N
            )
  class(L) <- c("hdrm")
  return(L)
}



#(x <- hdrm_test_internal(data = Y, group = group, hypothesis = "sub", alpha = 0.05, B = "500*N"))

### bei load_all() wird iwie nicht die richtige print methode aufgerufen...

#' @method print hdrm
#' @export
print.hdrm <- function(x,...){
  cat("\n")
  cat("          Multi Group Repeated Measure
        \nAnalysis of", x$dim$N, "individuals", paste0("(", x$removed.cases, " removed)"), "in", x$groups$a, "groups", "and", paste0(x$dim$d), "dimensions:",
      "\nW =", x$statisitc, " f =", round(x$f,4), " p.value =", round(x$p.value, 4),
      "\nNull-Hypothesis:", ifelse(is.list(x$hypothesis), "custom", paste0("No effect in ",x$hypothesis,"plot-factor")),
      "\nConvergence parameter \u03c4 =", round(x$tau,4))
  cat("\n")
}
#' @method summary hdrm
#' @export
summary.hdrm <- function(object,...){
  print(object,...)
}
