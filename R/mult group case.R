
#' @export
hdrm_test <- function(data, formula,...){
  UseMethod("hdrm_test")
}

#' @method hdrm_test default
#' @export
hdrm_test.default <- function(data, formula, hypothesis, TW, TS, alpha, B, subsample){
  stop("Your data must be either a list, a data frame or a matrix")
}


#' @method hdrm_test list
#' @export
hdrm_test.list <- function(data, hypothesis, TW, TS, alpha = 0.05, B = "500*N", subsample = FALSE){
  ...
}

#' @method hdrm_test matrix
#' @export
hdrm_test.matrix <- function(data, group, hypothesis, TW, TS, alpha = 0.05, B = "500*N", subsample = FALSE){
  ...
}

#' @method hdrm_test data.frame
#' @export
hdrm_test.data.frame <- function(data, formula, hypothesis, TW, TS, alpha = 0.05, B = "100*N", subsample = FALSE){
  ...
}



#### Was soll eingegeben werden??
## X = matrix mit dimensionen c(d,N) -> Individuen in Spalten, Dimension in Zeilen
## group = vector der Länge N mit der Gruppenzuweisung
## hypothesis = character %in% c("time","group","interaction") oder Matrix mit Dimension c(d,d)
#' @export
hdrm_test_internal <- function(data, group, hypothesis = c("whole","sub","interaction"), TW, TS, alpha = 0.05, B = "500*N", subsample = subsample){

  #browser()

  # N,d,a bestimmen
  N_with_NA <- dim(data)[2]
  group <- group[complete.cases(t(data))] # hier mit factor arbeiten

  data <- na.omit(data)
  N <- dim(data)[2]
  d <- dim(data)[1]
  a <- length(table(group))
  n <- as.integer(table(group))

  ### Stopkriterien
  # mind 2 Gruppen
  # mind 6 Individuen
  stopifnot(a >= 2, all(n >= 6))


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

  if(is.character(hypothesis)){ # hypothesis gegeben und TW,TS missing
    if(hypothesis[1] == "whole"){
      TW <- diag(a) - matrix(1,a,a)/a
      TS <- matrix(1,d,d)/d
    }
    if(hypothesis[1] == "sub"){
      TW <- matrix(1,a,a)/a
      TS <- diag(d) - matrix(1,d,d)/d
    }
    if(hypothesis[1] == "interaction"){
      TW <- diag(a) - matrix(1,a,a)/a
      TS <- diag(d) - matrix(1,d,d)/d
    }
    TM <- kronecker(TW, TS) # TM = T
  }


  ### Spurenschätzer -> hier weiß ich nicht, ob alles richtig ist...

  # X vorbereiten
  X <- TS %*% data

  #
  spur1 <- spur3 <- spur4 <- 0
  spur2 <- matrix(0, a, a)
  # spur1
  for (i in 1:a) {
    spur1 <- spur1 + A1(X = X[, group == i], n = n[i])
  }

  # spur2
  for (i in 1:(a-1)) {
    for(r in (i+1):a){
      spur2[i,r] <- A2(X = X[, group == i], Y = X[, group == r])
    }
  }

  # spur3
  for (i in 1:a) {
    spur3 <- spur3 + A3(X = X[, group == i])
  }


  ########## Hier stimm iwas mit der Funktion C5star nicht -> muss verändert werden -> beide Versionen klappen nicht :(
  # spur4
  #spur4 <- C5stern(X = X, w = eval(parse(text = B)), a = a, n = n)
  ## jetzt erstmal im alten stil -> das heißt: X in Liste umwandeln
  X_list <- vector(mode = "list")
  for(i in 1:a){
    X_list <- c(X_list, list(X[, group == i]))
  }
  spur4 <- C5stern_alt(X = X_list, w = eval(parse(text = B)), N = N, a = a, n = n)

  # Erwartungswert
  EW <- sum((N/n) * diag(TW)) * spur1 # falls spur1 nicht vektorisiert

  # Varianz
  hilf <- 0
  for (i in 1:(a-1)) {
    for(r in (i+1):a){
      hilf = hilf + ( (N^2 / (n[i]*n[r])) * TW[i,r]^2 * spur2[i,r])
    }
  }
  Var <- (sum((N/n)^2 * TW^2) * spur3 + 2*hilf) # das ist eigentlich nicht Var, sondern A4!!!
  rm(hilf)

  ### Teststatistik
  X_quer <- NULL
  for (i in 1:a) {
    X_quer <- c(X_quer, colMeans(data[, group == i]))
  }
  Q <- N * sum(T %*% X_quer^2)
  W <- as.numeric((Q - EW) / sqrt(2*Var))

  # Test
  f <- as.numeric(Var^3 / spur4^2)
  crit.value <- (qchisq(1-alpha, df = f) - f) / sqrt(2*f)
  pWert_Kf <- function(p, param, statistic) abs( ((qchisq(1-p, param) - param)/sqrt(2*param)) - statistic)
  # hier wird auf[0, 1] das minimum der funktion pWert_Kf gesucht
  # tol = .Machine$double.eps für maximale Accuracy -> sonst kommt ab einer bestimmten Extremität immer der gleiche Wert raus
  p.value <- optimise(pWert_Kf, interval = c(0,1), param = f, statistic = W, tol = .Machine$double.eps)$minimum

  #browser()

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
      "\nW =", x$statisitc, " f =", x$f, " p.value =", x$p.value,
      "\nNull-Hypothesis:", ifelse(is.list(x$hypothesis), "custom", paste0("No effect in ",x$hypothesis,"plot-factor")),
      "\nConvergence parameter \u03c4 =", x$tau)
  cat("\n")
}
#' @method summary hdrm
#' @export
summary.hdrm <- function(x,...){
  print(x,...)
}
