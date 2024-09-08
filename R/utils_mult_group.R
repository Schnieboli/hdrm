# Da es zwei interne Funktionen gibt, ist hier zu besseren übersicht eine funktion,
# die die hypothesenmatrizen extrahiert -> so muss nicht mehr in beiden
# internen Funktionen etwas geändert werden...

### get_hypothesis_mult - erstellt die Hypothesenmatrix fuer den mehrfaktoriellen Fall
#
# Input:
#  hypothesis - entweder ein character oder eine Liste mit Komponenten TW und TS
#  a - eine Zahl, Anzahl der Gruppen
#  d - eine Zahl, Anzahl der Dimensionen
#
# Output: eine Liste mit komponenten TW und TS
#  TW - Matrix TW
#  TS - Matrix TS

#' @keywords internal
get_hypothesis_mult <- function(hypothesis, a, d){
  stopifnot(is.list(hypothesis) | is.character(hypothesis))
  # Fall: hypothesis ist eine Liste
  if(is.list(hypothesis)){
    # ueberpruefen, ob es die gesuchten listenelemente gibt
    if(is.null(hypothesis$TW)) stop("no list entry with name 'TW' in hypothesis", call. = FALSE)
    if(is.null(hypothesis$TS)) stop("no list entry with name 'TS' in hypothesis", call. = FALSE)
    # Liste "auspacken"
    TW <- hypothesis$TW
    TS <- hypothesis$TS
    ## ueberpruefen, ob Matrizen richtige dimensionen haben
    stopifnot(is.matrix(TW), is.matrix(TS), is.numeric(TW), is.numeric(TS))
    stopifnot(dim(TW) == c(a,a), dim(TS) == c(d,d))
    # idempotenz
    if(!all.equal(TW^2,TW) | !all.equal(TS^2,TS) | any(TW != t(TW)) | any(TS != t(TS))) stop("TW and TS must be idempotent")
  }
  # Fall: hypothesis ist character
  if(is.character(hypothesis)){
    if(length(hypothesis) > 1) warning("length of hypothesis is > 1, only first element used.", call. = FALSE)
    if(!all(hypothesis %in% c("whole","sub","interaction"))) stop("hypothesis must be one of c('whole', 'sub', 'interaction') or a list with components TW and TS")
    if(hypothesis[1] == "whole"){
      # Matrizen erstellen
      TW <- diag(a) - matrix(1/a,a,a)
      TS <- matrix(1,d,d)/d
    }
    if(hypothesis[1] == "sub"){
      # Matrizen erstellen
      TW <- matrix(1,a,a)/a
      TS <- diag(d) - matrix(1/d,d,d)
    }
    if(hypothesis[1] == "interaction"){
      # Matrizen erstellen
      TW <- diag(a) - matrix(1,a,a)/a
      TS <- diag(d) - matrix(1,d,d)/d
    }
  }
  return(list(TW = TW, TS = TS))
}



################################################################################
######################### Trace estimators ####################################
################################################################################


# Exact -------------------------------------------------------------------

#' @keywords internal
A1 <- function(X){# braucht Individuen in Spalten und dimension in Zeilen
  Result <- 0
  n <- ncol(X)
  for (i in 1:(n - 1)){
    for (j in (i + 1):n){
      a1 <- X[, i] - X[, j]
      Result <- Result + sum(a1^2)
    }
  }
  return(Result/(2*choose(n,2)))
}


#' @keywords internal
A2 <- function(X, Y){# braucht Individuen in Zeilen und dimension in Spalten
  nX <- ncol(X) # Anzahl Individuen in Gruppe
  nY <- ncol(Y)
  PX <- diag(1, nX, nX) - matrix(data = 1 / nX, # Zentrierungsmatrix für Gruppe X
                                 nrow = nX,
                                 ncol = nX)
  PY <- diag(1, nY, nY) - matrix(data = 1 / nY, # Zentrierungsmatrix für Gruppe Y
                                 nrow = nY,
                                 ncol = nY)
  MX <- tcrossprod(PX, X)
  MY <- tcrossprod(Y, PY)
  MXY <-  MX %*% MY # Seite 38 im Paper
  EBSchaetzer = sum(MXY^2)/((nX-1)*(nY-1))
  return(EBSchaetzer)
}


#' @keywords internal
A3 <- function(X){
  nX = dim(X)[2]
  Part1 = 0
  Part2 = 0
  Part3 = 0
  Part4 = 0
  Part5 = 0
  Part6 = 0
  Part7 = 0

  for (l2 in 1:nX){
    a22 = sum(X[,l2]^2)
    Part7 = Part7 + a22
    for (l1 in 1:nX){
      a12 = sum(X[, l1] * X[, l2])
      Part1 = Part1 + a12^2 * (l1 != l2)
      for (l3 in 1:nX){
        a23 = sum(X[, l2] * X[, l3])
        a13 = sum(X[, l1] * X[, l3])
        Part5 = Part5 + a12 * a23 * (l1 != l2)
        Part2 = Part2 + a12 * a13 * (l1 != l2) * (l1 != l3) * (l2 != l3)
        Part3 = Part3 + a13 * (a23 + a12) * (l1 != l3) * (l2 != l3)
        Part4 = Part4 + a13 * a22 * (l1 != l2) * (l1 != l3) * (l2 != l3)
      }
    }
  }

  Part1 = Part1 * (nX - 2) * (nX - 3)
  Part2 = Part2 * (2 * nX - 5)
  a8 = rowMeans(X)
  Part6 = nX ^ 2 * sum(a8 * a8)

  PSSchaetzer = (Part1 - Part2 - Part3 - Part4 - Part5 + Part6 * (Part6 - Part7)) / (nX * (nX - 1) * (nX - 2) * (nX - 3))
  return(as.numeric(PSSchaetzer))
}

# bootstrap ------------------------------------------------------

A1star <- function(X, B){
  res <- 0
  n <- ncol(X)
  d <- nrow(X)
  sigma <- matrix(0, d, 2)

  for (b in 1:B) {
    sigma <- X[, sample.int(n, 2)]
    res <- res + sum((sigma[, 1] - sigma[, 2])^2)
  }
  return(res/(2*B))
}


A2star <- function(X, Y, B){
  res <- 0
  nX <- ncol(X)
  nY <- ncol(Y)
  d <- nrow(X)
  sigma <- matrix(0, d, 4)

  for (i in 1:B) {
    sigma[, c(1,2)] <- X[, sample.int(nX, 2)]
    sigma[, c(3,4)] <- Y[, sample.int(nY, 2)]
    res <- res + sum((sigma[, 1] - sigma[, 2]) * (sigma[, 3] - sigma[, 4]))^2
  }
  return(res/(4*B))
}


A3star <- function(X, B){
  res <- 0
  n <- ncol(X)
  d <- nrow(X)
  sigma <- matrix(0, d, 4)

  for (b in 1:B) {
    sigma <- X[, sample.int(n, 4)]
    res <- res + sum((sigma[, 1] - sigma[, 2]) * (sigma[, 3] - sigma[, 4]))^2
  }
  return(res/(4*B))
}


## old version of C5 without the "shortcut"
# C5star_alt <- function(X, group, TM, B){ # X mit Individuen in Spalten und Dimensionen in Zeilen
#
#   stopifnot(length(group) == ncol(X))
#   a <- length(table(group))
#   d <- nrow(X)
#   N <- ncol(X)
#   n <- as.integer(table(group))
#   for (i  in 1:a) {
#     X[, group == i] <- X[, group == i] * sqrt(N/n[i]) # Vorfaktor jetzt schon dranmultiplizieren
#   }
#
#   Rout <- numeric(1)
#   Z12 <- Z34 <- Z56 <- numeric(d*a)
#   for (b in 1:B) {
#     # Zufallsindividuen auswählen
#     sigma = matrix(0, d, 6*a)
#     for(i in 1:a){
#       # Indizes zeihen
#       sigma[,(i-1)*6 + (1:6)] <- X[, group == i][, sample.int(n[i], 6)]
#       ## Z_ij erstellen
#       Z12[1:d + d*(i-1)] <- (sigma[, 1 + 6*(i-1)] - sigma[, 2 + 6*(i-1)])
#       Z34[1:d + d*(i-1)] <- (sigma[, 3 + 6*(i-1)] - sigma[, 4 + 6*(i-1)])
#       Z56[1:d + d*(i-1)] <- (sigma[, 5 + 6*(i-1)] - sigma[, 6 + 6*(i-1)])
#     }
#     mprod1 <- crossprod(TM, Z34)
#     mprod2 <- crossprod(TM, Z56)
#     mprod3 <- crossprod(TM, Z12)
#     Rout <- Rout + (sum(Z12 * (mprod1)) * sum(Z34 * (mprod2)) * sum(Z56 * (mprod3)))
#   }
#   return(Rout/(8*B))
# }


#' @keywords internal
C5star <- function(X, group, TW, TS, B){ # X mit Individuen in Spalten und Dimensionen in Zeilen

  stopifnot(length(group) == ncol(X))
  # X sortieren nach Gruppe
  X <- X[, order(group)]
  group <- sort(group)

  a <- length(table(group))
  d <- nrow(X)
  N <- ncol(X)
  n <- as.integer(table(group))
  stopifnot(all(n >= 6))

  ind <- cumsum(c(1, n)) # zum indexen der Gruppen, ausnutzen, dass X sortiert ist!

  Y <- matrix(0, a*d, N)
  for(i in 1:a){ # Teil von X der i-ten Gruppe mit dem entsprechenden Teil der Hypothesenmatrix und den Vorfaktoren multiplizieren
    Y[, ind[i]:(ind[i+1]-1)] <- kronecker(TW[, i], TS%*%(X[, ind[i]:(ind[i+1]-1)] * sqrt(N/n[i])))
  }

  Rout <- numeric(1)
  for (b in 1:B) {
    # Vektoren Initialisieren/zurücksetzen
    Z12 <- numeric(d*a)
    Z34 <- numeric(d*a)
    Z56 <- numeric(d*a)
    # Zufallsindividuen auswählen
    sigma = matrix(0.0, a*d, 6*a)
    for(i in 1:a){
      # Indizes ziehen
      sigma[,(i-1)*6 + (1:6)] <- Y[ , sample(ind[i]:(ind[i+1]-1), 6)]
      ## Z_ij erstellen
      Z12 <- Z12 + (sigma[, 1 + 6*(i-1)] - sigma[, 2 + 6*(i-1)])
      Z34 <- Z34 + (sigma[, 3 + 6*(i-1)] - sigma[, 4 + 6*(i-1)])
      Z56 <- Z56 + (sigma[, 5 + 6*(i-1)] - sigma[, 6 + 6*(i-1)])
    }
    Rout <- Rout + (sum(Z12 * Z34) * sum(Z34 * Z56) * sum(Z56 * Z12))
  }
  return(Rout/(8*B))
}


#' @keywords internal
C5star_cpp <- function(X, group, TW, TS, B){ # ist zu groß; ich glaube es liegt daran, dass die interne Funktion nicht die Gruppen beachtet (obwohl sie das mMn sollte). Denn wenn ich bei der R-Funktion die Gruppenzuweisung weglasse, dann kommen Ergebnisse in ähnlicher Größenordnung :(
  # browser()
  stopifnot(length(group) == ncol(X))
  # X sortieren nach Gruppe
  X <- X[, order(group)]
  group <- sort(group)

  a <- length(table(group))
  d <- nrow(X)
  N <- ncol(X)
  n <- as.integer(table(group))
  stopifnot(all(n >= 6))

  ind <- cumsum(c(1, n)) # zum indexen der Gruppen, ausnutzen, dass X sortiert ist!

  Y <- matrix(0, a*d, N)
  for(i in 1:a){ # Teil von X der i-ten Gruppe mit dem entsprechenden Teil der Hypothesenmatrix und den Vorfaktoren multiplizieren
    Y[, ind[i]:(ind[i+1]-1)] <- kronecker(TW[, i], TS%*%(X[, ind[i]:(ind[i+1]-1)] * sqrt(N/n[i])))
  }


  return(C5star_cpp_internal(X = Y, group = group, B = B, n = n))
}
