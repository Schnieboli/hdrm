# Da es zwei interne Funktionen gibt, ist hier zu besseren uebersicht eine funktion,
# die die hypothesenmatrizen extrahiert -> so muss nicht mehr in beiden
# internen Funktionen etwas geaendert werden...

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
    if(any(dim(TW) != c(a,a))) stop(paste0("TW should be a quadratic matrix with ", a, " rows"))
    if(any(dim(TS) != c(d,d))) stop(paste0("TS should be a quadratic matrix with ", d, " rows"))
    # Symmetrie
    if(!isSymmetric.matrix(TW) | !isSymmetric.matrix(TS)) stop("TW and TS must be symmetric")
    ############# Annahmen TW
    # Symmetrie und Idempotenz pruefen
    # Symmetrie
    if((mean(t(TW) - TW) >= sqrt(.Machine$double.eps))) warning(paste0("TW is not symmetric (mean difference = ", mean(t(TW) - TW),"). This will likely have a big influence on the test result!"))
    # Idempotenz
    # checke ob all.equal TRUE ist -> wenn ja, dann milde warnung
    if((mean(TW%*%TW - TW) >= sqrt(.Machine$double.eps))) warning(paste0("TW is not idempotent (mean difference = ", mean(TW%*%TW - TW),"). This will likely have a big influence on the test result!"))

    ############# Annahmen TS
    # Symmetrie und Idempotenz pruefen
    # Symmetrie
    if((mean(t(TS) - TS) >= sqrt(.Machine$double.eps))) warning(paste0("TS is not symmetric (mean difference = ", mean(t(TS) - TS),"). This will likely have a big influence on the test result!"))
    # Idempotenz
    # checke ob all.equal TRUE ist -> wenn ja, dann milde warnung
    if((mean(TS%*%TS - TS) >= sqrt(.Machine$double.eps))) warning(paste0("TS is not idempotent (mean difference = ", mean(TS%*%TS - TS),"). This will likely have a big influence on the test result!"))

  }

  # Fall: hypothesis ist character
  if(is.character(hypothesis)){
    if(length(hypothesis) > 1) warning("length of hypothesis is > 1, only first element used.", call. = FALSE)
    if(!(hypothesis[1] %in% c("whole","sub","interaction", "identical", "flat"))) stop("hypothesis must be one of c('whole', 'sub', 'interaction', 'identical', 'flat') or a list with components TW and TS")
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
    if(hypothesis[1] == "identical"){
      # Matrizen erstellen
      TW <- diag(a) - matrix(1,a,a)/a
      TS <- diag(d)
    }
    if(hypothesis[1] == "flat"){
      # Matrizen erstellen
      TW <- diag(a)
      TS <- diag(d) - matrix(1,d,d)/d
    }
  }

  ### ueberpruefen, ob die matrizen numeric sind oder NAs haben
  # sollte an dieser stelle alles unmoeglich sein, aber just in case
  if(!is.matrix(TW)) stop("TW must be a matrix")
  if(!is.numeric(TW)) stop("TW must be numeric")
  if(any(is.na(TW))) stop("TW must not contain NAs")
  if(!is.matrix(TS)) stop("TS must be a matrix")
  if(!is.numeric(TS)) stop("TS must be numeric")
  if(any(is.na(TS))) stop("TS must not contain NAs")

  # Ausgabe
  return(list(TW = TW, TS = TS))
}


# check_criteria_grouped - ueberprueft, ob alles so ist, wie es sein soll
#
# Input:    X - Datenmatrix X
#       group - Vektor group
#  hypothesis - hypothesis
#        reps - reps
#   bootstrap - bootstrap
#
# Output: bricht den Vorgang mit Fehlermeldung ab
#
# Motivation: auf diese Weise sind alle Kriterien an einer Stelle und falls ein
# neues Kriterium auftaucht, muss es nicht in beiden Methoden geaendert werden...
#
#' @keywords internal
check_criteria_grouped <- function(X, group, hypothesis, reps, subsampling){

  n <- as.integer(table(group))
  a <- length(n)
  d <- nrow(X)
  N <- ncol(X)

  if(a < 2) stop("At least two groups needed. Did you mean to call hdrm_single_...?", call. = FALSE)
  if(d < 2) stop("At least 2 dimensions are needed", call. = FALSE)
  if(any(n < 6)) stop("all group sizes must be >= 6", call. = FALSE)
  if(length(group) != N) stop("length(group) must be ncol(data)", call. = FALSE)
  if(!is.character(hypothesis) & !is.list(hypothesis)) stop("hypothesis must be a character or a list", call. = FALSE)
  if(length(subsampling) > 1) warning("length(subsampling) > 1. Only first element used", call. = FALSE)
  if(!(subsampling[1] %in% c(1, 0, TRUE, FALSE, T, F))) stop("subsampling must be logical", call. = FALSE)
  if(reps < 0) stop("subsamples must be > 0", call. = FALSE)
  if(!is.numeric(reps)) stop("subsamples could not be converted to a number", call. = FALSE)

}


################################################################################
#######################   Trace estimators   ###################################
################################################################################

## all trace estimators expect subject in columns and factor levels in rows
## estimators that are not used anywhere in the package are commented

# Exact -------------------------------------------------------------------

# #' @keywords internal
# A1 <- function(X){
#   Result <- 0
#   n <- ncol(X)
#   for (i in 1:(n - 1)){
#     for (j in (i + 1):n){
#       Result <- Result + crossprod(X[, i] - X[, j])
#     }
#   }
#   return(as.numeric(Result)/(n*(n-1)))
# }


#' @keywords internal
A2 <- function(X, Y){
  nX <- ncol(X)
  nY <- ncol(Y)
  PX <- diag(1, nX, nX) - matrix(data = 1 / nX, # Zentrierungsmatrix fuer Gruppe X
                                 nrow = nX,
                                 ncol = nX)
  PY <- diag(1, nY, nY) - matrix(data = 1 / nY, # Zentrierungsmatrix fuer Gruppe Y
                                 nrow = nY,
                                 ncol = nY)
  MX <- tcrossprod(PX, X)
  MY <- tcrossprod(Y, PY)
  MXY <-  MX %*% MY # Seite 38 im Paper
  EBSchaetzer = sum(MXY^2)/((nX-1)*(nY-1))
  return(EBSchaetzer)
}


# #' @keywords internal
# A3 <- function(X){
#   nX <- ncol(X)
#   Part1 <- 0
#   Part2 <- 0
#   Part3 <- 0
#   Part4 <- 0
#   Part5 <- 0
#   Part6 <- 0
#   Part7 <- 0
#
#   for (l2 in 1:nX){
#     a22 = as.numeric(crossprod(X[,l2]))
#     Part7 = Part7 + a22
#     for (l1 in 1:nX){
#       a12 = crossprod(X[, l1], X[, l2])
#       Part1 = Part1 + as.numeric(a12^2) * (l1 != l2)
#       for (l3 in 1:nX){
#         a23 = as.numeric(crossprod(X[, l2], X[, l3]))
#         a13 = as.numeric(crossprod(X[, l1], X[, l3]))
#         Part5 = Part5 + a12 * a23 * (l1 != l2)
#         Part2 = Part2 + a12 * a13 * (l1 != l2) * (l1 != l3) * (l2 != l3)
#         Part3 = Part3 + a13 * (a23 + a12) * (l1 != l3) * (l2 != l3)
#         Part4 = Part4 + a13 * a22 * (l1 != l2) * (l1 != l3) * (l2 != l3)
#       }
#     }
#   }
#
#   Part1 = Part1 * (nX - 2) * (nX - 3)
#   Part2 = Part2 * (2 * nX - 5)
#   a8 = rowMeans(X)
#   Part6 = nX ^ 2 * crossprod(a8)
#
#   PSSchaetzer = (Part1 - Part2 - Part3 - Part4 - Part5 + Part6 * (Part6 - Part7)) / (nX * (nX - 1) * (nX - 2) * (nX - 3))
#   return(as.numeric(PSSchaetzer))
# }


# bootstrap ------------------------------------------------------

# A1star <- function(X, B){
#   res <- 0
#   n <- ncol(X)
#   d <- nrow(X)
#   sigma <- matrix(0, d, 2)
#
#   for (b in 1:B) {
#     sigma <- X[, sample.int(n, 2)]
#     res <- res + sum((sigma[, 1] - sigma[, 2])^2)
#   }
#   return(res/(2*B))
# }


# A2star <- function(X, Y, B){
#   res <- 0
#   nX <- ncol(X)
#   nY <- ncol(Y)
#   d <- nrow(X)
#   sigma <- matrix(0, d, 4)
#
#   for (i in 1:B) {
#     sigma[, c(1,2)] <- X[, sample.int(nX, 2)]
#     sigma[, c(3,4)] <- Y[, sample.int(nY, 2)]
#     res <- res + sum((sigma[, 1] - sigma[, 2]) * (sigma[, 3] - sigma[, 4]))^2
#   }
#   return(res/(4*B))
# }


# A3star <- function(X, B){
#   res <- 0
#   n <- ncol(X)
#   d <- nrow(X)
#   sigma <- matrix(0, d, 4)
#
#   for (b in 1:B) {
#     sigma <- X[, sample.int(n, 4)]
#     res <- res + sum((sigma[, 1] - sigma[, 2]) * (sigma[, 3] - sigma[, 4]))^2
#   }
#   return(res/(4*B))
# }


## old version of C5 without the "shortcut"
# C5star_alt <- function(X, group, TM, B){
#
#   stopifnot(length(group) == ncol(X))
#   a <- length(table(group))
#   d <- nrow(X)
#   N <- ncol(X)
#   n <- as.integer(table(group))
#   for (i  in 1:a) {
#     X[, group == i] <- X[, group == i] * sqrt(N/n[i])
#   }
#
#   Rout <- numeric(1)
#   Z12 <- Z34 <- Z56 <- numeric(d*a)
#   for (b in 1:B) {
#     # Zufallsindividuen auswaehlen
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


# #' @keywords internal
# C5star <- function(X, group, TW, TS, B){ # X mit Individuen in Spalten und Dimensionen in Zeilen
#
#   stopifnot(length(group) == ncol(X))
#
#   a <- length(table(group))
#   d <- nrow(X)
#   N <- ncol(X)
#   n <- as.integer(table(group))
#   stopifnot(all(n >= 6))
#
#   ind <- cumsum(c(1, n)) # zum indexen der Gruppen, ausnutzen, dass X sortiert ist!
#
#   Y <- matrix(0, a*d, N)
#   for(i in 1:a){ # Teil von X der i-ten Gruppe mit dem entsprechenden Teil der Hypothesenmatrix und den Vorfaktoren multiplizieren
#     Y[, ind[i]:(ind[i+1]-1)] <- kronecker(TW[, i], TS%*%(X[, ind[i]:(ind[i+1]-1)] * sqrt(N/n[i])))
#   }
#
#   Rout <- numeric(1)
#   for (b in 1:B) {
#     # Vektoren Initialisieren/zuruecksetzen
#     Z12 <- numeric(d*a)
#     Z34 <- numeric(d*a)
#     Z56 <- numeric(d*a)
#     # Zufallsindividuen auswaehlen
#     sigma = matrix(0.0, a*d, 6*a)
#     for(i in 1:a){
#       # Indizes ziehen
#       sigma[,(i-1)*6 + (1:6)] <- Y[ , sample(ind[i]:(ind[i+1]-1), 6)]
#       ## Z_ij erstellen
#       Z12 <- Z12 + (sigma[, 1 + 6*(i-1)] - sigma[, 2 + 6*(i-1)])
#       Z34 <- Z34 + (sigma[, 3 + 6*(i-1)] - sigma[, 4 + 6*(i-1)])
#       Z56 <- Z56 + (sigma[, 5 + 6*(i-1)] - sigma[, 6 + 6*(i-1)])
#     }
#     Rout <- Rout + (sum(Z12 * Z34) * sum(Z34 * Z56) * sum(Z56 * Z12))
#   }
#   return(Rout/(8*B))
# }


#' @keywords internal
C5star_cpp <- function(X, group, TW, TS, B){ # ist zu gross; ich glaube es liegt daran, dass die interne Funktion nicht die Gruppen beachtet (obwohl sie das mMn sollte). Denn wenn ich bei der R-Funktion die Gruppenzuweisung weglasse, dann kommen Ergebnisse in aehnlicher Groessenordnung :(

    stopifnot(length(group) == ncol(X))

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
  # cpp-teil aufrufen
  return(C5star_cpp_internal(X = Y, group = group, B = B, n = n))
}
