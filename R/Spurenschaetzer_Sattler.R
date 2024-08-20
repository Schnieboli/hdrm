

# A1 ----------------------------------------------------------------------


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


  # A2 ----------------------------------------------------------------------

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



# A3 ----------------------------------------------------------------------
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

# C5star ------------------------------------------------------

#' @keywords internal
C5star <- function(X, group, TM, B){ # X mit Individuen in Spalten und Dimensionen in Zeilen
  stopifnot(length(group) == ncol(X))
  a <- length(table(group))
  d <- nrow(X)
  N <- ncol(X)
  n <- as.integer(table(group))
  for (i  in 1:a) {
    X[, group == i] <- X[, group == i] * sqrt(N/n[i]) # Vorfaktor jetzt schon dranmultiplizieren
  }

  Rout <- numeric(1)
  Z12 <- Z34 <- Z56 <- numeric(d*a)
  for (b in 1:B) {
    # Zufallsindividuen auswählen
    sigma = matrix(0, d, 6*a)
    for(i in 1:a){
      # Indizes zeihen
      sigma[,(i-1)*6 + (1:6)] <- X[, group == i][, sample.int(n[i], 6)]
      ## Z_ij erstellen
      Z12[1:d + d*(i-1)] <- (sigma[, 1 + 6*(i-1)] - sigma[, 2 + 6*(i-1)])
      Z34[1:d + d*(i-1)] <- (sigma[, 3 + 6*(i-1)] - sigma[, 4 + 6*(i-1)])
      Z56[1:d + d*(i-1)] <- (sigma[, 5 + 6*(i-1)] - sigma[, 6 + 6*(i-1)])
    }
    mprod1 <- crossprod(TM, Z34)
    mprod2 <- crossprod(TM, Z56)
    mprod3 <- crossprod(TM, Z12)
    Rout <- Rout + (sum(Z12 * (mprod1)) * sum(Z34 * (mprod2)) * sum(Z56 * (mprod3)))
  }
  return(Rout/(8*B))
}


C5star_list <- function(X, group, TM, B){ # X mit Individuen in Spalten und Dimensionen in Zeilen
  #browser()
  stopifnot(length(group) == ncol(X))
  a <- length(table(group))
  d <- nrow(X)
  N <- ncol(X)
  n <- as.integer(table(group))
  stopifnot(all(n >= 6))
  for (i  in 1:a) {
    X[, group == i] <- X[, group == i] * sqrt(N/n[i]) # Vorfaktor jetzt schon dranmultiplizieren
  }

  X_list <- vector(mode = "list", length = a)
  for(i in 1:a){
    X_list[[i]] <- TM[, (1:d) + d*(i-1)] %*% X[, group == i]
  }


  Rout <- numeric(1)
  for (b in 1:B) {
    # Vektoren Initialisieren
    Z12 <- numeric(d*a)
    Z34 <- numeric(d*a)
    Z56 <- numeric(d*a)
    # Zufallsindividuen auswählen
    sigma = matrix(0, a*d, 6*a)
    for(i in 1:a){
      # Indizes zeihen
      sigma[,(i-1)*6 + (1:6)] <- X_list[[i]][, sample.int(n[i], 6)]
      ## Z_ij erstellen
      Z12 <- Z12 + (sigma[, 1 + 6*(i-1)] - sigma[, 2 + 6*(i-1)])
      Z34 <- Z34 + (sigma[, 3 + 6*(i-1)] - sigma[, 4 + 6*(i-1)])
      Z56 <- Z56 + (sigma[, 5 + 6*(i-1)] - sigma[, 6 + 6*(i-1)])
    }
    Rout <- Rout + (sum(Z12 * Z34) * sum(Z34 * Z56) * sum(Z56 * Z12))
  }
  return(Rout/(8*B))
}
