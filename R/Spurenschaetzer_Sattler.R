

# A1 ----------------------------------------------------------------------


#' @keywords internal
A1 <- function(X, n){# braucht Individuen in Zeilen und dimension in Spalten
  Result <- 0
  n <- dim(X)[1]
  for (k in 1:(n - 1))
  {
    for (l in (k + 1):n)
    {
      a1 <- (X[l, ] - X[k, ])
      Result <- Result + sum(a1 * a1)
    }
  }
  return(Result / (2 * choose(n, 2)))
}

# Idee: ziehe 6 zufällige Vektoren aus jeder Gruppe
A1_star <- function(B){
  for (i in 1:B) {
    #sample()
  }
}


# A2 ----------------------------------------------------------------------

#' @keywords internal
A2 <- function(X, Y){# braucht Individuen in Zeilen und dimension in Spalten
  nX <- dim(X)[1] # Anzahl Individuen in Gruppe
  nY <- dim(Y)[1]
  PX <- diag(1, nX, nX) - matrix(data = 1 / nX, # Zentrierungsmatrix für Gruppe X
                                 nrow = nX,
                                 ncol = nX)
  PY <- diag(1, nY, nY) - matrix(data = 1 / nY, # Zentrierungsmatrix für Gruppe Y
                                 nrow = nY,
                                 ncol = nY)
  MXY <- PX %*% (X) %*% t(Y) %*% t(PY) # Seite 38 im Paper
  EBSchaetzer = matrix(data = 1, nrow = 1, ncol = nX) %*%
    (MXY * MXY) %*% matrix(data = 1, nrow = nY, ncol = 1) / ((nX - 1) * (nY - 1))
  # (Edgar Brunner Schätzer)
  return(as.vector(EBSchaetzer))
}



# A3 ----------------------------------------------------------------------

# ## unterschätzt total, keine ahnung, wo der fehler liegt :(
# A3_Nils <- function(X){ # vielleicht mal ausprobieren mit einem datenssatz, wo man genau wüsste, wo im endeffekt der fehler liegt
#   # -> nur 1en oder so
#   nX <- dim(X)[1]
#
#   S <- 0
#   comb <- combn(nX, 2)
#   # Schleife über die Spalten von comb -> iwas funktioniert nicht
#   for(i in choose(nX, 2)){ # wahrscheinlich etwas falsch mit der Reihenfolge der Schleifen
#     for (k2 in 1:(nX-1)) {
#       for(k1 in k2:nX){
#         if((comb[1,i] != k1) & (comb[2,i] != k1) & (comb[1,i] != k2) & (comb[2,i] != k2)){
#           S <- S + sum((X[comb[1,i],] - X[comb[2,i],]) * (X[k2,] - X[k1,]))^2
#         }
#       }
#     }
#   }
#   return(S / (24*choose(nX,4)))
# }
#
# # gleiches ergebnis wie bei a3_nils -> warum???
# A3_Nils2 <- function(X){
#   nX <- dim(X)[1]
#   S <- 0
#   comb <- combn(nX, 2)
#
#   for(i in choose(nX, 2)){
#     l1 <- comb[1,i]
#     l2 <- comb[2,i]
#     for (k2 in 1:(nX-1)) {
#       for(k1 in k2:nX){
#         if(k2 != l1 & k2 != l2 & k1 != l1 & k1 != l2) S <- S + sum((X[l1, ] - X[l2, ]) * (X[k1, ] - X[k2, ]))^2
#       }
#     }
#   }
#   return(S / (24*choose(nX,4)))
# }
#
# A3_schlecht <- function(X){
#   nX <- dim(X)[1]
#   S <- 0
#   comb <- combn(nX, 2)
#
#   for (i in 1:choose(nX,2)) {
#     l1 <- comb[1,i]
#     l2 <- comb[2,i]
#     for(k2 in 1:(nX-1)){
#       for(k1 in k2:nX){
#         if(k2 != l1 & k2 != l2 & k1 != l1 & k1 != l2) S <- S + sum((X[l1, ] - X[l2, ]) * (X[k1, ] - X[k2, ]))^2
#       }
#     }
#   }
#   return(S / (24 * choose(nX,4)))
# }

#' @keywords internal
A3 <- function(X){ # braucht Individuen in Zeilen und dimension in Spalten
  nX = dim(X)[1]
  Part1 = 0
  Part2 = 0
  Part3 = 0
  Part4 = 0
  Part5 = 0
  Part6 = 0
  Part7 = 0

  # browser()

  for (l2 in 1:(nX)){
    a22 = t(X[l2, ]) %*% X[l2, ]
    Part7 = Part7 + a22
    for (l1 in 1:nX){
      a12 = (t(X[l1, ]) %*% X[l2, ])
      Part1 = Part1 + (a12) ^ 2 * (l1 != l2)
      for (l3 in 1:nX){
        a23 = (t(X[l2, ]) %*% X[l3, ])
        a13 = (t(X[l1, ]) %*% X[l3, ])
        Part5 = Part5 + a12 * a23 * (l1 != l2)
        Part2 = Part2 + a12 * a13 * (l1 != l2) * (l1 != l3) * (l2 != l3)
        Part3 = Part3 + a13 * (a23 + a12) * (l1 != l3) * (l2 != l3)
        Part4 = Part4 + a13 * a22 * (l1 != l2) * (l1 != l3) * (l2 != l3)
      }
    }
  }
  Part1 = Part1 * (nX - 2) * (nX - 3)
  Part2 = Part2 * (2 * nX - 5)
  a8 = colMeans(X)
  Part6 = nX ^ 2 * sum(a8 * a8)

  PSSchaetzer = (Part1 - Part2 - Part3 - Part4 - Part5 + Part6 * (Part6 - Part7)) / (nX * (nX - 1) * (nX - 2) * (nX - 3))
  return(PSSchaetzer)
}



# C5star=0
# for (l in 1:(w))
# {
#   i=sample(seq(1,n1,by=1),6,replace=FALSE)
#   j=sample(seq(1,n2,by=1),6,replace=FALSE)
#   a9=(sqrt(1/0.4)*(X1T[i[1],1:d]-X1T[i[2],1:d])+sqrt(1/0.6)*(X2T[j[1],1:d]-X2T[j[2],1:d]))
#   a10=(sqrt(1/0.4)*(X1T[i[3],1:d]-X1T[i[4],1:d])+sqrt(1/0.6)*(X2T[j[3],1:d]-X2T[j[4],1:d]))
#   a11=(sqrt(1/0.4)*(X1T[i[5],1:d]-X1T[i[6],1:d])+sqrt(1/0.6)*(X2T[j[5],1:d]-X2T[j[6],1:d]))
#   C5star=C5star+sum(a9*a10)*sum(a10*a11)*sum(a11*a9)
# }

# Eigener Versuch C5 ------------------------------------------------------

#' @keywords internal
C5stern_alt <- function(X, w, N, a, n){ # die funktion ist für eine liste geschrieben!!! -> in der neuen funktion hdrm_test ist X jedoch theoretisch eine matrix...!!!
  ind <- matrix(0,6,a)
  C5 <- 0
  for (p in 1:w){
    for(i in 1:a){ # geht safe auch ohne schleife
      ind[,i] <- sample(1:n[i], 6)
    }
    Z12 <- Z34 <- Z56 <- numeric(0)
    for (i in 1:a) {
      Z12 <- c(Z12, (X[[i]][ind[1,i],] - X[[i]][ind[2,i],]) * sqrt(N / n[i]))
      Z34 <- c(Z34, (X[[i]][ind[3,i],] - X[[i]][ind[4,i],]) * sqrt(N / n[i]))
      Z56 <- c(Z56, (X[[i]][ind[5,i],] - X[[i]][ind[6,i],]) * sqrt(N / n[i]))
    }
    C5 <- C5 + (sum(Z12 * Z34) + sum(Z34 * Z56) + sum(Z56*Z12))
  }

  return(C5/(8*w))
}

# C5stern_neu <- function(X, w, N, a, n){ # die funktion ist für eine liste geschrieben!!! -> in der neuen funktion hdrm_test ist X jedoch theoretisch eine matrix...!!!
#   ind <- matrix(0,6,a)
#   C5 <- 0
#   for (p in 1:w){
#     for(i in 1:a){ # geht safe auch ohne schleife
#       ind[,i] <- sample(1:n[i], 6)
#     }
#     Z12 <- Z34 <- Z56 <- numeric(0)
#     for (i in 1:a) {
#       Z12 <- c(Z12, (X[[i]][,ind[1,i]] - X[[i]][,ind[2,i]]) * sqrt(N / n[i]))
#       Z34 <- c(Z34, (X[[i]][,ind[3,i]] - X[[i]][,ind[4,i]]) * sqrt(N / n[i]))
#       Z56 <- c(Z56, (X[[i]][,ind[5,i]] - X[[i]][,ind[6,i]]) * sqrt(N / n[i]))
#     }
#     C5 <- C5 + (sum(Z12 * Z34) + sum(Z34 * Z56) + sum(Z56*Z12))
#   }
#
#   return(C5/(8*w))
# }
#
# Daten
# C5 <- 0
# ind <- matrix(0,6,a)
# L1 <- L2 <- L3 <- numeric(a)
# group <- rep(1:3, each = 7)
# X_groups <- split(X, group)
# X_groups <- lapply(X_groups, as.matrix)
# N <- 21
#
# for (w in 1:10){
#   for(i in 1:a){ # geht safe auch ohne schleife
#     ind[,i] <- sample(1:n_i[i], 6)
#   }
#   Z12 <- Z34 <- Z56 <- numeric(0)
#   for (i in 1:a) {
#     Z12 <- c(Z12, X_groups[[i]][ind[1,i],] - X_groups[[i]][ind[2,i],] * sqrt(N / n_i[i]))
#     Z34 <- c(Z34, X_groups[[i]][ind[3,i],] - X_groups[[i]][ind[4,i],] * sqrt(N / n_i[i]))
#     Z56 <- c(Z56, X_groups[[i]][ind[5,i],] - X_groups[[i]][ind[6,i],] * sqrt(N / n_i[i]))
#
#
#   }
#
#   C5 <- C5 + (sum(Z12 * Z34) + sum(Z34 * Z56) + sum(Z56*Z12))
# }
