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
