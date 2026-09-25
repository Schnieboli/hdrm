
# Unequal Covariance ------------------------------------------------------

A1 <- function(X_TS_list, subsampling, B){
  a <- length(X_TS_list)
  if(subsampling){
    A1 <- sapply(X_TS_list, A1star_i_cpp, B = B)
  }else{
    A1 <- sapply(X_TS_list, A1_i_cpp)
  }
  if (any(is.na(A1)) || any(!is.finite(A1)) || any(A1 < 0)) {
    stop("Internal error: please contact 'hdrm' package maintainer",
         call. = FALSE)
  }
  A1
}


A4 <- function(X_TS_list, TW, subsampling, B){
  a <- length(X_TS_list)
  n <- sapply(X_TS_list, ncol)
  N <- sum(n)

  if(subsampling){
    part1 <- sum(sapply(X_TS_list, A3star_i_cpp, B = B) * (N/n)^2 * diag(TW)^2)
  }else{
    part1 <- sum(sapply(X_TS_list, A3_i_cpp) * (N/n)^2 * diag(TW)^2)
  }

  part2 <- 0
  for(i in 1:(a-1)){
    for(r in (i+1):a){
      if(subsampling){
        part2 = part2 +
          (N^2 / (n[i]*n[r])) *
          TW[i,r]^2 *
          A2star_ir_cpp(X_TS_list[[i]], X_TS_list[[r]], B = B)
      }else{
        part2 = part2 +
          (N^2 / (n[i]*n[r])) *
          TW[i,r]^2 *
          A2_ir_cpp(X_TS_list[[i]], X_TS_list[[r]])
      }
    }
  }
  A4 <- part1 + 2 * part2

  if (is.na(A4) || !is.finite(A4) || A4 <= 0) {
    stop("Internal error: please contact 'hdrm' package maintainer",
         call. = FALSE)
  }

  A4
}


C5star <- function(data_list, TW, TS, B){
  a <- length(data_list)
  n <- sapply(data_list, ncol)
  N <- sum(n)

  Y_list = list()
  for(i in 1:a){
    Y_list[[i]] <- kronecker(TW[, i], TS %*% data_list[[i]]) * sqrt(N / n[i])
  }

  joint_B <- a * B

  C5star_cpp(Y_list, B = joint_B)
}


# Equal covariance --------------------------------------------------------

A1_eq <- function(X_TS_list){
  n_i <- sapply(X_TS_list, ncol)
  A1_body <- sum(sapply(X_TS_list, A1_i_eq_cpp))
  A1 <- A1_body / sum(n_i*(n_i-1))
  if (is.na(A1) || !is.finite(A1) || A1 < 0)
    stop("Internal error: please contact 'hdrm' package maintainer", call. = FALSE)
  A1
}


A2_eq <- function(X_TS_list){
  n_i <- sapply(X_TS_list, ncol)
  A2_body <- sum(sapply(X_TS_list, A2_i_eq_cpp))
  A2 <- A2_body / (24 * sum(choose(n_i, 4)))
  if (is.na(A2) || !is.finite(A2) || A2 < 0)
    stop("Internal error: please contact 'hdrm' package maintainer", call. = FALSE)
  A2
}


C1star_eq <- function(X_TS_list, B){
  n_i <- sapply(X_TS_list, ncol)
  B_i <- allocate_C1_subsamples(group_sizes = n_i, B = B)
  C1_body <- 0
  for(i in seq_along(X_TS_list)){
    C1_body <- C1_body + C1star_i_eq_cpp(X_TS_list[[i]], B_i[i])
  }
  C1 <- C1_body / sum(8 * B_i)
  if (is.na(C1) || !is.finite(C1) || C1 < 0)
    stop("Internal error: please contact 'hdrm' package maintainer", call. = FALSE)
  C1
}
