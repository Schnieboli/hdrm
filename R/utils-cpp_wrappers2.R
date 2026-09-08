Exp_Q <- function(X_list, TW, subsampling, B){
  a <- length(X_list)
  if(subsampling){
    A1 <- sapply(X_list, A1star_i_cpp, B = B)
  }else{
    A1 <- sapply(X_list, A1_i_cpp)
  }
  n <- sapply(X_list, ncol)
  N <- sum(n)
  sum(N/n * diag(TW) * A1)
}

A4 <- function(X_list, TW, subsampling, B){
  a <- length(X_list)
  n <- sapply(X_list, ncol)
  N <- sum(n)
  
  if(subsampling){
    part1 <- sum(sapply(X_list, A3star_i_cpp, B = B) * (N/n)^2 * diag(TW)^2)
  }else{
    part1 <- sum(sapply(X_list, A3_i_cpp) * (N/n)^2 * diag(TW)^2)
  }
  
  part2 <- 0
  for(i in 1:(a-1)){
    for(r in (i+1):a){
      if(subsampling){
        part2 = part2 +
          (N^2 / (n[i]*n[r])) *
          TW[i,r]^2 *
          A2star_ir_cpp(X_list[[i]], X_list[[r]], B = B)
      }else{
        part2 = part2 +
          (N^2 / (n[i]*n[r])) *
          TW[i,r]^2 *
          A2_ir_cpp(X_list[[i]], X_list[[r]])
      }
    }
  }
  part1 + 2 * part2
}
