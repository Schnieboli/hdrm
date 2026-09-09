A1_eq <- function(X_list){
  n_i <- sapply(X_list, ncol)
  A1_body <- sum(sapply(X_list, A1_i_eq_cpp))
  A1_body / sum(n_i*(n_i-1))
}


A2_eq <- function(X_list){
  n_i <- sapply(X_list, ncol)
  A2_body <- sum(sapply(X_list, A2_i_eq_cpp))
  A2_body / (24 * sum(choose(n_i, 4)))
}


C1star_eq <- function(X_list, B){
  n_i <- sapply(X_list, ncol)
  B_i <- allocate_C1_subsamples(group_sizes = n_i, B = B)
  C1_body <- 0
  for(i in seq_along(X_list)){
    C1_body <- C1_body + C1star_i_eq_cpp(X_list[[i]], B_i[i])
  }
  C1_body / sum(8 * B_i)
}
