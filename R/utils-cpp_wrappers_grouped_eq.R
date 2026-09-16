A1_eq <- function(X_TS_list){
  n_i <- sapply(X_TS_list, ncol)
  A1_body <- sum(sapply(X_TS_list, A1_i_eq_cpp))
  A1 <- A1_body / sum(n_i*(n_i-1))
  if (is.na(A1) || !is.finite(A1) || A1 < 0)
    stop("The test is degenerate for these data.", call. = FALSE)
  A1
}


A2_eq <- function(X_TS_list){
  n_i <- sapply(X_TS_list, ncol)
  A2_body <- sum(sapply(X_TS_list, A2_i_eq_cpp))
  A2 <- A2_body / (24 * sum(choose(n_i, 4)))
  if (is.na(A2) || !is.finite(A2) || A2 < 0)
    stop("The test is degenerate for these data.", call. = FALSE)
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
    stop("The test is degenerate for these data.", call. = FALSE)
  C1
}
