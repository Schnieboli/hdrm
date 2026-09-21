hdrm_grouped_internal <- function(X_list, H, cov.equal, subsampling, B, seed){
  
  # Determine the number of samples (N), dimensions (d), groups (a), and group sizes (n)
  a <- length(X_list)
  n <- sapply(X_list, ncol)
  d <- nrow(X_list[[1]])
  N <- sum(n)
  
  ## multiply data with TSalt
  X_TS_list <- lapply(X_list, function(x) H$TSalt %*% x)
  
  # Temporarily set the seed and restore the previous RNG state on exit
  if (!is.null(seed)) {
    withr::local_seed(seed)
  }
  
  if(cov.equal){
    EW <- sum((N / n) * diag(H$TW)) * A1(X_TS_list)
    second_order <- A2_eq(X_TS_list)
    third_order <- C1star_eq(X_TS_list, B = B)
    Var <- 2 * second_order * sum((H$TW)^2 * (N^2 / outer(n, n)))
    eta_Na <- compute_eta_Na(TW = H$TW, group_sizes = n)
  } else{
    EW <- sum((N / n) * diag(H$TW)) * A1_eq(X_TS_list)
    second_order <- A4(X_TS_list, TW = H$TW, subsampling = subsampling, B = B)
    third_order <- C5star(X_list, TW = H$TWalt, TS = H$TSalt, B = B)
    Var <- 2 * second_order
    eta_Na <- 1
  }
  
  ### Calculate test statistic
  X_bar <- c(sapply(X_list, rowMeans))
  TM_Xbar <- (H$TMalt %*% X_bar)
  QN <- N * sum(TM_Xbar * TM_Xbar)
  W <- compute_grouped_statistic(QN = QN, EW = EW, variance = Var)
  f <- compute_grouped_df(second_order = second_order, third_order = third_order, design_factor = eta_Na)
  p.value <- compute_grouped_p_value(statistic = W, degrees_of_freedom = f)
  
  ## Output
  list(
    f = f,
    statistic = W,
    H = list(TW = H$TW, TS = H$TS),  # Hypothesis matrices (TW and TS)
    p.value = p.value,
    B = B
  )
}
