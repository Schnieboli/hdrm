#'@title Intern function to conduct the grouped test
#'
#'@description If all parameter and options are suitable  when execute the
#'  function hdrm_grouped and after preparing the arguments, this function
#'  conducts the test.
#'@param data a list of matrices, each matrix representing a group with subjects
#'  in columns and observations in rows
#'@param H a list containing TW, TS, TM, TWalt, TSalt and TMalt
#'@param subsampling logical value specifying whether the subsampling versions
#'  of all trace estimators should be used
#'@param B a character or numeric value determining the base subsampling budget.
#'  The grouped third-trace estimator uses `a * B` draws.
#'@param seed optional value used to set the random seed for reproducible
#'  computations.
#'@noRd
hdrm_grouped_internal <- function(X_list, H, subsampling, B, seed){

  # Determine the number of samples (N), dimensions (d), groups (a), and group sizes (n)
  a <- length(X_list)
  n <- sapply(X_list, ncol)
  d <- nrow(X_list[[1]])
  N <- sum(N)

  ## multiply data with TSalt
  X_TS_list <- lapply(X_list, function(x) H$TSalt %*% x)
  
  }

  EW <- Exp_Q(X_TS_list, TW = H$TW, subsampling = subsampling, B = B)
  
  A4 <- A4(X_TS_list, TW = H$TW, subsampling = subsampling, B = B)
  C5 <- C5star(data_list, TW = H$TWalt, TS = H$TSalt, B = B)

  ### Calculate test statistic
  X_bar <- c(sapply(X_list, rowMeans))

  # Calculate expectation values (EW), variances (Var), and test statistic components (QN and W)
  Var <- 2 * A4
  TMbar <- (H$TMalt %*% X_bar)
  QN <- N * sum(TMbar * TMbar)
  W <- compute_grouped_statistic(QN = QN, EW = EW, variance = Var)

  # Calculate f and the p-value based on the test statistic
  f <- compute_grouped_df(second_order = A4, third_order = C5)
  p.value <- compute_grouped_p_value(statistic = W, degrees_of_freedom = f)

  ## Output
  list(
    f = f,
    statistic = W,
    tau = 1 / f,  # Inverse of f
    H = list(H$TW, H$TS),  # Hypothesis matrices (TW and TS)
    p.value = p.value,
    dim = list(d = d, N = N),  # Dimensions of the input data
    groups = list(a = a, table = table(group))  # Grouping information
    dim = list(d = d, N = N)  # Dimensions of the input data
  )
}
