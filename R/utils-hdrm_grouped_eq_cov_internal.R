#'@title Internal function to conduct the grouped test with equal covariance
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
hdrm_grouped_eq_cov_internal <- function(X_list, H, B, seed){
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
  
  # Calculate the first- and second-order trace estimators
  A1 <- A1_eq(X_TS_list)
  A2 <- A2_eq(X_TS_list)
  
  ### Compute the test statistic
  X_bar <- c(sapply(X_list, rowMeans))
  
  # Calculate the expectation values (EW), variance (Var), and the test statistic components
  EW <- sum((N / n) * diag(H$TW)) * A1  # Expectation values
  Var <- 2 * A2 * sum((H$TW)^2 * (N^2 / outer(n, n)))  # Variance calculation
  
  # Calculate C1 only after confirming that the second-order estimate is valid
  C1 <- C1star_eq(X_TS_list, B = B)
  
  TMbar <- H$TMalt %*% X_bar  # Calculate the transformed means using the alternative matrices
  QN <- N * sum(TMbar * TMbar)  # Sum of squared values for QN
  W <- compute_grouped_statistic(QN = QN, EW = EW, variance = Var)
  
  # Calculate the known whole-plot factor eta_{N,a}
  eta_Na <- compute_eta_Na(TW = H$TW, group_sizes = n)
  
  # Calculate f and the p-value using the chi-square distribution
  f <- compute_grouped_df(
    second_order = A2,
    third_order = C1,
    design_factor = eta_Na
  )
  p.value <- compute_grouped_p_value(statistic = W, degrees_of_freedom = f)
  
  ## Return the results as a list
    list(
      f = f,  # The scale factor f
      statistic = W,  # The test statistic W
      tau = 1 / f,  # The tau value (inverse of f)
      H = list(H$TW, H$TS),  # The hypothesis matrices (TW and TS)
      p.value = p.value,  # The computed p-value
      dim = list(d = d, N = N)  # Dimensions of the data (d: number of dimensions, N: number of samples)
    )
}
