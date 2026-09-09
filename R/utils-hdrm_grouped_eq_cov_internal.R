#' @keywords internal
hdrm_grouped_eq_cov_internal <- function(data, group, hypothesis = c("whole", "sub", "interaction"), AM, B, seed){
  
  
  # Temporarily set the seed and restore the previous RNG state on exit
  if (!is.null(seed)) {
    withr::local_seed(seed)
  }
  
  # Determine key variables: N (number of samples), d (dimension), a (number of groups), and n (group sizes)
  N <- ncol(data)  # Number of samples (columns in the data)
  d <- nrow(data)  # Number of dimensions (rows in the data)
  a <- length(table(group))  # Number of groups
  n <- as.integer(table(group))  # Size of each group
  
  ## write each group matrix to list
  data_list <- vector("list", a)
  for(i in 1:a){
    data_list[[i]] = data[, group == i,drop=FALSE]
  }
  
  # Determine the hypothesis matrices based on the hypothesis parameter
  H <- get_hypothesis_mult(hypothesis, AM, a, d)
  
  # Prepare the X_TS matrix by multiplying TSalt with the data matrix
  X_TS <- H$TSalt %*% data
  
  
  X_TS_list <- lapply(data_list, function(X) H$TSalt %*% X)
  # Calculate the first- and second-order trace estimators
  A1 <- A1_eq(X_TS_list)
  A2 <- A2_eq(X_TS_list)
  
  if (
    anyNA(c(A1, A2)) ||
    any(!is.finite(c(A1, A2))) ||
    A1 < 0 ||
    A2 < 0
  ) {
    stop(
      "The equal-covariance trace estimators must be finite and non-negative.",
      call. = FALSE
    )
  }
  
  ### Compute the test statistic
  X_bar <- numeric(a * d)  # Initialize a vector for the group-wise means
  for (i in 1:a) {
    # For each group, calculate the row means across dimensions (using rowMeans as the data is transposed)
    X_bar[1:d + ((i - 1) * d)] <- rowMeans(data[, group == i,drop=FALSE])
  }
  
  # Calculate the expectation values (EW), variance (Var), and the test statistic components
  EW <- sum((N / n) * diag(H$TW)) * A1  # Expectation values
  # tmp = 0
  # for (i in 1:a) {
  #   for(r in 1:a){
  #     tmp = tmp + (H$TW[i, r]^2 * (N^2 / (n[i] * n[r])))  # Accumulate variance terms
  #   }
  # }
  Var <- 2 * A2 * sum((H$TW)^2 * (N^2 / outer(n, n)))  # Variance calculation
  
  if (
    length(Var) != 1L ||
    is.na(Var) ||
    !is.finite(Var) ||
    Var <= 0
  ) {
    stop(
      "The estimated variance of the test statistic must be finite and positive.",
      call. = FALSE
    )
  }
  
  # Calculate C1 only after confirming that the second-order estimate is valid
  C1 <- C1star_eq(X_TS_list, B = B)
  
  TMbar <- H$TMalt %*% X_bar  # Calculate the transformed means using the alternative matrices
  QN <- N * sum(TMbar * TMbar)  # Sum of squared values for QN
  W <- compute_grouped_statistic(
    QN = QN,
    EW = EW,
    variance = Var
  )
  
  # Calculate the known whole-plot factor eta_{N,a}
  eta_Na <- compute_eta_Na(
    TW = H$TW,
    group_sizes = n
  )
  
  # Calculate f and the p-value using the chi-square distribution
  f <- compute_grouped_df(
    second_order = A2,
    third_order = C1,
    design_factor = eta_Na
  )
  p.value <- compute_grouped_p_value(
    statistic = W,
    degrees_of_freedom = f
  )
  
  ## Return the results as a list
  return(
    list(
      f = f,  # The scale factor f
      statistic = W,  # The test statistic W
      tau = 1 / f,  # The tau value (inverse of f)
      H = list(H$TW, H$TS),  # The hypothesis matrices (TW and TS)
      hypothesis = ifelse(is.character(hypothesis), hypothesis[1], "custom"),  # Description of the hypothesis
      p.value = p.value,  # The computed p-value
      dim = list(d = d, N = N),  # Dimensions of the data (d: number of dimensions, N: number of samples)
      groups = list(a = a, table = table(group))  # Group information (a: number of groups, group sizes)
    )
  )
}

