#' @title Intern function to conduct the grouped test
#'
#' @description If all parameter and options are suitable  when execute the
#' function hdrm_grouped and after preparing the arguments, this function
#' conducts the test.
#' @param data the data for which the test is applied. Subjects are represented
#' by columns of the numeric matrix.
#' @param group a vector specifying the group allocation of the subjects.
#' @param hypothesis either one of `"whole"`, `"sub"`, `"interaction"`,
#' `"identical"`, or `"flat"`, or a named list containing the quadratic
#' matrices `TW` and `TS` (see Details).
#' @param AM binary variable specifying whether an alternative hypothesis
#' matrix based on \insertCite{Sattler2025;textual}{hdrm} should be used.
#' This matrix has fewer rows but does not affect the resulting test statistic.
#' @param subsampling logical value specifying whether the subsampling versions
#' of all trace estimators should be used (see Details).
#' @param B a character or numeric value determining the base subsampling
#' budget. The grouped third-trace estimator uses `a * B` draws.
#' @param seed optional value used to set the random seed for reproducible
#' computations.
#'@return a named list of class "hdrm_grouped" with the components
#'@returns \item{data}{the input data used.}
#'@returns \item{f}{the degrees of freedom \eqn{f}.}
#'@returns \item{tau}{the convergence parameter \eqn{\tau}.}
#'@returns \item{H}{a named list with components `TW` and `TS` that give the
#'  components of the hypothesis matrix.}
#'@returns \item{hypothesis}{a character. Will be "custom" if `hypothesis` is a
#'  list, otherwise `hypothesis[1]`.}
#'@returns \item{p.value}{the \eqn{p}-value of the test statistic.}
#'@returns \item{dim}{a named list with with number of factor levels \eqn{d} and
#'  number of subjects \eqn{N} of `data`.}
#'@returns \item{groups}{a named list with components number of groups `a` and
#'  distribution of groups `table`.}
#'@returns \item{removed.cases}{number of incomplete subjects removed.}
#'@returns \item{subsamples}{evaluated base subsampling budget `B`.}
#' @noRd
hdrm_grouped_internal <- function(data, group, hypothesis = c("whole", "sub", "interaction"), AM, subsampling, B, seed){

  # Temporarily set the seed and restore the previous RNG state on exit
  if (!is.null(seed)) {
    withr::local_seed(seed)
  }

  # Determine the number of samples (N), dimensions (d), groups (a), and group sizes (n)
  N <- ncol(data)
  d <- nrow(data)
  a <- length(table(group))
  n <- as.integer(table(group))  # Number of samples in each group

  # Get the hypothesis matrices based on the provided hypothesis
  H <- get_hypothesis_mult(hypothesis, a, d)
  TW <- H$TW
  TS <- H$TS

  # Alternative matrices for the setting when AM = TRUE
  TWalt <- TW
  TSalt <- TS
  if(AM == 1) {
    TWalt <- MSrootcompact(TW)   # Apply a transformation to TW if AM = 1
    TSalt <- MSrootcompact(TS)  # Apply a transformation to TS if AM = 1
  }


  # Kronecker product of hypothesis matrices
  TM <- kronecker(TW, TS)
  TMalt <- kronecker(TWalt, TSalt)

  # Prepare the transformed data matrix (X_TS)
  X_TS <- TSalt %*% data

  # Initialize vectors for estimators
  A1 <- A3 <- numeric(a)
  A2 <- matrix(0, a, a)
  C5 <- numeric(1)

  # Estimate A1 and A3 using the appropriate method based on subsampling
  for (i in 1:a) {
    if(subsampling){
      # Use bootstrap sampling if subsampling is true
      A1[i] <- A1star_cpp(X = X_TS[, group == i,drop=FALSE], B)
      A3[i] <- A3star_cpp(X = X_TS[, group == i,drop=FALSE], B)
    } else {
      # Use the original method without subsampling
      A1[i] <- A1_cpp(mat = X_TS[, group == i,drop=FALSE])
      A3[i] <- A3_cpp(mat = X_TS[, group == i,drop=FALSE], Part6 = sum(rowMeans(X_TS[, group == i,drop=FALSE])^2))
    }
  }

  # Estimate A2 for pairwise group comparisons
  for (i in 1:(a-1)) {
    for(r in (i+1):a){
      if(subsampling){
        A2[i, r] <- A2star_cpp(X = X_TS[, group == i,drop=FALSE], Y = X_TS[, group == r,drop=FALSE], B)
      } else {
        A2[i, r] <- A2(X = X_TS[, group == i,drop=FALSE], Y = X_TS[, group == r,drop=FALSE])
      }
    }
  }

  trace_estimates <- c(
    A1,
    A2[upper.tri(A2)],
    A3
  )

  if (
    anyNA(trace_estimates) ||
    any(!is.finite(trace_estimates)) ||
    any(trace_estimates < 0)
  ) {
    stop(
      "The grouped trace estimators must be finite and non-negative.",
      call. = FALSE
    )
  }

  # Estimate A4 using A1, A2, and A3
  temp1 <- temp2 <- 0
  for (i in 1:a) {
    temp1 <- temp1 + ((N/n[i])^2 * TW[i, i]^2 * A3[i])
  }
  for (i in 1:(a-1)) {
    for(r in (i+1):a){
      temp2 <- temp2 + ( (N^2 / (n[i]*n[r])) * TW[i, r]^2 * A2[i, r])
    }
  }
  A4 <- temp1 + 2 * temp2  # Combine terms to get A4

  if (
    length(A4) != 1L ||
    is.na(A4) ||
    !is.finite(A4) ||
    A4 <= 0
  ) {
    stop(
      "The estimated variance of the test statistic must be finite and positive.",
      call. = FALSE
    )
  }

  # Calculate C5 only after confirming that the second-order estimate is valid
  C5 <- C5star_cpp(
    X = data,
    group = group,
    TW = TWalt,
    TS = TSalt,
    B = B
  )

  ### Calculate test statistic
  X_bar <- numeric(a * d)
  for (i in 1:a) {
    # Calculate row means for each group (since 'data' is not transposed)
    X_bar[1:d + ((i - 1) * d)] <- rowMeans(data[, group == i,drop=FALSE])
  }

  # Calculate expectation values (EW), variances (Var), and test statistic components (QN and W)
  EW <- sum((N/n) * diag(TW) * A1)
  Var <- 2 * A4
  TMbar <- (TMalt %*% X_bar)
  QN <- N * sum(TMbar * TMbar)
  W <- compute_grouped_statistic(
    QN = QN,
    EW = EW,
    variance = Var
  )

  # Calculate f and the p-value based on the test statistic
  f <- compute_grouped_df(
    second_order = A4,
    third_order = C5
  )
  p.value <- compute_grouped_p_value(
    statistic = W,
    degrees_of_freedom = f
  )

  ## Output
  L <- list(
    f = f,
    statistic = W,
    tau = 1 / f,  # Inverse of f
    H = H,  # Hypothesis matrices (TW and TS)
    hypothesis = ifelse(is.character(hypothesis), hypothesis[1], "custom"),  # Description of the hypothesis
    p.value = p.value,
    dim = list(d = d, N = N),  # Dimensions of the input data
    groups = list(a = a, table = table(group))  # Grouping information
  )

  # Assign class 'hdrm' to the result
  class(L) <- c("hdrm")

  return(L)  # Return the results as a list
}





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

  # Determine the hypothesis matrices based on the hypothesis parameter
  H <- get_hypothesis_mult(hypothesis, a, d)
  TW <- H$TW  # Matrix TW from the hypothesis
  TS <- H$TS  # Matrix TS from the hypothesis

  # Modify TW and TS if AM is set to 1 (apply root compact transformation)
  TWalt <- TW
  TSalt <- TS
  if(AM == 1){
    TWalt <- MSrootcompact(TW)
    TSalt <- MSrootcompact(TS)
  }

  # Create the Kronecker product of TW and TS, and their alternative versions if AM = 1
  TM <- kronecker(TW, TS)
  TMalt <- kronecker(TWalt, TSalt)

  # Prepare the X_TS matrix by multiplying TSalt with the data matrix
  X_TS <- TSalt %*% data

  # Calculate the first- and second-order trace estimators
  A1 <- make_A1_eq(X = X_TS, group = group)
  A2 <- make_A2_eq(X = X_TS, group = group)

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
  EW <- sum((N / n) * diag(TW)) * A1  # Expectation values
  tmp = 0
  for (i in 1:a) {
    for(r in 1:a){
      tmp = tmp + (TW[i, r]^2 * (N^2 / (n[i] * n[r])))  # Accumulate variance terms
    }
  }
  Var <- 2 * A2 * tmp  # Variance calculation

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
  C1 <- make_C1_star_eq(
    X = X_TS,
    group = group,
    B = B
  )

  TMbar <- TMalt %*% X_bar  # Calculate the transformed means using the alternative matrices
  QN <- N * sum(TMbar * TMbar)  # Sum of squared values for QN
  W <- compute_grouped_statistic(
    QN = QN,
    EW = EW,
    variance = Var
  )

  # Calculate the known whole-plot factor eta_{N,a}
  eta_Na <- compute_eta_Na(
    TW = TW,
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
      H = H,  # The hypothesis matrices (TW and TS)
      hypothesis = ifelse(is.character(hypothesis), hypothesis[1], "custom"),  # Description of the hypothesis
      p.value = p.value,  # The computed p-value
      dim = list(d = d, N = N),  # Dimensions of the data (d: number of dimensions, N: number of samples)
      groups = list(a = a, table = table(group))  # Group information (a: number of groups, group sizes)
    )
  )
}

