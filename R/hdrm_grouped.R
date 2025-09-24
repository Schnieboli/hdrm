#' @title Internal function to conduct the grouped test
#'
#' @description If all parameter and options are suitable when executing the function
#' hdrm_grouped and after preparing the arguments, this function conducts the test.
#' @param data a matrix where cols represent subjects and rows represent factor levels
#'
#' @param hypothesis either one of "whole", "sub" "interaction", "identical" or
#'  "flat" or a named list with quadratic matrices `TW` and `TS`
#' @param AM binary variable, specifying whether an alternative hypothesis
#'  matrix based on Sattler and Rosenbaum (2025) should be used, which has
#'  fewer rows but does not influence any values. The usage of this matrix is
#'  the predefined setting.
#' @param group parameter vector, which gives the allocation of the measurements
#' to the single groups
#' @param subject optional parameter vector, which gives the allocation to the
#' subjects of the measurement.
#' @param B a `string` specifying a function of the number of subjects \eqn{N}.
#'  Determines the number of subsamples used by the subsampling trace
#'  estimators.
#' @param subsampling `logical`. Specifying whether the subsampling versions for
#'  all trace estimators should be used (see details).
#' @param seed an optional natural number as seed, if it should be set for
#' reproducibility.
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
#'@returns \item{subsamples}{number of subsamples used for subsampling
#'  estimators.}
#' @keywords internal
hdrm_grouped_internal <- function(data, group, hypothesis = c("whole", "sub", "interaction"), AM, subsampling, B, seed){

  # Set seed if provided
  if(!is.null(seed)){
    if(exists(".Random.seed")){
      old_seed <- .Random.seed
      on.exit({ .Random.seed <<- old_seed })  # Ensure seed is restored after function exit
    }
    if(!exists(".Random.seed")){
      on.exit({ set.seed(NULL) })  # Default seed restoration
    }
    set.seed(seed)  # Set the random seed
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
  if(AM) {
    TWalt <- MSrootcompact(TW)   # Apply a transformation to TW if AM is TRUE
    TSalt <- MSrootcompact(TS)  # Apply a transformation to TS if AM is TRUE
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
        A2[i, r] <- A2star_cpp(X = X_TS[, group == i, drop = FALSE], Y = X_TS[, group == r, drop = FALSE], B)
      } else {
        A2[i, r] <- A2(X = X_TS[, group == i, drop = FALSE], Y = X_TS[, group == r, drop = FALSE])
      }
    }
  }

  # Calculate C5 estimator using the alternative matrices
  C5 <- C5star_cpp(X = data, group = group, TW = TWalt, TS = TSalt, B = B)
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

  ### Calculate test statistic
  X_bar <- numeric(a * d)
  for (i in 1:a) {
    # Calculate row means for each group (since 'data' is not transposed)
    X_bar[1:d + ((i - 1) * d)] <- rowMeans(data[, group == i, drop = FALSE])
  }

  # Calculate expectation values (EW), variances (Var), and test statistic components (QN and W)
  EW <- sum((N/n) * diag(TW) * A1)
  Var <- 2 * A4
  TMbar <- (TMalt %*% X_bar)
  QN <- N * sum(TMbar * TMbar)
  W <- as.numeric((QN - EW) / sqrt(Var))

  # Calculate f and the p-value based on the test statistic
  f <- max(1, as.numeric(A4^3 / C5^2))  # Calculate f-statistic
  p.value <- max(1 - stats::pchisq(W * sqrt(2 * f) + f, df = f), .Machine$double.eps)

  # Return the results as a list
  return(
    list(
      f = f,
      statistic = W,
      tau = 1 / f,
      # Inverse of f
      H = H,
      # Hypothesis matrices (TW and TS)
      hypothesis = ifelse(is.character(hypothesis), hypothesis[1], "custom"),
      # Description of the hypothesis
      p.value = p.value,
      dim = list(d = d, N = N),
      # Dimensions of the input data
      groups = list(a = a, table = table(group))  # Grouping information
    )
  )
}




#' @title Internal function to conduct the grouped test for equal covariances
#'
#' @description If all parameter and options are suitable  when executing the function
#' hdrm_grouped and after preparing the arguments, this function conducts the test,
#' if equal covariance matrices are assumed.
#' @param data the data for which the test is applied, where two formats are
#' possible:
#' -A matrix where subjects represented by columns.
#' -A value column from a data frame, where in this case an additional parameter
#' specifying the subjects is necessary.
#'@param hypothesis either one of "whole", "sub" "interaction", "identical" or
#'  "flat" or a named list with quadratic matrices `TW` and `TS` (see details).
#'@param AM binary variable, specifying whether an alternative hypothesis
#'  matrix based on
#'  Sattler and Rosenbaum (2025) should be used, which has
#'  fewer rows but does not influence any values. The usage of this matrix is
#'  the predefined setting.
#'@param group parameter vector, which gives the allocation of the measurements
#' to the single groups
#'@param subject optional parameter vector, which gives the allocation to the
#' subjects of the measurement.
#'@param B a `string` specifying a function of the number of subjects \eqn{N}.
#'  Determines the number of subsamples used by the subsampling trace
#'  estimators.
#'@param subsampling `logical`. Specifying whether the subsampling versions for
#'  all trace estimators should be used (see details).
#'@param seed an optional natural number as seed, if it should be set for
#'reproducibility.
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
#'@returns \item{subsamples}{number of subsamples used for subsampling
#'  estimators.}
#' @keywords internal
#' @keywords internal
hdrm_grouped_eq_cov_internal <- function(data, group, hypothesis = c("whole", "sub", "interaction"), AM, B, seed){

  # Set the random seed for reproducibility, if provided
  if(!is.null(seed)){
    if(exists(".Random.seed")){
      old_seed <- .Random.seed
      on.exit({ .Random.seed <<- old_seed })  # Restore the original seed upon exit
    }
    if(!exists(".Random.seed")){
      on.exit({ set.seed(NULL) })  # If no previous seed, reset the seed to NULL
    }
    set.seed(seed)  # Set the provided seed
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

  # Calculate estimators A1, A2, and C1 using the provided functions
  A1 <- make_A1_eq(X = X_TS, group = group)
  A2 <- make_A2_eq(X = X_TS, group = group)
  C1 <- make_C1_star_eq(X = X_TS, group = group, B = B)

  ### Compute the test statistic
  X_bar <- numeric(a * d)  # Initialize a vector for the group-wise means
  for (i in 1:a) {
    # For each group, calculate the row means across dimensions (using rowMeans as the data is transposed)
    X_bar[1:d + ((i - 1) * d)] <- rowMeans(data[, group == i, drop = FALSE])
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
  TMbar <- TMalt %*% X_bar  # Calculate the transformed means using the alternative matrices
  QN <- N * sum(TMbar * TMbar)  # Sum of squared values for QN
  W <- as.numeric((QN - EW) / sqrt(Var))  # Compute the test statistic W

  # Calculate f (scale factor) and the p-value using the chi-square distribution
  f <- max(1, as.numeric(A2^3 / C1^2))  # Scale factor f, ensuring it's at least 1
  p.value <- max(1 - stats::pchisq(W * sqrt(2 * f) + f, df = f), .Machine$double.eps)  # p-value using chi-square distribution

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

