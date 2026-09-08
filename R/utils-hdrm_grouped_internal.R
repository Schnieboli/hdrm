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
  
  ## write each group matrix to list
  X_TS_list <- vector("list", a)
  for(i in 1:a){
    X_TS_list[[i]] = X_TS[, group == i,drop=FALSE]
  }

  EW <- Exp_Q(X_TS_list, TW = TW, subsampling = subsampling, B = B)
  
  # if (
  #   anyNA(EW) ||
  #   any(!is.finite(EW)) ||
  #   any(EW < 0)
  # ) {
  #   stop(
  #     "The grouped trace estimators must be finite and non-negative.",
  #     call. = FALSE
  #   )
  # }
  
  A4 <- A4(X_TS_list, TW = TW, subsampling = subsampling, B = B)
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
