#' @title Build the hypothesis matrix for a multifactorial setting
#' @description Since there are two internal functions, for better clarity we create a function
#' that extracts the hypothesis matrices. This way, changes don't need to be made
#' in both internal functions.
#' @param hypothesis - either a character or a list with components TW and TS
#' @param a - natural number, representing the number of groups
#' @param d - natural number, representing the dimension of the observation vectors
#
#' @returns Returns a list  with the components
#' @return \item{TW}{wholeplot hypothesis matrix TW}
#' @return \item{TS}{subplot hypothesis matrix TS}
#' @keywords internal
get_hypothesis_mult <- function(hypothesis, a, d){

  # Check if 'hypothesis' is either a list or a character
  stopifnot(is.list(hypothesis) | is.character(hypothesis))

  # Case 1: Hypothesis is a list
  if(is.list(hypothesis)){
    # Check for required elements (TW and TS) in the list
    if(is.null(hypothesis$TW)) stop("No list entry with name 'TW' in hypothesis", call. = FALSE)
    if(is.null(hypothesis$TS)) stop("No list entry with name 'TS' in hypothesis", call. = FALSE)

    # Extract the matrices from the list
    TW <- hypothesis$TW
    TS <- hypothesis$TS

    # Check the dimensions of TW and TS
    if(any(dim(TW) != c(a, a))) stop(paste0("TW should be a square matrix with ", a, " rows"))
    if(any(dim(TS) != c(d, d))) stop(paste0("TS should be a square matrix with ", d, " rows"))

    # Check if TW and TS are symmetric matrices
    if(!isSymmetric.matrix(TW) | !isSymmetric.matrix(TS)) stop("TW and TS must be symmetric")

    ############# Check properties of TW
    # Check if TW is symmetric
    if((mean(t(TW) - TW) >= sqrt(.Machine$double.eps)))
      warning(paste0("TW is not symmetric (mean difference = ", mean(t(TW) - TW), "). This could significantly affect the test result!"))

    # Check if TW is idempotent
    if((mean(TW %*% TW - TW) >= sqrt(.Machine$double.eps)))
      warning(paste0("TW is not idempotent (mean difference = ", mean(TW %*% TW - TW), "). This could significantly affect the test result!"))

    ############# Check properties of TS
    # Check if TS is symmetric
    if((mean(t(TS) - TS) >= sqrt(.Machine$double.eps)))
      warning(paste0("TS is not symmetric (mean difference = ", mean(t(TS) - TS), "). This could significantly affect the test result!"))

    # Check if TS is idempotent
    if((mean(TS %*% TS - TS) >= sqrt(.Machine$double.eps)))
      warning(paste0("TS is not idempotent (mean difference = ", mean(TS %*% TS - TS), "). This could significantly affect the test result!"))
  }

  # Case 2: Hypothesis is a character
  if(is.character(hypothesis)){
    # Ensure that only one hypothesis is provided
    if(length(hypothesis) > 1) warning("Length of hypothesis is > 1, only first element will be used.", call. = FALSE)

    # Check if the provided hypothesis is valid
    if(!(hypothesis[1] %in% c("whole", "sub", "interaction", "identical", "flat")))
      stop("hypothesis must be one of c('whole', 'sub', 'interaction', 'identical', 'flat') or a list with components TW and TS")

    # Build hypothesis matrices based on the provided hypothesis
    if(hypothesis[1] == "whole"){
      TW <- diag(a) - matrix(1/a, a, a)
      TS <- matrix(1, d, d) / d
    }
    if(hypothesis[1] == "sub"){
      TW <- matrix(1, a, a) / a
      TS <- diag(d) - matrix(1/d, d, d)
    }
    if(hypothesis[1] == "interaction"){
      TW <- diag(a) - matrix(1, a, a) / a
      TS <- diag(d) - matrix(1/d, d, d)
    }
    if(hypothesis[1] == "identical"){
      TW <- diag(a) - matrix(1/a, a, a)
      TS <- diag(d)
    }
    if(hypothesis[1] == "flat"){
      TW <- diag(a)
      TS <- diag(d) - matrix(1/d, d, d)
    }
  }

  ### Check if matrices are numeric and do not contain NAs
  if(!is.matrix(TW)) stop("TW must be a matrix")
  if(!is.numeric(TW)) stop("TW must be numeric")
  if(any(is.na(TW))) stop("TW must not contain NAs")

  if(!is.matrix(TS)) stop("TS must be a matrix")
  if(!is.numeric(TS)) stop("TS must be numeric")
  if(any(is.na(TS))) stop("TS must not contain NAs")

  # Return the hypothesis matrices as a list
  return(list(TW = TW, TS = TS))
}

#' @title Checks whether everything is as expected for the grouped test
#' @description This function ensures that all criteria are checked in one place,
# making it easier to maintain. If new criteria arise, they only need to be
# added in one location, rather than in multiple places.
#' @param X data matrix which is investigated
#' @param group vector containing the allocation to the
#' subjects of the measurement
#' @param hypothesis a hypothesis given through a character to specify a predefined
#' standard hypothesis or by a list
#' @param subsampling a logical number specifying whether subsampling versions of
#' estimators should be used if existing
#' @param reps a `string` determines the number of subsamples used by the
#' subsampling trace estimators
#' @return an error message if any condition is violated
#
#' @keywords internal
check_criteria_grouped <- function(X, group, hypothesis, reps, subsampling){

  n <- as.integer(table(group))  # Group sizes
  a <- length(n)  # Number of groups
  d <- nrow(X)  # Number of dimensions
  N <- ncol(X)  # Number of samples

  # Check if at least two groups are present
  if(a < 2) stop("At least two groups are needed. Did you mean to call hdrm_single_...?", call. = FALSE)

  # Check if at least two dimensions are present
  if(d < 2) stop("At least 2 dimensions are needed", call. = FALSE)

  # Ensure all group sizes are >= 6
  if(any(n < 6)) stop("All group sizes must be >= 6", call. = FALSE)

  # Check if the length of the group vector matches the number of columns in X
  if(length(group) != N) stop("Length(group) must be equal to ncol(data)", call. = FALSE)

  # Ensure that 'hypothesis' is either a character or a list
  if(!is.character(hypothesis) & !is.list(hypothesis)) stop("hypothesis must be a character or a list", call. = FALSE)

  # Ensure that subsampling is a logical value and has only one element
  if(length(subsampling) > 1) warning("Length(subsampling) > 1. Only the first element will be used", call. = FALSE)
  if(!(subsampling[1] %in% c(1, 0, TRUE, FALSE, T, F))) stop("subsampling must be logical", call. = FALSE)

  # Ensure the number of repetitions is a positive integer
  if(reps < 0) stop("Repetitions must be >= 0", call. = FALSE)
  if(!is.numeric(reps)) stop("Repetitions should be a number", call. = FALSE)

}



################################################################################
#######################   Trace estimators   ###################################
################################################################################

## all trace estimators expect subject in columns and factor levels in rows


#' @keywords internal
A2 <- function(X, Y){
  nX <- ncol(X)
  nY <- ncol(Y)
  PX <- diag(1, nX, nX) - matrix(data = 1 / nX, # Zentrierungsmatrix fuer Gruppe X
                                 nrow = nX,
                                 ncol = nX)
  PY <- diag(1, nY, nY) - matrix(data = 1 / nY, # Zentrierungsmatrix fuer Gruppe Y
                                 nrow = nY,
                                 ncol = nY)
  MX <- tcrossprod(PX, X)
  MY <- tcrossprod(Y, PY)
  MXY <-  MX %*% MY # Seite 38 im Paper
  EBSchaetzer = sum(MXY^2)/((nX-1)*(nY-1))
  return(EBSchaetzer)
}




#' @keywords internal
C5star_cpp <- function(X, group, TW, TS, B){ # ist zu gross; ich glaube es liegt daran, dass die interne Funktion nicht die Gruppen beachtet (obwohl sie das mMn sollte). Denn wenn ich bei der R-Funktion die Gruppenzuweisung weglasse, dann kommen Ergebnisse in aehnlicher Groessenordnung :(

    stopifnot(length(group) == ncol(X))

  a <- length(table(group))
  d <- nrow(X)
  N <- ncol(X)
  n <- as.integer(table(group))
  stopifnot(all(n >= 6))

  ind <- cumsum(c(1, n)) # zum indexen der Gruppen, ausnutzen, dass X sortiert ist!

  Y <- matrix(0, nrow(TW)*nrow(TS), N)
  for(i in 1:a){ # Teil von X der i-ten Gruppe mit dem entsprechenden Teil der Hypothesenmatrix und den Vorfaktoren multiplizieren
    Y[, ind[i]:(ind[i+1]-1)] <- kronecker(TW[, i], TS%*%(X[, ind[i]:(ind[i+1]-1)] * sqrt(N/n[i])))
  }
  # cpp-teil aufrufen
  return(C5star_cpp_internal(X = Y, group = group, B = B, n = n))
}

## make_A1_eq - function that combines the calculators A1_eq for all groups.
##
## Input:
##      X - a numeric matrix, where cols represent the individuals of one group
##          and rows represent dimensions
##  group - an integer vector, specifies the group membership
##
## Output: numeric; gives the estimator A1 for equal covariance matrices
make_A1_eq <- function(X, group){
  ns <- unname(table(group))
  a <- length(ns)
  out <- 0.0
  prefactor <- 0
  for(i in 1:a){
    n <- sum
    out <- out + A1_eq_cpp(X[, group == i,drop=FALSE])
    prefactor <- prefactor + (ns[i]* (ns[i] - 1))
  }
  return(out/(prefactor))
}

## make_A2_eq - function that combines the calculators A2_eq for all groups.
##
## Input:
##      X - a numeric matrix, where cols represent the individuals of one group
##          and rows represent dimensions
##  group - an integer vector, specifies the group membership
##
## Output: numeric; gives the estimator A2 for equal covariances
make_A2_eq <- function(X, group){
  ns <- unname(table(group))
  a <- length(ns)
  out <- 0.0
  prefactor <- 0
  for(i in 1:a){
    n <- sum
    out <- out + A2_eq_cpp(X[, group == i,drop=FALSE])
    prefactor <- prefactor + (24* choose(ns[i], 4))
  }
  return(out/(prefactor))
}
## make_C1_eq - function that combines the calculations C1_eq for all groups.
##
## Input:
##      X - a numeric matrix, where cols represent the individuals of one group
##          and rows represent dimensions
##  group - an integer vector, specifies the group membership
##
## Output: numeric; gives the estimator A2 for equal covariances
make_C1_eq <- function(X, group){
  ns <- unname(table(group))
  a <- length(ns)
  out <- 0.0
  prefactor <- 0
  for(i in 1:a){
    n <- sum
    out <- out + C1_eq_cpp(X[, group == i,drop=FALSE])
    prefactor <- prefactor + (choose(ns[i], 6))
  }
  return(out/(prefactor * 5760)) ## 5760 = 6! * 8
}


make_C1_star_eq <- function(X, B, group){
  ns <- unname(table(group))
  a <- length(ns)
  out <- 0.0
  prefactor <- 0
  for(i in 1:a){
    n <- sum
    out <- out + C1_star_eq_cpp(X[, group == i,drop=FALSE], B)
  }
  return(out/(8 * a * B)) ## 8 * a*B considering B is equal for all groups
}

#' @title Checks whether everything is as expected for the single test
#' @description This function ensures that all criteria are checked in one place,
# making it easier to maintain. If new criteria arise, they only need to be
# added in one location, rather than in multiple places.
#' @param X data matrix which is investigated
#' @param hypothesis a hypothesis given through a character to specify a predefined
#' standard hypothesis or by a list
#' @return an error message if any condition is violated
#' @keywords internal
check_criteria_single <- function(X, hypothesis){

  d <- nrow(X)
  N <- ncol(X)

  if(d < 2) stop("at least 2 dimensions are needed", call. = FALSE)
  if(!is.character(hypothesis) & !is.matrix(hypothesis)) stop("hypothesis must be a character or a matrix", call. = FALSE)
  if(is.vector(hypothesis) & length(hypothesis) > 1) warning("length of hypothesis is > 1, only first element used.", call. = FALSE)

}


#' @title Compact root of a matrix
#'
#' @description calculates the compact root of a matrix, as it was introduced in
#' Sattler and Rosenbaum (2025).
#'
#' @param H a matrix which compact root should be determined
#' @return a compact root of the original matrix
#' @keywords internal
MSrootcompact<- function(H){
  if(length(H) == 1){
    MSroot <- matrix(sqrt(H), 1, 1)
    }else {
    r <- qr(H)$rank
    SVD <- svd(H)
    MSroot <-  diag(sqrt(SVD$d[1:r]), r, r) %*% t((SVD$u)[, 1:r])
  }
  return(MSroot)
}
