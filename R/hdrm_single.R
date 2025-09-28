#' @title Test for one group high dimensional repeated measures
#'
#' @description This function implements the methods outlined in
#'   Pauly et al. (2015).
#'
#' @param data the data for which the test is applied, where two formats are
#' possible:
#' \itemize{
#' \item a matrix where subjects represented by columns.
#' \item a value column from a data frame, where in this case an additional parameter
#' specifying the subjects is necessary.
#' }

#' @param subject optional parameter vector, which gives the allocation of the
#' measurements to the subjects
#' @param hypothesis either "flat" or a quadratic numeric matrix (see Details).
#' @param AM `logical`. Specifying whether an alternative hypothesis matrix
#'  based on Sattler and Rosenbaum (2025) should be used (default: `TRUE`),
#'  which has fewer rows but does not influence any values.
#' @param ... further arguments are currently ignored.
#' @details
#' If `data` is a matrix, the parameter `subject` has no effect. \cr
#' If `data` is a vector, `subject` allocates the measurements to the subjects.
#' `subject` then must be of length `length(data)`.
#'
#' The function can deal with missing values only in the measurements.
#' The test is then performed without the affected subjects. Missing values in
#' `subject` will result in an error.
#'
#' The test outlined in Pauly et al. (2015) is performed for
#' the hypothesis given by `hypothesis`. The `hypothesis` can either be given as
#' a `character` or a `list`. The only legal character is "flat" (default).
#' "flat" stands for \eqn{H_0}: time profile is flat. Other
#' characters will result in an error.
#'
#' Alternatively, `hypothesis` can be the quadratic hypothesis matrix \eqn{T} with
#' the number of rows equal to the number of rows of `data`. \eqn{T} must be
#' a projection matrix, meaning symmetrical and \eqn{T^2 = T}. A matrix that does not
#' match those criteria will result in an error.
#'
#'
#' @returns Returns a list with class "hdrm_single" with the components
#' @return \item{data}{the input data used.}
#' @return \item{f}{the degrees of freedom \eqn{f}.}
#' @return \item{statistic}{the test statistic \eqn{W}.}
#' @return \item{tau}{the convergence parameter \eqn{\tau}.}
#' @return \item{H}{the hypothesis-matrix used.}
#' @return \item{hypothesis}{a character. Will be 'custom' if `hypothesis` was
#'   given as a matrix, otherwise `hypothesis[1]`.}
#' @return \item{p.value}{the p-value of the test statistic.}
#' @return \item{dim}{a named vector giving the dimensions \eqn{d \times N} of `data`.}
#' @return \item{removed.cases}{number of subjects removed for having missing values.}
#'
#' @examples
#' ## Load data set EEG (data frame)
#' data(EEG)
#' # ?EEG
#'
#'
#' ## call test
#' hdrm_single(data = EEG$value,
#'             ## test if time profiles are flat
#'             hypothesis = "flat",
#'             ## if data is a vector, subject must be specified
#'             subject = EEG$subject
#' )
#'
#'
#' ## hypothesis as list
#' hypothesis_matrix <- diag(40) - matrix(1/40, 40, 40)
#'
#' hdrm_single(data = EEG$value,
#'             ## test if time profiles are flat
#'             hypothesis = hypothesis_matrix,
#'             subject = EEG$subject
#' )
#'
#' ## Load data set birthrates (matrix)
#' data(birthrates)
#' # ?birthrates
#'
#' ## call test
#' hdrm_single(data = birthrates,
#'             ## test if time profiles are flat
#'             hypothesis = "flat"
#'             ## if data is a matrix, subject does not need to be specified
#' )
#'
#'
#' ## hypothesis as list
#' hypothesis_matrix <- diag(34) - matrix(1/34, 34, 34)
#'
#' hdrm_single(data = birthrates,
#'             ## equivalent to hypothesis = "flat"
#'             hypothesis = hypothesis_matrix
#' )
#'
#' @references {Pauly, M., Ellenberger, D., & Brunner, E. (2015). Analysis of high-dimensional one group repeated measures designs. Statistics, 49(6), 1243–1261. doi:10.1080/02331888.2015.1050022
#' }
#' @references {Sattler, P., & Rosenbaum, M. (2025). Choice of the hypothesis matrix for using the Anova-type-statistic. Statistics & Probability Letters, 219(110356), 110356. doi:10.1016/j.spl.2025.110356
#' }
#'
#' @export
hdrm_single <- function(data, hypothesis = "flat", subject = NULL, AM = TRUE, ...) {

  # Input validation: Ensure that 'data' is either a vector or a matrix
  if(!is.vector(data) & !is.matrix(data)) stop("data must be a data column or a matrix")

  # Ensure that all entries in 'data' are numeric
  if(!all(sapply(data, is.numeric))) stop("data must contain only numeric entries")

  # Check that 'subject' does not contain NA if 'data' is a vector
  if(is.vector(data) & sum(is.na(subject)) > 0) stop("subject contains NA")

  # Ensure that the length of 'subject' matches the length of 'data' when 'data' is a vector
  if(is.vector(data) & length(subject) != length(data)) stop("length of data vector does not coincide with subject")

  # Warn if 'subject' is provided when 'data' is a matrix (since matrix data doesn't require 'subject')
  if(is.matrix(data) & !is.null(subject)) warning("for data in a matrix format, no additional subject parameter is necessary")

  # Case 1: If 'data' is a matrix
  if(is.matrix(data)) {
    X <- data
    N <- ncol(X)  # Number of columns (subjects)
    d <- nrow(X)  # Number of rows (dimensions)
  }

  # Case 2: If 'data' is a vector and 'subject' is provided
  if(is.vector(data) & length(subject) == length(data) & !any(is.na(subject))) {

    # Create a data frame with 'data' and 'subject' columns
    dframe <- data.frame(value = data, subject = subject)

    # Ensure that all subjects have the same number of dimensions
    if(length(unique(table(dframe$subject))) != 1) stop("All subjects must have the same number of dimensions")

    # Get the number of subjects (N) and number of dimensions (d)
    N <- nlevels(as.factor(dframe$subject))
    d <- unique(table(dframe$subject))

    # Build the matrix X by adding data for each subject
    X <- matrix(NA, d, N)
    subjects <- levels(as.factor(subject))
    for (i in 1:N) {
      X[, i] <- dframe$value[dframe$subject == subjects[i]]
    }
  }

  # Check that the data meets the necessary criteria for the test (specific to the 'hdrm_single' method)
  check_criteria_single(X, hypothesis = hypothesis)

  # Handle missing values by removing rows with NA values from the matrix X
  N_with_NA <- ncol(X)
  X <- t(stats::na.omit(t(X)))  # Transpose, remove missing rows, and transpose back
  N <- ncol(X)  # Update the number of subjects after removing missing values
  if(N_with_NA > N) warning("Subjects with missing values dropped", call. = FALSE)

  # Prepare the hypothesis matrix (TM) based on the hypothesis input
  TM <- NA
  if(is.matrix(hypothesis)) {
    TM <- hypothesis  # If hypothesis is already a matrix, use it directly
  }

  # If hypothesis is a character string, specifically "flat", build the hypothesis matrix
  if(is.character(hypothesis)) {
    if(hypothesis[1] == "flat") {
      TM <- diag(d) - matrix(1/d, d, d)  # Flat hypothesis: diagonal minus equal off-diagonal elements
    } else {
      stop("The only legal character is 'flat'. Other hypotheses can be specified by a matrix.")
    }
  }

  # Check that the hypothesis matrix has correct dimensions (square matrix with dimension d)
  if(!all(dim(TM) == c(d, d))) stop("hypothesis must be a quadratic matrix with the number of rows equal to the number of dimension d.")

  # Check that hypothesis matrix is valid (numeric and no NAs)
  if(any(is.na(TM)) | !is.numeric(TM)) stop("Please specify valid hypothesis.")

  # Check symmetry of TM
  if((mean(t(TM) - TM) >= sqrt(.Machine$double.eps))) {
    warning(paste0("TM is not symmetric (mean difference = ", mean(t(TM) - TM), "). This will likely have a big influence on the test result!"))
  }

  # Check idempotence of TM (TM * TM should equal TM)
  if((mean(TM %*% TM - TM) >= sqrt(.Machine$double.eps))) {
    warning(paste0("TM is not idempotent (mean difference = ", mean(TM %*% TM - TM), "). This will likely have a big influence on the test result!"))
  }

  # Compute the test statistic Q
  if(AM) {
    TM <- MSrootcompact(TM)  # If 'AM' is set to TRUE, apply the MSrootcompact transformation to TM
  }
  XT <- TM %*% X  # Multiply TM with the data matrix X
  Xquer <- rowMeans(XT)  # Compute the mean of each row (dimension)
  Qn <- N * sum(Xquer^2)  # Compute the Q statistic (sum of squared row means)

  # Calculate the trace values used in the test statistic W
  traceNormal <- B0_cpp(XT)  # Call to an external C++ function to compute the trace
  traceSquare <- B2_cpp(XT)  # Another trace calculation
  traceCubic <- B3_cpp(XT)  # Another trace calculation

  # Compute the test statistic W
  W <- (Qn - traceNormal) / sqrt(2 * traceSquare)

  # Estimate the parameter f
  f <- max(1, traceSquare^3 / traceCubic^2)

  # Compute the p-value for the test
  p.value <- max(1 - stats::pchisq(W * sqrt(2 * f) + f, df = f), .Machine$double.eps)

  # Prepare the output
  out <- list(data = X,
            f = f,
            statistic = W,
            tau = 1/f,
            H = TM,
            hypothesis = ifelse(is.character(hypothesis), hypothesis[1], "custom"),
            p.value = p.value,
            dim = c(d = d, N = N),
            removed.cases = N_with_NA - N
  )

# Set the class of the output to 'hdrm_single' for print function
  class(out) <- "hdrm_single"

  # Return the result
  return(out)
}






# Print Function --------------------------------------------------------

#' @method print hdrm_single
#' @export
print.hdrm_single <- function(x, digits = 4,...){
  ## prepare p-value
  p <- round(x$p.value, digits)
  if(p <= 0){
    p = paste0("< ", 10^(-digits))
  }else{
    p = paste0("= ", p)
  }

  # print-output
  cat("\n")
  cat("           One Group Repeated Measure
       \nAnalysis of", x$dim[2], "subjects in", paste0(x$dim[1]), "dimensions:",
      "\nW =", round(x$statistic, digits), " f =", round(x$f, digits), " p.value", p,
      "\nHypothesis type:", x$hypothesis,
      "\nConvergence parameter \u03c4 =", round(x$tau, digits))
  cat("\n")
}
