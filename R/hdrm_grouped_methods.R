#' @title Test for multiple group high dimensional repeated measures
#' @description This function implements the methods outlined in
#'  Sattler and Pauly (2018).
#' @param data the data for which the test is applied, where two formats are
#' possible:
#' \itemize{
#'    \item A matrix where subjects are represented by columns
#'    \item A value column from a data frame. In this case an additional parameter
#' specifying the subjects and their groups is necessary.
#'    }
#' @param hypothesis either one of "whole", "sub" "interaction", "identical" or
#'  "flat" or a named list with quadratic matrices `TW` and `TS` (see details).
#' @param AM `logical`. If `TRUE` (default), an alternative hypothesis
#'  matrix based on Sattler and Rosenbaum (2025) is used, which has
#'  fewer rows but does not influence any values.
#' @param group parameter vector, which allocates the measurements to the single
#' groups (see Details).
#' @param subject optional parameter vector of `length(data)`, which gives the allocation of the
#'  measurements to the subjects (see Details).
#' @param B a `string` of the pattern "1000*N" specifying a function of the
#'  number of subjects \eqn{N}.
#'  Determines the number of subsamples used by the subsampling trace estimators.
#' @param cov.equal `logical`. If `TRUE` (default = `FALSE`), the test for equal group covariances is performed as outlined in Sattler (2021).
#' @param subsampling `logical`. If `TRUE` (default = `FALSE`), the subsampling versions for
#'  all trace estimators are used (see details).
#'@param seed an optional natural number as seed, if it should be set for reproducibility. Default
#'  is `NULL`, which means no seed is set.
#'@param ... further arguments. Currently ignored.
#'
#'@details
#' If `data` is a matrix, the parameter `group` allocates the columns to
#' the groups. Therefore `group` must be of length `ncol(data)`. \cr
#'  If `data` is a vector, `subject` allocates the measurements to the subjects
#'  and `group` then allocates the subjects to a group. `subject` and `group`
#'  then must be of length `length(data)`.
#'
#' The function can deal with missing values only in the measurements.
#'  The test is then performed without the affected subjects. Missing values in
#'  `group` or `subject` will result in an error.
#'
#'The test outlined in Sattler and Pauly (2018) is performed for
#'the hypothesis \eqn{(\bm T_W \otimes \bm T_S)\bm\mu = \bm 0}, with \eqn{\bm
#'T_W}, \eqn{\bm T_S} given by `hypothesis`. The `hypothesis` can either be
#'given as a `character` or a `list`. Legal characters are "whole" (default),
#'"sub", "interaction", "identical" and "flat". The first three are given by (as
#'outlined in Sattler and Pauly (2018)):
#' \itemize{
#'  \item \eqn{H_0^{\overline W}\text{: } \left(\bm P_a \otimes \frac{1}{d}\bm J_d \right)\bm\mu=\bm 0}
#' (no main effect in whole plot factor group)
#'  \item  \eqn{H_0^{\overline S}\text{: } \left(\frac{1}{a}\bm J_a \otimes \bm P_d \right)\bm\mu=\bm 0}
#' (no main effect in sub plot factor dimension)
#'  \item \eqn{H_0^{WS}\text{: } \left(\bm P_a \otimes \bm P_d \right)\bm\mu=\bm 0}
#' (no interaction effect of whole plot and sub plot factor).
#'}
#'Additionally, "identical" and "flat" specify the hypotheses
#' \itemize{
#' \item  \eqn{H_0^{W}\text{: } \left(\bm P_a \otimes \bm I_d \right)\bm\mu=\bm 0}
#' (vectors of expected values are equal between all groups)
#' \item \eqn{H_0^{S}\text{: } \left(\bm I_a \otimes \bm P_d \right)\bm\mu=\bm 0}
#' (flat profile in all groups).
#' }
#'Other characters will result in an error.
#'
#'Alternatively, `hypothesis` can be a named list with quadratic
#' matrices `TW` for the wholeplot-part and `TS` for the subplot-part of \eqn{\bm T
#' = \bm T_W \otimes \bm T_S}. If `hypothesis` is a list, `TW` and `TS` must be
#'projection matrices, meaning symmetrical and \eqn{\bm T_W^2 = \bm T_W},
#'\eqn{\bm T_S^2 = \bm T_S}. Differences up to the tolerance in `all.equal` are
#'ignored, while bigger differences will result in a warning. Also, the number of
#'rows and columns of `TW` must be equal to the number of groups and the number
#'of rows and columns of `TS` must be equal to the number of factor levels.
#'Lists that do not match those criteria will result in an error.
#'
#'Note that for `subsampling = FALSE` results are still seed dependent,
#'because the computational heaviest trace estimator is always calculated using
#'subsamples. Also, depending on the choice of `B`, the non-subsampling versions
#'might not be faster for small data.
#'
#'That also means, that even for `subsampling = FALSE` the `p.value` depends
#'heavily on the choice of `B`, as it depends on the degrees of freedom `f`,
#'which are always estimated by the subsampling version. The test statistic
#'\eqn{W} is only seed dependent for `subsampling = TRUE`
#'
#'The number of subsamples `B` can also be a numeric value. However, this is not
#'advised, as `B` should be a function of \eqn{N} that goes to \eqn{\infty} for
#'\eqn{N\to\infty}.
#'
#'
#'@return a named list of class "hdrm_grouped" with the components
#'@returns \item{data}{the input data used.}
#'@returns \item{f}{the degrees of freedom \eqn{f}.}
#'@returns \item{tau}{the convergence parameter \eqn{\tau}.}
#'@returns \item{H}{a named list with components `TW` and `TS` that give the
#'  components of the hypothesis matrix.}
#'@returns \item{hypothesis}{a character. Will be 'custom' if `hypothesis` is a
#'  list, otherwise `hypothesis[1]`.}
#'@returns \item{p.value}{the \eqn{p}-value of the test statistic.}
#'@returns \item{dim}{a named list with with number of factor levels \eqn{d} and
#'  number of subjects \eqn{N} of `data`.}
#'@returns \item{groups}{a named list with components number of groups `a` and
#'  distribution of groups `table`.}
#'@returns \item{removed.cases}{number of incomplete subjects removed.}
#'@returns \item{subsamples}{number of subsamples used for subsampling
#'  estimators.}
#'@returns \item{cov.equal}{`logical`. Indicates whether the test for equal covariances was used.}
#'
#' @references {
#'  Sattler, P., & Pauly, M. (2018). Inference for high-dimensional split-plot-designs: A unified approach for small to large numbers of factor levels. Electronic Journal of Statistics, 12(2), 2743–2805. doi:10.1214/18-ejs1465
#' }
#' @references {
#'  Sattler, P. (2021). A comprehensive treatment of quadratic-form-based inference in repeated measures designs under diverse asymptotics. Electronic Journal of Statistics, 15(1). doi:10.1214/21-ejs1865
#' }
#'@references {
#'  Sattler, P., & Rosenbaum, M. (2025). Choice of the hypothesis matrix for using the Anova-type-statistic. Statistics & Probability Letters, 219(110356), 110356. doi:10.1016/j.spl.2025.110356
#'  }
#'
#'@examples
#'## load data set EEG (data frame)
#' data(EEG)
#'
#' # call test
#' hdrm_grouped(data = EEG$value,
#'              # test for no group effect
#'              hypothesis = "whole",
#'              group = EEG$group,
#'              ## if data is a vector, subject has to be specified
#'              subject = EEG$subject,
#'              subsampling = FALSE,
#'              B = "100*N"
#' )
#'
#'
#' # test using all subsampling version of the trace estimators
#' hdrm_grouped(data = EEG$value,
#'              # test for no time effect
#'              hypothesis = "sub",
#'              group = EEG$group,
#'              subject = EEG$subject,
#'              subsampling = TRUE,
#'              # B can also be given as a number
#'              B = 10000
#' )
#'
#'
#' # hypothesis as list
#' hypothesis_list <- list(TW = matrix(1/4, 4, 4),
#'                         TS = diag(40) - matrix(1/40, 40, 40)
#' ) # equivalent to hypothesis = "sub"
#'
#' # call test
#' hdrm_grouped(data = EEG$value,
#'              hypothesis = hypothesis_list,
#'              group = EEG$group,
#'              subject = EEG$subject,
#'              subsampling = FALSE,
#'              B = "100*N"
#' )
#' rm(hypothesis_list)
#'
#'
#' ## load data set birthrates (matrix)
#' data(birthrates)
#' # ?birthrates
#'
#' ## divide states into (former) east and west with Berlin as east
#' group <- factor(c(1,1,2,2,1,1,1,2,1,1,1,1,2,2,1,2), labels = c("west", "east"))
#'
#' ## call test for interaction effect
#' hdrm_grouped(data = birthrates,
#'              ## test for interaction effect between group and dimension
#'              hypothesis = "interaction",
#'              group = group,
#'              ## if data is a matrix, subject does not have to be specified
#'              B = "100*N"
#' )
#'
#'
#' ## test using all subsampling version of the trace estimators
#' hdrm_grouped(data = birthrates,
#'              ## test for group effect between group and dimension
#'              hypothesis = "whole",
#'              group = group,
#'              ## use the subsampling version of all trace estimators
#'              subsampling = TRUE,
#'              ## B can also be given as a number
#'              B = 10000
#' )
#'
#'
#' ## hypothesis as list, equivalent to hypothesis = "sub"
#' hypothesis_list <- list(TW = matrix(1/2, 2, 2),
#'                         TS = diag(34) - matrix(1/34, 34, 34)
#' )
#'
#' ## call test with hypothesis as list
#' hdrm_grouped(birthrates,
#'              ## equivalent to hypothesis = "sub"
#'              hypothesis = hypothesis_list,
#'              group = group,
#'              subsampling = TRUE,
#'              B = "100*N"
#' )
#'
#'@export
hdrm_grouped <- function(data, hypothesis = c("whole", "sub", "interaction", "identical", "flat"), AM = TRUE, group, subject = NULL, cov.equal = FALSE, subsampling = FALSE, B = "1000*N", seed = NULL, ...) {

  # Check if data is either a vector or matrix
  if(!is.vector(data) & !is.matrix(data)) stop("data must be a vector or matrix")

  # Check if 'group' is a vector or factor
  if(!is.vector(group) & !is.factor(group)) stop("group must be a vector containing group levels.")

  # Check if 'subject' is a vector or factor and not NULL
  if(!is.vector(subject) & !is.factor(subject) & !is.null(subject)) stop("subject must be a vector.")

  # If 'data' is a vector, check that subject has no missing values (NA)
  if(is.vector(data) & sum(is.na(subject)) > 0) stop("subject contains NA")

  # Check that the length of 'subject' and 'data' match if 'data' is a vector
  if(is.vector(data) & length(subject) != length(data)) stop("length of data vector does not coincide with subject")

  # Check that the length of 'group' and 'data' match if 'data' is a vector
  if(is.vector(data) & length(group) != length(data)) stop("length of data vector does not coincide with group")

  # If 'data' is a matrix, check that the number of rows matches the length of 'group'
  if(is.matrix(data) & max(!is.null(ncol(data)), ncol(data)) != length(group)) stop("number of rows of matrix does not coincide with group")

  # Warning if 'subject' is provided when 'data' is a matrix
  if(is.matrix(data) & !is.null(subject)) warning("for data in a matrix format, no additional subject parameter is necessary")

  # Check that 'group' does not contain missing values
  if(any(is.na(group))) stop("group must not contain missing values")


  if(is.vector(data)) {  # If 'data' is a vector

    # Create a data frame with 'data', 'subject', and 'group'
    dframe <- data.frame(value = data,
                         subject = subject,
                         whole = group
    )

    # Ensure that there are exactly three columns in the data frame
    if(length(names(dframe)) != 3) stop("could not find at least one column")

    # Ensure that 'value' column contains numeric data
    stopifnot(is.numeric(dframe$value))

    # Convert 'subject' and 'whole' to factors and remove unnecessary levels
    dframe$subject <- droplevels(as.factor(dframe$subject))
    dframe$whole <- droplevels(as.factor(dframe$whole))

    # Order the data frame by 'subject'
    dframe <- dframe[order(dframe$subject), ]

    ## Extract number of levels (N) and number of groups (a)
    N_with_NA <- nlevels(dframe$subject)
    a <- nlevels(dframe$whole)

    ## Remove missing values
    for (i in levels(dframe$subject)) {
      if(any(is.na(dframe$value[dframe$subject == i]))) {
        dframe$value[dframe$subject == i] <- NA
      }
    }

    # Filter out rows with missing values (NA) and remove unnecessary levels
    dframe <- dframe[stats::complete.cases(dframe), ]
    dframe <- droplevels(dframe)

    # Split the data frame by group ('whole')
    L <- split(dframe, dframe$whole)

    # Initialize an array to store dimensions for each group
    dimensions = rep(NA, a)

    # Check that all groups have the same number of unique subjects
    for (i in 1:a) {
      dimension <- unique(table(L[[1]]$subject))
      dimension <- dimension[dimension != 0]
      if(length(dimension) > 1) stop("all dimensions in the single groups have to be equal")
      dimensions[i] <- dimension
    }

    # Ensure that all groups have the same dimensions
    if(length(unique(dimensions)) > 1) stop("all groups must have the same dimensions")

    d = dimension[1]
    Nv = rep(0, a)

    # Get the number of subjects in each group
    for (i in 1:a) {
      Nv[i] <- length(unique(L[[i]]$subject))
    }

    # Calculate the total number of subjects
    N <- sum(Nv)

    # Initialize matrices for building the data
    X <- matrix(NA, d, 0)
    M <- matrix(NA, d, N)

    group <- numeric(0)

    # Fill the matrix with data for each group
    for(j in 1:a) {
      temp = droplevels(L[[j]])
      M <- matrix(0, d, nlevels(temp$subject))
      group <- c(group, rep(j, nlevels(temp$subject)))
      k = 1
      for (i in levels(temp$subject)) {
        M[, k] <- temp$value[temp$subject == i]
        k <- k + 1
      }
      X <- cbind(X, M)
    }

    # Store the processed data in the output
    out <- list(data = dframe)
  }

  if(is.matrix(data)) {  # If 'data' is a matrix

    # Remove missing values from the data
    N_with_NA <- ncol(data)
    group <- group[stats::complete.cases(t(data))]
    group <- droplevels(as.factor(group))
    X <- t(stats::na.omit(t(data)))

    # Get the number of subjects (N) and the number of features (d)
    N <- ncol(X)
    d <- nrow(X)



    # Store the processed data in the output
    out <- list(data = data)
  }

  # Warning if there were subjects with missing values
  if(N_with_NA > N) warning("Subjects with missing values dropped", call. = FALSE)

  # Validate the input for B (number of iterations for subsampling)
  if(!is.vector(B)) stop("B must be a vector of class character or numeric")
  if(!is.numeric(B) & !is.character(B)) stop("B must either be numeric or a character")
  if(length(B) > 1) {
    B <- B[1]
    warning("length(B) > 1. Only first element used")
  }

  # Evaluate B to get the number of repetitions for the subsampling process
  reps <- eval(parse(text = B[1]))
  reps <- ceiling(reps)

  # Check the grouping criteria
  check_criteria_grouped(X = X, group = group, hypothesis = hypothesis, reps = reps, subsampling = subsampling)

  subsampling <- subsampling[1]


  ### Output

  # If covariance is assumed to be equal across groups
  if(cov.equal) {
    # Call the internal function with equal covariance
    out <- c(out, hdrm_grouped_eq_cov_internal(
      data = X[, order(group)],
      group = sort(as.integer(group)),
      hypothesis = hypothesis,
      AM = AM,
      B = reps,
      subsampling =as.logical(subsampling),
      seed = seed
    ))
  }

  # If covariance is not assumed to be equal across groups
  if(!cov.equal) {
    # Call the internal function with unequal covariance
    out <- c(out, hdrm_grouped_internal(
      data = X[, order(group)],
      group = sort(as.integer(group)),
      hypothesis = hypothesis,
      AM = AM,
      subsampling = as.logical(subsampling),
      B = reps,
      seed = seed
    ))
  }

    # Add further output to the result
  out$groups$table <- table(group)
  out$removed.cases <- N_with_NA - N
  out$subsamples <- reps
  out$cov.equal <- cov.equal
  class(out) <- "hdrm_grouped"
  return(out)

}




# Print generics ----------------------------------------------------------

#' @method print hdrm_grouped
#' @export
print.hdrm_grouped <- function(x, digits = 4,...){
  ## p-Wert vorbereiten
  p <- round(x$p.value, digits)
  if(p <= 0){
    p = paste0("< ", 10^(-digits))
  }else{
    p = paste0("= ", p)
  }

  # print-output
  cat("\n")
  cat("          Multi Group Repeated Measure
      \nAnalysis of", x$dim$N, "individuals in", x$groups$a, "groups", "and", paste0(x$dim$d), "dimensions:",
      "\nW =", round(x$statistic, digits), " f =", round(x$f,digits), " p.value", p,
      "\nHypothesis type:", ifelse(is.list(x$hypothesis), "custom", x$hypothesis),
      "\nConvergence parameter \u03c4 =", round(x$tau,digits))
  cat("\n")
}
