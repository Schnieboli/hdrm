#' Test for multiple group high dimensional repeated measures
#'
#' @param data either a matrix, a list or a data.frame in longtable format containing the repeated measures data.
#' @param group if data is a matrix: a factor of length N that specifies groups.
#' @param formula a formula object that specifies the corresponding columns if data is a data.frame. Must be of the form `value ~ subject + group + dimension`
#' @param hypothesis either one of `whole`, `sub`  or `interaction` or a named list specifying \eqn{T_W} and \eqn{T_S}.
#' @param alpha the \eqn{\alpha}-level used for calculating the test statistic \eqn{W}.
#' @param B a string specifying a function of the number of subjects \eqn{N}. Determinates how many subsamples are used for estimating C5.
#' @param ... further arguments. Currently ignored
#'
#' @details
#' If `data` is a matrix or a list containing multiple matrices of the same amount of rows then it is assumed, that that individuals are represented by columns and fator levels by rows.
#'
#' If `data` is a data.frame then formula must be a formula object of the form `value ~ subject + group + dimension` specifying the corresponding columns of data.frame. Each subject must be represented by a unique ID.
#' Note, that specifying an interaction effect in  `formula` will give an error. For testing interaction effects, refer to `hypotheis`.
#'
#' `hypothesis` can either be given as "whole" (default), "sub" or "interaction". "whole" and "sub" stand for \eqn{H_0}: no differences in wholeplot/subplot-factor, while "interaction" means \eqn{H_0}: no interaction between subplot- and wholeplot-factor.
#' Alternatively `hypothesis` can be a named list with components `TW` for the wholeplot-part and `TS` for the subplot-part of \eqn{T = T_W \otimes T_S}.
#'
#' Note, that `p.value` and `critical.value` may depend heavily on the choice of `B`, as both depend on the degrees of freedom `f`, which are estimated by C5.
#'
#'
#' @return a named list of class "hdrm_test" with the components
#' @returns \item{data}{the initial input data converted to a matrix.}
#' @returns \item{f}{the degrees of freedom \eqn{f} of the critical value.}
#' @returns \item{statistic}{the test statistic \eqn{W} depending on \eqn{\alpha}.}
#' @returns \item{tau}{the convergence parameter \eqn{\tau}}
#' @returns \item{p.value}{the \eqn{p}-value of the test statistic}
#' @returns \item{dim}{a named list with components d and N of the input data}
#' @returns \item{groups}{a named list with components a and table giving the number and distribution of groups}
#' @returns \item{removed.cases}{number of incomplete subjects removed.}
#'
#' @references Paavo Sattler, Markus Pauly. "Inference for high-dimensional split-plot-designs: A unified approach for small to large numbers of factor levels." Electronic Journal of Statistics, 12(2) 2743-2805 2018. https://doi.org/10.1214/18-EJS1465
#'
#' @export
hdrm_test <- function(data, hypothesis = c("whole","sub","interaction"), group, value, subject, factor, bootstrap = FALSE, B = "500*N",...){
  UseMethod("hdrm_test")
}

#' @method hdrm_test default
#' @export
hdrm_test.default <- function(data, hypothesis = c("whole","sub","interaction"), group, value, subject, factor, bootstrap = FALSE, B = "500*N",...){
  stop("Your data must be either a list, a data frame or a matrix")
}


#' @method hdrm_test matrix
#' @export
hdrm_test.matrix <- function(data, hypothesis = c("whole","sub","interaction"), group, bootstrap = FALSE, B = "500*N",...){

  # df erstellen für output
  N <- ncol(data)
  d <- nrow(data)
  df <- data.frame(value = as.vector(data))
  df$subject <- as.factor(rep(1:N, each = d))
  df$whole <- as.factor(rep(group, each = d))
  df$sub <- as.factor(rep(1:d, N))

  ## Voraussetzungen überprüfen
  if(any(table(group) < 6)) stop("all group sizes must be >= 6")
  if(length(table(group)) < 2) stop("at least two groups needed")
  if(length(group) != N) stop("length of group must be equal to ncol(data")
  if(d < 2) stop("at least two dimensions are needed")
  if(!is.character(hypothesis) | is.list(hypothesis)) stop("hypothesis must be a character or a list")

  if(!(as.logical(bootstrap) %in% c(TRUE, FALSE))) stop("bootstrap must be logical")

  # B extrahieren
  if(!(is.numeric(B) | is.character(B))) stop("B must either be numeric or a character")
  reps <- eval(parse(text = B))[1]
  reps <- ceiling(reps)
  stopifnot(is.numeric(reps), reps > 0)

  out <- list(data = df)
  if(!bootstrap[1]) out <- c(out, hdrm_test_internal(data = data[, order(group)], group = sort(as.integer(group)), hypothesis = hypothesis, B = reps))
  if(bootstrap[1]) out <- c(out, hdrm_test_internal_bootstrap(data = data, group = as.integer(group), hypothesis = hypothesis, B = reps))
  out$groups$table <- table(df$whole)/(nlevels(df$sub))
  out$subsamples = reps
  class(out) <- "hdrm"
  return(out)
}


#' @method hdrm_test data.frame
#' @export
hdrm_test.data.frame <- function(data, hypothesis = c("whole","sub","interaction"), group, value, subject, factor, bootstrap = FALSE, B = "500*N",...){


  df <- data.frame(value = data[[value]],
                   subject = data[[subject]],
                   whole = data[[group]],
                   sub = data[[factor]]
                   )


  if(length(names(df)) != 4) stop("could not find at least one column")
  stopifnot(is.numeric(df$value))

  df$subject <- droplevels(as.factor(df$subject))
  df$whole <- droplevels(as.factor(df$whole))
  df$sub <- droplevels(as.factor(df$sub))
  df <- df[order(df$subject, df$sub), ]

  ## Dimensionen
  N <- nlevels(df$subject)
  a <- nlevels(df$whole)
  if(a < 2) stop("at least two groups needed")
  d <- nlevels(df$sub)
  if(d < 2) stop("at least 2 dimensions are needed")
  n <- as.integer(table(df$group))/d
  if(any(n < 6)) stop("all group sizes must be >= 6")
  if(!is.character(hypothesis) | is.list(hypothesis)) stop("hypothesis must be a character or a list")
  # alle Individuen haben gleich viele dimensionen
  if(length(unique(table(df$subject))) != 1) stop("All subjects must have the same number of dimensions")
  # ??
  if(length(unique(table(df$sub))) != 1) stop("Unequal distribution of dimensions to subjects")
  if(!(as.logical(bootstrap) %in% c(TRUE, FALSE))) stop("bootstrap must be logical")


  # B extrahieren
  if(!(is.numeric(B) | is.character(B))) stop("B must either be numeric or a character")
  reps <- eval(parse(text = B))[1]
  reps <- ceiling(reps)
  stopifnot(is.numeric(reps), reps > 0)

  # Matrix erstellen
  X <- matrix(NA, d, 0)
  M <- matrix(NA, d, N)
  L <- split(df, df$whole)
  group <- numeric(0)

  for(j in 1:a){
    temp = droplevels(L[[j]])
    M <- matrix(0, d, nlevels(temp$subject))
    group <- c(group, rep(j, nlevels(temp$subject)))
    k = 1
    for (i in levels(temp$subject)) {
      M[,k] <- temp$value[temp$subject == i]
      k <- k + 1
    }
    X <- cbind(X,M)
  }

  ### Output
  out <- list(data = df)
  if(!bootstrap[1]) out <- c(out, hdrm_test_internal(data = X[, order(group)], group = sort(as.integer(group)), hypothesis = hypothesis, B = reps))
  if(bootstrap[1]) out <- c(out, hdrm_test_internal_bootstrap(data = X[, order(group)], group = sort(as.integer(group)), hypothesis = hypothesis, B = reps))
  out$groups$table <- table((df$whole))/d
  out$subsamples = reps
  class(out) <- "hdrm"
  return(out)
}



# Print generics ----------------------------------------------------------

#' @method print hdrm
#' @export
print.hdrm <- function(x,...){
  cat("\n")
  cat("          Multi Group Repeated Measure
        \nAnalysis of", x$dim$N, "individuals in", x$groups$a, "groups", "and", paste0(x$dim$d), "dimensions:",
      "\nW =", x$statisitc, " f =", round(x$f,4), " p.value =", round(x$p.value, 4),
      "\nNull-Hypothesis:", ifelse(is.list(x$hypothesis), "custom", paste0("No effect in ",x$hypothesis,"plot-factor")),
      "\nConvergence parameter \u03c4 =", round(x$tau,4))
  cat("\n")
}
#' @method summary hdrm
#' @export
summary.hdrm <- function(object,...){
  print(object,...)
}
