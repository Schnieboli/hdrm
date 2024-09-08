#' Test for multiple group high dimensional repeated measures.
#'
#' @description See [hdrm_grouped.matrix] and [hdrm_grouped.data.frame] for data in  wide- and longtable versions, respectively.
#'
#' @param data either a matrix in widetable or a data.frame in longtable format
#' @param hypothesis either one of "whole", "sub" and "interaction" or a named list with elements TW and TS.
#' @param group if `data` is a matrix: a factor of length N that specifies the groups. If `data` is a `data.frame`: name or number of corresponding column.
#' @param value if `data` is a `data.frame`: name or number of value column.
#' @param subject if `data` is a `data.frame`: name or number of subject column.
#' @param factor if `data` is a `data.frame`: name or number of factor column.
#' @param bootstrap a logical indicating whether the bootstrap version of the trace estimators should be used.
#' @param B a string specifying a function of the number of subjects \eqn{N}. Determinates how many subsamples are used for the bootstrap version of the trace estimators.
#' @param ... further arguments. Currently ignored
#'
#' @details
#' `hypothesis` can either be given as a character or a list. Legal characters are "whole" (default), "sub" or "interaction". "whole" and "sub" stand for \eqn{H_0}: no differences in wholeplot/subplot-factor, while "interaction" means \eqn{H_0}: no interaction between subplot- and wholeplot-factor.
#' Alternatively `hypothesis` can be a named list with components `TW` for the wholeplot-part and `TS` for the subplot-part of \eqn{T = T_W \otimes T_S}.
#'
#' Note that even for `bootstrap = FALSE` the bootstrap version for the trace estimator C5 is used.
#' Therefore the `p.value` and the parameter `f` are seed dependent and may vary a lot for small B's.
#' For \eqn{p}-values near \eqn{alpha} it is advised to run the test again with `B = "5000*N"` at least.
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
hdrm_grouped <- function(data, hypothesis = c("whole","sub","interaction"), group, value, subject, factor, bootstrap = FALSE, B = "500*N",...){
  UseMethod("hdrm_grouped")
}

#' @method hdrm_grouped default
#' @export
hdrm_grouped.default <- function(data, hypothesis = c("whole","sub","interaction"), group, value, subject, factor, bootstrap = FALSE, B = "500*N",...){
  stop("Your data must be either a data frame or a matrix")
}


#' Test for multiple group high dimensional repeated measures.
#'
#' @description This is the widetable version. For longtable version see [hdrm_grouped.data.frame].
#'
#' @param data a matrix in widetable format with subjects in columns and factor levels in rows.
#' @param hypothesis either one of "whole", "sub" and "interaction" or a named list with elements TW and TS.
#' @param B a string specifying a function of the number of subjects \eqn{N}. Determinates how many subsamples are used for estimating C5.
#' @param ... further arguments. Currently ignored
#'
#' @details
#' If `data` is a matrix it is assumed, that individuals are represented by columns and factor levels by rows.
#'
#'
#' `hypothesis` can either be given as a character or a list. Legal characters are "whole" (default), "sub" or "interaction". "whole" and "sub" stand for \eqn{H_0}: no differences in wholeplot/subplot-factor, while "interaction" means \eqn{H_0}: no interaction between subplot- and wholeplot-factor.
#' Alternatively `hypothesis` can be a named list with components `TW` for the wholeplot-part and `TS` for the subplot-part of \eqn{T = T_W \otimes T_S}.
#'
#' Note, that the `p.value` may depend heavily on the choice of `B`, as both depend on the degrees of freedom `f`, which are estimated by C5.
#'
#'
#' @return a named list of class "hdrm_grouped" with the components
#' @returns \item{data}{the initial input data.}
#' @returns \item{f}{the degrees of freedom \eqn{f} of the critical value.}
#' @returns \item{statistic}{the test statistic \eqn{W} depending on \eqn{\alpha}.}
#' @returns \item{tau}{the convergence parameter \eqn{\tau}}
#' @returns \item{p.value}{the \eqn{p}-value of the test statistic}
#' @returns \item{dim}{a named list with components d and N of data}
#' @returns \item{groups}{a named list with components a and table giving the number and distribution of groups.}
#' @returns \item{removed.cases}{number of subjects that werde removed for missing values.}

#' @method hdrm_grouped matrix
#' @export
hdrm_grouped.matrix <- function(data, hypothesis = c("whole","sub","interaction"), group, bootstrap = FALSE, B = "500*N",...){

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
  if(!bootstrap[1]) out <- c(out, hdrm_grouped_internal(data = data[, order(group)], group = sort(as.integer(group)), hypothesis = hypothesis, B = reps))
  if(bootstrap[1]) out <- c(out, hdrm_grouped_internal_bootstrap(data = data, group = as.integer(group), hypothesis = hypothesis, B = reps))
  out$groups$table <- table(df$whole)/(nlevels(df$sub))
  out$subsamples = reps
  class(out) <- "hdrm_grouped"
  return(out)
}

#' Test for multiple group high dimensional repeated measures
#'
#' @description This is the longtable version. For widetable version see [hdrm_grouped.data.frame].
#'
#' @param data a data.frame in longtable format.
#' @param hypothesis either one of "whole", "sub" and "interaction" or a named list with elements TW and TS.
#' @param group name or number of corresponding column.
#' @param value name or number of value column.
#' @param subject name or number of subject column.
#' @param factor name or number of factor column.
#' @param B a string specifying a function of the number of subjects \eqn{N}. Determinates how many subsamples are used for estimating C5.
#' @param ... further arguments. Currently ignored
#'
#' @details
#' `hypothesis` can either be given as a character or a list. Legal characters are "whole" (default), "sub" or "interaction". "whole" and "sub" stand for \eqn{H_0}: no differences in wholeplot/subplot-factor, while "interaction" means \eqn{H_0}: no interaction between subplot- and wholeplot-factor.
#' Alternatively `hypothesis` can be a named list with components `TW` for the wholeplot-part and `TS` for the subplot-part of \eqn{T = T_W \otimes T_S}.
#'
#' Note, that the `p.value` may depend heavily on the choice of `B`, as both depend on the degrees of freedom `f`, which are estimated by C5.
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

#' @method hdrm_grouped data.frame
#' @export
hdrm_grouped.data.frame <- function(data, hypothesis = c("whole","sub","interaction"), group, value, subject, factor, bootstrap = FALSE, B = "500*N",...){


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
  ## für eine dimension gibt es eine andere anzahl an subjects
  if(length(unique(table(df$subject))) != 1) stop("All subjects must have the same number of dimensions")
  ## das ist für den fall, dass imdf eine Zeile "fehlt"!!!
  if(length(unique(table(df$sub))) != 1) stop("Unequal distribution of dimensions to subjects")
  if(!(as.logical(bootstrap) %in% c(TRUE, FALSE))) stop("bootstrap must be logical")


  # B extrahieren
  if(!(is.numeric(B) | is.character(B))) stop("B must either be numeric or a character")
  reps <- eval(parse(text = B))[1]
  reps <- ceiling(reps)
  stopifnot(is.numeric(reps), reps > 0)


  print(paste0("Performing test with group = ", group, ", value = ", value, ", subject = ", subject, "and factor = ", factor))

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
  if(!bootstrap[1]) out <- c(out, hdrm_grouped_internal(data = X[, order(group)], group = sort(as.integer(group)), hypothesis = hypothesis, B = reps))
  if(bootstrap[1]) out <- c(out, hdrm_grouped_internal_bootstrap(data = X[, order(group)], group = sort(as.integer(group)), hypothesis = hypothesis, B = reps))
  out$groups$table <- table((df$whole))/d
  out$subsamples = reps
  class(out) <- "hdrm_grouped"
  return(out)
}



# Print generics ----------------------------------------------------------

#' @method print hdrm_grouped
#' @export
print.hdrm_grouped <- function(x,...){
  cat("\n")
  cat("          Multi Group Repeated Measure
        \nAnalysis of", x$dim$N, "individuals in", x$groups$a, "groups", "and", paste0(x$dim$d), "dimensions:",
      "\nW =", x$statisitc, " f =", round(x$f,4), " p.value =", round(x$p.value, 4),
      "\nNull-Hypothesis:", ifelse(is.list(x$hypothesis), "custom", paste0("No effect in ",x$hypothesis,"plot-factor")),
      "\nConvergence parameter \u03c4 =", round(x$tau,4))
  cat("\n")
}

