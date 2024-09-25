#' Test for multiple group high dimensional repeated measures
#'
#' @description This function implements the methods outlined
#'   \insertCite{Sattler2018;textual}{hdrm} for data in a widetable format. For data in
#'   a longtable format see [hdrm_grouped_longtable].
#'
#' @param data a data frame with subjects represented by columns and factor
#' levels represented by rows.
#' @param hypothesis either one of "whole", "sub" and "interaction" or a named
#'   list with quadratic matrices `TW` and `TS`.
#' @param group a `factor` specifying the groups.
#' @param B a `string` specifying a function of the number of subjects \eqn{N}.
#'   Determines the number of subsamples used by the  subsampling trace estimators.
#' @param subsampling `logical`. Specifying whether the subsampling versions for all
#'   trace estimators should be used (see details).
#' @param ... further arguments. Currently ignored
#'
#' @details
#' If `data` contains missing values, affected subjects will be dropped with a
#' warning.
#'
#' The test outlined in \insertCite{Sattler2018;textual}{hdrm} is performed for
#' the hypothesis given by `hypothesis`. The `hypothesis` can either be given as
#' a `character` or a `list`. Legal characters are "whole" (default), "sub" or
#' "interaction". "whole" and "sub" stand for \eqn{H_0}: no differences in
#' wholeplot/subplot-factor, while "interaction" means \eqn{H_0}: no interaction
#' between subplot- and wholeplot-factor. Other characters will result in an
#' error.
#'
#' Alternatively, `hypothesis` can be a named list with quadratic
#' matrices `TW` for the wholeplot-part and `TS` for the subplot-part of \eqn{T
#' = T_W \otimes T_S}. If `hypothesis` is a list, `TW` and `TS` must be
#' projection matrices, meaning symmetrical and \eqn{T_W^2 = T_W}, \eqn{T_S^2 = T_S}.
#' Also, the number of rows and columns of `TW` must be equal to the number of
#' groups and the number of rows and columns of `TS` must be equal to the number
#' of factor levels. Lists that do not match those criteria will result in an
#' error.
#'
#' Note that for `subsampling = FALSE` your results are still seed dependent, because
#' the computational heaviest trace estimator is still calculated using
#' subsamples Also, the non subsampling versions might not be faster for small
#' data, depending on the choice of `B`.
#'
#' That also means, that even for `subsampling = FALSE` the `p.value` depends
#' heavily on the choice of `B`, as it depends on the estimated degrees of
#' freedom `f`, which are always estimated by the subsampling version. For
#' `subsampling = TRUE`
#' the test statistic \eqn{W} is also seed dependent.
#'
#' The number of subsamples `B` can also be a numeric value. However, this is
#' not advised, as `B` should be a function of \eqn{N} that goes to \eqn{\infty}
#' for \eqn{N\to\infty}.
#'
#'
#' @return a named list of class "hdrm_grouped" with the components
#' @returns \item{data}{the input data used.}
#' @returns \item{f}{the estimated degrees of freedom \eqn{f}.}
#' @returns \item{tau}{the convergence parameter \eqn{\tau}}
#' @returns \item{H}{a named list with components `TW` and `TS` that give the
#'   components of the hypothesis matrix}
#' @returns \item{hypothesis}{a character. Will be "custom" if `hypothesis` is a
#'   list, otherwise `hypothesis[1]`.}
#' @returns \item{p.value}{the \eqn{p}-value of the test statistic.}
#' @returns \item{dim}{a named list with with number of factor levels \eqn{d}
#'   and number of subjects \eqn{N} of `data`}
#' @returns \item{groups}{a named list with components number of groups `a`
#'   and distribution of groups `table`.}
#' @returns \item{removed.cases}{number of incomplete subjects removed.}
#' @returns \item{subsamples}{number of subsamples used for subsampling estimators.}
#'
#' \insertNoCite{Sattler2018}{hdrm}
#' @references \insertAllCited
#'
#' @example examples/examples_hdrm_grouped_widetable.R
#'
#' @export
hdrm_grouped_widetable <- function(data, hypothesis = c("whole","sub","interaction", "identical", "flat"), group, subsampling = FALSE, B = "500*N",...){

  # data muss data.frame sein
  if(!is.data.frame(data)) stop("data must be a a data.frame with numeric entries")
  # group muss vektor sein
  if(!is.vector(group) & !is.factor(group)) stop("group must be a character vector or a factor.")

  # fehlende werte entfernen
  M <- as.matrix(data)
  N_with_NA <- ncol(M)
  group <- group[stats::complete.cases(t(M))]
  group <- droplevels(as.factor(group))
  M <- t(stats::na.omit(t(M)))

  # Dimensionen extrahieren
  N <- ncol(M)
  d <- nrow(M)

  # Warnung fuer NAs
  if(N_with_NA > N) warning("Subjects with missing values dropped", call. = FALSE)

  # B extrahieren
  if(!is.vector(B)) stop("B must be a vector of class character or numeric")
  if(!is.numeric(B) & !is.character(B)) stop("B must either be numeric or a character")
  if(length(B) > 1){
    B <- B[1]
    warning("length(B) > 1. Only first elemt used")
  }
  # soll ganze zahl groesser 0 sein
  reps <- eval(parse(text = B[1]))
  reps <- ceiling(reps)

  # Voraussetzungen ueberpruefen
  check_criteria_grouped(X = M, group = group, hypothesis = hypothesis, reps = reps, subsampling = subsampling)
  subsampling <- subsampling[1]

  # output erstellen
  out <- list(data = data)
  #funktionsaufruf
  out <- c(out, hdrm_grouped_internal(
    data = M[, order(group)],
    group = sort(as.integer(group)),
    hypothesis = hypothesis,
    B = reps,
    subsampling = subsampling
  ))
  # weiteren output hinzufügen
  out$groups$table <- table(group)
  out$removed.cases <- N_with_NA - N
  out$subsamples = reps

  class(out) <- "hdrm_grouped"
  return(out)
}

#' Test for multiple group high dimensional repeated measures
#'
#' @description This function implements the methods outlined in
#'   \insertCite{Sattler2018;textual}{hdrm} for data in a longtable format. For data in
#'   a widetable format see [hdrm_grouped_widetable]
#' @param data a data.frame in longtable format.
#' @param hypothesis either one of "whole", "sub" and "interaction" or a named
#'   list with quadratic matrices `TW` and `TS`.
#' @param group name or index of group column.
#' @param value name or index of value column.
#' @param subject name or index of subject column.
#' @param dimension name or index of dimension column.
#' @param B a `string` specifying a function of the number of subjects \eqn{N}.
#'   Determines the number of subsamples used by the  subsampling trace estimators.
#' @param subsampling `logical`. Specifying whether the subsampling versions for all
#'   trace estimators should be used (see details).
#' @param ... further arguments. Currently ignored
#'
#' @details
#' The function can deal with missing values only in `value`. The test is then
#' performed without the affected subjects. Missing values in any of the other
#' columns of `data` will result in an error.
#'
#' The `data` is converted to a matrix and the test outlined in
#' \insertCite{Sattler2018;textual}{hdrm} is performed for the hypothesis given by
#' `hypothesis`. The `hypothesis` can either be given as a `character` or a `list`.
#' Legal characters are "whole" (default), "sub" or "interaction". "whole" and
#' "sub" stand for \eqn{H_0}: no differences in wholeplot/subplot-factor, while
#' "interaction" means \eqn{H_0}: no interaction between subplot- and
#' wholeplot-factor. Other characters will result in an error.
#'
#' Alternatively, `hypothesis` can be a named list with quadratic
#' matrices `TW` for the wholeplot-part and `TS` for the subplot-part of \eqn{T
#' = T_W \otimes T_S}. If `hypothesis` is a list, `TW` and `TS` must be
#' projection matrices, meaning symmetrical and \eqn{T_W^2 = T_W}, \eqn{T_S^2 = T_S}.
#' Also, the number of rows and columns of `TW` must be equal to the number of
#' groups and the number of rows and columns of `TS` must be equal to the number
#' of factor levels. Lists that do not match those criteria will cause an
#' error.
#'
#' Note that for `subsampling = FALSE` your results are still seed dependent, because
#' the computational heaviest trace estimator is still calculated using
#' subsamples. Also, depending on the choice of `B`, the non subsampling versions might not be faster for small
#' data.
#'
#' That also means, that even for `subsampling = FALSE` the `p.value` depends
#' heavily on the choice of `B`, as it depends on the degrees of freedom `f`,
#' which are always estimated by the subsampling version. For `subsampling = TRUE`
#' the test statistic \eqn{W} is also seed dependent.
#'
#' The number of subsamples `B` can also be a numeric value. However, this is
#' not advised, as `B` should be a function of \eqn{N} that goes to \eqn{\infty}
#' for \eqn{N\to\infty}.
#'
#'
#' @return a named list of class "hdrm_grouped" with the components
#' @returns \item{data}{the input data used.}
#' @returns \item{f}{the degrees of freedom \eqn{f}.}
#' @returns \item{tau}{the convergence parameter \eqn{\tau}}
#' @returns \item{H}{a named list with components `TW` and `TS` that give the
#'   components of the hypothesis matrix}
#' @returns \item{hypothesis}{a character. Will be "custom" if `hypothesis` is a
#'   list, otherwise `hypothesis[1]`.}
#' @returns \item{p.value}{the \eqn{p}-value of the test statistic.}
#' @returns \item{dim}{a named list with with number of factor levels \eqn{d}
#'   and number of subjects \eqn{N} of `data`}
#' @returns \item{groups}{a named list with components number of groups `a`
#'   and distribution of groups `table`.}
#' @returns \item{removed.cases}{number of incomplete subjects removed.}
#' @returns \item{subsamples}{number of subsamples used for subsampling estimators.}
#'
#' \insertNoCite{Sattler2018}{hdrm}
#' @references \insertAllCited
#'
#' @example examples/examples_hdrm_grouped_longtable.R
#'
#' @export
hdrm_grouped_longtable <- function(data, hypothesis = c("whole","sub","interaction"), group, value, subject, dimension, subsampling = FALSE, B = "500*N",...){

  # data muss data.frame sein
  if(!is.data.frame(data)) stop("data must be a data.frame")
  # muss alles laenge 1 sein
  if(length(group) != 1) stop("group must be of length 1")
  if(length(value) != 1) stop("value must be of length 1")
  if(length(subject) != 1) stop("subject must be of length 1")
  if(length(dimension) != 1) stop("dimension must be of length 1")

  # relevante spalten extrahieren
  df <- data.frame(value = data[[value]],
                   subject = data[[subject]],
                   whole = data[[group]],
                   sub = data[[dimension]]
                   )

  # wenn eine spalte nicht gefunden werden konnte, sit sie NULL und df zu kurz
  if(length(names(df)) != 4) stop("could not find at least one column")
  stopifnot(is.numeric(df$value))

  # sicher stellen, dass alles ein faktor ist und ueberfuessige levels entfernen
  df$subject <- droplevels(as.factor(df$subject))
  df$whole <- droplevels(as.factor(df$whole))
  df$sub <- droplevels(as.factor(df$sub))
  df <- df[order(df$subject, df$sub), ]

  # keine NAs in group, value, factor erlaubt
  if(any(is.na(df$whole))) stop("data$group must not contain missing values")
  if(any(is.na(df$sub))) stop("data$factor must not contain missing values")
  if(any(is.na(df$subject))) stop("data$subject must not contain missing values")

  ## Dimensionen extrahieren
  N <- nlevels(df$subject)
  a <- nlevels(df$whole)
  d <- nlevels(df$sub)

  # erste Voraussetzungen ueberpruefen
  ## fuer eine dimension gibt es eine andere anzahl an subjects
  if(length(unique(table(df$subject))) != 1) stop("All subjects must have the same number of dimensions")
  ## das ist fuer den fall, dass mind eine Zeile "fehlt"!!!
  if(length(unique(table(df$sub))) != 1) stop("Unequal distribution of dimensions to subjects")


  # B extrahieren
  if(!is.vector(B)) stop("B must be a vector of class character or numeric")
  if(!is.numeric(B) & !is.character(B)) stop("B must either be numeric or a character")
  if(length(B) > 1){
    B <- B[1]
    warning("length(B) > 1. Only first elemt used")
  }
  # reps soll ganze zahl groesser 0 sein
  reps <- eval(parse(text = B[1]))
  reps <- ceiling(reps)


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
      M[, k] <- temp$value[temp$subject == i]
      k <- k + 1
    }
    X <- cbind(X,M)
  }
  # X ist nach Gruppen sortiert, da gruppenweise aufgebaut


  # fehlende werte entfernen
  N_with_NA <- ncol(X)
  group <- group[stats::complete.cases(t(X))]
  group <- droplevels(as.factor(group))
  X <- t(stats::na.omit(t(X)))
  # neues N
  N <- ncol(X)

  # Voraussetzungen erneut ueberpruefen
  check_criteria_grouped(X = X, group = group, hypothesis = hypothesis, reps = reps, subsampling = subsampling)
  subsampling <- subsampling[1]
  # Warnung fuer fehlende Werte
  if(N_with_NA > N) warning("Subjects with missing values dropped", call. = FALSE)

  ### Output
  out <- list(data = df)
  out <- c(out, hdrm_grouped_internal(
    data = X[, order(group)],
    group = sort(as.integer(group)),
    hypothesis = hypothesis,
    B = reps,
    subsampling = as.logical(subsampling)
  ))
  # noch output hinzufuegen
  out$groups$table <- table(group)
  out$removed.cases <- N_with_NA - N
  out$subsamples <- reps
  class(out) <- "hdrm_grouped"
  return(out)
}



# Print generics ----------------------------------------------------------

#' @method print hdrm_grouped
#' @export
print.hdrm_grouped <- function(x, digits = 4,...){
  # print-output
  cat("\n")
  cat("          Multi Group Repeated Measure
      \nAnalysis of", x$dim$N, "individuals in", x$groups$a, "groups", "and", paste0(x$dim$d), "dimensions:",
      "\nW =", round(x$statisitc, digits), " f =", round(x$f,digits), " p.value =", round(x$p.value, digits),
      "\nHypothesis type:", ifelse(is.list(x$hypothesis), "custom", x$hypothesis),
      "\nConvergence parameter \u03c4 =", round(x$tau,digits))
  cat("\n")
}

