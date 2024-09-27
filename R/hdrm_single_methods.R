#' Test for one group high dimensional repeated measures
#'
#' @description This function implements the methods outlined in
#'   \insertCite{Pauly2015;textual}{hdrm} for data in a widetable format. For
#'   data in a longtable format see [hdrm_single_longtable].
#'
#' @param data a data.frame with subjects represented by columns and factor levels
#'   represented by rows.
#' @param hypothesis either "equal" or a quadratic numeric matrix.
#' @param ... further arguments are currently ignored.
#' @details
#' If `data` contains missing values, affected subjects will be dropped with a
#' warning.
#'
#' The test outlined in \insertCite{Pauly2015;textual}{hdrm} is performed for
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
#' @return \item{data}{the initial input data.}
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
#' @example examples/examples_hdrm_single_widetable.R
#'
#' @references \insertAllCited
#'
#' @export
hdrm_single_widetable <- function(data, hypothesis = "flat",...){

  # data matrix?
  if(!is.data.frame(data) | !all(sapply(data, is.numeric))) stop("data must be a data.frame with numeric entries")
  M <- as.matrix(data)
  #Anforderungen ueberpruefen
  check_criteria_single(X = M, hypothesis = hypothesis)

  # Output erstellen
  out <- list(data = data)
  # Funktionsaufruf
  out <- c(out, hdrm1_internal(X = M, hypothesis = hypothesis,...))
  class(out) <- "hdrm_single"
  return(out)
}



#' Test for one group high dimensional repeated measures
#'
#' @description This function implements the methods outlined in
#'   \insertCite{Pauly2015;textual}{hdrm} for data in a longtable format. For
#'   data in a widetable format see [hdrm_single_widetable].
#'
#' @param data a data.frame in longtable format.
#' @param hypothesis either "equal" or a quadratic numeric matrix.
#' @param value name or index of value column.
#' @param subject name or index of subject column.
#' @param dimension name or index of dimension column.
#' @param ... further arguments are currently ignored.
#' @details
#' The function can deal with missing values only in the `value` column.
#' Affected subjects are dropped with a warning. Missing values in any other
#' column of data will result in an error.
#'
#' The test outlined in \insertCite{Pauly2015;textual}{hdrm} is performed for
#' the hypothesis given by `hypothesis`. The `hypothesis` can either be given as
#' a `character` or a `list`. The only legal character is "flat" (default).
#' "flat" stands for \eqn{H_0}: time profile is flat. Other
#' characters will result in an error.
#'
#' Alternatively `hypothesis` can be a quadratic hypothesis matrix \eqn{T} with
#' the number of rows equal to the number of levels in the  of `data`. \eqn{T} must be
#' a projection matrix, meaning symmetrical and \eqn{T^2 = T}. A matrix that does not
#' match those criteria will result in an error.
#'
#'
#' @returns Returns a list with class "hdrm_single" with the components
#' @return \item{data}{the initial input data.}
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
#' @example examples/examples_hdrm_single_longtable.R
#'
#' @references \insertAllCited
#'
#' @export
hdrm_single_longtable <- function(data, hypothesis = "flat", value, subject, dimension,...){

  if(!is.data.frame(data)) stop("data must be a data.frame")
  if(length(value) != 1) stop("value must be of length 1")
  if(length(subject) != 1) stop("subject must be of length 1")
  if(length(dimension) != 1) stop("dimension must be of length 1")

  ## relevante spalten extrahieren
  # auf diese weise koennen value, subject, factor sowohl die namen als auch die nummer der spalte sein
  df <- data.frame(value = data[[value]],
                   subject = data[[subject]],
                   factor = data[[dimension]])
  # falls eine der spalten nicht gefunden werden kann ist diese NULL und df kuerzer als 3:
  if(ncol(df) != 3) stop("could not find at least one of the columns")

  # nicht benutzte Faktorlevel loeschen, sicherstellen, dass subject und factor
  # faktoren sind und value numerisch ist
  df$subject <- droplevels(as.factor(df$subject))
  df$factor <- droplevels(as.factor(df$factor))
  if(!is.numeric(df$value)) stop("value must be numeric")

  # Überpruefen, ob alle dimensionen gleich sind
  if(length(unique(table(df$subject))) != 1) stop("All subjects must have the same number of dimensions")
  if(length(unique(table(df$factor))) != 1) stop("Unequal distribution of dimensions to subjects")
  d <- nlevels(df$factor)
  N_with_NA <- nlevels(df$subject)


  ## fehlende Werte entfernen
  for (i in levels(df$subject)) {
    if(any(is.na(df$value[df$subject == i]))){
      df$value[df$subject == i] <- NA
    }
  }
  df <- df[complete.cases(df),]
  df <- droplevels(df)
  N <- nlevels(df$subject)

  ### Matrix bauen
  k <- 1
  M <- matrix(NA, d, N)
  for (i in levels(df$subject)) { # baue Matrix spaltenweise/fuege in jedem Schritt ein Individuum in die i-te Spalte ein
    M[ ,k] <- df$value[df$subject == i][order(df$factor[df$subject == i])]
    k <- k+1
  }

  #Anforderungen ueberpruefen
  check_criteria_single(X = data, hypothesis = hypothesis)


  # Warnung fuer fehlende Werte -> erst hier, da vorher Abbruch aus anderen Gruenden moeglich
  if(N_with_NA > N) warning("Subjects with missing values dropped", call. = FALSE)

  # output erstellen
  out <- list(data = df)
  ## Funktionsaufruf
  out <- c(out, hdrm1_internal(X = M, hypothesis = hypothesis, alpha = 0.05))
  out$removed.cases <- N_with_NA - N
  class(out) <- "hdrm_single"
  return(out)
}



# Print Funktionen --------------------------------------------------------

#' @method print hdrm_single
#' @export
print.hdrm_single <- function(x, digits = 4,...){
  # print-output
  cat("\n")
  cat("           One Group Repeated Measure
       \nAnalysis of", x$dim[2], "subjects in", paste0(x$dim[1]), "dimensions:",
      "\nW =", round(x$statisitc, digits), " f =", round(x$f, digits), " p.value =", round(x$p.value, digits),
      "\nHypothesis type:", x$hypothesis,
      "\nConvergence parameter \u03c4 =", round(x$tau, digits))
  cat("\n")
}
