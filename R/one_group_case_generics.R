#' Test for one group high dimensional repeated measures
#'
#' @description See [hdrm1.matrix] and [hdrm1.data.frame] for data in  wide- and longtable versions, respectively.
#'
#' @param data either a matrix with subjects in columns and factor levels in rows or a data.frame of longtable format
#' @param hypothesis a character or a quadratic matrix
#' @param value if `data` is a `data.frame`: name or number of value column.
#' @param subject if `data` is a `data.frame`: name or number of subject column.
#' @param factor if `data` is a `data.frame`: name or number of factor column.
#' @param ... further arguments are currently ignored.
#' @details
#' If `data` is a matrix or a list of multiple matrices with the same amount of rows then it is assumed that individuals are represented by columns and factor levels by rows.
#'
#' If `data` is a data.frame, `formula` must be a formula object of the form `value ~ subject + factor` that specifies the corresponding columns of the `data.frame`.
#' Each subject and factor level must be represented by a unique ID.
#'
#' `hypothesis`` can either be given as one of "flat" or "equal" for testing whether the time profile is flat: \eqn{H_0: \mu_1 = \ldots = \mu_d} or equal: \eqn{H_0: \mu = 0}. Alternatively 'hypothesis' can also be a custom contrast matrix \eqn{T \in \mathbb{R}^{d\times d}}
#'
#' Subjects with missing values will be disregarded.
#'
#' @returns Returns a list with class "hdrm1" with the components
#' @return \item{data}{the initial input data. If a data.frame was given it will have been transformed to a matrix.}
#' @return \item{f}{the degrees of freedom \eqn{f} of the distribution of the test statistic}
#' @return \item{statistic}{the test statistic \eqn{W}}
#' @return \item{tau}{the convergence parameter \eqn{\tau}}
#' @return \item{H0}{the Nullhypothesis tested. Will be 'custom' if `hypothesis` was given as a matrix.}
#' @return \item{H}{the hypothesis-matrix used.}
#' @return \item{p.value}{the p-value of the test statistic}
#' @return \item{critical.value}{the critical value depending on \eqn{\alpha}}
#' @return \item{dim}{a vector of length 2, giving the dimensions \eqn{d \times N} of the data}
#' @return \item{removed.cases}{number of subjects removed for having missing values.}
#'
#' \insertNoCite{Pauly2015}{hdrm}
#' @references \insertAllCited
#'
#' @examples
#' # see help pages for the respective generics
#'
#' @export
hdrm1 <- function(data, hypothesis = c("flat", "equal"),  value, subject, factor,...){
  UseMethod("hdrm1")
}

#' @method hdrm1 default
#' @export
hdrm1.default <- function(data, hypothesis = c("flat", "equal"),  value, subject, factor,...){
  stop("Your data needs to be either a matrix or a data.frame")
}

#' Test for one group high dimensional repeated measures
#'
#' @description
#' S3 generic methods for performing a one sample test for one group repeated measure data
#' @param data either a matrix with subjects in columns and factor levels in rows or a data.frame of longtable format
#' @param hypothesis a character or a quadratic matrix
#' @param ... further arguments are currently ignored.
#' @details
#' If `data` is a matrix or a list of multiple matrices with the same amount of rows then it is assumed that individuals are represented by columns and factor levels by rows.
#'
#' If `data` is a data.frame, `formula` must be a formula object of the form `value ~ subject + factor` that specifies the corresponding columns of the `data.frame`.
#' Each subject and factor level must be represented by a unique ID.
#'
#' `hypothesis`` can either be given as one of "flat" or "equal" for testing whether the time profile is flat: \eqn{H_0: \mu_1 = \ldots = \mu_d} or equal: \eqn{H_0: \mu = 0}. Alternatively 'hypothesis' can also be a custom contrast matrix \eqn{T \in \mathbb{R}^{d\times d}}
#'
#' Subjects with missing values will be disregarded.
#'
#' @returns Returns a list with class "hdrm1" with the components
#' @return \item{data}{the initial input data. If a data.frame was given it will have been transformed to a matrix.}
#' @return \item{f}{the degrees of freedom \eqn{f} of the distribution of the test statistic}
#' @return \item{statistic}{the test statistic \eqn{W}}
#' @return \item{tau}{the convergence parameter \eqn{\tau}}
#' @return \item{H0}{the Nullhypothesis tested. Will be 'custom' if `hypothesis` was given as a matrix.}
#' @return \item{H}{the hypothesis-matrix used.}
#' @return \item{p.value}{the p-value of the test statistic}
#' @return \item{critical.value}{the critical value depending on \eqn{\alpha}}
#' @return \item{dim}{a vector of length 2, giving the dimensions \eqn{d \times N} of the data}
#' @return \item{removed.cases}{number of subjects removed for having missing values.}
#'
#' \insertNoCite{Pauly2015}{hdrm}
#' @references \insertAllCited
#'
#' @method hdrm1 matrix
#' @export
hdrm1.matrix <- function(data, hypothesis = c("flat","equal"),...){ # na.action muss raus -> klappt nicht so, wie ich mir das vorgestellt habe und ist ja auch egal, weil ja eh nur na.omit funktionieren soll

  # data.frame vorbereiten (zur Ausgabe)
  d <- nrow(data)
  N <- ncol(data)
  df <- data.frame(value = as.vector(data))
  df$subject <- rep(1:N, each = d)
  df$factor <- rep(1:d, N)
  df <- df[order(df$subject),]

  out <- data.frame(data = df)
  out <- c(out, hdrm1_internal(X = data, hypothesis = hypothesis))
  class(out) <- "hdrm1"
  return(out)
}



#' Test for one group high dimensional repeated measures
#'
#' @description
#' S3 generic methods for performing a one sample test for one group repeated measure data
#' @param data either a matrix with subjects in columns and factor levels in rows or a data.frame of longtable format
#' @param hypothesis a character or a quadratic matrix
#' @param value if `data` is a `data.frame`: name or number of value column.
#' @param subject if `data` is a `data.frame`: name or number of subject column.
#' @param factor if `data` is a `data.frame`: name or number of factor column.
#' @param ... further arguments are currently ignored.
#' @details
#' If `data` is a matrix or a list of multiple matrices with the same amount of rows then it is assumed that individuals are represented by columns and factor levels by rows.
#'
#' If `data` is a data.frame, `formula` must be a formula object of the form `value ~ subject + factor` that specifies the corresponding columns of the `data.frame`.
#' Each subject and factor level must be represented by a unique ID.
#'
#' `hypothesis`` can either be given as one of "flat" or "equal" for testing whether the time profile is flat: \eqn{H_0: \mu_1 = \ldots = \mu_d} or equal: \eqn{H_0: \mu = 0}. Alternatively 'hypothesis' can also be a custom contrast matrix \eqn{T \in \mathbb{R}^{d\times d}}
#'
#' Subjects with missing values will be disregarded.
#'
#' @returns Returns a list with class "hdrm1" with the components
#' @return \item{data}{the initial input data. If a data.frame was given it will have been transformed to a matrix.}
#' @return \item{f}{the degrees of freedom \eqn{f} of the distribution of the test statistic}
#' @return \item{statistic}{the test statistic \eqn{W}}
#' @return \item{tau}{the convergence parameter \eqn{\tau}}
#' @return \item{H0}{the Nullhypothesis tested. Will be 'custom' if `hypothesis` was given as a matrix.}
#' @return \item{H}{the hypothesis-matrix used.}
#' @return \item{p.value}{the p-value of the test statistic}
#' @return \item{critical.value}{the critical value depending on \eqn{\alpha}}
#' @return \item{dim}{a vector of length 2, giving the dimensions \eqn{d \times N} of the data}
#' @return \item{removed.cases}{number of subjects removed for having missing values.}
#'
#' \insertNoCite{Pauly2015}{hdrm}
#' @references \insertAllCited
#'
#' @method hdrm1 data.frame
#' @export
hdrm1.data.frame <- function(data, hypothesis = c("flat", "equal"), value, subject, factor,...){

  # df enthält nur die relevanten Spalten
  df <- data.frame(value = data[[value]], # auf diese weise können value, subject, factor sowohl die namen als auch die nummer der spalte sein
                   subject = data[[subject]],
                   factor = data[[factor]])
  # falls eine der spalten nicht gefunden werden kann ist diese NULL und df kürzer als 3:
  if(ncol(df) != 3) stop("could not find at least one of the columns")

  # nicht benutzte Faktorlevel löschen
  df$subject <- droplevels(as.factor(df$subject))
  df$factor <- droplevels(as.factor(df$factor))
  stopifnot(is.numeric(df$value))

  # Überprüfen, ob alle dimensionen gleich sind
  if(length(unique(table(df$subject))) != 1) stop("All subjects must have the same number of dimensions")
  if(length(unique(table(df$factor))) != 1) stop("Your dimensions dont make sense") # der fehler sollte anders heißen
  d <- nlevels(df$factor)
  N <- nlevels(df$subject)

  ### Matrix bauen
  M <- matrix(NA, d, N)
  for (i in 1:N) { # baue Matrix spaltenweise/füge in jedem Schritt ein Individuum in die i-te Spalte ein
    M[ ,i] <- df$value[df$subject == i][order(df$factor[df$subject == i])]
  }

  ## Funktionsaufruf
  out <- list(data = df)
  out <- c(out, hdrm1_internal(X = M, hypothesis = hypothesis, alpha = 0.05))
  class(out) <- "hdrm1"
  return(out)
}



# Print Funktionen --------------------------------------------------------

#' @method print hdrm1
#' @export
print.hdrm1 <- function(x,...){
  cat("\n")
  cat("          One Group Repeated Measure
        \nAnalysis of", x$dim[2], "subjects in", paste0(x$dim[1]), "dimensions:",
      "\nW =", x$statisitc, " f =", x$f, " p.value =", round(x$p.value, 5),
      "\nNull-Hypothesis:", x$H,
      "\nConvergence parameter \u03c4 =", x$tau)
  cat("\n")
}
