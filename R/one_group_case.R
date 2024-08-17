#' Test for one group high dimensional repeated measures
#'
#' @description
#' S3 generic methods for performing a one sample test for one group repeated measure data
#' @param data either a matrix with subjects in columns and factor levels in rows or a data.frame of longtable format
#' @param formula a formula object that specifies the corresponding cols if data is a data.frame. Must be of the form `value ~ subject + factor`
#' @param hypothesis One of 'equal' or 'flat' for H0: equal/flat time profile. Can also be a matrix of dimensions \eqn{d\times d}
#' @param alpha alpha level used for calculating the critical value
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
#' @aliases hdrm1.matrix
#' @aliases hdrm1.data.frame
#' @aliases hdrm1.list
#' @usage hdrm1.matrix(data, hypothesis, alpha = 0.05,...)
#' @usage hdrm1.data.frame(data, formula, hypothesis, alpha = 0.05,...)
#' @usage hdrm1.list(data, hypothesis, alpha, ...)
#' @example examples_one_group_case.R
#'
#' @references Pauly, M., Ellenberger, D., & Brunner, E. (2015). Analysis of high-dimensional one group repeated measures designs. Statistics, 49(6), 1243–1261. https://doi.org/10.1080/02331888.2015.1050022
#' @export
hdrm1 <- function(data, formula, hypothesis = c("flat", "equal"), alpha = 0.05,...){
  UseMethod("hdrm1")
}

#' @method hdrm1 default
#' @export
hdrm1.default <- function(data, formula, hypothesis = c("flat","equal"), alpha=0.05,...){
  stop("Your data needs to be either a matrix or a data.frame")
}


#' @method hdrm1 matrix
#' @export
hdrm1.matrix <- function(data, formula, hypothesis = c("flat","equal"), alpha = 0.05,...){ # na.action muss raus -> klappt nicht so, wie ich mir das vorgestellt habe und ist ja auch egal, weil ja eh nur na.omit funktionieren soll

  return(hdrm1_internal(X = data, hypothesis = hypothesis, alpha = alpha))
}

#' @method hdrm1 list
#' @export
hdrm1.list <- function(data, formula, hypothesis = c("flat","equal"), alpha = 0.05,...){
  dims <- sapply(data, dim)
  stopifnot(dim(t(dims)) == c(1,2)) # checken, ob alle matrizen die gleiche dimension haben
  stopifnot(all(sapply(data, is.matrix)) | all(sapply(data, is.data.frame)))

  mat <- matrix(0, dims[1,2],0)
  for (i in 1:length) {
    mat <- cbind(mat, data[[i]])
  }

  return(hdrm1_internal(X = mat, hypothesis = hypothesis, alpha = alpha))
}




#' @method hdrm1 data.frame
#' @export
hdrm1.data.frame <- function(data, formula, hypothesis = c("flat", "equal"), alpha = 0.05,...){

  # df enthält nur die relevanten Spalten
  df <- stats::model.frame(formula = formula, data = data)
  # df soll genau 3 Spalten haben
  stopifnot(ncol(df) == 3)
  names(df) <- c("value", "subject", "factor")
  # nicht benutzte Faktorlevel löschen
  df$subject <- droplevels(as.factor(df$subject))
  df$factor <- droplevels(as.factor(df$factor))

  # Überprüfen, ob alle dimensionen gleich sind
  if(length(unique(table(df$subject))) != 1) stop("All subjects must have the same number of dimensions")
  if(length(unique(table(df$factor))) != 1) stop("Your dimensions dont make sense")
  d <- unique(table(df$subject))
  N <- unique(table(df$factor))

  ### Matrix bauen
  M <- matrix(NA, d, N)
  for (i in 1:N) { # baue Matrix spaltenweise/füge in jedem Schritt ein Individuum in die i-te Spalte ein
    M[ ,i] <- df$value[df$subject == i][order(df$factor[df$subject == i])]
  }

  ## Funktionsaufruf
  hdrm1_internal(X = M, hypothesis = hypothesis, alpha = alpha)
}



#' interne Funktion für die Berechnung
#' @keywords internal
hdrm1_internal <- function(X, hypothesis, alpha, na.action = "na.omit"){

  # Matrix X kommt eingegeben als: dim(X) = c(d,N)
  ## Fehlende Werte
  stopifnot(na.action %in% c("na.omit"))
  N_with_NA <- dim(X)[2]
  X <- t(stats::na.omit(X)) # -> hier muss dann iwie mit attr gearbeitet werden -> auf diese weise kann man dann auch sehen, welche Individuen entfernt wurden


  ## Schätzer definieren
  N <- dim(X)[1]
  d <- dim(X)[2]
  stopifnot(is.matrix(X), is.numeric(X), N >= 3)


  ### Hypothesenmatrizen
  TM <- NA
  if(is.matrix(hypothesis) & all(dim(hypothesis) == c(d,d))) TM <- hypothesis
  if(hypothesis[1] == "equal") TM <- diag(d)
  if(hypothesis[1] == "flat") TM <- diag(d) -  matrix(1/d, d, d)
  # wenn keiner der oberen fälle zutrifft oder ein NA in TM ist, dann breche ab
  if(any(is.na(TM))) stop("Please specify valid hypothesis.")


  ### Teststatistik Q
  XT <- X %*% TM
  Xquer <- colMeans(XT)
  Qn = N * sum(Xquer * Xquer)


  ### Schätzer berechnen
  spurNormal <- B0(XT)/N
  spurQuadrat <- B2_cpp(XT)/(N*(N-1))
  spurHoch3 <- B3_cpp(XT)/choose(N,3)

  ### Teststatistik W
  W <- (Qn - spurNormal) / sqrt(2*spurQuadrat)

  ### Verteilungsparameter f schätzen
  f <- spurQuadrat^3 / spurHoch3^2

  ### Kritischer Wert und p-Werte
  critical.value <- (stats::qchisq(1- alpha, df = f) - f) / sqrt(2*f)

  # pWert_Kf gibt die für ein alpha und param = f den Abstand des Quantils zur Teststatistik aus
  pWert_Kf <- function(p, param, statistic) abs( ((stats::qchisq(1-p, param) - param)/sqrt(2*param)) - statistic )
  # hier wird auf[0, 1] das minimum der funktion pWert_Kf gesucht
  # tol = .Machine$double.eps für maximale Accuracy -> sonst kommt ab einer bestimmten Extremität immer der gleiche Wert raus
  p.value <- stats::optimise(pWert_Kf, interval = c(0,1), param = f, statistic = W, tol = .Machine$double.eps)$minimum
#
  ### von Paavo -> bin mir nicht sicher, ob das richtig ist, deswegen lass ich mal das alte...
  ## p-Wert
  #p.value = 1 - dchisq(W * sqrt(2 * f) + f, df = f)


  ### Ausgabe
  if(is.matrix(hypothesis)) H <- "custom"
  else H <- paste0(hypothesis[1], " time profile")

  L <- list(data = t(X), # X ohne NAs -> muss transponiert sein, da X am Anfang transponiert wurde
            f = f,
            statisitc = W,
            tau = 1/f,
            H0 = H,
            hypothesis.matrix = TM,
            p.value = p.value,
            critical.value = critical.value,
            dim = c(d = d, N = N),
            na.action = na.action,
            removed.cases = N_with_NA - N
  )
  class(L) <- c("hdrm1")
  return(L)

}


# generische print-Funktion -----------------------------------------------

#' @method print hdrm1
#' @export
print.hdrm1 <- function(x,...){
  cat("\n")
  cat("          One Group Repeated Measure
        \nAnalysis of", x$dim[2], "individuals", paste0("(", x$removed.cases, " removed)"), "and", paste0(x$dim[1]), "dimensions:",
      "\nW =", x$statisitc, " f =", x$f, " p.value =", round(x$p.value, 5),
      "\nNull-Hypothesis:", x$H,
      "\nConvergence parameter \u03c4 =", x$tau)
  cat("\n")
}

#' @method summary hdrm1
#' @export
summary.hdrm1 <- function(object,...){
  print(object,...)
}



