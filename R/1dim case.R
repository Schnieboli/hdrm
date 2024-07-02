#' Test for one dimensional repeated measures
#'
#' @description
#' S3 generic methods for performing a one sample test for high dimensional repeated measures
#' @param data either a matrix or a data.frame with N subjects in rows and d factor levels in columns
#' @param hypothesis either one of c("equal","flat") or a matrix of dimensions c(d,d)
#' @param alpha alpha level used for calculating the critical value
#' @return \item{data}{the initial input data}
#' @return \item{f}{the degrees of freedom f of the Distribution of the test statistic}
#' @return \item{statistic}{the value of the test statistic}
#' @return \item{tau}{the convergence parameter tau}
#' @return \item{hypothesi}{the type of hypothesis}
#' @return \item{critical.value}{the critical value depending on alpha}
#' @return \item{p.value}{the $p$-value of the test statistic}
#' @return \item{dim}{a vector of length 2, giving the dimensions c(d,N) of the data}
#' @export
hdrm1 <- function(data, hypothesis = c("equal","flat"), alpha = 0.05){
  UseMethod("hdrm1")
}

#' @method hdrm1 default
#' @export
hdrm1.default <- function(data,...){
  stop("Your data needs to be either a matrix or a data.frame")
}

#' @method hdrm1 matrix
#' @export
hdrm1.matrix <- function(data, hypothesis, alpha){
  X <- as.matrix(data)
  return(hdrm1_internal(X = t(X), hypothesis = hypothesis, alpha = alpha))
}

#' @keywords internal
hdrm1_internal <- function(X, hypothesis, alpha){

  ## Schätzer definieren
  N <- dim(X)[1]
  d <- dim(X)[2]
  stopifnot(is.matrix(X), is.numeric(X), N >= 3)



  ### Hypothesenmatrizen
  if(is.matrix(hypothesis) & all(dim(hypothesis) == c(d,d))) TM <- hypothesis
  if(hypothesis[1] == "equal") TM <- diag(d)
  if(hypothesis[1] == "flat") TM <- diag(d) -  matrix(1/d, d, d)


  ### Teststatistik Q
  XT <- X %*% TM
  Xquer <- colMeans(XT)
  Qn = N * sum(Xquer * Xquer)


  ### Schätzer berechnen
  spurNormal <- B0(X = XT, N = N)
  spurQuadrat <- B2(X = XT, N = N)
  spurHoch3 <- B3(X = XT, N = N)

  ### Teststatistik W
  W <- (Qn - spurNormal) / sqrt(2*spurQuadrat)

  ### Verteilungsparameter f schätzer
  f <- spurQuadrat^3 / spurHoch3^2

  ### Kritischer Wert und p-Werte
  critical.value <- (qchisq(1- alpha, df = f) - f) / sqrt(2*f)

  # pWert_Kf gibt die für ein alpha und param = f den Abstand des Quantils zur Teststatistik aus
  pWert_Kf <- function(p, param, statistic) abs( ((qchisq(1-p, param) - param)/sqrt(2*param)) - statistic )
  # hier wird auf[0, 1] das minimum der funktion pWert_Kf gesucht
  # tol = .Machine$double.eps für maximale Accuracy -> sonst kommt ab einer bestimmten Extremität immer der gleiche Wert raus
  p.value <- optimise(pWert_Kf, interval = c(0,1), param = f, statistic = W, tol = .Machine$double.eps)$minimum

  ### Ausgabe
  if(is.matrix(hypothesis)) H <- "custom"
  else H <- paste0(hypothesis[1], " time profile")

  L <- list(data = X,
            f = f,
            statisitc = W,
            tau = 1/f,
            hypothesis = H,
            p.value = p.value,
            critical.value = critical.value,
            dim = c(d = d, N = N)
  )
  class(L) <- c("hdrm1","list")
  return(L)

}


# generische print-Funktion -----------------------------------------------

#' @method print hdrm1
#' @export
print.hdrm1 <- function(X,...){
  cat("\n")
  cat("       One Dimensional Repeated Measure
      \nW =", X$statisitc, " f =", X$f, " p.value =", X$p.value,
      "\nNull-Hypothesis:", X$hypothesis,
      "\nConvergence parameter \u03c4 =", X$tau,
      "\n")
}
