#' Test for one dimensional repeated measures
#'
#' @description
#' S3 generic methods for performing a one sample test for high dimensional repeated measures
#' @param data either a matrix or a data.frame with N subjects in rows and d factor levels in columns
#' @param hypothesis either one of c("equal","flat") or a matrix of dimensions c(d,d)
#' @param alpha alpha level used for calculating the critical value
#' @param na.action determines how to deal with individuals with missing values. only `copmplete.obs` is supported right now.
#' @return \item{data}{the initial input data}
#' @return \item{f}{the degrees of freedom f of the Distribution of the test statistic}
#' @return \item{statistic}{the value of the test statistic}
#' @return \item{tau}{the convergence parameter tau}
#' @return \item{hypothesi}{the type of hypothesis}
#' @return \item{critical.value}{the critical value depending on alpha}
#' @return \item{p.value}{the $p$-value of the test statistic}
#' @return \item{dim}{a vector of length 2, giving the dimensions c(d,N) of the data}
#' @return \item{na.action}{na.action}
#' @return \item{removed.cases}{number of cases removed for having missing values}
#' @export
hdrm1 <- function(data, formula,...){
  UseMethod("hdrm1")
}

#' @method hdrm1 default
#' @export
hdrm1.default <- function(data,...){
  stop("Your data needs to be either a matrix or a data.frame")
}

#' @method hdrm1 matrix
#' @export
hdrm1.matrix <- function(data, hypothesis = c("flat","equal"), alpha = 0.05, na.action = "na.omit"){ # na.action muss raus -> klappt nicht so, wie ich mir das vorgestellt habe und ist ja auch egal, weil ja eh nur na.omit funktionieren soll
  ## Mat wird so eingegeben, dass Individuen in Spalten und Dimension in Zeilen ist -> für na.action besser
  return(hdrm1_internal(X = data, hypothesis = hypothesis, alpha = alpha, na.action = na.action))
}


#' @method hdrm1 data.frame
#' @export
hdrm1.data.frame <- function(formula, data, hypothesis = c("flat","equal"), alpha = 0.05, na.action = "na.omit"){

  #### Fall: nur ein data.frame gegeben, kein formelobjekt
  if(missing(formula) || (is.data.frame(formula) & missing(data))){

    ## falls nur ein df als formula eingegeben wird, dann checkt die funktion, dass es sich um data handeln soll
    if(!missing(formula) && (is.data.frame(formula) & missing(data))) data <- formula

    if(all(c("subject","value","factor") %in% names(data))){
      value <- data$value
      stopifnot(is.numeric(value)) # value muss numeric sein
      subject <- droplevels(as.factor(data$subject))
      factor <- droplevels(as.factor(data$factor))
    } else stop("data must either be a data.frame with the columns 'value', 'subject' and 'factor' or the columns have to be specified by formula")
  }

  #### Fall: formula und data gegeben
  ### Formelobjekt der Form value ~ subject + dimension
  if(!missing(formula) & (rlang::is_formula(formula) & is.data.frame(data))){
    # genau 2 unabhängige variablen
    if(length(attr(terms(formula), "term.labels")) != 2) stop("formula must be of the form 'value ~ subject + factor'")

    # aus Formel extrahieren
    value <- data[[rlang::f_lhs(formula)]]
    subject <- droplevels(as.factor(data[[attr(terms(formula), "term.labels")[1]]]))
    factor <- droplevels(as.factor(data[[attr(terms(formula), "term.labels")[2]]]))

  }

  ## N und d definieren
  N <- nlevels(subject)
  d <- nlevels(factor)

  ## Kontrolle, ob alle Individuen gleich viele Dimensionen haben
  if(!all(table(factor) == table(factor)[1])) stop("All subjects need to have the same number of dimensions")

  ### Matrix bauen
  Mat <- matrix(NA, N, d)
  for (i in 1:N) {
    Mat[i,] <- value[subject == i][order(factor[subject == i])]
  }
  ## Mat wird so eingegeben, dass Individuen in Spalten und Dimension in Zeilen ist -> für na.action besser
  return(hdrm1_internal(t(Mat), hypothesis = hypothesis, alpha = alpha, na.action = na.action))
}


#' @keywords internal
hdrm1_internal <- function(X, hypothesis, alpha, na.action = na.action){

  # Matrix X kommt eingegeben als: dim(X) = c(d,N)
  ## Fehlende Werte
  stopifnot(na.action %in% c("na.omit","na.exclude","na.fail"))
  N_with_NA <- dim(X)[2]
  X <- t(na.omit(X)) # -> hier muss dann iwie mit attr gearbeitet werden -> auf diese weise kann man dann auch sehen, welche Individuen


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
  spurNormal <- B0(X = XT, N = N)/N
  spurQuadrat <- B2(X = XT, N = N)/(N*(N-1))
  spurHoch3 <- B3(X = XT, N = N)/choose(N,3)

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

  ### von Paavo -> bin mir nicht sicher, ob das richtig ist, deswegen lass ich mal das alte...
  ## p-Wert
  #p.value = 1 - dchisq(W * sqrt(2 * f) + f, df = f)


  ### Ausgabe
  if(is.matrix(hypothesis)) H <- "custom"
  else H <- paste0(hypothesis[1], " time profile")

  L <- list(data = X, # hier vielleicht X ohne NAs?
            f = f,
            statisitc = W,
            tau = 1/f,
            hypothesis = H,
            hypothesis.matrix = TM,
            p.value = p.value,
            critical.value = critical.value,
            dim = c(d = d, N = N),
            na.action = na.action,
            removed.cases = N_with_NA - N
  )
  class(L) <- c("hdrm1","list")
  return(L)

}


# generische print-Funktion -----------------------------------------------

#' @method print hdrm1
#' @export
print.hdrm1 <- function(x,...){
  cat("\n")
  cat("       One Dimensional Repeated Measure
        \nAnalysis of", x$dim[2], "individuals", paste0("(", x$removed.cases, " removed)"), "and", paste0(x$dim[1]), "dimensions:",
      "\nW =", x$statisitc, " f =", x$f, " p.value =", x$p.value,
      "\nNull-Hypothesis:", x$hypothesis,
      "\nConvergence parameter \u03c4 =", x$tau)
  cat("\n")
}

#' @method summary hdrm1
#' @export
summary.hdrm1 <- function(X,...){
  print(X,...)
}



