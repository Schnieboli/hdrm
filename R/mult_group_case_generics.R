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
#' @returns \item{critical.valie}{the critical value based on \eqn{\alpha}}
#' @returns \item{dim}{a named list with components d and N of the input data}
#' @returns \item{groups}{a named list with components a and table giving the number and distribution of groups}
#' @returns \item{removed.cases}{number of incomplete subjects removed.}
#'
#' @references Paavo Sattler, Markus Pauly. "Inference for high-dimensional split-plot-designs: A unified approach for small to large numbers of factor levels." Electronic Journal of Statistics, 12(2) 2743-2805 2018. https://doi.org/10.1214/18-EJS1465
#'
#' @export
hdrm_test <- function(data, group, formula, hypothesis = c("whole","sub","interaction"), alpha = 0.05, B = "500*N",...){
  UseMethod("hdrm_test")
}

#' @method hdrm_test default
#' @export
hdrm_test.default <- function(data, hypothesis = c("whole","sub","interaction"), alpha = 0.05, B = "500*N", ...){
  stop("Your data must be either a list, a data frame or a matrix")
}


#' @method hdrm_test list
#' @export
hdrm_test.list <- function(data, hypothesis = c("whole","sub","interaction"), alpha = 0.05, B = "500*N",...){
  #browser()
  # Was muss sein?
  # in jedem Eintrag eine numerische Matrix
  stopifnot(all(sapply(data, is.matrix)), all(sapply(data, is.numeric)))
  # midestens 2 Einträge
  stopifnot(length(data) >= 2)
  # in jeder Matrix gleich viele Zeilen
  stopifnot(length(unique(sapply(data, function(X)dim(X)[1]))) == 1)
  # in jeder Matrix mind. 6 Spalten
  stopifnot(all(sapply(data, function(X) dim(X)[2]) >= 6))

  ### Matrix bauen
  # group erstellen
  if(is.null(names(data))){
    group <- factor(1:length(data))
  }else{
    group <- as.factor(names(data))
  }

  # Matrix bauen
  X <- data[[1]]
  for (i in 2:length(data)) { # so ist das halt sehr Laufzeitunfreundlich (da immer neuer Speicher für X gesucht werden muss), deswegen mal schauen, ob man das so lassen kann...
    X <- cbind(X, data[[i]])
  }

  ## Aufruf der internen Funktion
  hdrm_test_internal(data = X, group = rep(group, each = dim(data[[1]])[2]), hypothesis = hypothesis, alpha = 0.05, B = "500*N",)
}

#' @method hdrm_test matrix
#' @export
hdrm_test.matrix <- function(data, group, hypothesis = c("whole","sub","interaction"), alpha = 0.05, B = "500*N", ...){
  hdrm_test_internal(data = data, group = group, hypothesis = hypothesis, alpha = alpha, B = B)
}

#' @method hdrm_test data.frame
#' @export
hdrm_test.data.frame <- function(data, formula, hypothesis = c("whole","sub","interaction"), TW, TS, alpha = 0.05, B = "100*N", ...){
  #browser()
  if(missing(formula)){
    # überprüfen, ob alle spalten gegeben sind
    stopifnot(all(c("value","whole","sub","subject") %in% names(data)))
    # Vektoren extrahieren
    df <- data.frame(value = data$value,
                     subject = as.factor(data$subject),
                     whole = as.factor(data$whole),
                     sub = as.factor(data$sub)
    )
  }else{ # formel soll gegeben sein als value ~ subject + whole + sub
    # df enthält nur die relevanten Spalten
    df <- stats::model.frame(formula = formula, data = data)
  }
  stopifnot(length(names(df)) == 4)

  ## Dimensionen
  N <- nlevels(df[,2])
  a <- nlevels(df[,3])
  d <- nlevels(df[,4])

  ### Bedingungen überprüfen
  # alle Individuen haben gleich viele dimensionen
  stopifnot(length(unique(table(df$sub))) == 1) # von Jessica geholfen
  # in jeder Gruppe mind 6 individuen
  # -> wird eig in der internen Funktion auch noch mal abgefragt...
  # mind 2 Gruppen
  stopifnot(a >= 2)

  # Matrix bauenn
  X <- matrix(NA, d, 0)
  M <- matrix(NA, d, N)
  L <- split(df, df[,3])
  group <- numeric(0)


  for(j in 1:a){
    temp = droplevels(L[[j]])
    M <- matrix(0, d, nlevels(temp[,2]))
    group <- c(group, rep(j, nlevels(temp[,2])))
    k = 1
    for (i in levels(temp[,2])) {
      M[,k] <- temp[,1][temp[,2] == i]
      k <- k + 1
    }
    X <- cbind(X,M)
  }


  ## Funktionsaufruf
  return(hdrm_test_internal(data = X, group = group, hypothesis = hypothesis, alpha = alpha, B = B))
}
