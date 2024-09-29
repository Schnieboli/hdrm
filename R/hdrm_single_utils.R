
# check_criteria_single - ueberprueft, ob alles so ist, wie es sein soll
#
# Input:    X - Datenmatrix X
#  hypothesis - hypothesis
#
# Output: bricht den Vorgang mit Fehlermeldung ab
#
# Motivation: auf diese Weise sind alle Kriterien an einer Stelle und falls ein
# neues Kriterium auftaucht, muss es nicht in beiden Methoden geaendert werden...
#' @keywords internal
check_criteria_single <- function(X, hypothesis){

  d <- nrow(X)
  N <- ncol(X)

  if(d < 2) stop("at least 2 dimensions are needed", call. = FALSE)
  if(!is.character(hypothesis) & !is.matrix(hypothesis)) stop("hypothesis must be a character or a matrix", call. = FALSE)
  if(is.vector(hypothesis) & length(hypothesis) > 1) warning("length of hypothesis is > 1, only first element used.", call. = FALSE)

}



# trace estimators --------------------------------------------------------


# #' @keywords internal
# B0 <- function(X){
#   return(sum(X^2) / ncol(X))
# }


# #' @keywords internal
# B2 <- function(X){
#   N <- ncol(X)
#   S <- (colSums(X^2))^2
#   M <- sum(X[,1]^2)^2
#   for (i in 2:N) {
#     Y <-  X[, c(i:N, 1:(i-1))] * X
#     S <- S + (colSums(Y))^2
#     M <- M + crossprod(X[, i])^2
#   }
#   return(sum(S) - M)
# }


# #' @keywords internal
# B3 <- function(X){
#   N <- ncol(X)
#   S <- 0
#   for (k in 1:(N-2)) {
#     for (l in (k+1):(N-1)) {
#       for (r in (l+1):N) {
#         S <- S + (crossprod(X[, k], X[, l]) * crossprod(X[, l], X[, r]) * crossprod(X[, r], X[, k]))
#       }
#     }
#   }
#   return(sum(S))
# }
