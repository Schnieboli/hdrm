# Da es zwei interne Funktionen gibt, ist hier zu besseren übersicht eine funktion,
# die die hypothesenmatrizen extrahiert -> so muss nicht mehr in beiden
# internen Funktionen etwas geändert werden...

### get_hypothesis_mult - erstellt die Hypothesenmatrix fuer den mehrfaktoriellen Fall
#
# Input:
#  hypothesis - entweder ein character oder eine Liste mit Komponenten TW und TS
#  a - eine Zahl, Anzahl der Gruppen
#  d - eine Zahl, Anzahl der Dimensionen
#
# Output: eine Liste mit komponenten TW und TS
#  TW - Matrix TW
#  TS - Matrix TS

#' @keywords internal
get_hypothesis_mult <- function(hypothesis, a, d){
  stopifnot(is.list(hypothesis) | is.character(hypothesis))
  # Fall: hypothesis ist eine Liste
  if(is.list(hypothesis)){
    # ueberpruefen, ob es die gesuchten listenelemente gibt
    if(is.null(hypothesis$TW)) stop("no list entry with name 'TW' in hypothesis", call. = FALSE)
    if(is.null(hypothesis$TS)) stop("no list entry with name 'TS' in hypothesis", call. = FALSE)
    # Liste "auspacken"
    TW <- hypothesis$TW
    TS <- hypothesis$TS
    ## ueberpruefen, ob Matrizen richtige dimensionen haben
    stopifnot(is.matrix(TW), is.matrix(TS), is.numeric(TW), is.numeric(TS))
    stopifnot(dim(TW) == c(a,a), dim(TS) == c(d,d))
    # idempotenz
    if(!all.equal(TW^2,TW) | !all.equal(TS^2,TS) | any(TW != t(TW)) | any(TS != t(TS))) stop("TW and TS must be idempotent")
  }
  # Fall: hypothesis ist character
  if(is.character(hypothesis)){
    if(length(hypothesis) > 1) warning("length of hypothesis is > 1, only first element used.", call. = FALSE)
    if(!all(hypothesis %in% c("whole","sub","interaction"))) stop("hypothesis must be one of c('whole', 'sub', 'interaction') or a list with components TW and TS")
    if(hypothesis[1] == "whole"){
      # Matrizen erstellen
      TW <- diag(a) - matrix(1/a,a,a)
      TS <- matrix(1,d,d)/d
    }
    if(hypothesis[1] == "sub"){
      # Matrizen erstellen
      TW <- matrix(1,a,a)/a
      TS <- diag(d) - matrix(1/d,d,d)
    }
    if(hypothesis[1] == "interaction"){
      # Matrizen erstellen
      TW <- diag(a) - matrix(1,a,a)/a
      TS <- diag(d) - matrix(1,d,d)/d
    }
  }
  return(list(TW = TW, TS = TS))
}
