#' @keywords internal
get_hypothesis_mult <- function(hypothesis, a, d){
  stopifnot(is.list(hypothesis) | is.character(hypothesis))

  if(is.list(hypothesis)){
    TW <- hypothesis$TW
    TS <- hypothesis$TS
    stopifnot(is.matrix(TW), is.matrix(TS), is.numeric(TW), is.numeric(TS))
    stopifnot(dim(TW) == c(a,a), dim(TS) == c(d,d))
    if(!all.equal(TW^2,TW) | !all.equal(TS^2,TS) | any(TW != t(TW)) | any(TS != t(TS))) stop("TW and TS must be idempotent")
  }

  if(is.character(hypothesis)){
    if(hypothesis[1] == "whole"){
      TW <- diag(a) - matrix(1/a,a,a)
      TS <- matrix(1,d,d)/d
    }
    if(hypothesis[1] == "sub"){
      TW <- matrix(1,a,a)/a
      TS <- diag(d) - matrix(1/d,d,d)
    }
    if(hypothesis[1] == "interaction"){
      TW <- diag(a) - matrix(1,a,a)/a
      TS <- diag(d) - matrix(1,d,d)/d
    }
  }

  return(list(TW = TW, TS = TS))
}
