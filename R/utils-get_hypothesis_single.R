get_hypothesis_single <- function(hypothesis, AM, d) {
  
  
  if(is.character(hypothesis) && length(hypothesis) == 1L && hypothesis == "flat"){
    ## checked, if hypothesis is legal character
    TM <- diag(d) - matrix(1 / d, nrow = d, ncol = d)
  } 
  else if(is.matrix(hypothesis) && is.numeric(hypothesis) && 
          !any(is.na(hypothesis)) && all(is.finite(hypothesis))){
    
    if(any(dim(hypothesis) != d)) ## check if symmetrical matrix
      stop(paste0("The hypothesis matrix must be a ", d, " by ", d, " matrix."),
           call. = FALSE)
    
    TM <- hypothesis
    
    ## check if idempotent
    tol <- sqrt(.Machine$double.eps)
    
    symmetry_error <- max(abs(TM - t(TM)))
    idempotence_error <- max(abs(TM %*% TM - TM))
    
    if (symmetry_error > tol)
      stop(paste0("The hypothesis matrix must be symmetric. Maximum deviation: ",
                  signif(symmetry_error, 4), "." ), call. = FALSE)
    
    if (idempotence_error > tol)
      stop(paste0("The hypothesis matrix must be idempotent. Maximum deviation: ",
                  signif(idempotence_error, 4), "."),call. = FALSE)
    
    if (qr(TM, tol = tol)$rank == 0L)
      stop("The hypothesis matrix must have positive rank.", call. = FALSE)
    
  }
  ## else: wrong input value
  else stop("'hypothesis' must be 'flat' or a numeric projection matrix containing only finite, non-missing values.",
            call. = FALSE)
  
  ## compute compact version of hypothesis matrix
  TMalt <- if (AM) {
    MSrootcompact(TM)
  } else {
    TM
  }
  list(TM = TM, TMalt = TMalt)
  
}
