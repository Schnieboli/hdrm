#' @keywords internal
hdrm1_internal <- function(X, hypothesis,...){

  ## Dimensionen definieren
  N <- ncol(X)
  d <- nrow(X)


  ### Hypothesenmatrizen
  TM <- NA
  # hypothesis is matrix
  if(is.matrix(hypothesis)) TM <- hypothesis
  #hypothesis ist character
  if(is.character(hypothesis)){
    if(hypothesis[1] == "flat") TM <- diag(d) -  matrix(1/d, d, d)
    else stop("The only legal character is 'flat'. Other hypotheses can be specified by a matrix.")
  }

  # wenn keiner der oberen faelle zutrifft oder ein NA in TM ist oder TM nicht
  # numeric ist, dann breche ab

  ## Dimensionen pruefen
  if(!all(dim(TM) == c(d,d))) stop("hypothesis must be a quadratic matrix with the number of rows equal to the number of dimension d.")
  if(any(is.na(TM)) | !is.numeric(TM)) stop("Please specify valid hypothesis.")
  # Symmetrie und Idempotenz pruefen
  # Symmetrie
  if((mean(t(TM) - TM) >= sqrt(.Machine$double.eps))) warning(paste0("TM is not symmetric (mean difference = ", mean(t(TM) - TM),"). This will likely have a big influence on the test result!"))
  # Idempootenz
  # checke ob all.equal TRUE ist -> wenn ja, dann milde warnung
  if((mean(TM%*%TM - TM) >= sqrt(.Machine$double.eps))) warning(paste0("TM is not idempotent (mean difference = ", mean(TM%*%TM - TM),"). This will likely have a big influence on the test result!"))

  ### Teststatistik Q
  XT <- TM %*% X
  Xquer <- rowMeans(XT)
  Qn = N * sum(Xquer^2)


  ### Schaetzer berechnen
  spurNormal <- B0_cpp(XT)
  spurQuadrat <- B2_cpp(XT)
  spurHoch3 <- B3_cpp(XT)/choose(N,3) # das muss so, weil ich nicht weiss, wie man in cpp bionom ausrechnet

  ### Teststatistik W
  W <- (Qn - spurNormal) / sqrt(2*spurQuadrat)

  ### Verteilungsparameter f schaetzen
  f <- max(1, spurQuadrat^3 / spurHoch3^2)

  ### p-Wert
  p.value <- 1 - stats::pchisq(W * sqrt(2 * f) + f, df = f)


  L <- list(f = f,
            statisitc = W,
            tau = 1/f,
            H = TM,
            hypothesis = ifelse(is.character(hypothesis), hypothesis[1], "custom"),
            p.value = p.value,
            dim = c(d = d, N = N)
  )
  return(L)

}
