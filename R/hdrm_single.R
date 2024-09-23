#' @keywords internal
hdrm1_internal <- function(X, hypothesis,...){

  # Matrix X kommt eingegeben als: dim(X) = c(d,N)
  ## Fehlende Werte
  N_with_NA <- ncol(X)
  # na.omit entfernt alle Zeilen mit fehlenden Werten, aber wir wollen alle Spalten entfernen
  # -> doppelt transponieren
  X <- t(stats::na.omit(t(X)))


  ## Dimensionen definieren
  N <- ncol(X)
  d <- nrow(X)
  # Warnung fuer fehlende Werte
  if(N_with_NA > N) warning("Subjects with missing values dropped", call. = FALSE)


  ### Hypothesenmatrizen
  TM <- NA
  # hypothesis is matrix
  if(is.matrix(hypothesis) & all(dim(hypothesis) == c(d,d))) TM <- hypothesis
  #hypothesis ist character
  if(is.character(hypothesis)){
    if(hypothesis[1] == "flat") TM <- diag(d) -  matrix(1/d, d, d)
    else stop("The only legal character is 'flat'. Other hypotheses can be specified by a matrix.")
  }

  # wenn keiner der oberen faelle zutrifft oder ein NA in TM ist oder TM nicht
  # numeric ist, dann breche ab
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
            dim = c(d = d, N = N),
            removed.cases = N_with_NA - N
  )
  return(L)

}
