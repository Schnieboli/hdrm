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
    if(hypothesis[1] == "equal") TM <- diag(d)
    if(hypothesis[1] == "flat") TM <- diag(d) -  matrix(1/d, d, d)
  }

  # wenn keiner der oberen faelle zutrifft oder ein NA in TM ist oder TM nicht
  # numeric ist, dann breche ab
  if(any(is.na(TM)) | !is.numeric(TM)) stop("Please specify valid hypothesis.")
  # Symmetrie und Idempotenz pruefen
  # Symmetrie
  if(all(TM != t(TM))) stop("hypothesis must be symmetric.")
  # Idempootenz
  if(!identical(TM%*%TM,TM)){ # wenn nicht identisch, dann checke ob all.equal TRUE ist
    # wenn ja, dann milde warnung
    if((mean(TM%*%TM - TM) < sqrt(.Machine$double.eps))) warning(paste0("TM is not exactly idempotent (mean difference = ", mean(TM%*%TM - TM),"). The effect on the test result is probably be negligible."), call. = FALSE)
    # wenn nein, dann starke warnung
    else warning(paste("TM is not idempotent (mean difference = ", mean(TM%*%TM - TM),"). This will affect the test result heavily!"))
  }

  ### Teststatistik Q
  XT <- TM %*% X
  Xquer <- rowMeans(XT)
  Qn = N * sum(Xquer^2)


  ### Schaetzer berechnen
  spurNormal <- B0(XT)
  spurQuadrat <- B2_cpp(XT)
  spurHoch3 <- B3_cpp(XT)/choose(N,3) # das muss so, weil ich nicht weiss, wie man in cpp bionom ausrechnet

  ### Teststatistik W
  W <- (Qn - spurNormal) / sqrt(2*spurQuadrat)

  ### Verteilungsparameter f schaetzen
  f <- spurQuadrat^3 / spurHoch3^2

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
