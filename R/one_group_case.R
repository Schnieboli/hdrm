
# interne Funktion für die Berechnung
#' @keywords internal
hdrm1_internal <- function(X, hypothesis, alpha, na.action = "na.omit"){

  # Matrix X kommt eingegeben als: dim(X) = c(d,N)
  ## Fehlende Werte
  stopifnot(na.action %in% c("na.omit"))
  N_with_NA <- ncol(X)
  # na.omit entfernt alle Zeilen mit fehlenden Werten, aber wir wollen alle Spalten entfernen
  # -> doppelt transponieren
  X <- t(stats::na.omit(t(X)))


  ## Schätzer definieren
  N <- ncol(X)
  d <- nrow(X)
  stopifnot(is.matrix(X), is.numeric(X), N >= 3)


  ### Hypothesenmatrizen
  TM <- NA
  if(is.matrix(hypothesis) & all(dim(hypothesis) == c(d,d))) TM <- hypothesis
  if(is.character(hypothesis)){
    if(hypothesis[1] == "equal") TM <- diag(d)
    if(hypothesis[1] == "flat") TM <- diag(d) -  matrix(1/d, d, d)
  }

  # wenn keiner der oberen fälle zutrifft oder ein NA in TM ist, dann breche ab
  if(any(is.na(TM))) stop("Please specify valid hypothesis.")


  ### Teststatistik Q
  XT <- TM %*% X
  Xquer <- rowMeans(XT)
  Qn = N * sum(Xquer * Xquer)


  ### Schätzer berechnen
  spurNormal <- B0(XT)
  spurQuadrat <- B2_cpp(XT)
  spurHoch3 <- B3_cpp(XT)/choose(N,3) # das muss so, weil ich nicht weiss, wie man in cpp bionom ausrechnet

  ### Teststatistik W
  W <- (Qn - spurNormal) / sqrt(2*spurQuadrat)

  ### Verteilungsparameter f schätzen
  f <- spurQuadrat^3 / spurHoch3^2

  ### p-Wert
  p.value <- 1 - stats::pchisq(W * sqrt(2 * f) + f, df = f)


  ### Ausgabe
  if(is.matrix(hypothesis)) H <- "custom"
  else H <- paste0(hypothesis[1], " time profile")

  L <- list(data = X, # X ohne NAs -> muss transponiert sein, da X am Anfang transponiert wurde
            f = f,
            statisitc = W,
            tau = 1/f,
            H0 = H,
            hypothesis.matrix = TM,
            p.value = p.value,
            dim = c(d = d, N = N),
            removed.cases = N_with_NA - N
  )
  return(L)

}
