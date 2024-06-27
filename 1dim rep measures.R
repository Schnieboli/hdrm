library(mvtnorm)


# Daten -------------------------------------------------------------------

N <- 20
means <- rep(1,N)
d <- 30
X <- t(rmvnorm(n = d, mean = means))
Y <- replicate(N, rnorm(d,sd = 0.1))


# Schätzer ----------------------------------------------------------------

# # hier wird nur das eine individuum eingegeben
# B0 <- function(X, N){
#   return(sum(X * X)/N)
# }
# 
# # hier wird die ganze matrix eigegeben
# B2 <- function(X){
#   S <- 0
#   for (k in 1:N) {
#     for(l in 1:N){
#       if(k != l) S <- S + (sum(X[k, ] * X[l, ])^2)
#     }
#   }
#   return(S / (n*(n-1)))
# }
# 
# 
# B3 <- function(X){
#   C <- combn(N, 2)
#   S <- 0
#   for(i in 1:choose(N,2)){
#     k <- C[1,i]
#     l <- C[2,i]
#     for (i in l:N) {
#       S <- S + (sum(X[k,] * X[l, ]) * sum(X[l, ] * X[r, ]) * sum(X[r, ] * X[k, ]))
#     }
#   }
#   return(S/choose(N, 3))
# }



# Funktion ----------------------------------------------------------------


hdrm1 <- function(X, hypothesis = c("equal", "flat"), alpha = 0.05){
  
  ## Schätzer definieren
  N <- dim(X)[1]
  d <- dim(X)[2]
  stopifnot(is.matrix(X), is.numeric(X), N >= 3)
  
  
  B0 <- function(X, N){
    S <- 0
    for (i in 1:N) {
      S <- S + sum(X[i, ]^2)
    }
    return(S/N)
  }
  # hier wird die ganze matrix eigegeben
  B2 <- function(X, N){
    S <- 0
    for (k in 1:N) {
      for(l in 1:N){
        if(k != l) S <- S + (sum(X[k, ] * X[l, ])^2)
      }
    }
    return(S / (N*(N-1)))
  }
  
### müsste eig richtig sein, aber ergibt was anderes als B3_2
### und B3_2 ist auf jeden Fall richtig
### hier müssen am ende von C spalten weggelassen werden, da sonst r>N wird
  # B3 <- function(X, N){
  #   C <- combn(N, 2)
  #   S <- 0
  #   browser()
  #   for(i in 1:choose(N,2)){
  #     k <- C[1,i]
  #     l <- C[2,i]
  #     for (r in (l+1):N) {
  #       S <- S + (sum(X[k,] * X[l, ]) * sum(X[l, ] * X[r, ]) * sum(X[r, ] * X[k, ]))
  #     }
  #   }
  #   return(S/choose(N, 3))
  # }
  
  B3_2 <- function(X, N){
    S <- 0
    #browser()
      for (k in 1:(N-2)) {
      for (l in (k+1):(N-1)) {
        for (r in (l+1):N) {
          S <- S + (sum(X[k, ] * X[l, ]) * sum(X[l, ] * X[r, ]) * sum(X[r, ] * X[k, ]))
        }
      }
    }
    return(S/choose(N,3))
  }
  
  #browser()
  
### Hypothesenmatrizen
  if(is.matrix(hypothesis) & all(dim(hypothesis) == c(d,d))) TM <- hypothesis
  if(hypothesis[1] == "equal") TM <- diag(d)
  if(hypothesis[1] == "flat") TM <- diag(d) -  matrix(1/d, d, d)
  
  
### Teststatistik Q
  XT <- X %*% TM
  Xquer <- colMeans(XT)
  Qn = N * sum(Xquer * Xquer)
  
  
### Schätzer berechnen
  spurNormal <- B0(X = XT, N = N)
  spurQuadrat <- B2(X = XT, N = N)
  spurHoch3 <- B3_2(X = XT, N = N)
  
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
  
### Ausgabe
  if(is.matrix(hypothesis)) H <- "custom"
    else H <- paste0(hypothesis[1], " time profile")
  
  L <- list(data = X,
            f = f,
            statisitc = W,
            tau = 1/f,
            hypothesis = H,
            p.value = p.value,
            critical.value = critical.value,
            dim = c(d = d, N = N)
            )
  class(L) <- c("hdrm1","list")
  return(L)
  
}

### es scheint so, als wären die schätzer für f un tau vertauscht???
### -> im paper steht tau in [0, 1]

hdrm1(X, hypothesis = "flat")
hdrm1(Z, hypothesis = "equal")

hdrm1(t(Y), hypothesis = "flat")
hdrm1(t(Y), hypothesis = "equal")

Erg <- hdrm1(X, hypothesis = "flat")

# generische print-Funktion -----------------------------------------------

print.hdrm1 <- function(X,...){
  cat("       One Dimensional Repeated Measure
      \nW = ", X$statisitc, " f = ", X$f, " p.value = ", X$p.value,
      "\nNull-Hypothesis: ", X$hypothesis,
      "\nConvergence parameter \u03c4 = ", X$tau,
      "\n")
}


# bisschen useless, weil man eig gar nix erkennt...
plot.hdrm1 <- function(X, lty = 1,...){
  ylims <- c(min(apply(X$data, 2, min)), max(apply(X$data, 2, max)))
  d <- X$dim[1]
  N <- X$dim[2]
  plot(0, col = 0, ylim = ylims, xlim = c(1,d))
  for (i in 1:N) {
    lines(X$data[i,],lty = lty,...)
  }
}

# Funktion für pWert ------------------------------------------------------

### der p-Wert ist das kleinste alpha, sodass die Teststatistik noch ablehnt
### bzw das quantil der größten teststatistik, für die noch verworfen wird?


## für Chi^2-Verteilung
pWert_Chisq <- function(p, df, statistic) abs(qchisq(1-p, df) - statistic)
optimise(pWert_Chisq, interval = c(0,1), df = 2, statistic = 5.99)$minimum


# für K_f-Verteilung
pWert_Kf <- function(p, param, statistic) abs( ((qchisq(1-p, param) - param)/sqrt(2*param)) - statistic )
optimise(pWert_Kf, interval = c(0,1), param = 2, statistic = 200, tol = .Machine$double.eps)$minimum




# Simulationsstudie, um zu überprüfen, ob ungefähr das richtige Er --------

Counter <- 0

for (i in 1:1000) {
  N <- 20
  means <- rep(1,N)
  d <- 30
  #X <- t(rmvnorm(n = d, mean = means))
  Y <- t(replicate(N, rnorm(d,sd = 1)))
  
  if(hdrm1(Y, hypothesis = "flat")$p.value < 0.05) Counter = Counter + 1
}

Counter / 1000


# bei X ergibt es sinn für "flat", aber für "equal" wirklich auch gar nicht :(
#   -> bei "equal" wird iwie immer abgelehnt -> ist da irgendwas mit der Matrix falsch?

# Bei Y mit sd = 0.01 ist der Fehler 1. Art allerdings bei "equal" und "flat" ungefähr 0.05
#   -> der Test ist wohl sehr sensibel bei "equal"?
# Bei sd = 1 allerdings auch. Dann liegt es wohl an der Kovarianz bei rmvnorm...











