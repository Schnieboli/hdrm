A1_eq <- function(X){
  n <- nrow(X)
  S <- 0.0
  for(j in 1:(n-1)){
    for(k in j:n){
      S <- S + crossprod(X[j, ], X[k, ])
    }
  }
  return(as.numeric(S))
}
