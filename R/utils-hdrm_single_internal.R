hdrm_single_internal <- function(X, H){
  
  d <- nrow(X)
  N <- ncol(X)
  
  X_TM <- H$TMalt %*% X
  Qn <- N * sum(rowMeans(X_TM)^2)
  
  traceNormal <- B0_cpp(X_TM)
  traceSquare <- B2_cpp(X_TM)
  traceCubic <- B3_cpp(X_TM)
  
  if (
    !is.finite(Qn) || !is.finite(traceNormal) || !is.finite(traceSquare) ||
    traceSquare <= 0 || !is.finite(traceCubic) || traceCubic == 0) {
    stop("The trace estimators are numerically degenerate for these data.",
         call. = FALSE)
  }
  
  W <- (Qn - traceNormal) / sqrt(2 * traceSquare)
  
  if (!is.finite(W)) {
    stop("The test statistic is numerically undefined for these data.",
         call. = FALSE)
  }
  
  f <- max(1, traceSquare^3 / traceCubic^2)
  
  p.value <- stats::pchisq(W * sqrt(2 * f) + f, df = f, lower.tail = FALSE)
  
  list(
    f = f,
    statistic = W,
    p.value = max(p.value, .Machine$double.eps),
    H = H$TM
  )
}