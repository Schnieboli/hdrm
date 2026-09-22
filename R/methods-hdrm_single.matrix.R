#' @export
hdrm_single.matrix <- function(
    data,
    hypothesis = "flat",
    AM = TRUE
) {
  ## checks for AM
  AM <- as.logical(AM)
  if (length(AM) != 1L || is.na(AM)) {
    stop("'AM' must be a single logical value.", call. = FALSE)
  }
  
  X <- t(data)
  if(any(dim(X) == 0) || any(is.na(X)) || any(!is.numeric(X)) || any(!is.finite(X))){
    stop("'data' must be a numeric matrix, containing only finite, non-missing values.")
  }
  d <- nrow(X)
  N <- ncol(X)
  
  ## mathematical checks
  if (d < 2L) stop("At least two repeated-measurement dimensions are required.", call. = FALSE)
  if (N < 3L) stop("At least three subjects are required.", call. = FALSE)
  
  test_result <- hdrm_single_internal(
    X = X, 
    H = get_hypothesis_single(hypothesis, AM, d)
  )
  
  out <- list(
    data = t(X),
    dim = c(N = N, d = d),
    H = list(
      "T" = test_result$H,
      label = ifelse(is.character(hypothesis), hypothesis, "custom")
    ),
    AM = AM,
    statistic = test_result$statistic,
    p.value = test_result$p.value,
    f = test_result$f,
    tau = 1/test_result$f
  )
  class(out) <- "hdrm_single"
  out
}
