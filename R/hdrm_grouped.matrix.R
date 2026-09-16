#' @export
hdrm_grouped.matrix <- function(data,
                         hypothesis = "whole",
                         group,
                         AM = TRUE,
                         cov.equal = FALSE,
                         subsampling = FALSE,
                         B = "1000*N",
                         seed = NULL) {
  ## checks for AM
  AM <- as.logical(AM)
  if (length(AM) != 1L || is.na(AM)) {
    stop("'AM' must be a single logical value.", call. = FALSE)
  }
  ## checks for cov.equal
  cov.equal <- as.logical(cov.equal)
  if (length(cov.equal) != 1L || is.na(cov.equal)) {
    stop("'cov.equal' must be a single logical value.", 
         call. = FALSE)
  }
  ## checks for subsampling
  subsampling <- as.logical(subsampling)
  if (length(subsampling) != 1L || is.na(subsampling)) {
    stop("'subsampling' must be a single logical value.",
         call. = FALSE)
  }
  
  ## checks for seed
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
        seed != floor(seed) || abs(seed) > .Machine$integer.max) {
      stop("'seed' must be NULL or a single finite integer-valued number.",
           call. = FALSE)
    }
    seed <- as.integer(seed)
  }
  
  ## initialize output object
  # Store the processed data in the public N x d orientation while the
  # internal calculations continue to use subjects in columns.
  out <- list(data = data)
  
  ## do all matrix related checks
  data <- t(data)
  if(any(dim(data) == 0)){
    stop("'data' must not be empty.")
  }
  if(any(is.na(data))){
    stop("'data' must not contain any missing values.", call. = FALSE)
  }
  # Check that 'group' is a one-dimensional atomic vector
  if (!is.atomic(group) || length(group) != ncol(data) || !is.null(dim(group))) {
    stop("'group' must be a one-dimensional vector or factor of length nrow(data).",
         call. = FALSE)
  }
  d <- nrow(data)
  N <- ncol(data)
  data_list <- lapply(split(seq_len(N), group), function(cols) data[, cols, drop = FALSE])
  
  
  ## check mathematical requirements
  a <- length(data_list)
  n <- sapply(data_list, ncol)
  if(a < 2) stop("there must be at least two groups", call. = FALSE)
  if(d < 2) stop("there must be at least two observations per subject", call. = FALSE)
  if(any(n < 6)) stop("there must be at least six subjects per group", call. = FALSE)
  
  # The grouped third-trace estimators use a * B draws. Validate the
  # effective budget before any stochastic estimator is evaluated.
  ### expand_subsample_budget(B = reps, multiplier = a) 
  # TODO kann das weg oder wurde hier vergessen, etwas zuzuweisen?
  
  out <- c(
    out,
    hdrm_grouped_internal(
      X_list = data_list,
      H = get_hypothesis_mult(hypothesis, AM, a, d),
      cov.equal = cov.equal,
      subsampling = subsampling,
      B = evaluate_subsample_budget(B = B, N = N),
      seed = seed
    )
  )
  
  # Add further output to the result
  out$groups = list(a = a, table = table(group) / d)  # Grouping information
  # Description of the hypothesis
  out$hypothesis = ifelse(is.character(hypothesis), hypothesis[1], "custom")
  class(out) <- "hdrm_grouped"
  return(out)
  
}
