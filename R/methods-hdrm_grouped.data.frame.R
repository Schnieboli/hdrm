#' @method hdrm_grouped data.frame
#' @export
hdrm_grouped.data.frame <- function(data,
                         hypothesis = "whole",
                         group,
                         cov.equal = FALSE,
                         subsampling = FALSE,
                         B = "1000*N",
                         AM = TRUE,
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
  
  
  if(is.null(data$value) || is.null(data$subject) || is.null(data$dimension)){
    stop("'data' must contain columns 'value', 'subject' and 'dimension'", 
         call. = FALSE)
  }
  
  if(nrow(data) < 1){
    stop("'data' must not be empty.", call. = FALSE)
  }
  if(any(is.na(data)) || any(!sapply(data, is.finite))){
    stop("'data' must only contain finite, non-missing values", call. = FALSE)
  }
  if(!is.numeric(data$value) || any(!is.finite(data$value))){
    stop("data$value must be numeric and finite", call. = FALSE)
  }
  
  data <- data.frame(value = data$value, 
                     subject = as.factor(data$subject), 
                     dimension = as.factor(data$dimension))
  data <- droplevels(data)
  
  if(!is.atomic(group) || !is.null(dim(group)) || length(group) != nrow(data)){
    stop("'group' must be a one-dimensional vector or factor of length nrow(data).", call. = FALSE)
  }
  
  if(any(is.na(group)) || any(!is.finite(group))){
    stop("'group' must only contain finite, non-missing values.", call. = FALSE)
  }
  
  if(any(table(data$subject, data$dimension) != 1))
    stop("each combination of subject and dimension must occur exactly once.", 
         call. = FALSE)
  
  data$group <- as.factor(group)
  data <- data[order(data$subject, data$dimension, data$group), ]
  ## reshape data to widetable format
  df_wide <- stats::reshape(
    data,
    idvar = c("subject", "group"),
    timevar = "dimension",
    direction = "wide"
  )
  ## split df_wide into groups, transform to matrix
  data_list <- lapply(split(df_wide[, -c(1, 2)], df_wide$group), 
                      function(x) unname(as.matrix(t(x))))
  N <- sum(sapply(data_list, ncol))
  d <- nrow(data_list[[1]])
  
  ## check mathematical requirements
  a <- length(data_list)
  n <- sapply(data_list, ncol)
  if(a < 2) stop("there must be at least two groups", call. = FALSE)
  if(d < 2) stop("there must be at least two observations per subject", call. = FALSE)
  if(any(n < 6)) stop("there must be at least six subjects per group", call. = FALSE)
  
  
  test_result <- hdrm_grouped_internal(
    X_list = data_list,
    H = get_hypothesis_grouped(hypothesis, AM, a, d),
    cov.equal = cov.equal,
    subsampling = subsampling,
    B = evaluate_subsample_budget(B = B, N = N, a = a),
    seed = seed
  )
  
  out <- list(
    data = t(do.call(cbind, data_list)),
    dim = c(N = N, a = a, d = d),
    group = rep(1:a, sapply(data_list, ncol)),
    H = list(
      TW = test_result$H$TW,
      TS = test_result$H$TS,
      label = ifelse(is.character(hypothesis), hypothesis, "custom")
    ),
    statistic = test_result$statistic,
    p.value = test_result$p.value,
    f = test_result$f,
    tau = 1/test_result$f,
    AM = AM,
    cov.equal = cov.equal,
    subsampling = subsampling,
    B = test_result$B,
    seed = seed
  )
  levels(out$group) <- names(data_list)
  class(out) <- "hdrm_grouped"
  out
  
}
