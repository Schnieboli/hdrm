#' @export
hdrm_single.data.frame <- function(
    data,
    hypothesis = "flat",
    AM = TRUE
) {
  ## checks for AM
  AM <- as.logical(AM)
  if (length(AM) != 1L || is.na(AM)) {
    stop("'AM' must be a single logical value.", call. = FALSE)
  }
  
  if(is.null(data$value) || is.null(data$subject) || is.null(data$dimension)){
    stop("'data' must contain columns 'value', 'subject' and 'dimension'", 
         call. = FALSE)
  }
  
  if(nrow(data) < 1){
    stop("'data' must not be empty.", call. = FALSE)
  }
  if(any(is.na(data)) || any(!sapply(data, is.finite))){
    stop("'data' must only contain finite, non-missing values.", call. = FALSE)
  }
  if(!is.numeric(data$value) || any(!is.finite(data$value))){
    stop("data$value must be numeric and finite", call. = FALSE)
  }
  
  data <- data.frame(value = data$value, 
                     subject = as.factor(data$subject), 
                     dimension = as.factor(data$dimension))
  data <- droplevels(data)
  

  if(any(table(data$subject, data$dimension) != 1))
    stop("each combination of subject and dimension must occur exactly once.", 
         call. = FALSE)
  
  data <- data[order(data$subject, data$dimension), ]
  ## reshape data to widetable format
  df_wide <- stats::reshape(
    data,
    idvar = c("subject"),
    timevar = "dimension",
    direction = "wide"
  )
  
  X <- t(as.matrix(df_wide[, -1]))
  dimnames(X) <- NULL
    
  
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
