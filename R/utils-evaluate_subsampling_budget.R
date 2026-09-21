# Evaluate the subsampling budget without evaluating arbitrary R code.
#
# Character expressions may contain only numeric constants, N, parentheses,
# and the arithmetic operators +, -, *, /, and ^.
evaluate_subsample_budget <- function(B, N, a) {
  invalid_expression <- function() {
    stop(
      paste0(
        "'B' must be an arithmetic expression involving only 'N', ",
        "numeric constants, parentheses, and the operators +, -, *, /, and ^."
      ),
      call. = FALSE
    )
  }
  
  evaluate_node <- function(node) {
    if (is.numeric(node) &&
        length(node) == 1L &&
        is.null(dim(node))) {
      return(as.numeric(node))
    }
    
    if (is.name(node)) {
      if (identical(as.character(node), "N")) {
        return(as.numeric(N))
      }
      
      invalid_expression()
    }
    
    if (!is.call(node)) {
      invalid_expression()
    }
    
    call_head <- node[[1L]]
    
    if (!is.name(call_head)) {
      invalid_expression()
    }
    
    operator <- as.character(call_head)
    
    if (length(operator) != 1L) {
      invalid_expression()
    }
    
    number_of_arguments <- length(node) - 1L
    
    if (identical(operator, "(") &&
        number_of_arguments == 1L) {
      return(evaluate_node(node[[2L]]))
    }
    
    if (operator %in% c("+", "-") &&
        number_of_arguments == 1L) {
      value <- evaluate_node(node[[2L]])
      
      return(if (operator == "+")
        value
        else-value)
    }
    
    if (!operator %in% c("+", "-", "*", "/", "^") ||
        number_of_arguments != 2L) {
      invalid_expression()
    }
    
    left <- evaluate_node(node[[2L]])
    right <- evaluate_node(node[[3L]])
    
    suppressWarnings(switch(
      operator,
      "+" = left + right,
      "-" = left - right,
      "*" = left * right,
      "/" = left / right,
      "^" = left^right
    ))
  }
  
  if (!is.numeric(N) ||
      length(N) != 1L ||
      is.na(N) ||
      !is.finite(N) ||
      N < 1 ||
      N != floor(N)) {
    stop("Internal error: 'N' must be a finite positive integer.",
         call. = FALSE)
  }
  
  if (!(is.numeric(B) || is.character(B)) ||
      length(B) != 1L ||
      !is.null(dim(B)) ||
      is.na(B)) {
    stop("'B' must be a single numeric or character value.", call. = FALSE)
  }
  
  if (is.numeric(B)) {
    value <- as.numeric(B)
  } else {
    parsed <- tryCatch(
      parse(text = B, keep.source = FALSE),
      error = function(e)
        NULL
    )
    
    if (is.null(parsed) ||
        length(parsed) != 1L) {
      invalid_expression()
    }
    
    value <- evaluate_node(parsed[[1L]])
  }
  
  value <- ceiling(value)
  
  if (!is.numeric(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !is.finite(value) ||
      value < 1 ||
      a * value > .Machine$integer.max) {
    stop(
      paste0(
        "'B' times the number of groups must evaluate to a single finite positive number not exceeding ",
        .Machine$integer.max,
        "."
      ),
      call. = FALSE
    )
  }
  
  as.integer(value)
}