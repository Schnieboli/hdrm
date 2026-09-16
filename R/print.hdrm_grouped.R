#' @method print hdrm_grouped
#' @export
print.hdrm_grouped <- function(x, digits = 4, ...) {
  # Format the p-value
  p <- round(x$p.value, digits)
  
  if (is.na(p)) {
    p_text <- "= NA"
  } else if (p <= 0) {
    p_text <- paste0("< ", 10^(-digits))
  } else {
    p_text <- paste0("= ", p)
  }
  
  hypothesis_text <- if (
    is.list(x$hypothesis)
  ) {
    "custom"
  } else {
    x$hypothesis
  }
  
  cat(
    "\n",
    "          Multi Group Repeated Measure\n",
    "Analysis of ", x$dim$N,
    " individuals in ", x$groups$a,
    " groups and ", x$dim$d,
    " dimensions:",
    "\nW = ", round(x$statistic, digits),
    "  f = ", round(x$f, digits),
    "  p.value ", p_text,
    "\nHypothesis type: ", hypothesis_text,
    "\nConvergence parameter \u03c4 = ", round(x$tau, digits),
    "\n",
    sep = ""
  )
  
  invisible(NULL)
}
