#' @method print hdrm_single
#' @export
print.hdrm_single <- function(x, digits = 4, ...) {
  p <- round(x$p.value, digits)
  
  if (is.na(p)) {
    p_text <- "= NA"
  } else if (p <= 0) {
    p_text <- paste0("< ", 10^(-digits))
  } else {
    p_text <- paste0("= ", p)
  }
  
  cat(
    "\n",
    "           One Group Repeated Measure\n",
    "Analysis of ", x$dim[["N"]],
    " subjects in ", x$dim[["d"]],
    " dimensions:",
    "\nW = ", round(x$statistic, digits),
    "  f = ", round(x$f, digits),
    "  p.value ", p_text,
    "\nHypothesis type: ", x$H$label,
    "\nConvergence parameter \u03c4 = ", round(x$tau, digits),
    "\n",
    sep = ""
  )
  
  invisible(NULL)
}