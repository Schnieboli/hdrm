#' Compute a compact matrix root
#'
#' Internal helper used for the compact representation of symmetric positive
#' semidefinite hypothesis matrices.
#'
#' @param H A finite numeric square matrix with positive rank.
#'
#' @returns A compact matrix root.
#'
#' @noRd
MSrootcompact <- function(H) {
  if (length(H) == 1L) {
    if (!is.numeric(H) ||
        is.na(H) ||
        !is.finite(H) ||
        H <= 0) {
      stop("'H' must have positive rank.", call. = FALSE)
    }
    return(matrix(sqrt(H), nrow = 1L, ncol = 1L))
  }
  
  if (!is.matrix(H) || !is.numeric(H) || anyNA(H) ||
      any(!is.finite(H)) || nrow(H) != ncol(H)) {
    stop("'H' must be a finite numeric square matrix.", call. = FALSE)
  }
  rank_H <- qr(H)$rank
  if (rank_H == 0L) stop("'H' must have positive rank.", call. = FALSE)

  decomposition <- svd(H)
  
  root <- diag(sqrt(decomposition$d[seq_len(rank_H)]), nrow = rank_H, ncol = rank_H) %*%
    t(decomposition$u[, seq_len(rank_H), drop = FALSE])
  
  root
}
