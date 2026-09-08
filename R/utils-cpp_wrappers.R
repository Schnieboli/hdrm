# All trace estimators expect subjects in columns and dimensions in rows.

#' @keywords internal
C5star_old <- function(X, group, TW, TS, B) {
  # if (length(group) != ncol(X)) {
  #   stop(
  #     "The length of 'group' must equal the number of columns of 'X'.",
  #     call. = FALSE
  #   )
  # }
  # 
  # if (is.unsorted(group)) {
  #   stop(
  #     "The columns of 'X' and the entries of 'group' must be ordered by group.",
  #     call. = FALSE
  #   )
  # }

  group_table <- table(group)
  a <- length(group_table)
  N <- ncol(X)
  n <- unname(as.integer(group_table))

  # if (any(n < 6L)) {
  #   stop(
  #     "All group sizes must be at least 6.",
  #     call. = FALSE
  #   )
  # }

  group_boundaries <- cumsum(c(1L, n))
  Y <- matrix(
    0,
    nrow = nrow(TW) * nrow(TS),
    ncol = N
  )

  for (i in seq_len(a)) {
    columns_i <- group_boundaries[i]:(group_boundaries[i + 1L] - 1L)

    Y[, columns_i] <- kronecker(
      TW[, i],
      TS %*% (
        X[, columns_i, drop = FALSE] *
          sqrt(N / n[i])
      )
    )
  }

  joint_B <- expand_subsample_budget(
    B = B,
    multiplier = a
  )

  C5star_cpp_internal(
    X = Y,
    group = group,
    B = joint_B,
    n = n
  )
}


#' @keywords internal
make_A1_eq <- function(X, group) {
  group_sizes <- unname(as.integer(table(group)))
  a <- length(group_sizes)
  total <- 0

  for (i in seq_len(a)) {
    total <- total + A1_i_eq_cpp(
      X[, group == i, drop = FALSE]
    )
  }

  denominator <- sum(
    group_sizes * (group_sizes - 1L)
  )

  total / denominator
}


#' @keywords internal
make_A2_eq <- function(X, group) {
  group_sizes <- unname(as.integer(table(group)))
  a <- length(group_sizes)
  total <- 0

  for (i in seq_len(a)) {
    total <- total + A2_i_eq_cpp(
      X[, group == i, drop = FALSE]
    )
  }

  denominator <- 24 * sum(
    choose(group_sizes, 4)
  )

  total / denominator
}


#' @keywords internal
make_C1_eq <- function(X, group) {
  group_sizes <- unname(as.integer(table(group)))
  a <- length(group_sizes)
  total <- 0

  for (i in seq_len(a)) {
    total <- total + C_i1_eq_cpp(
      X[, group == i, drop = FALSE]
    )
  }

  denominator <- 5760 * sum(
    choose(group_sizes, 6)
  )

  total / denominator
}


#' @keywords internal
make_C1_star_eq <- function(X, B, group) {
  group_sizes <- unname(as.integer(table(group)))
  a <- length(group_sizes)

  group_subsamples <- allocate_C1_subsamples(
    group_sizes = group_sizes,
    B = B
  )

  total <- 0

  for (i in seq_len(a)) {
    total <- total + C1star_i_eq_cpp(
      X[, group == i, drop = FALSE],
      group_subsamples[i]
    )
  }

  total / (8 * sum(group_subsamples))
}
