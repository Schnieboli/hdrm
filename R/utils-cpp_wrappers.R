# All trace estimators expect subjects in columns and dimensions in rows.


#' @keywords internal
A2 <- function(X, Y) {
  nX <- ncol(X)
  nY <- ncol(Y)

  PX <- diag(nX) - matrix(
    1 / nX,
    nrow = nX,
    ncol = nX
  )
  PY <- diag(nY) - matrix(
    1 / nY,
    nrow = nY,
    ncol = nY
  )

  MX <- tcrossprod(PX, X)
  MY <- tcrossprod(Y, PY)
  MXY <- MX %*% MY

  sum(MXY^2) / ((nX - 1) * (nY - 1))
}


#' @keywords internal
C5star_cpp <- function(X, group, TW, TS, B) {
  if (length(group) != ncol(X)) {
    stop(
      "The length of 'group' must equal the number of columns of 'X'.",
      call. = FALSE
    )
  }

  if (is.unsorted(group)) {
    stop(
      "The columns of 'X' and the entries of 'group' must be ordered by group.",
      call. = FALSE
    )
  }

  group_table <- table(group)
  a <- length(group_table)
  N <- ncol(X)
  n <- unname(as.integer(group_table))

  if (any(n < 6L)) {
    stop(
      "All group sizes must be at least 6.",
      call. = FALSE
    )
  }

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
A1_eq <- function(X) {
  n <- ncol(X)
  total <- 0

  for (l2 in seq_len(n - 1L)) {
    for (l1 in seq.int(l2 + 1L, n)) {
      total <- total + sum(
        (X[, l1] - X[, l2])^2
      )
    }
  }

  total
}


#' @keywords internal
make_A1_eq <- function(X, group) {
  group_sizes <- unname(as.integer(table(group)))
  a <- length(group_sizes)
  total <- 0

  for (i in seq_len(a)) {
    total <- total + A1_eq_cpp(
      X[, group == i, drop = FALSE]
    )
  }

  denominator <- sum(
    group_sizes * (group_sizes - 1L)
  )

  total / denominator
}


#' @keywords internal
A2_eq <- function(X) {
  n <- ncol(X)
  total <- 0

  for (l2 in seq_len(n - 1L)) {
    for (l1 in seq.int(l2 + 1L, n)) {
      for (k2 in seq_len(n - 1L)) {
        if (k2 != l1 && k2 != l2) {
          for (k1 in seq.int(k2 + 1L, n)) {
            if (k1 != l1 && k1 != l2) {
              total <- total +
                crossprod(
                  X[, l1] - X[, l2],
                  X[, k1] - X[, k2]
                )^2
            }
          }
        }
      }
    }
  }

  as.numeric(total)
}


#' @keywords internal
make_A2_eq <- function(X, group) {
  group_sizes <- unname(as.integer(table(group)))
  a <- length(group_sizes)
  total <- 0

  for (i in seq_len(a)) {
    total <- total + A2_eq_cpp(
      X[, group == i, drop = FALSE]
    )
  }

  denominator <- 24 * sum(
    choose(group_sizes, 4)
  )

  total / denominator
}


#' @keywords internal
C1_eq <- function(X) {
  n <- ncol(X)
  total <- 0

  for (l1 in seq_len(n)) {
    for (l2 in seq_len(n)) {
      if (l1 != l2) {
        for (l3 in seq_len(n)) {
          if (l2 != l3) {
            for (l4 in seq_len(n)) {
              if (l3 != l4) {
                for (l5 in seq_len(n)) {
                  if (l5 != l4) {
                    for (l6 in seq_len(n)) {
                      if (
                        length(
                          unique(c(l1, l2, l3, l4, l5, l6))
                        ) == 6L
                      ) {
                        total <- total +
                          crossprod(
                            X[, l1] - X[, l2],
                            X[, l3] - X[, l4]
                          ) *
                          crossprod(
                            X[, l3] - X[, l4],
                            X[, l5] - X[, l6]
                          ) *
                          crossprod(
                            X[, l5] - X[, l6],
                            X[, l1] - X[, l2]
                          )
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  as.numeric(total)
}


#' @keywords internal
make_C1_eq <- function(X, group) {
  group_sizes <- unname(as.integer(table(group)))
  a <- length(group_sizes)
  total <- 0

  for (i in seq_len(a)) {
    total <- total + C1_eq_cpp(
      X[, group == i, drop = FALSE]
    )
  }

  denominator <- 5760 * sum(
    choose(group_sizes, 6)
  )

  total / denominator
}


#' @keywords internal
C1_star_eq <- function(X, B) {
  n <- ncol(X)
  total <- 0

  for (i in seq_len(B)) {
    indices <- sample.int(
      n = n,
      size = 6L,
      replace = FALSE
    )

    total <- total +
      crossprod(
        X[, indices[1L]] - X[, indices[2L]],
        X[, indices[3L]] - X[, indices[4L]]
      ) *
      crossprod(
        X[, indices[3L]] - X[, indices[4L]],
        X[, indices[5L]] - X[, indices[6L]]
      ) *
      crossprod(
        X[, indices[5L]] - X[, indices[6L]],
        X[, indices[1L]] - X[, indices[2L]]
      )
  }

  as.numeric(total)
}


#' @title Allocate the equal-covariance third-trace budget across groups
#'
#' @description
#' Treats `B` as a base budget per group and allocates the resulting exact total
#' budget \eqn{aB} approximately proportionally to the numbers of available
#' six-subject subsets, \eqn{\binom{n_i}{6}}. Every group receives at least one
#' subsample. The remaining budget is distributed by the largest-remainder
#' method, with ties resolved by group order.
#'
#' @param group_sizes Integer group sizes, all at least six.
#' @param B A positive integer giving the base subsampling budget.
#'
#' @returns An integer vector of group-specific subsample counts whose sum is
#' exactly \eqn{aB}.
#'
#' @noRd
allocate_C1_subsamples <- function(group_sizes, B) {
  if (
    !is.numeric(group_sizes) ||
    length(group_sizes) < 1L ||
    anyNA(group_sizes) ||
    any(!is.finite(group_sizes)) ||
    any(group_sizes != floor(group_sizes)) ||
    any(group_sizes < 6L)
  ) {
    stop(
      "'group_sizes' must contain finite integers of at least 6.",
      call. = FALSE
    )
  }

  group_sizes <- as.integer(group_sizes)
  a <- length(group_sizes)

  if (
    !is.numeric(B) ||
    length(B) != 1L ||
    is.na(B) ||
    !is.finite(B) ||
    B != floor(B) ||
    B < 1
  ) {
    stop(
      "'B' must be a finite positive integer.",
      call. = FALSE
    )
  }

  total_budget <- expand_subsample_budget(
    B = B,
    multiplier = a
  )

  combination_counts <- choose(group_sizes, 6)

  if (
    any(!is.finite(combination_counts)) ||
    sum(combination_counts) <= 0
  ) {
    stop(
      "The six-subject combination weights could not be computed.",
      call. = FALSE
    )
  }

  remaining_budget <- total_budget - a

  if (remaining_budget == 0L) {
    return(rep.int(1L, a))
  }

  raw_allocation <- remaining_budget *
    combination_counts /
    sum(combination_counts)

  integer_allocation <- floor(raw_allocation)
  allocation <- 1L + as.integer(integer_allocation)

  remainder <- total_budget - sum(allocation)

  if (remainder > 0L) {
    fractional_parts <- raw_allocation - integer_allocation

    priority <- order(
      -fractional_parts,
      seq_along(fractional_parts)
    )

    selected_groups <- priority[seq_len(remainder)]
    allocation[selected_groups] <-
      allocation[selected_groups] + 1L
  }

  if (
    any(allocation < 1L) ||
    sum(allocation) != total_budget
  ) {
    stop(
      "Internal error while allocating the subsampling budget.",
      call. = FALSE
    )
  }

  allocation
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
    total <- total + C1_star_eq_cpp(
      X[, group == i, drop = FALSE],
      group_subsamples[i]
    )
  }

  total / (8 * sum(group_subsamples))
}


