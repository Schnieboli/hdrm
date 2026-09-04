#' @title Build the hypothesis matrices for a multifactorial setting
#'
#' @description
#' Extracts or constructs the whole-plot and subplot hypothesis matrices used
#' by the grouped procedures.
#'
#' @param hypothesis Either one of `"whole"`, `"sub"`, `"interaction"`,
#' `"identical"`, or `"flat"`, or a named list containing `TW` and `TS`.
#' @param a A positive integer giving the number of groups.
#' @param d A positive integer giving the repeated-measurement dimension.
#'
#' @returns A named list with components `TW` and `TS`.
#'
#' @noRd
get_hypothesis_mult <- function(hypothesis, a, d) {
  if (!is.list(hypothesis) && !is.character(hypothesis)) {
    stop(
      "'hypothesis' must be a character value or a list.",
      call. = FALSE
    )
  }

  if (is.list(hypothesis)) {
    if (is.null(hypothesis$TW)) {
      stop(
        "No list entry named 'TW' was found in 'hypothesis'.",
        call. = FALSE
      )
    }

    if (is.null(hypothesis$TS)) {
      stop(
        "No list entry named 'TS' was found in 'hypothesis'.",
        call. = FALSE
      )
    }

    TW <- hypothesis$TW
    TS <- hypothesis$TS

    if (!is.matrix(TW) || !is.numeric(TW)) {
      stop(
        "'TW' must be a numeric matrix.",
        call. = FALSE
      )
    }

    if (!is.matrix(TS) || !is.numeric(TS)) {
      stop(
        "'TS' must be a numeric matrix.",
        call. = FALSE
      )
    }

    if (anyNA(TW) || any(!is.finite(TW))) {
      stop(
        "'TW' must contain only finite, non-missing values.",
        call. = FALSE
      )
    }

    if (anyNA(TS) || any(!is.finite(TS))) {
      stop(
        "'TS' must contain only finite, non-missing values.",
        call. = FALSE
      )
    }

    if (
      length(dim(TW)) != 2L ||
      any(dim(TW) != c(a, a))
    ) {
      stop(
        paste0(
          "'TW' must be a ",
          a,
          " by ",
          a,
          " matrix."
        ),
        call. = FALSE
      )
    }

    if (
      length(dim(TS)) != 2L ||
      any(dim(TS) != c(d, d))
    ) {
      stop(
        paste0(
          "'TS' must be a ",
          d,
          " by ",
          d,
          " matrix."
        ),
        call. = FALSE
      )
    }

    tol <- sqrt(.Machine$double.eps)

    TW_symmetry_error <- max(abs(TW - t(TW)))
    TW_idempotence_error <- max(abs(TW %*% TW - TW))
    TS_symmetry_error <- max(abs(TS - t(TS)))
    TS_idempotence_error <- max(abs(TS %*% TS - TS))

    if (TW_symmetry_error > tol) {
      stop(
        paste0(
          "'TW' must be symmetric. Maximum deviation: ",
          signif(TW_symmetry_error, 4),
          "."
        ),
        call. = FALSE
      )
    }

    if (TW_idempotence_error > tol) {
      stop(
        paste0(
          "'TW' must be idempotent. Maximum deviation: ",
          signif(TW_idempotence_error, 4),
          "."
        ),
        call. = FALSE
      )
    }

    if (TS_symmetry_error > tol) {
      stop(
        paste0(
          "'TS' must be symmetric. Maximum deviation: ",
          signif(TS_symmetry_error, 4),
          "."
        ),
        call. = FALSE
      )
    }

    if (TS_idempotence_error > tol) {
      stop(
        paste0(
          "'TS' must be idempotent. Maximum deviation: ",
          signif(TS_idempotence_error, 4),
          "."
        ),
        call. = FALSE
      )
    }

    if (qr(TW, tol = tol)$rank == 0L) {
      stop(
        "'TW' must have positive rank.",
        call. = FALSE
      )
    }

    if (qr(TS, tol = tol)$rank == 0L) {
      stop(
        "'TS' must have positive rank.",
        call. = FALSE
      )
    }
  } else {
    if (
      length(hypothesis) != 1L ||
      is.na(hypothesis)
    ) {
      stop(
        "'hypothesis' must be a single non-missing character value.",
        call. = FALSE
      )
    }

    if (
      !hypothesis %in%
      c("whole", "sub", "interaction", "identical", "flat")
    ) {
      stop(
        paste0(
          "'hypothesis' must be one of 'whole', 'sub', 'interaction', ",
          "'identical', or 'flat', or a list containing 'TW' and 'TS'."
        ),
        call. = FALSE
      )
    }

    P_a <- diag(a) - matrix(1 / a, nrow = a, ncol = a)
    P_d <- diag(d) - matrix(1 / d, nrow = d, ncol = d)
    J_a <- matrix(1 / a, nrow = a, ncol = a)
    J_d <- matrix(1 / d, nrow = d, ncol = d)

    if (hypothesis == "whole") {
      TW <- P_a
      TS <- J_d
    } else if (hypothesis == "sub") {
      TW <- J_a
      TS <- P_d
    } else if (hypothesis == "interaction") {
      TW <- P_a
      TS <- P_d
    } else if (hypothesis == "identical") {
      TW <- P_a
      TS <- diag(d)
    } else {
      TW <- diag(a)
      TS <- P_d
    }
  }

  list(TW = TW, TS = TS)
}


#' @title Check the prepared inputs for the grouped test
#'
#' @description
#' Checks the assumptions required by the internal grouped procedures after the
#' public wrapper has prepared the data.
#'
#' @param X A numeric data matrix with subjects in columns.
#' @param group An integer vector containing consecutive group codes starting
#' at one.
#' @param hypothesis A predefined hypothesis name or a list containing `TW`
#' and `TS`.
#' @param reps A positive integer giving the requested number of subsamples.
#' @param subsampling A single non-missing logical value.
#'
#' @returns `NULL` invisibly. An error is raised when a condition is violated.
#'
#' @noRd
check_criteria_grouped <- function(
    X,
    group,
    hypothesis,
    reps,
    subsampling
) {
  if (!is.matrix(X) || !is.numeric(X)) {
    stop(
      "'X' must be a numeric matrix.",
      call. = FALSE
    )
  }

  if (anyNA(X) || any(!is.finite(X))) {
    stop(
      "'X' must contain only finite, non-missing values.",
      call. = FALSE
    )
  }

  d <- nrow(X)
  N <- ncol(X)

  if (
    !is.null(dim(group)) ||
    length(group) != N ||
    anyNA(group)
  ) {
    stop(
      "'group' must contain one non-missing group label per column of 'X'.",
      call. = FALSE
    )
  }

  if (is.factor(group)) {
    group_codes <- as.integer(group)
  } else {
    if (
      !is.numeric(group) ||
      any(!is.finite(group)) ||
      any(group != floor(group))
    ) {
      stop(
        "'group' must be a factor or an integer-valued vector.",
        call. = FALSE
      )
    }

    group_codes <- as.integer(group)
  }

  observed_groups <- sort(unique(group_codes))
  a <- length(observed_groups)

  if (!identical(observed_groups, seq_len(a))) {
    stop(
      "Internal group codes must be consecutive integers starting at one.",
      call. = FALSE
    )
  }

  n <- tabulate(group_codes, nbins = a)

  if (a < 2L) {
    stop(
      "At least two groups are required.",
      call. = FALSE
    )
  }

  if (d < 2L) {
    stop(
      "At least two repeated-measurement dimensions are required.",
      call. = FALSE
    )
  }

  if (any(n < 6L)) {
    stop(
      "All group sizes must be at least 6.",
      call. = FALSE
    )
  }

  if (!is.character(hypothesis) && !is.list(hypothesis)) {
    stop(
      "'hypothesis' must be a character value or a list.",
      call. = FALSE
    )
  }

  if (
    !is.logical(subsampling) ||
    length(subsampling) != 1L ||
    is.na(subsampling)
  ) {
    stop(
      "'subsampling' must be a single non-missing logical value.",
      call. = FALSE
    )
  }

  if (
    !is.numeric(reps) ||
    length(reps) != 1L ||
    is.na(reps) ||
    !is.finite(reps) ||
    reps < 1 ||
    reps != floor(reps) ||
    reps > .Machine$integer.max
  ) {
    stop(
      "The number of subsamples must be a finite positive integer.",
      call. = FALSE
    )
  }

  invisible(NULL)
}


################################################################################
# Whole-plot design factor
################################################################################

# Internal helper for the known factor eta_{N,a} in the equal-covariance
# procedure. The calculation is kept separate so that it can be tested
# independently from the stochastic trace estimators.
compute_eta_Na <- function(TW, group_sizes) {
  if (
    !is.matrix(TW) ||
    !is.numeric(TW) ||
    nrow(TW) != ncol(TW) ||
    anyNA(TW) ||
    any(!is.finite(TW))
  ) {
    stop(
      "'TW' must be a finite numeric square matrix.",
      call. = FALSE
    )
  }

  a <- nrow(TW)

  if (
    !is.numeric(group_sizes) ||
    !is.null(dim(group_sizes)) ||
    length(group_sizes) != a ||
    anyNA(group_sizes) ||
    any(!is.finite(group_sizes)) ||
    any(group_sizes <= 0) ||
    any(group_sizes != floor(group_sizes))
  ) {
    stop(
      paste0(
        "'group_sizes' must contain one positive integer for each ",
        "row of 'TW'."
      ),
      call. = FALSE
    )
  }

  symmetry_error <- max(abs(TW - t(TW)))

  if (symmetry_error > sqrt(.Machine$double.eps)) {
    stop(
      "'TW' must be symmetric.",
      call. = FALSE
    )
  }

  N <- sum(group_sizes)
  D_N <- diag(
    N / group_sizes,
    nrow = a,
    ncol = a
  )
  DNTW <- D_N %*% TW

  trace2_whole <- sum(
    diag(DNTW %*% DNTW)
  )
  trace3_whole <- sum(
    diag(DNTW %*% DNTW %*% DNTW)
  )

  if (
    !is.finite(trace2_whole) ||
    !is.finite(trace3_whole) ||
    trace2_whole <= 0 ||
    trace3_whole == 0
  ) {
    stop(
      "The whole-plot trace factor is degenerate.",
      call. = FALSE
    )
  }

  eta_Na <- trace2_whole^3 / trace3_whole^2

  if (
    !is.finite(eta_Na) ||
    eta_Na <= 0
  ) {
    stop(
      "The whole-plot trace factor is not finite and positive.",
      call. = FALSE
    )
  }

  as.numeric(eta_Na)
}


################################################################################
# Numerical safeguards for grouped procedures
################################################################################

#' @keywords internal
compute_grouped_statistic <- function(QN, EW, variance) {
  quantities <- c(
    QN = QN,
    EW = EW,
    variance = variance
  )

  if (
    any(lengths(list(QN, EW, variance)) != 1L) ||
    !is.numeric(quantities) ||
    anyNA(quantities) ||
    any(!is.finite(quantities))
  ) {
    stop(
      paste0(
        "The quadratic form, its estimated expectation, and its estimated ",
        "variance must be finite numeric scalars."
      ),
      call. = FALSE
    )
  }

  if (variance <= 0) {
    stop(
      "The estimated variance of the test statistic must be positive.",
      call. = FALSE
    )
  }

  statistic <- as.numeric(
    (QN - EW) / sqrt(variance)
  )

  if (
    length(statistic) != 1L ||
    is.na(statistic) ||
    !is.finite(statistic)
  ) {
    stop(
      "The standardized test statistic is not finite.",
      call. = FALSE
    )
  }

  statistic
}


#' @keywords internal
compute_grouped_df <- function(
    second_order,
    third_order,
    design_factor = 1
) {
  quantities <- c(
    second_order = second_order,
    third_order = third_order,
    design_factor = design_factor
  )

  if (
    any(lengths(list(second_order, third_order, design_factor)) != 1L) ||
    !is.numeric(quantities) ||
    anyNA(quantities) ||
    any(!is.finite(quantities))
  ) {
    stop(
      paste0(
        "The second-order trace estimate, third-order trace estimate, and ",
        "design factor must be finite numeric scalars."
      ),
      call. = FALSE
    )
  }

  if (second_order <= 0) {
    stop(
      "The second-order trace estimate must be positive.",
      call. = FALSE
    )
  }

  if (third_order == 0) {
    stop(
      paste0(
        "The third-order trace estimate is zero. Increase 'B' or check ",
        "whether the data are degenerate."
      ),
      call. = FALSE
    )
  }

  if (design_factor <= 0) {
    stop(
      "The design factor must be positive.",
      call. = FALSE
    )
  }

  degrees_of_freedom <- as.numeric(
    (second_order^3 / third_order^2) * design_factor
  )

  if (
    length(degrees_of_freedom) != 1L ||
    is.na(degrees_of_freedom) ||
    !is.finite(degrees_of_freedom) ||
    degrees_of_freedom <= 0
  ) {
    stop(
      paste0(
        "The estimated degrees-of-freedom parameter is not finite and ",
        "positive."
      ),
      call. = FALSE
    )
  }

  max(
    1,
    degrees_of_freedom
  )
}


#' @keywords internal
compute_grouped_p_value <- function(statistic, degrees_of_freedom) {
  quantities <- c(
    statistic = statistic,
    degrees_of_freedom = degrees_of_freedom
  )

  if (
    any(lengths(list(statistic, degrees_of_freedom)) != 1L) ||
    !is.numeric(quantities) ||
    anyNA(quantities) ||
    any(!is.finite(quantities)) ||
    degrees_of_freedom <= 0
  ) {
    stop(
      paste0(
        "'statistic' must be finite and 'degrees_of_freedom' must be ",
        "finite and positive."
      ),
      call. = FALSE
    )
  }

  p_value <- max(
    stats::pchisq(
      statistic * sqrt(2 * degrees_of_freedom) +
        degrees_of_freedom,
      df = degrees_of_freedom,
      lower.tail = FALSE
    ),
    .Machine$double.eps
  )

  if (
    length(p_value) != 1L ||
    is.na(p_value) ||
    !is.finite(p_value) ||
    p_value < 0 ||
    p_value > 1
  ) {
    stop(
      "The p-value calculation did not return a finite value in [0, 1].",
      call. = FALSE
    )
  }

  as.numeric(p_value)
}


################################################################################
# Subsampling-budget helpers
################################################################################

#' Expand a base subsampling budget by an integer multiplier
#'
#' Internal helper used when a grouped third-trace estimator deliberately uses
#' more than `B` draws. It also protects the C++ interface from integer
#' overflow.
#'
#' @param B A positive integer base budget.
#' @param multiplier A positive integer multiplier.
#'
#' @returns The integer product `B * multiplier`.
#'
#' @noRd
expand_subsample_budget <- function(B, multiplier) {
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

  if (
    !is.numeric(multiplier) ||
    length(multiplier) != 1L ||
    is.na(multiplier) ||
    !is.finite(multiplier) ||
    multiplier != floor(multiplier) ||
    multiplier < 1
  ) {
    stop(
      "'multiplier' must be a finite positive integer.",
      call. = FALSE
    )
  }

  effective_budget <- as.double(B) * as.double(multiplier)

  if (
    !is.finite(effective_budget) ||
    effective_budget > .Machine$integer.max
  ) {
    stop(
      paste0(
        "The effective subsampling budget must not exceed ",
        .Machine$integer.max,
        "."
      ),
      call. = FALSE
    )
  }

  as.integer(effective_budget)
}


################################################################################
# Trace estimators
################################################################################

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
    if (
      !is.numeric(H) ||
      is.na(H) ||
      !is.finite(H) ||
      H <= 0
    ) {
      stop(
        "'H' must have positive rank.",
        call. = FALSE
      )
    }

    return(
      matrix(
        sqrt(H),
        nrow = 1L,
        ncol = 1L
      )
    )
  }

  if (
    !is.matrix(H) ||
    !is.numeric(H) ||
    anyNA(H) ||
    any(!is.finite(H)) ||
    nrow(H) != ncol(H)
  ) {
    stop(
      "'H' must be a finite numeric square matrix.",
      call. = FALSE
    )
  }

  rank_H <- qr(H)$rank

  if (rank_H == 0L) {
    stop(
      "'H' must have positive rank.",
      call. = FALSE
    )
  }

  decomposition <- svd(H)

  root <- diag(
    sqrt(decomposition$d[seq_len(rank_H)]),
    nrow = rank_H,
    ncol = rank_H
  ) %*%
    t(
      decomposition$u[
        ,
        seq_len(rank_H),
        drop = FALSE
      ]
    )

  root
}
