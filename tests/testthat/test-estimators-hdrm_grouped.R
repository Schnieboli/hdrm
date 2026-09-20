test_that("A2 agrees with independent covariance-trace references", {
  A2_centered_R <- function(X, Y) {
    X_centered <- sweep(X,
                        MARGIN = 1L,
                        STATS = rowMeans(X),
                        FUN = "-")
    Y_centered <- sweep(Y,
                        MARGIN = 1L,
                        STATS = rowMeans(Y),
                        FUN = "-")
    
    cross_inner_products <- t(X_centered) %*% Y_centered
    
    sum(cross_inner_products^2) /
      ((ncol(X) - 1L) * (ncol(Y) - 1L))
  }
  
  A2_covariance_R <- function(X, Y) {
    covariance_X <- stats::cov(t(X))
    covariance_Y <- stats::cov(t(Y))
    
    sum(covariance_X * covariance_Y)
  }
  
  X_fixed <- matrix(c(1, 2, 3, 0, 2, 4), nrow = 2L, ncol = 3L)
  
  Y_fixed <- matrix(c(0, 1, 2, 3, 4, -1, 1, 2), nrow = 2L, ncol = 4L)
  
  # Independently calculated value:
  # tr(cov(t(X_fixed)) %*% cov(t(Y_fixed))) = 71 / 4.
  expect_equal(hdrm:::A2_ir_cpp(X_fixed, Y_fixed), 17.75, tolerance = 1e-12)
  
  expect_equal(hdrm:::A2_ir_cpp(X_fixed, Y_fixed),
               A2_centered_R(X_fixed, Y_fixed),
               tolerance = 1e-12)
  
  expect_equal(
    hdrm:::A2_ir_cpp(X_fixed, Y_fixed),
    A2_covariance_R(X_fixed, Y_fixed),
    tolerance = 1e-12
  )
  
  set.seed(123)
  
  X_random <- matrix(rnorm(4L * 5L), nrow = 4L, ncol = 5L)
  
  Y_random <- matrix(rnorm(4L * 8L), nrow = 4L, ncol = 8L)
  
  expect_equal(
    hdrm:::A2_ir_cpp(X_random, Y_random),
    A2_centered_R(X_random, Y_random),
    tolerance = 1e-12
  )
  
  expect_equal(
    hdrm:::A2_ir_cpp(X_random, Y_random),
    A2_covariance_R(X_random, Y_random),
    tolerance = 1e-12
  )
  
  X_one_dimension <- matrix(c(-2, -0.5, 1, 3, 4.5), nrow = 1L)
  
  Y_one_dimension <- matrix(c(-3, 0, 2, 7), nrow = 1L)
  
  expect_equal(
    hdrm:::A2_ir_cpp(X_one_dimension, Y_one_dimension),
    stats::var(as.numeric(X_one_dimension)) *
      stats::var(as.numeric(Y_one_dimension)),
    tolerance = 1e-12
  )
  
  # The trace product is symmetric in the two samples.
  expect_equal(
    hdrm:::A2_ir_cpp(X_random, Y_random),
    hdrm:::A2_ir_cpp(Y_random, X_random),
    tolerance = 1e-12
  )
  
  # Centering inside the estimator makes it invariant to group-specific
  # translations of every dimension.
  X_translated <- X_random + c(10, -4, 2.5, 100)
  
  Y_translated <- Y_random + c(-7, 3, 0.25, -20)
  
  expect_equal(
    hdrm:::A2_ir_cpp(X_translated, Y_translated),
    hdrm:::A2_ir_cpp(X_random, Y_random),
    tolerance = 1e-12
  )
})

# test_that("C5star_cpp wrapper and internal transformation agree", {
#
#   X <- matrix(
#     c(
#       -2, 1,
#       -1, 3,
#       0, -2,
#       2, 0,
#       3, 2,
#       5, -1,
#       1, 4,
#       -3, 2,
#       2, 5,
#       4, -2,
#       6, 1,
#       0, 3
#     ),
#     nrow = 2L,
#     ncol = 12L
#   )
#
#   group <- as.numeric(
#     rep(
#       1:2,
#       each = 6L
#     )
#   )
#
#   TW <- diag(2L)
#   TS <- diag(2L)
#   B <- 250L
#
#   group_sizes <- unname(
#     as.integer(
#       table(group)
#     )
#   )
#   N <- ncol(X)
#   group_boundaries <- cumsum(
#     c(
#       1L,
#       group_sizes
#     )
#   )
#
#   Y <- matrix(
#     0,
#     nrow = nrow(TW) * nrow(TS),
#     ncol = N
#   )
#
#   for (i in seq_len(length(group_sizes))) {
#     columns_i <- group_boundaries[[i]]:(
#       group_boundaries[[i + 1L]] - 1L
#     )
#
#     Y[, columns_i] <- kronecker(
#       TW[, i],
#       TS %*% (
#         X[, columns_i, drop = FALSE] *
#           sqrt(N / group_sizes[[i]])
#       )
#     )
#   }
#
#   withr::local_seed(2718)
#
#   wrapper_result <- hdrm:::C5star_cpp(
#     X = X,
#     group = group,
#     TW = TW,
#     TS = TS,
#     B = B
#   )
#
#   set.seed(2718)
#
#   internal_result <- hdrm:::C5star_cpp_internal(
#     X = Y,
#     group = group,
#     B = length(group_sizes) * B,
#     n = group_sizes
#   )
#
#   expect_identical(
#     wrapper_result,
#     internal_result
#   )
# })


test_that("exact grouped C++ estimators agree with R references", {
  A1_R <- function(X) {
    n <- ncol(X)
    out <- 0
    
    for (i in seq_len(n - 1L)) {
      for (j in seq.int(i + 1L, n)) {
        out <- out + sum((X[, i] - X[, j])^2)
      }
    }
    
    out / (n * (n - 1))
  }
  
  A3_R <- function(X) {
    n <- ncol(X)
    out <- 0
    number_of_terms <- 0L
    
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        for (k in seq_len(n)) {
          for (l in seq_len(n)) {
            indices <- c(i, j, k, l)
            
            if (length(unique(indices)) == 4L) {
              inner_product <- sum((X[, i] - X[, j]) *
                                     (X[, k] - X[, l]))
              
              out <- out + inner_product^2 / 4
              number_of_terms <- number_of_terms + 1L
            }
          }
        }
      }
    }
    
    out / number_of_terms
  }
  
  X_fixed <- matrix(1:12, nrow = 3L, ncol = 4L)
  
  set.seed(123)
  X_random <- matrix(rnorm(3L * 6L), nrow = 3L, ncol = 6L)
  
  X_one_dimension <- matrix(c(-2, -0.5, 1, 3, 4.5), nrow = 1L)
  
  # Independently calculated values for the fixed matrix
  expect_equal(hdrm:::A1_i_cpp(X_fixed), 45, tolerance = 1e-12)
  
  expect_equal(hdrm:::A3_i_cpp(X_fixed), 1579.5, tolerance = 1e-12)
  
  matrices <- list(fixed = X_fixed,
                   random = X_random,
                   one_dimension = X_one_dimension)
  
  for (current_matrix in matrices) {
    expect_equal(hdrm:::A1_i_cpp(current_matrix),
                 A1_R(current_matrix),
                 tolerance = 1e-12)
    
    expect_equal(hdrm:::A3_i_cpp(current_matrix),
                 A3_R(current_matrix),
                 tolerance = 1e-10)
  }
})


test_that("grouped C++ subsampling estimators have the intended kernels", {
  X_pair <- matrix(c(1, 3, 5, 7), nrow = 2L, ncol = 2L)
  
  Y_pair <- matrix(c(2, -1, 8, 5), nrow = 2L, ncol = 2L)
  
  B_pair <- 25L
  
  # With exactly two columns, every draw selects the same pair up to order.
  expect_equal(hdrm:::A1star_i_cpp(X_pair, B_pair), sum((X_pair[, 1L] - X_pair[, 2L])^2) / 2, tolerance = 1e-12)
  
  # Independent reversals of the two pairs only change the sign of the
  # inner product, which disappears after squaring.
  expect_equal(hdrm:::A2star_ir_cpp(X_pair, Y_pair, B_pair),
               sum((X_pair[, 1L] - X_pair[, 2L]) *
                     (Y_pair[, 1L] - Y_pair[, 2L]))^2 / 4,
               tolerance = 1e-12)
  
  A3_kernel_values <- function(X) {
    n <- ncol(X)
    values <- numeric(0)
    
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        for (k in seq_len(n)) {
          for (l in seq_len(n)) {
            indices <- c(i, j, k, l)
            
            if (length(unique(indices)) == 4L) {
              inner_product <- sum((X[, i] - X[, j]) *
                                     (X[, k] - X[, l]))
              
              values <- c(values, inner_product^2 / 4)
            }
          }
        }
      }
    }
    
    values
  }
  
  X_subsampling <- matrix(c(-2, 1, 0, -1, 3, 2, 0, -2, 4, 2, 0, -1, 3, 2, 1, 5, -1, 3),
                          nrow = 3L,
                          ncol = 6L)
  
  kernel_values <- A3_kernel_values(X_subsampling)
  expected_value <- mean(kernel_values)
  
  B_subsampling <- 20000L
  standard_error <- sqrt(stats::var(kernel_values) / B_subsampling)
  
  withr::local_seed(3141)
  
  estimate_1 <- hdrm:::A3star_i_cpp(X_subsampling, B_subsampling)
  
  set.seed(3141)
  
  estimate_2 <- hdrm:::A3star_i_cpp(X_subsampling, B_subsampling)
  
  expect_identical(estimate_1, estimate_2)
  
  expect_true(is.finite(estimate_1))
  expect_gte(estimate_1, 0)
  
  # Ten Monte Carlo standard errors leave a wide margin while still detecting
  # missing factors such as the division by four.
  expect_lte(abs(estimate_1 - expected_value), 10 * standard_error + 1e-12)
})


test_that("equal-covariance exact C++ estimators agree with R references", {
  A1_eq_R_raw <- function(X) {
    n <- ncol(X)
    out <- 0
    
    for (l2 in seq_len(n - 1L)) {
      for (l1 in seq.int(l2 + 1L, n)) {
        difference <- X[, l1] - X[, l2]
        out <- out + sum(difference^2)
      }
    }
    
    out
  }
  
  A2_eq_R_raw <- function(X) {
    n <- ncol(X)
    out <- 0
    
    for (l2 in seq_len(n - 1L)) {
      for (l1 in seq.int(l2 + 1L, n)) {
        remaining <- setdiff(seq_len(n), c(l1, l2))
        
        for (position_2 in seq_len(length(remaining) - 1L)) {
          for (position_1 in seq.int(position_2 + 1L, length(remaining))) {
            k2 <- remaining[[position_2]]
            k1 <- remaining[[position_1]]
            
            inner_product <- sum((X[, l1] - X[, l2]) *
                                   (X[, k1] - X[, k2]))
            
            out <- out + inner_product^2
          }
        }
      }
    }
    
    out
  }
  
  all_permutations <- function(values) {
    if (length(values) == 1L) {
      return(matrix(values, nrow = 1L))
    }
    
    do.call(rbind, lapply(seq_along(values), function(index) {
      cbind(values[[index]], all_permutations(values[-index]))
    }))
  }
  
  C1_kernel_values <- function(X) {
    permutations <- all_permutations(seq_len(ncol(X)))
    
    apply(permutations, 1L, function(index) {
      difference_12 <- X[, index[[1L]]] - X[, index[[2L]]]
      difference_34 <- X[, index[[3L]]] - X[, index[[4L]]]
      difference_56 <- X[, index[[5L]]] - X[, index[[6L]]]
      
      inner_product_12_34 <- sum(difference_12 * difference_34)
      inner_product_34_56 <- sum(difference_34 * difference_56)
      inner_product_56_12 <- sum(difference_56 * difference_12)
      
      inner_product_12_34 *
        inner_product_34_56 *
        inner_product_56_12
    })
  }
  
  C1_eq_R_raw <- function(X) {
    sum(C1_kernel_values(X))
  }
  
  X_fixed <- matrix(c(0, 1, 1, 0, 2, 1, 0, 2, 1, 3, 3, 1),
                    nrow = 2L,
                    ncol = 6L)
  
  set.seed(123)
  
  X_random <- matrix(rnorm(3L * 6L), nrow = 3L, ncol = 6L)
  
  # Independently calculated values for the fixed matrix
  expect_equal(hdrm:::A1_i_eq_cpp(X_fixed), 73, tolerance = 1e-12)
  
  expect_equal(hdrm:::A2_i_eq_cpp(X_fixed), 844, tolerance = 1e-12)
  
  expect_equal(hdrm:::C1_i_eq_cpp(X_fixed), 12960, tolerance = 1e-12)
  
  matrices <- list(fixed = X_fixed, random = X_random)
  
  for (current_matrix in matrices) {
    expect_equal(hdrm:::A1_i_eq_cpp(current_matrix),
                 A1_eq_R_raw(current_matrix),
                 tolerance = 1e-12)
    
    expect_equal(hdrm:::A2_i_eq_cpp(current_matrix),
                 A2_eq_R_raw(current_matrix),
                 tolerance = 1e-12)
    
    expect_equal(hdrm:::C1_i_eq_cpp(current_matrix),
                 C1_eq_R_raw(current_matrix),
                 tolerance = 1e-10)
  }
})


test_that("equal-covariance C1 subsampling has the intended kernel", {
  all_permutations <- function(values) {
    if (length(values) == 1L) {
      return(matrix(values, nrow = 1L))
    }
    
    do.call(rbind, lapply(seq_along(values), function(index) {
      cbind(values[[index]], all_permutations(values[-index]))
    }))
  }
  
  C1_kernel_values <- function(X) {
    permutations <- all_permutations(seq_len(ncol(X)))
    
    apply(permutations, 1L, function(index) {
      difference_12 <- X[, index[[1L]]] - X[, index[[2L]]]
      difference_34 <- X[, index[[3L]]] - X[, index[[4L]]]
      difference_56 <- X[, index[[5L]]] - X[, index[[6L]]]
      
      inner_product_12_34 <- sum(difference_12 * difference_34)
      inner_product_34_56 <- sum(difference_34 * difference_56)
      inner_product_56_12 <- sum(difference_56 * difference_12)
      
      inner_product_12_34 *
        inner_product_34_56 *
        inner_product_56_12
    })
  }
  
  X_fixed <- matrix(c(0, 1, 1, 0, 2, 1, 0, 2, 1, 3, 3, 1),
                    nrow = 2L,
                    ncol = 6L)
  
  kernel_values <- C1_kernel_values(X_fixed)
  expected_kernel_mean <- mean(kernel_values)
  
  expect_equal(expected_kernel_mean, 18, tolerance = 1e-12)
  
  B_subsampling <- 20000L
  standard_error <- sqrt(stats::var(kernel_values) / B_subsampling)
  
  withr::local_seed(3141)
  
  raw_estimate_1 <- hdrm:::C1star_i_eq_cpp(X_fixed, B_subsampling)
  
  set.seed(3141)
  
  raw_estimate_2 <- hdrm:::C1star_i_eq_cpp(X_fixed, B_subsampling)
  
  expect_identical(raw_estimate_1, raw_estimate_2)
  
  estimate_per_draw <- raw_estimate_1 / B_subsampling
  
  expect_true(is.finite(estimate_per_draw))
  
  # Ten Monte Carlo standard errors provide a wide stability margin while
  # still detecting incorrect kernels or substantial scaling errors.
  expect_lte(abs(estimate_per_draw - expected_kernel_mean),
             10 * standard_error + 1e-12)
  
  # The R wrapper normalises the raw C++ sum by 8 * B.
  expect_lte(abs(raw_estimate_1 / (8 * B_subsampling) -
                   expected_kernel_mean / 8),
             10 * standard_error / 8 + 1e-12)
})


test_that("C5star_cpp_internal respects groups and matches its R kernel", {
  all_permutations <- function(values) {
    if (length(values) == 1L) {
      return(matrix(values, nrow = 1L))
    }
    
    do.call(rbind, lapply(seq_along(values), function(index) {
      cbind(values[[index]], all_permutations(values[-index]))
    }))
  }
  
  second_group_values <- c(-3, -1, 0, 2, 4, 7)
  
  permutations <- all_permutations(seq_along(second_group_values))
  
  kernel_values <- apply(permutations, 1L, function(index) {
    Z12 <- second_group_values[index[[1L]]] -
      second_group_values[index[[2L]]]
    Z34 <- second_group_values[index[[3L]]] -
      second_group_values[index[[4L]]]
    Z56 <- second_group_values[index[[5L]]] -
      second_group_values[index[[6L]]]
    
    (Z12 * Z34) *
      (Z34 * Z56) *
      (Z56 * Z12) / 8
  })
  
  reference_mean <- mean(kernel_values)
  
  # Exact mean over all 6! ordered permutations
  expect_equal(reference_mean, 63907 / 60, tolerance = 1e-12)
  
  # The first group is identically zero. A nonzero result can therefore only
  # arise when the C++ routine applies the correct offset to the second group.
  X_internal <- list(matrix(0, ncol = 6), matrix(second_group_values, ncol = 6))
  
  B_internal <- 20000L
  standard_error <- sqrt(stats::var(kernel_values) / B_internal)
  
  withr::local_seed(3141)
  
  estimate_1 <- hdrm:::C5star_cpp(X = X_internal, B = B_internal)
  
  set.seed(3141)
  
  estimate_2 <- hdrm:::C5star_cpp(X = X_internal, B = B_internal)
  
  expect_identical(estimate_1, estimate_2)
  
  expect_true(is.finite(estimate_1))
  
  expect_gt(estimate_1, 0)
  
  # Ten Monte Carlo standard errors give a wide stability margin while still
  # detecting an omitted group offset or a missing factor of one eighth.
  expect_lte(abs(estimate_1 - reference_mean), 10 * standard_error + 1e-12)
})
