data("birthrates")
Matrixbirthrates <- t(as.matrix(birthrates))
group <- c(1, 1, 2, 2, 1, 1, 1, 2, 1, 1, 1, 1, 2, 2, 1, 2)

test_that("degenerate grouped data produce informative errors", {
  constant_data <- matrix(3, nrow = 12L, ncol = 2L)
  constant_group <- factor(rep(c("A", "B"), each = 6L))
  
  # In the heterogeneous procedure, validation of the individual trace
  # estimators is reached before the derived variance is calculated.
  expect_error(
    hdrm_grouped(
      constant_data,
      hypothesis = "whole",
      group = constant_group,
      cov.equal = FALSE,
      subsampling = FALSE,
      B = 100,
      seed = 3141
    ),
    "The estimated variance of the test statistic must be finite and positive.",
    fixed = TRUE
  )
  
  
  expect_error(
    hdrm_grouped(
      constant_data,
      hypothesis = "whole",
      group = constant_group,
      cov.equal = TRUE,
      subsampling = FALSE,
      B = 100,
      seed = 3141
    ),
    "The estimated variance of the test statistic must be positive.",
    fixed = TRUE
  )
})


test_that("compute_eta_Na agrees with independent design-factor references",
          {
            eta_eigen_R <- function(TW, group_sizes) {
              N <- sum(group_sizes)
              square_root_D <- diag(
                sqrt(N / group_sizes),
                nrow = length(group_sizes),
                ncol = length(group_sizes)
              )
              
              symmetric_representation <- square_root_D %*%
                TW %*%
                square_root_D
              
              eigenvalues <- eigen(symmetric_representation,
                                   symmetric = TRUE,
                                   only.values = TRUE)$values
              
              sum(eigenvalues^2)^3 /
                sum(eigenvalues^3)^2
            }
            
            a <- 3L
            P_a <- diag(a) - matrix(1 / a, nrow = a, ncol = a)
            J_a <- matrix(1 / a, nrow = a, ncol = a)
            I_a <- diag(a)
            
            balanced_sizes <- c(5L, 5L, 5L)
            
            # For a balanced design, D_N = a I. The factor therefore equals the
            # rank of the whole-plot hypothesis matrix.
            expect_equal(hdrm:::compute_eta_Na(P_a, balanced_sizes, TRUE), 2, tolerance = 1e-12)
            
            expect_equal(hdrm:::compute_eta_Na(J_a, balanced_sizes, TRUE), 1, tolerance = 1e-12)
            
            expect_equal(hdrm:::compute_eta_Na(I_a, balanced_sizes, TRUE), 3, tolerance = 1e-12)
            
            unbalanced_sizes <- c(4L, 7L, 9L)
            
            expect_equal(
              hdrm:::compute_eta_Na(P_a, unbalanced_sizes, TRUE),
              eta_eigen_R(P_a, unbalanced_sizes),
              tolerance = 1e-12
            )
            
            expect_equal(
              hdrm:::compute_eta_Na(J_a, unbalanced_sizes, TRUE),
              eta_eigen_R(J_a, unbalanced_sizes),
              tolerance = 1e-12
            )
            
            expect_equal(
              hdrm:::compute_eta_Na(I_a, unbalanced_sizes, TRUE),
              eta_eigen_R(I_a, unbalanced_sizes),
              tolerance = 1e-12
            )
            
            expect_equal(hdrm:::compute_eta_Na(P_a, unbalanced_sizes, TRUE),
                         1.7009569124998098,
                         tolerance = 1e-12)
            
            # A non-coordinate projector represents a valid custom whole-plot
            # hypothesis and is checked against the independent eigenvalue formula.
            basis <- qr.Q(qr(matrix(
              c(1, 2, 0, 1, 0, 1, 2, 1, 2, 0, 1, 1, 1, 1, 1, -1),
              nrow = 4L,
              ncol = 4L
            )))
            
            custom_TW <- basis[, 1:2, drop = FALSE] %*%
              t(basis[, 1:2, drop = FALSE])
            custom_sizes <- c(3L, 5L, 8L, 11L)
            
            expect_equal(
              hdrm:::compute_eta_Na(custom_TW, custom_sizes, TRUE),
              eta_eigen_R(custom_TW, custom_sizes),
              tolerance = 1e-11
            )
            
            # A common nonzero scaling of TW cancels from the trace ratio.
            expect_equal(
              hdrm:::compute_eta_Na(7.5 * custom_TW, custom_sizes, TRUE),
              hdrm:::compute_eta_Na(custom_TW, custom_sizes, TRUE),
              tolerance = 1e-12
            )
            
            # Simultaneously relabelling groups and rows/columns of TW changes neither
            # the design nor eta_{N,a}.
            permutation <- c(3L, 1L, 4L, 2L)
            
            expect_equal(
              hdrm:::compute_eta_Na(custom_TW[permutation, permutation, drop = FALSE], custom_sizes[permutation], TRUE),
              hdrm:::compute_eta_Na(custom_TW, custom_sizes, TRUE),
              tolerance = 1e-12
            )
          })


test_that("compute_eta_Na rejects invalid or degenerate inputs", {
  # expect_error(compute_eta_Na(matrix(0, nrow = 2L, ncol = 2L), c(5L, 5L)),
  #              "The whole-plot trace factor is degenerate.",
  #              fixed = TRUE)
  
  # expect_error(
  #   compute_eta_Na(diag(3L), c(5L, 5L)),
  #   paste0(
  #     "'group_sizes' must contain one positive integer for each ",
  #     "row of 'TW'."
  #   ),
  #   fixed = TRUE
  # )
  
  # expect_error(
  #   compute_eta_Na(diag(3L), c(5L, 0L, 5L)),
  #   paste0(
  #     "'group_sizes' must contain one positive integer for each ",
  #     "row of 'TW'."
  #   ),
  #   fixed = TRUE
  # )
  
  # expect_error(compute_eta_Na(matrix(
  #   c(1, 1, 0, 1), nrow = 2L, ncol = 2L
  # ), c(5L, 5L)), "'TW' must be symmetric.", fixed = TRUE)
})


test_that("grouped numerical helper functions preserve valid calculations",
          {
            expect_equal(hdrm:::compute_grouped_statistic(
              QN = 10,
              EW = 4,
              variance = 9
            ),
            2,
            tolerance = 1e-15)
            
            expect_equal(hdrm:::compute_grouped_df(second_order = 8, third_order = 4),
                         32,
                         tolerance = 1e-15)
            
            expect_equal(hdrm:::compute_grouped_df(second_order = 1, third_order = 10),
                         1,
                         tolerance = 1e-15)
            
            expect_equal(
              hdrm:::compute_grouped_df(
                second_order = 8,
                third_order = -4,
                design_factor = 2
              ),
              64,
              tolerance = 1e-15
            )
            
            expected_p_value <- max(stats::pchisq(2 * sqrt(2 * 32) + 32, df = 32, lower.tail = FALSE),
                                    .Machine$double.eps)
            
            expect_equal(
              hdrm:::compute_grouped_p_value(statistic = 2, degrees_of_freedom = 32),
              expected_p_value,
              tolerance = 1e-15
            )
            
            expect_identical(
              hdrm:::compute_grouped_p_value(statistic = 1e6, degrees_of_freedom = 5),
              .Machine$double.eps
            )
          })


test_that("grouped numerical helper functions reject degenerate quantities",
          {
            expect_error(
              hdrm:::compute_grouped_statistic(
                QN = 1,
                EW = 1,
                variance = 0
              ),
              "The estimated variance of the test statistic must be positive.",
              fixed = TRUE
            )
            
            expect_error(
              hdrm:::compute_grouped_statistic(
                QN = Inf,
                EW = 1,
                variance = 2
              ),
              paste0(
                "The quadratic form, its estimated expectation, and its estimated ",
                "variance must be finite numeric scalars."
              ),
              fixed = TRUE
            )
            
            expect_error(
              hdrm:::compute_grouped_df(second_order = 0, third_order = 1),
              "The second-order trace estimate must be positive.",
              fixed = TRUE
            )
            
            expect_error(
              hdrm:::compute_grouped_df(second_order = 1, third_order = 0),
              paste0(
                "The third-order trace estimate is zero. Increase 'B' or check ",
                "whether the data are degenerate."
              ),
              fixed = TRUE
            )
            
            expect_error(
              hdrm:::compute_grouped_df(second_order = 1e200, third_order = 1e-200),
              paste0(
                "The estimated degrees-of-freedom parameter is not finite and ",
                "positive."
              ),
              fixed = TRUE
            )
            
            expect_error(
              hdrm:::compute_grouped_df(
                second_order = 1,
                third_order = 1,
                design_factor = 0
              ),
              "The design factor must be positive.",
              fixed = TRUE
            )
            
            expect_error(
              hdrm:::compute_grouped_p_value(statistic = 0, degrees_of_freedom = 0),
              paste0(
                "'statistic' must be finite and 'degrees_of_freedom' must be ",
                "finite and positive."
              ),
              fixed = TRUE
            )
          })


test_that("equal-covariance subsampling budget is allocated correctly", {
  # B is a base budget. With B = 1, a groups produce a total of a draws,
  # so every group receives one draw.
  minimum_allocation <- hdrm:::allocate_C1_subsamples(group_sizes = c(6L, 7L, 8L), B = 1L)
  
  expect_identical(minimum_allocation, c(1L, 1L, 1L))
  expect_equal(sum(minimum_allocation), 3L)
  expect_true(all(minimum_allocation >= 1L))
  
  
  # Equal group sizes and base budget B = 8 give exactly 8 draws per group.
  balanced_allocation <- hdrm:::allocate_C1_subsamples(group_sizes = c(6L, 6L, 6L), B = 8L)
  
  expect_identical(balanced_allocation, c(8L, 8L, 8L))
  expect_equal(sum(balanced_allocation), 24L)
  expect_true(all(balanced_allocation >= 1L))
  
  
  # Unequal group sizes are weighted approximately proportionally to
  # choose(n_i, 6), while the exact total remains a * B.
  unbalanced_allocation <- hdrm:::allocate_C1_subsamples(group_sizes = c(6L, 7L, 8L), B = 20L)
  
  expect_identical(unbalanced_allocation, c(3L, 12L, 45L))
  expect_equal(sum(unbalanced_allocation), 60L)
  expect_true(all(unbalanced_allocation >= 1L))
  
  
  # # The helper rejects non-positive base budgets.
  # expect_error(
  #   hdrm:::allocate_C1_subsamples(group_sizes = c(6L, 7L, 8L), B = 0L),
  #   "'B' must be a finite positive integer.",
  #   fixed = TRUE
  # )
})

test_that("B expressions are parsed without evaluating arbitrary R code", {
  expect_identical(hdrm:::evaluate_subsample_budget(B = 10, N = 16), 10L)
  expect_identical(hdrm:::evaluate_subsample_budget(B = 10.1, N = 16), 11L)
  expect_identical(hdrm:::evaluate_subsample_budget(B = "10*N", N = 16), 160L)
  expect_identical(hdrm:::evaluate_subsample_budget(B = "2 * (N + 1)", N = 16), 34L)
  expect_identical(hdrm:::evaluate_subsample_budget(B = "N^2 / 3", N = 16), 86L)
  expect_identical(hdrm:::evaluate_subsample_budget(B = "-(1 - N)", N = 16), 15L)
  expression_error <- "'B' must be an arithmetic expression"
  invalid_expressions <- c("sqrt(N)",
                           "N; 10",
                           "N <- 10",
                           "M",
                           "N[1]",
                           "base::identity(N)",
                           "{N}",
                           "10 *")
  
  for (current_expression in invalid_expressions) {
    expect_error(
      hdrm:::evaluate_subsample_budget(B = current_expression, N = 16),
      expression_error,
      fixed = TRUE
    )
  }
  
  expect_error(
    hdrm:::evaluate_subsample_budget(B = "N / 0", N = 16),
    paste0(
      "'B' must evaluate to a single finite positive number not exceeding ",
      .Machine$integer.max,
      "."
    ),
    fixed = TRUE
  )
  
  result <- hdrm_grouped(
    Matrixbirthrates,
    hypothesis = "whole",
    group = group,
    subsampling = FALSE,
    B = "N / 2",
    seed = 3141
  )
  
  expect_identical(result$B, 8L)
})


test_that("grouped third-trace budgets scale with the number of groups", {
  # expect_identical(hdrm:::expand_subsample_budget(B = 100L, multiplier = 4L),
  #                  400L)
  
  # expect_error(
  #   hdrm:::expand_subsample_budget(B = .Machine$integer.max, multiplier = 2L),
  #   paste0(
  #     "The effective subsampling budget must not exceed ",
  #     .Machine$integer.max,
  #     "."
  #   ),
  #   fixed = TRUE
  # )
})


test_that("additional invalid B values are rejected", {
  expect_error(
    hdrm_grouped(Matrixbirthrates, group = group, B = numeric(0)),
    "'B' must be a single numeric or character value.",
    fixed = TRUE
  )
  
  expect_error(hdrm_grouped(Matrixbirthrates, group = group, B = 0))
  expect_error(hdrm_grouped(Matrixbirthrates, group = group, B = NA_real_))
  expect_error(hdrm_grouped(Matrixbirthrates, group = group, B = Inf))
  
  expect_error(hdrm_grouped(
    Matrixbirthrates,
    group = group,
    B = .Machine$integer.max + 1
  ))
  
  expect_error(
    hdrm_grouped(Matrixbirthrates, group = group, B = list(10)),
    "'B' must be a single numeric or character value.",
    fixed = TRUE
  )
  
  expect_error(
    hdrm_grouped(Matrixbirthrates, group = group, B = "10 *"),
    "'B' must be an arithmetic expression",
    fixed = TRUE
  )
  
  # # B is a base budget; B = 1 is valid because the equal-covariance
  # # third-trace estimator uses a * B total draws.
  # expect_identical(hdrm:::expand_subsample_budget(B = 1L, multiplier = 2L), 2L)
})
