test_that("single-group C++ trace estimators agree with R references", {
  B0_R <- function(X) {
    N <- ncol(X)
    out <- 0
    
    for (i in seq_len(N)) {
      out <- out + sum(X[, i]^2)
    }
    
    out / N
  }
  
  B2_R <- function(X) {
    N <- ncol(X)
    out <- 0
    
    for (i in seq_len(N)) {
      for (j in seq_len(N)) {
        if (i != j) {
          inner_product <- sum(X[, i] * X[, j])
          
          out <- out + inner_product^2
        }
      }
    }
    
    out / (N * (N - 1))
  }
  
  B3_R <- function(X) {
    triples <- combn(seq_len(ncol(X)), 3L)
    
    sum(apply(triples, 2L, function(index) {
      i <- index[[1L]]
      j <- index[[2L]]
      k <- index[[3L]]
      
      sum(X[, i] * X[, j]) *
        sum(X[, j] * X[, k]) *
        sum(X[, k] * X[, i])
    })) / choose(ncol(X), 3)
  }
  
  X_fixed <- matrix(1:8, nrow = 2L, ncol = 4L)
  
  set.seed(123)
  
  X_random <- matrix(rnorm(5L * 7L), nrow = 5L, ncol = 7L)
  
  X_one_dimension <- matrix(c(-2, -0.5, 1, 3, 4.5), nrow = 1L)
  
  matrices <- list(fixed = X_fixed,
                   random = X_random,
                   one_dimension = X_one_dimension)
  
  # Independently known values for the fixed matrix
  expect_equal(hdrm:::B0_cpp(X_fixed), 51, tolerance = 1e-12)
  
  expect_equal(hdrm:::B2_cpp(X_fixed), 6079 / 3, tolerance = 1e-12)
  
  for (current_matrix in matrices) {
    expect_equal(hdrm:::B0_cpp(current_matrix),
                 B0_R(current_matrix),
                 tolerance = 1e-12)
    
    expect_equal(hdrm:::B2_cpp(current_matrix),
                 B2_R(current_matrix),
                 tolerance = 1e-12)
    
    cpp_B3 <- hdrm:::B3_cpp(current_matrix)
    r_B3 <- B3_R(current_matrix)
    
    expect_equal(cpp_B3, r_B3, tolerance = 1e-12)
  }
})
