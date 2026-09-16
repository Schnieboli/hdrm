#Tests for vector format

# Loading the dataset
data(EEG)

L <- list(1:160, 1:40)

df <- data.frame(value = EEG$value, subject = EEG$subject, dimension = EEG$dimension)

test_that("perfect case works", {
  # hypothesis = flat
  expect_no_condition(hdrm_single(
    df,
    hypothesis = "flat"
  ))
  
  
  # Multiple hypothesis values are rejected
  expect_error(
    hdrm_single(
      df,
      hypothesis = c("flat", "SUB")
    ),
    "'hypothesis' must be 'flat' or a numeric projection matrix containing only finite, non-missing values.",
    fixed = TRUE
  )
  
  # legal list to hypothesis
  expect_no_condition(hdrm_single(
    df,
    hypothesis = diag(40)
  ))
  
  
  # unknown argument is rejected
  expect_error(
    hdrm_single(
      df,
      hypothesis = "flat",
      subsampling = FALSE
    )
  )
})


test_that("wrong input: data", {
  # data as list
  expect_error(hdrm_single(L, hypothesis = "flat"))
})


test_that("missing values", {
  # NA in value
  df2 <- df
  df2$value[123] <- NA
  expect_error(hdrm_single(df2, hypothesis = "flat"),
               "'data' must not contain any missing values.")
  
  
  # NA in subject
  df2 <- df
  df2$subject[145] <- NA
  expect_error(hdrm_single(df2, hypothesis = "flat"))
  
  
  # NA in dimension
  df2 <- df
  df2$dimension[145] <- NA
  expect_error(hdrm_single(df2, hypothesis = "flat"))
  
  # different length of subject and data
  
  expect_error(hdrm_single(
    df[-1, ],
    hypothesis = "flat"
  ), "each combination of subject and dimension must occur exactly one.")
  
  # nonexistant column
  df2 <- data.frame(df$value, df$subject)
  expect_error(hdrm_single(df2, hypothesis = "flat"),
               "'data' must contain columns 'value', 'subject' and 'dimension'")
  
})


test_that("wrong input: hypothesis", {
  # a number
  expect_error(hdrm_single(
    df,
    hypothesis = 1,
  ))
  
  # illegal character
  expect_error(hdrm_single(
    df,
    hypothesis = c("flart")
  ))
  
  
  # list TS missing
  expect_error(hdrm_single(
    df,
    hypothesis = diag(1:41)
  ))
  
  
  
  # wrong dimension of hypothesis
  expect_error(hdrm_single(
    df,
    hypothesis = matrix(0, 41, 39)
  ))
  
  # TW not symmetrical
  expect_error(hdrm_single(
    df,
    hypothesis = list(TW = diag(4) + c(0, 1), TS = diag(40))
  ))
  
  # list
  expect_error(hdrm_single(
    df,
    hypothesis = list(diag(40))
  ))
  
})


test_that("hdrm_single test statistics", {
  # hypothesis = flat
  expect_equal(
    hdrm_single(
      df,
      hypothesis = "flat"
    )$statistic,
    110.236926
  )
})


test_that("hdrm_single  p.values", {
  # hypothesis = flat
  expect_equal(
    hdrm_single(
      df,
      hypothesis = "flat"
    )$p.value,
    2.220446e-16
  )
})

test_that("hdrm_single f", {
  # hypothesis = flat
  expect_equal(
    hdrm_single(
      df,
      hypothesis = "flat"
    )$f,
    1.0021019186641513
  )
})


#Tests for matrix format


data("birthrates")
Matrixbirthrates = t(as.matrix(birthrates))
L2 <- list(1:160, 1:40)


test_that("perfect case works", {
  # hypothesis = flat
  expect_no_condition(hdrm_single(Matrixbirthrates, hypothesis = "flat"))
  
  
  # Multiple hypothesis values are rejected
  expect_error(
    hdrm_single(Matrixbirthrates, hypothesis = c("flat", "SUB")),
    "'hypothesis' must be 'flat' or a numeric projection matrix containing only finite, non-missing values.",
    fixed = TRUE
  )
  
  # legal list to hypothesis
  expect_no_condition(hdrm_single(Matrixbirthrates, hypothesis = diag(34)))
  
  
  # unknown argument is rejected
  expect_error(hdrm_single(
    Matrixbirthrates,
    hypothesis = "flat",
    subsampling = FALSE
  ))
})

test_that("wrong input: data", {
  # data as list
  expect_error(
    hdrm_single(
      L2,
      hypothesis = "flat",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )
})


test_that("missing values", {
  # NA in value
  M <- Matrixbirthrates
  M[4, 7] <- NA
  expect_error(hdrm_single(data = M, hypothesis = "flat"),
               "'data' must be a numeric matrix, containing only finite, non-missing values")
})


test_that("wrong input: hypothesis", {
  # a number
  expect_error(hdrm_single(Matrixbirthrates, hypothesis = 1))
  
  # illegal character
  expect_error(hdrm_single(Matrixbirthrates, hypothesis = c("flart")))
  
  
  # wrong dimension of hypothesis
  expect_error(hdrm_single(Matrixbirthrates, hypothesis = matrix(0, 41, 39)))
  
  # TW not symmetrical
  expect_error(hdrm_single(Matrixbirthrates, hypothesis = diag(34) + c(0, 1)))
  
  # list
  expect_error(hdrm_single(Matrixbirthrates, hypothesis = list(diag(40))))
  
})

test_that("wrong input: AM", {
  # AM must not be missing
  expect_error(
    hdrm_single(
      df,
      hypothesis = "flat",
      AM = NA
    ),
    "'AM' must be a single logical value.",
    fixed = TRUE
  )
  
  # AM must have length one
  expect_error(
    hdrm_single(
      df,
      hypothesis = "flat",
      AM = c(0, 1)
    ),
    "'AM' must be a single logical value.",
    fixed = TRUE
  )
})

test_that("hdrm_single test statistics", {
  # hypothesis = flat
  expect_equal(hdrm_single(Matrixbirthrates, hypothesis = "flat")$statistic,
               7.481676)
})


test_that("hdrm_single p.values", {
  # hypothesis = flat
  expect_equal(
    hdrm_single(Matrixbirthrates, hypothesis = "flat")$p.value,
    0.00037711148523941827
    ,
    tolerance = 0.0000002
  )
})


test_that("hdrm_single f", {
  # hypothesis = flat
  expect_equal(hdrm_single(Matrixbirthrates, hypothesis = "flat")$f,
               1.4360121127976566)
})

# Additional tests for the revised public interface and preprocessing ---------

test_that("data and subject inputs are validated explicitly", {
  df2 <- data.frame(value = numeric(0), subject = integer(0), dimension = integer(0))
  expect_error(hdrm_single(df2, "flat"),
               "'data' must not be empty.",
               fixed = TRUE)
  
  expect_error(hdrm_single(matrix(
    numeric(0), nrow = 0L, ncol = 0L)
    ), "'data' must be a numeric matrix, containing only finite, non-missing values", fixed = TRUE)
})


test_that("wide matrix input uses subjects in rows and dimensions in columns",
          {
            result <- hdrm_single(Matrixbirthrates, hypothesis = "flat")
            
            expect_equal(result$dim, c(
              d = ncol(Matrixbirthrates),
              N = nrow(Matrixbirthrates)
            ))
            expect_equal(result$data, Matrixbirthrates)
          })

print("auskommentiert: within-subject measurement order is preserved")
#### Grund: diese eigenschaft ist nicht unnötig, aber auch nicht unbedingt
#### notwendig -> mal sehen, ob man dan och was dribbeln kann
# test_that("within-subject measurement order is preserved", {
#   subject_blocks <- split(seq_along(df$value), as.character(df$subject))
#
#   original_subject_order <- unique(as.character(df$subject))
#   reversed_subject_order <- rev(original_subject_order)
#
#   permutation <- unlist(subject_blocks[reversed_subject_order], use.names = FALSE)
#
#   df2 <- df[permutation, ]
#
#   result_original <- hdrm_single(df, hypothesis = "flat",)
#
#   result_permuted <- hdrm_single(df2, hypothesis = "flat",)
#
#   expected_subject_order <- match(reversed_subject_order, original_subject_order)
#
#   expect_equal(result_permuted$data, result_original$data[expected_subject_order, , drop = FALSE])
#   expect_equal(result_permuted$statistic, result_original$statistic)
#   expect_equal(result_permuted$p.value, result_original$p.value)
#   expect_equal(result_permuted$f, result_original$f)
# })


test_that("minimum dimensions and complete subject counts are enforced", {
  expect_error(
    hdrm_single(matrix(
      seq_len(6), nrow = 6L, ncol = 1L
    )),
    "At least two repeated-measurement dimensions are required.",
    fixed = TRUE
  )
  
  expect_error(hdrm_single(matrix(
    seq_len(4), nrow = 2L, ncol = 2L
  )),
  "At least three subjects are required.",
  fixed = TRUE)
})


test_that("non-finite data and degenerate hypotheses are rejected", {
  data_with_inf <- Matrixbirthrates
  data_with_inf[1L, 1L] <- Inf
  
  expect_error(
    hdrm_single(data_with_inf, hypothesis = "flat"),
    "'data' must be a numeric matrix, containing only finite, non-missing values.",
    fixed = TRUE
  )
  
  expect_error(
    hdrm_single(Matrixbirthrates, hypothesis = matrix(
      0,
      nrow = ncol(Matrixbirthrates),
      ncol = ncol(Matrixbirthrates)
    )),
    "The hypothesis matrix must have positive rank.",
    fixed = TRUE
  )
  
  hypothesis_with_inf <- diag(ncol(Matrixbirthrates))
  hypothesis_with_inf[1L, 1L] <- Inf
  
  expect_error(
    hdrm_single(Matrixbirthrates, hypothesis = hypothesis_with_inf),
    "'hypothesis' must be 'flat' or a numeric projection matrix containing only finite, non-missing values",
    fixed = TRUE
  )
})


test_that("AM representations agree for the one-group method", {
  result_am0 <- hdrm_single(
    df,
    hypothesis = "flat",
    AM = FALSE
  )
  
  result_am1 <- hdrm_single(
    df,
    hypothesis = "flat",
    AM = TRUE
  )
  
  expect_equal(result_am0$statistic, result_am1$statistic, tolerance = 1e-12)
  expect_equal(result_am0$p.value, result_am1$p.value, tolerance = 1e-12)
  expect_equal(result_am0$f, result_am1$f, tolerance = 1e-12)
})


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
