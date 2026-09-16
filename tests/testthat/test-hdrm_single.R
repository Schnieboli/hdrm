#Tests for vector format

# Loading the dataset
data(EEG)

L <- list(1:160, 1:40)


test_that("perfect case works", {
  # hypothesis = flat
  expect_no_condition(hdrm_single(
    EEG$value,
    hypothesis = "flat",
    subject = EEG$subject
  ))
  
  
  # Multiple hypothesis values are rejected
  expect_error(
    hdrm_single(
      EEG$value,
      hypothesis = c("flat", "SUB"),
      subject = EEG$subject
    ),
    "'hypothesis' must be a single non-missing character value.",
    fixed = TRUE
  )
  
  # legal list to hypothesis
  expect_no_condition(hdrm_single(
    EEG$value,
    hypothesis = diag(40),
    subject = EEG$subject
  ))
  
  
  # unknown argument is rejected
  expect_error(
    hdrm_single(
      EEG$value,
      hypothesis = "flat",
      subject = EEG$subject,
      subsampling = FALSE
    )
  )
})


test_that("wrong input: data", {
  # data as list
  expect_error(hdrm_single(L, hypothesis = "flat", subject = EEG$subject))
})


test_that("missing values", {
  # NA in value
  vector <- EEG$value
  vector[123] <- NA
  expect_warning(hdrm_single(vector, hypothesis = "flat", subject = EEG$subject))
  
  
  # NA in subject
  subjectvector <- EEG$subject
  subjectvector[145] <- NA
  expect_error(hdrm_single(EEG$value, hypothesis = "flat", subject = subjectvector))
  
  
  # different length of subject and data
  expect_error(hdrm_single(
    EEG$value[-1],
    hypothesis = "flat",
    subject = EEG$subject
  ))
  
  # nonexistant column
  expect_error(hdrm_single(EEG$value, hypothesis = "flat", subject = "nonexistant"))
  
})


test_that("wrong input: hypothesis", {
  # a number
  expect_error(hdrm_single(
    EEG$value,
    hypothesis = 1,
    subject = EEG$subject
  ))
  
  # illegal character
  expect_error(hdrm_single(
    EEG$value,
    hypothesis = c("flart"),
    subject = EEG$subject
  ))
  
  
  # list TS missing
  expect_error(hdrm_single(
    EEG$value,
    hypothesis = diag(1:41),
    subject = EEG$subject
  ))
  
  
  
  # wrong dimension of hypothesis
  expect_error(hdrm_single(
    EEG$value,
    hypothesis = matrix(0, 41, 39),
    subject = EEG$subject
  ))
  
  # TW not symmetrical
  expect_error(hdrm_single(
    EEG$value,
    hypothesis = list(TW = diag(4) + c(0, 1), TS = diag(40)),
    subject = EEG$subject
  ))
  
  # list
  expect_error(hdrm_single(
    EEG$value,
    hypothesis = list(diag(40)),
    subject = EEG$subject
  ))
  
})





test_that("hdrm_single test statistics", {
  # hypothesis = flat
  expect_equal(
    hdrm_single(
      EEG$value,
      hypothesis = "flat",
      subject = EEG$subject
    )$statistic,
    110.236926
  )
})


test_that("hdrm_single  p.values", {
  # hypothesis = flat
  expect_equal(
    hdrm_single(
      EEG$value,
      hypothesis = "flat",
      subject = EEG$subject
    )$p.value,
    2.220446e-16
  )
})

test_that("hdrm_single f", {
  # hypothesis = flat
  expect_equal(
    hdrm_single(
      EEG$value,
      hypothesis = "flat",
      subject = EEG$subject
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
    "'hypothesis' must be a single non-missing character value.",
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
  
  # data as matrix with additional subject values
  expect_warning(hdrm_single(
    Matrixbirthrates,
    hypothesis = "flat",
    subject = seq_len(nrow(Matrixbirthrates))
  ))
  # data as df
  expect_error(hdrm_single(data = EEG, hypothesis = "flat"))
})


test_that("missing values", {
  # NA in value
  M <- Matrixbirthrates
  M[4, 7] <- NA
  expect_warning(hdrm_single(data = M, hypothesis = "flat"))
  
  expect_warning(expect_true(all(
    dim(hdrm_single(
      data = M, hypothesis = "flat"
    )$data) == dim(M) - c(1, 0)
  )))
  
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
  # AM must be binary
  expect_error(
    hdrm_single(
      EEG$value,
      hypothesis = "flat",
      subject = EEG$subject,
      AM = 2
    ),
    "'AM' must be a single logical value or 0/1.",
    fixed = TRUE
  )
  
  # AM must not be missing
  expect_error(
    hdrm_single(
      EEG$value,
      hypothesis = "flat",
      subject = EEG$subject,
      AM = NA
    ),
    "'AM' must be a single logical value or 0/1.",
    fixed = TRUE
  )
  
  # AM must have length one
  expect_error(
    hdrm_single(
      EEG$value,
      hypothesis = "flat",
      subject = EEG$subject,
      AM = c(0, 1)
    ),
    "'AM' must be a single logical value or 0/1.",
    fixed = TRUE
  )
  
  # character values are rejected
  expect_error(
    hdrm_single(
      EEG$value,
      hypothesis = "flat",
      subject = EEG$subject,
      AM = "TRUE"
    ),
    "'AM' must be a single logical value or 0/1.",
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
  expect_error(hdrm_single(numeric(0), subject = character(0)),
               "'data' must not be empty.",
               fixed = TRUE)
  
  expect_error(hdrm_single(matrix(
    numeric(0), nrow = 0L, ncol = 0L
  )), "'data' must not be empty.", fixed = TRUE)
  
  expect_error(hdrm_single(EEG$value),
               "'subject' must be provided when 'data' is a vector.",
               fixed = TRUE)
  
  expect_error(
    hdrm_single(EEG$value, subject = EEG$subject[-1L]),
    "The lengths of 'data' and 'subject' must be equal.",
    fixed = TRUE
  )
  
  expect_error(
    hdrm_single(EEG$value, subject = as.list(EEG$subject)),
    "'subject' must be a one-dimensional atomic vector or factor.",
    fixed = TRUE
  )
  
  expect_warning(
    hdrm_single(Matrixbirthrates, subject = seq_len(nrow(Matrixbirthrates))),
    "'subject' is ignored when 'data' is a matrix.",
    fixed = TRUE
  )
})


test_that("named numeric vectors are accepted", {
  named_values <- EEG$value
  names(named_values) <- seq_along(named_values)
  
  expect_no_condition(hdrm_single(
    named_values,
    hypothesis = "flat",
    subject = EEG$subject
  ))
})


test_that("vector preprocessing removes complete subject blocks", {
  data_with_na <- EEG$value
  removed_subject <- EEG$subject[[1L]]
  removed_index <- which(EEG$subject == removed_subject)[[1L]]
  data_with_na[removed_index] <- NA_real_
  
  result <- NULL
  
  expect_warning(
    result <- hdrm_single(
      data_with_na,
      hypothesis = "flat",
      subject = EEG$subject
    ),
    "Subjects with missing values dropped",
    fixed = TRUE
  )
  
  expect_equal(result$removed.cases, 1L)
  expect_equal(nrow(result$data), length(unique(EEG$subject)) - 1L)
  expect_false(anyNA(result$data))
})


test_that("matrix preprocessing is reflected in the returned object", {
  data_with_na <- Matrixbirthrates
  data_with_na[4L, 7L] <- NA_real_
  
  result <- NULL
  
  expect_warning(
    result <- hdrm_single(data_with_na, hypothesis = "flat"),
    "Subjects with missing values dropped",
    fixed = TRUE
  )
  
  complete_subjects <- stats::complete.cases(data_with_na)
  
  expect_equal(result$data, data_with_na[complete_subjects, , drop = FALSE])
  expect_equal(result$removed.cases, 1L)
  expect_equal(result$dim[["N"]], sum(complete_subjects))
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

test_that("within-subject measurement order is preserved", {
  subject_blocks <- split(seq_along(EEG$value), as.character(EEG$subject))
  
  original_subject_order <- unique(as.character(EEG$subject))
  reversed_subject_order <- rev(original_subject_order)
  
  permutation <- unlist(subject_blocks[reversed_subject_order], use.names = FALSE)
  
  result_original <- hdrm_single(EEG$value,
                                 hypothesis = "flat",
                                 subject = EEG$subject)
  
  result_permuted <- hdrm_single(EEG$value[permutation],
                                 hypothesis = "flat",
                                 subject = EEG$subject[permutation])
  
  expected_subject_order <- match(reversed_subject_order, original_subject_order)
  
  expect_equal(result_permuted$data, result_original$data[expected_subject_order, , drop = FALSE])
  expect_equal(result_permuted$statistic, result_original$statistic)
  expect_equal(result_permuted$p.value, result_original$p.value)
  expect_equal(result_permuted$f, result_original$f)
})


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
  "At least three complete subjects are required.",
  fixed = TRUE)
})


test_that("non-finite data and degenerate hypotheses are rejected", {
  data_with_inf <- Matrixbirthrates
  data_with_inf[1L, 1L] <- Inf
  
  expect_error(
    hdrm_single(data_with_inf, hypothesis = "flat"),
    "'X' must contain only finite, non-missing values.",
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
    paste0(
      "The hypothesis matrix must contain only finite, ",
      "non-missing values."
    ),
    fixed = TRUE
  )
})


test_that("AM representations agree for the one-group method", {
  result_am0 <- hdrm_single(
    EEG$value,
    hypothesis = "flat",
    subject = EEG$subject,
    AM = FALSE
  )
  
  result_am1 <- hdrm_single(
    EEG$value,
    hypothesis = "flat",
    subject = EEG$subject,
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
