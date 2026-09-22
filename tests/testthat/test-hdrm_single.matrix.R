data("EEG")
EEGmatrix <- hdrm_grouped(EEG, group = EEG$group, B = 2)$data


legal <- list()
legal$hypotheses <- list("flat",
                         diag(40),
                         diag(40) - matrix(1/40, 40, 40)
)
legal$AM <- c(TRUE, FALSE)


test_that("perfect input does not produce any conditions and results can be reproduced from output", {
  
  # legal list to hypothesis
  for (hypothesis in legal$hypotheses) {
    for (AM in legal$AM)
      expect_no_condition({
        initial <- hdrm_single(
          EEGmatrix,
          hypothesis = hypothesis,
          AM = AM
        )
        
        initial$H$label <- "custom"
        
        secondary <- hdrm_single(
          initial$data,
          hypothesis = initial$H$T,
          AM = initial$AM
        )
      })
    expect_identical(initial, secondary) 
  }
})


test_that("illegal 'data' input", {
  ## test missing values, non-finite values and non-numeric values
  mat_NA <- mat_Inf <- mat_char <- EEGmatrix
  mat_NA[1,1] <- NA
  mat_Inf[1,1] <- Inf
  mat_char <- matrix(as.character(mat_char), nrow = nrow(EEGmatrix))
  
  for(M in list(mat_NA, mat_Inf, mat_char))
    expect_error(
      hdrm_single(M),
      "'data' must be a numeric matrix, containing only finite, non-missing values."
    )
  
  expect_error(
    hdrm_single(matrix(numeric(0), nrow = 0, ncol = 0))
    )
  
})

test_that("wrong input: hypothesis", {
  TM_idem <- diag(40)
  TM_idem[c(1,40),40] <- c(1,0)
  TM_symm <- 2*diag(34)
  illegal_hyp <- list("flart",             ## illegal char
                      c("flat, flat"),     ## two legal chars
                      1,                   ## numeric value
                      list(TM = diag(40)), ## list input
                      TM_idem,             ## TM not symmetrical
                      TM_symm              ## TM not idempotent
  )
  
  for(h in illegal_hyp){
    expect_error(
      hdrm_single(
        EEGmatrix,
        hypothesis = h
      )
    )
  }
})


test_that("hdrm_single test statistics", {
  res_flat <- hdrm_single(EEGmatrix, hypothesis = "flat")
  
  expect_equal(res_flat$statistic, 110.236926, tolerance = 1e-6)
  expect_equal(res_flat$p.value, 2.220446e-16)
  expect_equal(res_flat$f, 1.002101918, tolerance = 1e-8)
})


test_that("AM representations agree for the one-group method", {
  hyp <- diag(40)
  hyp[1,1] <- 0
  result_am0 <- hdrm_single(
    EEGmatrix,
    hypothesis = hyp,
    AM = FALSE
  )
  result_am1 <- hdrm_single(
    EEGmatrix,
    hypothesis = hyp,
    AM = TRUE
  )
  expect_equal(result_am0$statistic, result_am1$statistic, tolerance = 1e-12)
  expect_equal(result_am0$p.value, result_am1$p.value, tolerance = 1e-12)
  expect_equal(result_am0$f, result_am1$f, tolerance = 1e-12)
  expect_identical(result_am0$H, result_am1$H, tolerance = 1e-12)
})
