# Loading the dataset
data(EEG)

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
          EEG,
          hypothesis = hypothesis,
          AM = AM
        )
        secondary <- hdrm_single(
          initial$data,
          hypothesis = initial$H,
          AM = initial$AM
        )
      })
    expect_identical(initial, secondary) 
  }
})


test_that("permutations result in the same test result",{
  ## test if permuted data yields identical output
  perm <- sample.int(nrow(EEG))
  
  res_sorted <- hdrm_single(EEG)
  res_perm <- hdrm_single(EEG[perm, ])
  expect_identical(res_sorted, res_perm)
})


test_that("illegal 'data' input", {
  ## test NAs and Inf in each column of data and non fitting levels
  for(col in c("value", "subject", "dimension")){
    for(illegal in c(NA, Inf, 2)){
      if(col == "value" & isTRUE(illegal == 2)) next;
      df <- EEG
      suppressWarnings(df[[col]][1] <- illegal)
      expect_error(
        hdrm_single(df)
      )
    }
  }
  ## missing colums
  for(col in c("value", "subject", "dimension")){
    df <- EEG
    df[[col]] <- NULL
    expect_error(
      hdrm_single(df), 
      "'data' must contain columns 'value', 'subject' and 'dimension'"
    )
  }
  expect_error(
    hdrm_grouped(
      data.frame(value = numeric(0), subject = numeric(0), dimension = numeric(0)),
      group = group_int,
      B = 10
    ), "'data' must not be empty.", fixed = TRUE)
  
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
      hdrm_single(EEG, hypothesis = h)
    )
  }
})


test_that("non continuous levels in data work", {
  # Non-continuous subject levels
  df <- EEG
  levels(df$subject) <- c(1, 3:161)
  levels(df$dimension) <- 2:41
  
  res_normal <- hdrm_single(df)
  expect_no_condition(
    res_non_cont <- hdrm_single(df)
  )
  expect_equal(res_normal$statistic, res_non_cont$statistic)
  expect_equal(res_normal$f, res_non_cont$f)
  expect_equal(res_normal$p.value, res_non_cont$p.value)
})



test_that("hdrm_single test statistics", {
  res_flat <- hdrm_single(EEG, hypothesis = "flat")
  
  expect_equal(res_flat$statistic, 110.236926, tolerance = 1e-6)
  expect_equal(res_flat$p.value, 2.220446e-16)
  expect_equal(res_flat$f, 1.002101918, tolerance = 1e-8)
})

test_that("AM representations agree for the one-group method", {
  hyp <- diag(40)
  hyp[1,1] <- 0
  result_am0 <- hdrm_single(EEG, hypothesis = hyp, AM = FALSE)
  result_am1 <- hdrm_single(EEG, hypothesis = hyp, AM = TRUE)
  
  expect_equal(result_am0$statistic, result_am1$statistic, tolerance = 1e-12)
  expect_equal(result_am0$p.value, result_am1$p.value, tolerance = 1e-12)
  expect_equal(result_am0$f, result_am1$f, tolerance = 1e-12)
  expect_identical(result_am0$H, result_am1$H, tolerance = 1e-12)
})
