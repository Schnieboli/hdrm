data(birthrates)
DFbirthrates <- data.frame(value = as.vector(as.matrix(birthrates)), 
                           subject = rep(1:16, each = 34),
                           dimension = rep(1:34, 16))
group_int <- rep(c(1, 1, 2, 2, 1, 1, 1, 2, 1, 1, 1, 1, 2, 2, 1, 2), each = 34)
group_fac <- factor(group_int, labels = c("west", "east"))


legal <- list()
legal$group <- list(int = group_int, fac = group_fac)
legal$hypotheses <- list("whole","sub","interaction","identical","flat",
                         "TW, TS, TR" = list(TW = diag(2), TS = diag(34), TR = diag(68)),
                         "TW, TS" = list(TW = diag(2), TS = diag(34)),
                         "TW, TS, TM" = list(TW = diag(2), TS = diag(34), TM = diag(68))
)
legal$subsampling <- legal$AM <- legal$cov.equal <- c(TRUE, FALSE)
legal$B <- list("2*N", 100)


test_that("perfect input does not produce any conditions and results can be reproduced from output", {
  
  # legal list to hypothesis
  for (hypothesis in legal$hypotheses) {
    for (subsampling in legal$subsampling)
      for (group in legal$group)
        for (AM in legal$AM)
          for (cov.equal in legal$cov.equal)
            for (B in legal$B){
              expect_no_condition({
                initial <- hdrm_grouped(
                  DFbirthrates,
                  hypothesis = hypothesis,
                  group = group,
                  subsampling = subsampling,
                  AM = AM,
                  cov.equal = cov.equal,
                  B = B,
                  seed = 3141
                )
                
                secondary <- hdrm_grouped(
                  initial$data,
                  hypothesis = initial$H,
                  group = initial$group,
                  subsampling = initial$subsampling,
                  AM = initial$AM,
                  cov.equal = initial$cov.equal,
                  B = initial$B,
                  seed = initial$seed
                )
              })
              expect_identical(initial, secondary)
            }
  }
})


test_that("permutations result in the same test result",{
  ## test if permuted data yields identical output
  perm <- sample.int(nrow(DFbirthrates))
  df <- DFbirthrates[perm, ]
  
  res_sorted <- hdrm_grouped(
    DFbirthrates,
    hypothesis = "interaction",
    group = group_int,
    subsampling = TRUE,
    AM = TRUE,
    cov.equal = FALSE,
    B = 100,
    seed = 3141
  )
  res_perm <- hdrm_grouped(
    df,
    hypothesis = "interaction",
    group = group_int[perm],
    subsampling = TRUE,
    AM = TRUE,
    cov.equal = FALSE,
    B = 100,
    seed = 3141
  )
  expect_identical(res_sorted, res_perm)
})


test_that("illegal 'data' input", {
  ## test NAs and Inf in each column of data and non fitting levels
  for(col in names(DFbirthrates)){
    for(illegal in c(NA, Inf, 2)){
      if(col == "value" & isTRUE(illegal == 2)) next;
      df <- DFbirthrates
      df[[col]][1] <- illegal
      expect_error(
        hdrm_grouped(
          df,
          hypothesis = "flat",
          group = group_int,
          B = 10
        )
      )
    }
  }
  ## missing colums
  for(col in names(DFbirthrates)){
    df <- DFbirthrates
    df[[col]] <- NULL
    expect_error(
      hdrm_grouped(
        df,
        hypothesis = "flat",
        group = group_int,
        B = 10
      ), "'data' must contain columns 'value', 'subject' and 'dimension'"
    )
  }
  
  expect_error(
    hdrm_grouped(
      data.frame(value = numeric(0), subject = numeric(0), dimension = numeric(0)),
      group = group_int,
      B = 10
    ), "'data' must not be empty.", fixed = TRUE)
  
})


test_that("illegal 'group' input",{
  
  ## test NA, Inf in group or group wrong length
  group_NA <- c(NA, group_int[-1])
  group_Inf <- c(Inf, group_int[-1])
  group_short <- group_int[-1]
  group_long <- c(group_int, 1)
  
  for(group in list(group_NA, group_Inf)){
    expect_error(
      hdrm_grouped(
        DFbirthrates,
        hypothesis = "sub",
        group = group,
        subsampling = FALSE,
        B = "10*N"
      ), "'group' must only contain finite, non-missing values."
    )
  }
  for(group in list(group_short, group_long, integer(0))){
    expect_error(
      hdrm_grouped(
        DFbirthrates,
        hypothesis = "sub",
        group = group,
        subsampling = FALSE,
        B = "10*N"
      ), "'group' must be a one-dimensional vector or factor of length nrow\\(data\\)\\."
    )
  }
})


test_that("wrong input: hypothesis", {
  TW_idem <- matrix(c(1,0,1,0), 2)
  TW_symm <- 2*diag(2)
  TW_rank0 <- matrix(0, 2, 2)
  TS_idem <- diag(34)
  TS_idem[c(1,34),34] <- c(1,0)
  TS_symm <- 2*diag(34)
  TS_rank0 <- matrix(0, 34, 34)
  illegal_hyp <- list("flart",                            ## illegal char
                      c("flat, sub"),                     ## two legal chars
                      1,                                  ## numeric value
                      diag(68),                           ## matrix
                      list(TW = diag(2),  TS = diag(33)), ## wrong dimension TS
                      list(TW = diag(3),  TS = diag(34)), ## wrong dimension TW
                      list(WT = diag(2),  TS = diag(34)), ## TW missing
                      list(TW = diag(2),  ST = diag(34)), ## TS missing
                      list(TW = TW_symm,  TS = diag(34)), ## TW not idempotent
                      list(TW = diag(2),  TS = TS_symm),  ## TS not idempotent
                      list(TW = TW_idem , TS = diag(34)), ## TW not symmetrical
                      list(TW = diag(2),  TS = TS_idem),  ## TW not symmetrical
                      list(TW = TW_rank0, TS = diag(34)), ## TW rank 0
                      list(TW = diag(2), TS = TS_rank0)   ## TS rank 0
  )
  
  for(h in illegal_hyp){
    expect_error(
      hdrm_grouped(
        DFbirthrates,
        hypothesis = h,
        group = group_int,
        subsampling = FALSE,
        B = "10*N"
      )
    )
  }
})


test_that("false input: B", {
  illegal_B <- list(
    c("10*N", "20*N"),
    c(100, 1000),
    c(100, "10*N"),
    -10,
    "-10*N",
    "10 *asdfghjkl",
    diag(3)
  )
  for(B in illegal_B){
    expect_error(
      hdrm_grouped(
        DFbirthrates,
        hypothesis = "sub",
        group = group_int,
        subsampling = FALSE,
        B = B
      )
    )
  }
})


test_that("non continuous levels in data and group work", {
  # Non-continuous subject levels
  df <- DFbirthrates
  df$subject <- df$subject + 1
  df$subject[df$subject == 2] <- 1
  df$dimension <- df$dimension + 1
  group_tmp <- group_int
  group_tmp[group_tmp == 4] <- 5
  
  res_normal <- hdrm_grouped(
    DFbirthrates,
    hypothesis = "whole",
    group = group_int,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  expect_no_condition(
    res_non_cont <- hdrm_grouped(
      df,
      hypothesis = "whole",
      group = group_tmp,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )
  )
  expect_equal(res_normal$statistic, res_non_cont$statistic)
  expect_equal(res_normal$f, res_non_cont$f)
  expect_equal(res_normal$p.value, res_non_cont$p.value)
})


test_that("hdrm_grouped statistics, p-values, and f-degrees of freedom", {
  
  # Helper: chechs if  statistic, p.value and f are equal to expected values
  test_metrics <- function(res, expected_stat, expected_p, expected_f, tolerance = 1e-6) {
    expect_equal(res$statistic, expected_stat, tolerance = tolerance)
    expect_equal(res$p.value,   expected_p,    tolerance = tolerance)
    expect_equal(res$f,         expected_f,    tolerance = tolerance)
  }

  res_whole <- hdrm_grouped(
    data = DFbirthrates,
    hypothesis = "whole",
    group = group_int,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  res_sub <- hdrm_grouped(
    data = DFbirthrates,
    hypothesis = "sub",
    group = group_int,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  res_interaction <- hdrm_grouped(
    data = DFbirthrates,
    hypothesis = "interaction",
    group = group_int,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  res_identical <- hdrm_grouped(
    data = DFbirthrates,
    hypothesis = "identical",
    group = group_int,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  res_flat <- hdrm_grouped(
    data = DFbirthrates,
    hypothesis = "flat",
    group = group_int,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  # Subsampling = TRUE
  res_subsampling <- hdrm_grouped(
    DFbirthrates,
    hypothesis = "whole",
    group = group_int,
    subsampling = TRUE,
    B = "10*N",
    seed = 3141
  )
  # AM = 0
  res_am0 <- hdrm_grouped(
    DFbirthrates,
    hypothesis = "whole",
    group = group_int,
    subsampling = FALSE,
    AM = 0,
    B = "10*N",
    seed = 3141
  )
  # cov.equal = TRUE, AM = 0
  res_cov_equal <- hdrm_grouped(
    DFbirthrates,
    hypothesis = "whole",
    group = group_int,
    cov.equal = TRUE,
    AM = 0,
    B = "10*N",
    seed = 3141
  )
  ## check if values are the expected values
  test_metrics(res_whole, 24.71615921, 1.482681e-12, 2.335495)
  test_metrics(res_sub, 321.685616, 2.220446e-16, 5.076035)
  test_metrics(res_interaction, 132.295574, 2.220446e-16, 6.601591)
  test_metrics(res_identical, 88.55248, 2.220446e-16, 5.509633)
  test_metrics(res_flat, 358.685785, 2.220446e-16, 67.916386)
  test_metrics(res_subsampling, 23.298551, 2.654115e-12, 3.244517)
  test_metrics(res_am0, 24.716159, 1.482681e-12, 2.335495)
  test_metrics(res_cov_equal, 11.078986, 2.155735e-07, 4.595081)
})


test_that("AM representations agree", {
  ## check if AM = TRUE/FALSE yields the same result for different input cases
  for(cov.equal in c(FALSE, TRUE)){
    for(subsampling in c(FALSE, TRUE)){
      result_am0 <- hdrm_grouped(
        DFbirthrates,
        hypothesis = "whole",
        AM = FALSE,
        group = group_int,
        cov.equal = cov.equal,
        subsampling = subsampling,
        B = "10*N",
        seed = 3141
      )
      result_am1 <- hdrm_grouped(
        DFbirthrates,
        hypothesis = "whole",
        AM = TRUE,
        group = group_fac,
        cov.equal = cov.equal,
        subsampling = subsampling,
        B = "10*N",
        seed = 3141
      )
      expect_equal(result_am0$statistic, result_am1$statistic, tolerance = 1e-6)
      expect_equal(result_am0$p.value, result_am1$p.value, tolerance = 1e-6)
      expect_equal(result_am0$f, result_am1$f, tolerance = 1e-6)
    }
  }
})
