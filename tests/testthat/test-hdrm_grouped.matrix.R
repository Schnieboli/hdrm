data(birthrates)
Matrixbirthrates = t(as.matrix(birthrates))
# M <- matrix(rnorm(1200), 40, 30)
group_int <- c(1, 1, 2, 2, 1, 1, 1, 2, 1, 1, 1, 1, 2, 2, 1, 2)
group_fac <- factor(group_int, labels = c("west", "east"))

legal <- list()
legal$group <- list(int = group_int, fac = group_fac)
legal$hypotheses <- list("whole","sub","interaction","identical","flat",
                         "TW, TS, TR" = list(TW = diag(2), TS = diag(34), TR = diag(120)),
                         "TW, TS" = list(TW = diag(2), TS = diag(34)),
                         "TW, TS, TM" = list(TW = diag(2), TS = diag(34), TM = diag(120))
)
legal$subsampling <- legal$AM <- legal$cov.equal <- c(TRUE, FALSE)
legal$B <- list("10*N", 100)



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
                  Matrixbirthrates,
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
                  seed = initial$seed[[1]]
                )
              })
              expect_identical(initial, secondary) 
            }
  }
})


test_that("illegal 'data' input", {
  ## test missing values, non-finite values and non-numeric values
  mat_NA <- mat_Inf <- mat_char <- Matrixbirthrates
  mat_NA[1,1] <- NA
  mat_Inf[1,1] <- Inf
  mat_char <- matrix(as.character(mat_char), nrow = nrow(Matrixbirthrates))
  
  for(M in list(mat_NA, mat_Inf, mat_char))
    expect_error(
      hdrm_grouped(
        M,
        hypothesis = "sub",
        group = group_int,
        subsampling = FALSE,
        B = "10*N"
      ),
      "'data' must be a finite numeric matrix without missing values."
    )
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
        Matrixbirthrates,
        hypothesis = "sub",
        group = group,
        subsampling = FALSE,
        B = "10*N"
      ), "'group' must only contain finite non-missing values."
    )
  }
  
  for(group in list(group_short, group_long)){
    expect_error(
      hdrm_grouped(
        Matrixbirthrates,
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
  TS_idem <- diag(34)
  TS_idem[c(1,34),34] <- c(1,0)
  TS_symm <- 2*diag(34)
  illegal_hyp <- list("flart",                            ## illegal char
                      c("flat, sub"),                     ## two legal chars
                      1,                                  ## numeric value
                      diag(68),                           ## matrix
                      list(TW = diag(2), TS = diag(33)),  ## wrong dimension TS
                      list(TW = diag(3), TS = diag(34)),  ## wrong dimension TW
                      list(WT = diag(2),  TS = diag(34)), ## TW missing
                      list(TW = diag(2), ST = diag(34)),  ## TS missing
                      list(TW = TW_symm, TS = diag(34)),  ## TW not idempotent
                      list(TW = diag(2), TS = TS_symm),   ## TS not idempotent
                      list(TW = TW_idem , TS = diag(34)), ## TW not symmetrical
                      list(TW = diag(2), TS = TS_idem)    ## TW not symmetrical
  )
  
  for(h in illegal_hyp){
    expect_error(
      hdrm_grouped(
        Matrixbirthrates,
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
    Matrixbirthrates,
    NA,
    Inf,
    2^34,
    "2^32 * N"
    
  )
  
  for(B in illegal_B){
    expect_error(
      hdrm_grouped(
        Matrixbirthrates,
        hypothesis = "sub",
        group = group_int,
        subsampling = FALSE,
        B = B
      )
    )
  }
})


test_that("hdrm_grouped statistics, p-values, and f-degrees of freedom", {
  
  # Helper: Vergleicht statistic, p.value und f in einem Rutsch
  test_metrics <- function(res, expected_stat, expected_p, expected_f, tolerance = 1e-6) {
    expect_equal(res$statistic, expected_stat, tolerance = tolerance)
    expect_equal(res$p.value,   expected_p,    tolerance = tolerance)
    expect_equal(res$f,         expected_f,    tolerance = tolerance)
  }
  
  ## define default arguments for readability
  res_whole <- hdrm_grouped(
    data = Matrixbirthrates,
    hypothesis = "whole",
    group = group_int,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  
  res_sub <- hdrm_grouped(
    data = Matrixbirthrates,
    hypothesis = "sub",
    group = group_int,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  
  res_interaction <- hdrm_grouped(
    data = Matrixbirthrates,
    hypothesis = "interaction",
    group = group_int,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  
  res_identical <- hdrm_grouped(
    data = Matrixbirthrates,
    hypothesis = "identical",
    group = group_int,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  
  res_flat <- hdrm_grouped(
    data = Matrixbirthrates,
    hypothesis = "flat",
    group = group_int,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  
  # Subsampling = TRUE
  res_subsampling <- hdrm_grouped(
    Matrixbirthrates,
    hypothesis = "whole",
    group = group_int,
    subsampling = TRUE,
    B = "10*N",
    seed = 3141
  )

  # AM = 0
  res_am0 <- hdrm_grouped(
    Matrixbirthrates,
    hypothesis = "whole",
    group = group_int,
    subsampling = FALSE,
    AM = 0,
    B = "10*N",
    seed = 3141
  )
  
  # cov.equal = TRUE, AM = 0
  res_cov_equal <- hdrm_grouped(
    Matrixbirthrates,
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

test_that("AM representations agree for the one-group method", {
  result_am0 <- hdrm_grouped(
    Matrixbirthrates,
    hypothesis = "whole",
    group = group_int,
    subsampling = FALSE,
    AM = 0,
    B = "10*N",
    seed = 3141
  )
  result_am1 <- hdrm_grouped(
    Matrixbirthrates,
    hypothesis = "whole",
    group = group_int,
    subsampling = FALSE,
    AM = 0,
    B = "10*N",
    seed = 3141
  )
  expect_equal(result_am0$statistic, result_am1$statistic, tolerance = 1e-12)
  expect_equal(result_am0$p.value, result_am1$p.value, tolerance = 1e-12)
  expect_equal(result_am0$f, result_am1$f, tolerance = 1e-12)
})
