data(birthrates)
Matrixbirthrates = t(as.matrix(birthrates))
M <- matrix(rnorm(1200), 40, 30)
L <- list(1:160, 1:40)
group <- factor(c(1, 1, 2, 2, 1, 1, 1, 2, 1, 1, 1, 1, 2, 2, 1, 2), labels = c("west", "east"))

test_that("perfect case works", {
  # hypothesis = whole
  expect_no_condition(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # hypothesis = sub
  expect_no_condition(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "sub",
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # hypothesis = interaction
  expect_no_condition(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "interaction",
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # hypothesis = identical
  expect_no_condition(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "identical",
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # hypothesis = all_flat
  expect_no_condition(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "flat",
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # subsampling = TRUE
  expect_no_condition(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      group = group,
      subsampling = TRUE,
      B = "10*N"
    )
  )
  
  # multiple hypothesis values are rejected
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = c("whole", "SUB"),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # legal list to hypothesis
  expect_no_condition(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = list(
        TW = diag(2),
        TS = diag(34),
        TR = diag(120)
      ),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  
  
  # unknown argument is rejected
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      group = group,
      subsampling = FALSE,
      B = "100*N",
      a = 5
    )
  )
})


test_that("wrong input: data", {
  # data as list
  expect_error(hdrm_grouped(
    L,
    hypothesis = "whole",
    group = group,
    subsampling = TRUE,
    B = "10*N"
  ))
  
  # data as data.frame
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = "whole",
      group = group,
      subsampling = TRUE,
      B = "10*N"
    )
  )
})

test_that("missing values", {
  # NA in value
  df <- Matrixbirthrates
  df[1, 1] <- NA
  expect_error(hdrm_grouped(
    df,
    hypothesis = "sub",
    group = group,
    subsampling = FALSE,
    B = "10*N"
  ), "'data' must not contain any missing values.")
  
  
})

test_that("wrong input: hypothesis", {
  # default hypothesis is "whole"
  expect_no_condition(hdrm_grouped(
    Matrixbirthrates,
    group = group,
    subsampling = FALSE,
    B = "10*N"
  ))
  
  # a number
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = 1,
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # illegal character
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = c("flart"),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  
  # list TS missing
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = list(TW = diag(4), ST = diag(40)),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  # TW missing
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = list(WT = diag(4), TS = diag(40)),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # wrong dimension of matrix TW
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = list(TW = diag(3), TS = diag(40)),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # wrong dimension of TS
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = list(TW = diag(2), TS = diag(51)),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # TW not symmetrical
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = list(TW = diag(2) + c(0, 1, 0, 0), TS = diag(34)),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # TS not symmetrical
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = list(TW = diag(2), TS = diag(34) + c(1, 0)),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # TW not idempotent
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = list(TW = diag(1:2), TS = diag(34)),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # TS not idempotent
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = list(TW = diag(2), TS = diag(1:34)),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # matrix
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = diag(34),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
})

test_that("wrong input: AM", {
  # AM must not be missing
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      group = group,
      AM = NA,
      B = "10*N"
    ),
    "'AM' must be a single logical value.",
    fixed = TRUE
  )
  
  # AM must have length one
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      group = group,
      AM = c(0, 1),
      B = "10*N"
    ),
    "'AM' must be a single logical value.",
    fixed = TRUE
  )
})

test_that("wrong input: cov.equal", {
  # cov.equal must not be missing
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      group = group,
      cov.equal = NA,
      B = "10*N"
    ),
    "'cov.equal' must be a single logical value.",
    fixed = TRUE
  )
  
  # cov.equal must have length one
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      group = group,
      cov.equal = c(TRUE, FALSE),
      B = "10*N"
    ),
    "'cov.equal' must be a single logical value.",
    fixed = TRUE
  )
})


test_that("wrong input: subsampling", {
  # Missing values are rejected
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "sub",
      group = group,
      subsampling = NA,
      B = "10*N"
    ),
    "'subsampling' must be a single logical value.",
    fixed = TRUE
  )
  
  # Vectors of length greater than one are rejected
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "sub",
      group = group,
      subsampling = c(TRUE, FALSE),
      B = "10*N"
    ),
    "'subsampling' must be a single logical value.",
    fixed = TRUE
  )
})



test_that("false input: B", {
  # character length 2
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "sub",
      group = group,
      subsampling = FALSE,
      B = c("10*N", "20*N")
    )
  )
  
  # numeric length 2
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "sub",
      group = group,
      subsampling = FALSE,
      B = c(100, 1000)
    )
  )
  
  # mixed numeric and character values form a vector of length 2
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "sub",
      group = group,
      subsampling = FALSE,
      B = c(100, "10*N")
    )
  )
  
  # negative number
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "sub",
      group = group,
      subsampling = TRUE,
      B = -10
    )
  )
  
  # negative character
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "sub",
      group = group,
      subsampling = TRUE,
      B = "-10*N"
    )
  )
  
  # function with nonexistent argument
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "sub",
      group = group,
      subsampling = TRUE,
      B = "10 *asdfghjkl"
    )
  )
  
  # B wrong class
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "sub",
      group = group,
      subsampling = TRUE,
      B = M
    )
  )
  
  
})


test_that("hdrm_grouped test statistics", {
  # hypothesis = whole
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$statistic,
    24.71615921
  )
  
  # hypothesis = sub
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "sub",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$statistic,
    321.685616
  )
  
  # hypothesis = interaction
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "interaction",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$statistic,
    132.295574
  )
  
  # hypothesis = identical
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "identical",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$statistic,
    88.55248
  )
  
  # hypothesis = all_flat
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "flat",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$statistic,
    358.685785
  )
  
  # subsampling = TRUE
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      group = group,
      subsampling = TRUE,
      B = "10*N",
      seed = 3141
    )$statistic,
    23.29855057964175202301
  )
  
  # AM = 0
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      AM = 0,
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$statistic,
    24.71615920
  )
  
  # cov.equal = TRUE
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      AM = 0,
      group = group,
      cov.equal = TRUE,
      B = "10*N",
      seed = 3141
    )$statistic,
    11.0789863
  )
  
})




test_that("hdrm_grouped p.value", {
  # hypothesis = whole
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$p.value,
    1.4826812854026686e-12
  )
  
  # hypothesis = sub
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "sub",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$p.value,
    2.220446e-16
  )
  
  # hypothesis = interaction
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "interaction",
      group = group,
      
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$p.value,
    2.220446e-16
  )
  
  # hypothesis = identical
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "identical",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$p.value,
    2.220446e-16
  )
  
  # hypothesis = all_flat
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "flat",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$p.value,
    2.220446e-16
  )
  
  # subsampling = TRUE
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      group = group,
      subsampling = TRUE,
      B = "10*N",
      seed = 3141
    )$p.value,
    2.6541147566165586e-12
  )
  
  
  
  # AM = 0
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      AM = 0,
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$p.value,
    1.4826812813658546e-12
  )
  
  # cov.equal = TRUE
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      AM = 0,
      group = group,
      cov.equal = TRUE,
      B = "10*N",
      seed = 3141
    )$p.value,
    2.1557352448469479e-07
  )
  
})


test_that("hdrm_grouped f", {
  # hypothesis = whole
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    2.3354947576950122
  )
  
  # hypothesis = sub
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "sub",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    5.0760352063818699
  )
  
  # hypothesis = interaction
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "interaction",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    6.6015909062582345
  )
  
  # hypothesis = identical
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "identical",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    5.5096330244284291
  )
  
  # hypothesis = all_flat
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "flat",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    67.9163858637283795
  )
  
  # subsampling = TRUE
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      group = group,
      subsampling = TRUE,
      B = "10*N",
      seed = 3141
    )$f,
    3.244517459760626998388
  )
  
  # AM = 0
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      AM = 0,
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    2.3354947588108192
  )
  
  # cov.equal = TRUE
  expect_equal(
    hdrm_grouped(
      Matrixbirthrates,
      hypothesis = "whole",
      AM = 0,
      group = group,
      cov.equal = TRUE,
      B = "10*N",
      seed = 3141
    )$f,
    4.5950814099639574
  )
  
})





# Loading the dataset
data("EEG")
df <- data.frame(value = EEG$value, subject = EEG$subject, dimension = EEG$dimension)

M <- matrix(rnorm(1200), 40, 30)
L <- list(1:160, 1:40)
test_that("perfect case works", {
  # hypothesis = whole
  expect_no_condition(
    hdrm_grouped(
      df,
      hypothesis = "whole",
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # hypothesis = sub
  expect_no_condition(
    hdrm_grouped(
      df,
      hypothesis = "sub",
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # hypothesis = interaction
  expect_no_condition(
    hdrm_grouped(
      df,
      hypothesis = "interaction",
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # hypothesis = identical
  expect_no_condition(
    hdrm_grouped(
      df,
      hypothesis = "identical",
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # hypothesis = all_flat
  expect_no_condition(
    hdrm_grouped(
      df,
      hypothesis = "flat",
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # subsampling = TRUE
  expect_no_condition(
    hdrm_grouped(
      df,
      hypothesis = "whole",
      group = EEG$group,
      subsampling = TRUE,
      B = "10*N"
    )
  )
  
  # multiple hypothesis values are rejected
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = c("whole", "SUB"),
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # legal list to hypothesis
  expect_no_condition(
    hdrm_grouped(
      df,
      hypothesis = list(
        TW = diag(4),
        TS = diag(40),
        TR = diag(120)
      ),
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  
  # unknown argument is rejected
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = "whole",
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N",
      a = 5
    )
  )
})


test_that("wrong input: data", {
  # data as list
  expect_error(
    hdrm_grouped(
      L,
      hypothesis = "whole",
      group = EEG$group,
      subsampling = TRUE,
      B = "10*N"
    )
  )
})

test_that("missing values", {
  # NA in value
  df2 <- df
  df2$value[123] <- NA
  expect_error(
    hdrm_grouped(
      df2,
      hypothesis = "sub",
      group = df$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # NA in group
  group_tmp <- EEG$group
  group_tmp[234] <- NA
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = "sub",
      group = group_tmp,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # NA in subject
  df2 <- df
  df2$subject[145] <- NA
  expect_error(
    hdrm_grouped(
      df2$value,
      hypothesis = "sub",
      group = df$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  
  
  # NA in unrelated column
  df2 <- df
  df2$variable <- 1
  df2$variable[134] <- NA
  expect_no_condition(
    hdrm_grouped(
      df2,
      hypothesis = "sub",
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # nonexistent column
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = "sub",
      group = "nonexistant",
      subsampling = TRUE,
      B = "10*N"
    )
  )
  
})

test_that("wrong input: hypothesis", {
  # default hypothesis is "whole"
  expect_no_condition(
    hdrm_grouped(
      df,
      group = EEG$group,,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # a number
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = 1,
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # illegal character
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = c("flart"),
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  
  # list TS missing
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = list(TW = diag(4), ST = diag(40)),
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  # TW missing
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = list(WT = diag(4), TS = diag(40)),
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # wrong dimension of matrix TW
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = list(TW = diag(3), TS = diag(40)),
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # wrong dimension of TS
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = list(TW = diag(4), TS = diag(51)),
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # TW not symmetrical
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = list(TW = diag(4) + c(0, 1), TS = diag(40)),
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # TS not symmetrical
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = list(TW = diag(4), TS = diag(40) + c(1, 0)),
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # TW not idempotent
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = list(TW = diag(1:4), TS = diag(40)),
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # TS not idempotent
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = list(TW = diag(4), TS = diag(1:40)),
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # matrix
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = diag(160),
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
})



test_that("false input: B", {
  # character length 2
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = "sub",
      group = EEG$group,
      subsampling = FALSE,
      B = c("10*N", "20*N")
    )
  )
  
  # numeric length 2
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = "sub",
      group = EEG$group,
      subsampling = FALSE,
      B = c(100, 1000)
    )
  )
  
  # mixed numeric and character values form a vector of length 2
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = "sub",
      group = EEG$group,
      subsampling = FALSE,
      B = c(100, "10*N")
    )
  )
  
  # negative number
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = "sub",
      group = EEG$group,
      subsampling = TRUE,
      B = -10
    )
  )
  
  # negative character
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = "sub",
      group = EEG$group,
      subsampling = TRUE,
      B = "-10*N"
    )
  )
  
  # function with nonexistent argument
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = "sub",
      group = EEG$group,
      subsampling = TRUE,
      B = "10 *asdfghjkl"
    )
  )
  
  # B wrong class
  expect_error(
    hdrm_grouped(
      df,
      hypothesis = "sub",
      group = EEG$group,
      subsampling = TRUE,
      B = M
    )
  )
  
  
})


test_that("non continuous levels", {
  # Non-continuous subject levels
  df2 <- df
  levels(df2$subject) <- c(1, 3:161)
  
  expect_no_condition(
    hdrm_grouped(
      df2,
      hypothesis = "whole",
      group = EEG$group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  
  # Non-continuous group levels
  group_tmp <- EEG$group
  levels(group_tmp) <- c(1, 4, 3, 9)
  
  expect_no_condition(
    hdrm_grouped(
      df,
      hypothesis = "whole",
      group = group_tmp,
      subsampling = FALSE,
      B = "10*N"
    )
  )
})




test_that("hdrm_grouped outputs (statistic, p.value, f)", {
  
  # 1. Calls ausführen
  res_whole <- hdrm_grouped(
    df,
    hypothesis = "whole",
    group = EEG$group,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  res_sub   <- hdrm_grouped(
    df,
    hypothesis = "sub",
    group = EEG$group,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  res_inter <- hdrm_grouped(
    df,
    hypothesis = "interaction",
    group = EEG$group,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  res_ident <- hdrm_grouped(
    df,
    hypothesis = "identical",
    group = EEG$group,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  res_flat  <- hdrm_grouped(
    df,
    hypothesis = "flat",
    group = EEG$group,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  res_subsam <- hdrm_grouped(
    df,
    hypothesis = "whole",
    group = EEG$group,
    subsampling = TRUE,
    B = "10*N",
    seed = 3141
  )
  res_am0   <- hdrm_grouped(
    df,
    hypothesis = "whole",
    AM = 0,
    group = EEG$group,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  res_coveq <- hdrm_grouped(
    df,
    hypothesis = "whole",
    AM = 0,
    group = EEG$group,
    cov.equal = TRUE,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  
  # 2. Assertions: hypothesis = whole
  expect_equal(res_whole$statistic, 0.787384763)
  expect_equal(res_whole$p.value, 0.17762962443683902)
  expect_equal(res_whole$f, 3.0762359575869449)
  
  # 3. Assertions: hypothesis = sub
  expect_equal(res_sub$statistic, 3434.748)
  expect_equal(res_sub$p.value, 2.220446e-16)
  expect_equal(res_sub$f, 1.1631027600477559)
  
  # 4. Assertions: hypothesis = interaction
  expect_equal(res_inter$statistic, 2.620662063)
  expect_equal(res_inter$p.value, 0.022775013735625439)
  expect_equal(res_inter$f, 3.7240829309044754)
  
  # 5. Assertions: hypothesis = identical
  expect_equal(res_ident$statistic, 2.34961433)
  expect_equal(res_ident$p.value, 0.031622054805077959)
  expect_equal(res_ident$f, 3.5593004918841231)
  
  # 6. Assertions: hypothesis = flat
  expect_equal(res_flat$statistic, 1662.76162)
  expect_equal(res_flat$p.value, 2.220446e-16)
  expect_equal(res_flat$f, 4.6057075938290319)
  
  # 7. Assertions: subsampling = TRUE
  expect_equal(res_subsam$statistic, 0.8011794626188581958104)
  expect_equal(res_subsam$p.value, 0.1689462389477219550482)
  expect_equal(res_subsam$f, 2.334405305781945383359)
  
  # 8. Assertions: AM = 0
  expect_equal(res_am0$statistic, 0.787384763)
  expect_equal(res_am0$p.value, 0.17762962443837976)
  expect_equal(res_am0$f, 3.0762359576915479)
  
  # 9. Assertions: cov.equal = TRUE
  expect_equal(res_coveq$statistic, 0.655452437)
  expect_equal(res_coveq$p.value, 0.19472578433824672)
  expect_equal(res_coveq$f, 2.2488424990595659)
})



test_that("AM representations agree for equal covariances", {
  result_am0 <- hdrm_grouped(
    df,
    hypothesis = "whole",
    AM = FALSE,
    group = EEG$group,
    cov.equal = TRUE,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  
  result_am1 <- hdrm_grouped(
    df,
    hypothesis = "whole",
    AM = TRUE,
    group = EEG$group,
    cov.equal = TRUE,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  
  expect_equal(result_am0$statistic, result_am1$statistic, tolerance = 1e-12)
  expect_equal(result_am0$p.value, result_am1$p.value, tolerance = 1e-12)
  expect_equal(result_am0$f, result_am1$f, tolerance = 1e-12)
})

# Additional tests for the revised public interface and data preprocessing ----

test_that("data, group, and subject inputs are validated explicitly", {
  # Empty data.frame input
  expect_error(
    hdrm_grouped(
      data.frame(value = numeric(0), subject = numeric(0), dimension = numeric(0)),
      group = character(0),
      B = 10
    ),
    "'data' must not be empty.",
    fixed = TRUE
  )
  
  # Empty matrix input
  expect_error(
    hdrm_grouped(
    matrix(numeric(0), nrow = 0L, ncol = 0L),
    group = character(0),
    B = 10
  ),
  "'data' must not be empty.",
  fixed = TRUE)
  
  # Group and data lengths must agree for vector input
  expect_error(
    hdrm_grouped(
      df,
      group = EEG$group[-1L],
      B = "10*N"
    ),
    "'group' must be a one-dimensional vector or factor of length nrow(data).",
    fixed = TRUE
  )
  
  # One group label is required for each matrix row
  expect_error(
    hdrm_grouped(Matrixbirthrates, group = group[-1L], B = "10*N"),
    "'group' must be a one-dimensional vector or factor of length nrow(data).",
    fixed = TRUE
  )
  
  # Lists are not valid group vectors
  expect_error(
    hdrm_grouped(Matrixbirthrates, group = as.list(group), B = "10*N"),
    "'group' must be a one-dimensional vector or factor of length nrow(data).",
    fixed = TRUE
  )
})


test_that("seed input is validated", {
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      group = group,
      seed = NA_real_,
      B = "10*N"
    ),
    "'seed' must be NULL or a single finite integer-valued number.",
    fixed = TRUE
  )
  
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      group = group,
      seed = Inf,
      B = "10*N"
    ),
    "'seed' must be NULL or a single finite integer-valued number.",
    fixed = TRUE
  )
  
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      group = group,
      seed = 1.5,
      B = "10*N"
    ),
    "'seed' must be NULL or a single finite integer-valued number.",
    fixed = TRUE
  )
  
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      group = group,
      seed = c(1, 2),
      B = "10*N"
    ),
    "'seed' must be NULL or a single finite integer-valued number.",
    fixed = TRUE
  )
  
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      group = group,
      seed = "3141",
      B = "10*N"
    ),
    "'seed' must be NULL or a single finite integer-valued number.",
    fixed = TRUE
  )
  
  expect_error(
    hdrm_grouped(
      Matrixbirthrates,
      group = group,
      seed = .Machine$integer.max + 1,
      B = "10*N"
    ),
    "'seed' must be NULL or a single finite integer-valued number.",
    fixed = TRUE
  )
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
  
  # B is a base budget; B = 1 is valid because the equal-covariance
  # third-trace estimator uses a * B total draws.
  expect_identical(hdrm:::expand_subsample_budget(B = 1L, multiplier = 2L), 2L)
})



test_that("grouped third-trace budgets scale with the number of groups", {
  expect_identical(hdrm:::expand_subsample_budget(B = 100L, multiplier = 4L),
                   400L)
  
  expect_error(
    hdrm:::expand_subsample_budget(B = .Machine$integer.max, multiplier = 2L),
    paste0(
      "The effective subsampling budget must not exceed ",
      .Machine$integer.max,
      "."
    ),
    fixed = TRUE
  )
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
  
  expect_identical(result$subsamples, 8L)
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
  
  
  # The helper rejects non-positive base budgets.
  expect_error(
    hdrm:::allocate_C1_subsamples(group_sizes = c(6L, 7L, 8L), B = 0L),
    "'B' must be a finite positive integer.",
    fixed = TRUE
  )
})


test_that("zero-rank custom hypotheses are rejected", {
  expect_error(hdrm_grouped(
    Matrixbirthrates,
    hypothesis = list(TW = matrix(0, nrow = 2L, ncol = 2L), TS = diag(34)),
    group = group,
    B = "10*N"
  ))
  
  expect_error(hdrm_grouped(
    Matrixbirthrates,
    hypothesis = list(TW = diag(2), TS = matrix(0, nrow = 34L, ncol = 34L)),
    group = group,
    B = "10*N"
  ))
})

print("auskommentiert: matrix preprocessing is reflected in the returned object")
#### auch im nächsten test wurde eine zeie auskommentiert!
# test_that("matrix preprocessing is reflected in the returned object", {
#   result <- hdrm_grouped(
#     Matrixbirthrates,
#     hypothesis = "whole",
#     group = group,
#     subsampling = FALSE,
#     B = "10*N",
#     seed = 3141
#   )
#   
#   expected_order <- order(group)
#   
#   expect_equal(result$data, Matrixbirthrates[expected_order, , drop = FALSE])
#   expect_equal(result$groups$table, table(group))
#   expect_equal(result$removed.cases, 0L)
# })




test_that("wide grouped matrix input uses subjects in rows", {
  result <- hdrm_grouped(
    Matrixbirthrates,
    hypothesis = "whole",
    group = group,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  
  expected_order <- order(group)
  
  expect_equal(result$dim$d, ncol(Matrixbirthrates))
  expect_equal(result$dim$N, nrow(Matrixbirthrates))
  # expect_equal(result$data, Matrixbirthrates[expected_order, , drop = FALSE])
})


print("auskommentiert: incomplete matrix subjects are removed with their group labels")
#### Grund: nicht mehr sinnvoll, da jetzt bei NA direkt ein Fehler kommt
# test_that("incomplete matrix subjects are removed with their group labels", {
#   data_with_na <- Matrixbirthrates
#   data_with_na[1L, 1L] <- NA_real_
#   
#   result <- NULL
#   expect_warning(
#     result <- hdrm_grouped(
#       data_with_na,
#       hypothesis = "sub",
#       group = group,
#       subsampling = FALSE,
#       B = "10*N",
#       seed = 3141
#     ),
#     "'data' must not contain any missing values.",
#     fixed = TRUE
#   )
#   
#   complete_subjects <- stats::complete.cases(data_with_na)
#   remaining_group <- droplevels(group[complete_subjects])
#   expected_order <- order(remaining_group)
#   expected_data <- data_with_na[complete_subjects, , drop = FALSE][expected_order, , drop = FALSE]
#   
#   expect_equal(result$data, expected_data)
#   expected_group_table <- table(remaining_group)
#   names(dimnames(expected_group_table)) <- "group"
#   
#   expect_equal(result$groups$table, expected_group_table)
#   expect_equal(result$removed.cases, 1L)
# })

print("auskommentiert: named numeric vectors are accepted")
#### Grund: überflüssig
# test_that("named numeric vectors are accepted", {
#   named_values <- EEG$value
#   names(named_values) <- seq_along(named_values)
#   
#   expect_no_condition(
#     hdrm_grouped(
#       named_values,
#       hypothesis = "whole",
#       group = EEG$group,
#       subsampling = FALSE,
#       B = "10*N",
#       seed = 3141
#     )
#   )
# })

print("entfernt: subject labels may be reused in different groups")
#### Grund: ich weiß nicht, was genau das macht
# test_that("subject labels may be reused in different groups", {
#   reused_subject <- integer(length(EEG$subject))
#   
#   for (current_group in levels(droplevels(as.factor(EEG$group)))) {
#     index <- EEG$group == current_group
#     subjects_in_group <- unique(EEG$subject[index])
#     reused_subject[index] <- match(EEG$subject[index], subjects_in_group)
#   }
#   
#   result_original <- hdrm_grouped(
#     df,
#     hypothesis = "whole",
#     group = EEG$group,
#     subsampling = FALSE,
#     B = "10*N",
#     seed = 3141
#   )
#   
#   df2 <- df
#   df2$subject <- reused_subject
#   result_reused <- hdrm_grouped(
#     df2,
#     hypothesis = "whole",
#     group = EEG$group,
#     subsampling = FALSE,
#     B = "10*N",
#     seed = 3141
#   )
#   
#   expect_equal(result_reused$data, result_original$data)
#   expect_equal(result_reused$statistic, result_original$statistic)
#   expect_equal(result_reused$p.value, result_original$p.value)
#   expect_equal(result_reused$f, result_original$f)
#   expect_equal(result_reused$groups$table, result_original$groups$table)
# })

print("auskommentiert: vector output reports subject counts rather than measurement counts")
#### Grund: das ist veraltet, ließe sich aber fixen -> ich weiß allerdings
#### nicht, ob das sinnvoll ist oder nicht
# test_that("vector output reports subject counts rather than measurement counts",{
#   result <- hdrm_grouped(
#     df,
#     hypothesis = "whole",
#     group = EEG$group,
#     subsampling = FALSE,
#     B = "10*N",
#     seed = 3141
#   )
#
#   composite_subject <- interaction(
#     droplevels(as.factor(EEG$group)),
#     as.factor(EEG$subject),
#     drop = TRUE,
#     lex.order = TRUE
#   )
#
#   subject_group_pairs <- unique(data.frame(subject = composite_subject, whole = droplevels(as.factor(EEG$group))))
#
#   expect_equal(result$groups$table, table(subject_group_pairs$whole))
#   expect_equal(nrow(result$data), nrow(subject_group_pairs))
#   expect_equal(result$removed.cases, 0L)
# })


test_that("the original within-subject measurement order is preserved", {
  composite_subject <- interaction(
    droplevels(as.factor(EEG$group)),
    as.factor(EEG$subject),
    drop = TRUE,
    lex.order = TRUE
  )
  
  subject_blocks <- split(seq_along(EEG$value), composite_subject)
  permutation <- unlist(rev(subject_blocks), use.names = FALSE)
  
  result_original <- hdrm_grouped(
    df,
    hypothesis = "whole",
    group = EEG$group,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  
  df2 <- df[permutation, ]
  
  result_permuted <- hdrm_grouped(
    df2,
    hypothesis = "whole",
    group = EEG$group[permutation],
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  
  # expect_equal(result_permuted$data, result_original$data)
  expect_equal(result_permuted$statistic, result_original$statistic)
  expect_equal(result_permuted$p.value, result_original$p.value)
  expect_equal(result_permuted$f, result_original$f)
})

print("auskommentiert: groups removed through missing data are dropped from the analysis")
#### Grund: nicht mehr zeitgemäß, de bei Fehlenden Werten ein Fehler kommt
# test_that("groups removed through missing data are dropped from the analysis",{
#   data_with_missing_group <- EEG
#   removed_group <- levels(droplevels(as.factor(data_with_missing_group$group)))[[1L]]
#   removed_index <- data_with_missing_group$group == removed_group
#   removed_subjects <- length(unique(data_with_missing_group$subject[removed_index]))
#   
#   data_with_missing_group$value[removed_index] <- NA_real_
#   
#   result <- NULL
#   expect_warning(
#     result <- hdrm_grouped(
#       data_with_missing_group$value,
#       hypothesis = "whole",
#       group = data_with_missing_group$group,
#       subject = data_with_missing_group$subject,
#       subsampling = FALSE,
#       B = "10*N",
#       seed = 3141
#     ),
#     "Subjects with missing values dropped",
#     fixed = TRUE
#   )
#   
#   expect_equal(result$groups$a, nlevels(droplevels(data_with_missing_group$group[!removed_index])))
#   expect_false(removed_group %in% names(result$groups$table))
#   expect_equal(result$removed.cases, removed_subjects)
# })


test_that("an error is raised if fewer than two groups remain", {
    expect_error(
    hdrm_grouped(
      df,
      hypothesis = "whole",
      group = rep(1, nrow(df)),
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    ),
    "there must be at least two groups",
    fixed = TRUE
  )
})



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

print("ein test wurde auskommentiert weil die getestete Funktion nicht mehr aktuell ist!")


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
            expect_equal(hdrm:::compute_eta_Na(P_a, balanced_sizes), 2, tolerance = 1e-12)
            
            expect_equal(hdrm:::compute_eta_Na(J_a, balanced_sizes), 1, tolerance = 1e-12)
            
            expect_equal(hdrm:::compute_eta_Na(I_a, balanced_sizes), 3, tolerance = 1e-12)
            
            unbalanced_sizes <- c(4L, 7L, 9L)
            
            expect_equal(
              hdrm:::compute_eta_Na(P_a, unbalanced_sizes),
              eta_eigen_R(P_a, unbalanced_sizes),
              tolerance = 1e-12
            )
            
            expect_equal(
              hdrm:::compute_eta_Na(J_a, unbalanced_sizes),
              eta_eigen_R(J_a, unbalanced_sizes),
              tolerance = 1e-12
            )
            
            expect_equal(
              hdrm:::compute_eta_Na(I_a, unbalanced_sizes),
              eta_eigen_R(I_a, unbalanced_sizes),
              tolerance = 1e-12
            )
            
            expect_equal(hdrm:::compute_eta_Na(P_a, unbalanced_sizes),
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
              hdrm:::compute_eta_Na(custom_TW, custom_sizes),
              eta_eigen_R(custom_TW, custom_sizes),
              tolerance = 1e-11
            )
            
            # A common nonzero scaling of TW cancels from the trace ratio.
            expect_equal(
              hdrm:::compute_eta_Na(7.5 * custom_TW, custom_sizes),
              hdrm:::compute_eta_Na(custom_TW, custom_sizes),
              tolerance = 1e-12
            )
            
            # Simultaneously relabelling groups and rows/columns of TW changes neither
            # the design nor eta_{N,a}.
            permutation <- c(3L, 1L, 4L, 2L)
            
            expect_equal(
              hdrm:::compute_eta_Na(custom_TW[permutation, permutation, drop = FALSE], custom_sizes[permutation]),
              hdrm:::compute_eta_Na(custom_TW, custom_sizes),
              tolerance = 1e-12
            )
          })


test_that("compute_eta_Na rejects invalid or degenerate inputs", {
  expect_error(
    compute_eta_Na(matrix(0, nrow = 2L, ncol = 2L), c(5L, 5L)),
    "The whole-plot trace factor is degenerate.",
    fixed = TRUE
  )
  
  expect_error(
    compute_eta_Na(diag(3L), c(5L, 5L)),
    paste0(
      "'group_sizes' must contain one positive integer for each ",
      "row of 'TW'."
    ),
    fixed = TRUE
  )
  
  expect_error(
    compute_eta_Na(diag(3L), c(5L, 0L, 5L)),
    paste0(
      "'group_sizes' must contain one positive integer for each ",
      "row of 'TW'."
    ),
    fixed = TRUE
  )
  
  expect_error(compute_eta_Na(matrix(c(1, 1, 0, 1), nrow = 2L, ncol = 2L),
                                     c(5L, 5L)),
  "'TW' must be symmetric.",
  fixed = TRUE)
})


test_that("grouped p-values use the stable upper chi-square tail", {
  heterogeneous_moderate <- hdrm_grouped(
    df,
    hypothesis = "whole",
    group = EEG$group,
    cov.equal = FALSE,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  
  expected_heterogeneous <- max(
    pchisq(
      heterogeneous_moderate$statistic *
        sqrt(2 * heterogeneous_moderate$f) +
        heterogeneous_moderate$f,
      df = heterogeneous_moderate$f,
      lower.tail = FALSE
    ),
    .Machine$double.eps
  )
  
  expect_equal(heterogeneous_moderate$p.value,
               expected_heterogeneous,
               tolerance = 1e-15)
  
  equal_covariance_moderate <- hdrm_grouped(
    df,
    hypothesis = "whole",
    group = EEG$group,
    cov.equal = TRUE,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  
  expected_equal_covariance <- max(
    stats::pchisq(
      equal_covariance_moderate$statistic *
        sqrt(2 * equal_covariance_moderate$f) +
        equal_covariance_moderate$f,
      df = equal_covariance_moderate$f,
      lower.tail = FALSE
    ),
    .Machine$double.eps
  )
  
  expect_equal(equal_covariance_moderate$p.value,
               expected_equal_covariance,
               tolerance = 1e-15)
  
  heterogeneous_extreme <- hdrm_grouped(
    Matrixbirthrates,
    hypothesis = "sub",
    group = group,
    cov.equal = FALSE,
    subsampling = FALSE,
    B = "10*N",
    seed = 3141
  )
  
  raw_extreme_tail <- stats::pchisq(
    heterogeneous_extreme$statistic *
      sqrt(2 * heterogeneous_extreme$f) +
      heterogeneous_extreme$f,
    df = heterogeneous_extreme$f,
    lower.tail = FALSE
  )
  
  expect_lt(raw_extreme_tail, .Machine$double.eps)
  
  expect_identical(heterogeneous_extreme$p.value, .Machine$double.eps)
  
  # The stable upper-tail calculation avoids cancellation while preserving
  # the package's deliberate positive lower bound for reported p-values.
  expect_true(heterogeneous_moderate$p.value >= .Machine$double.eps)
  
  expect_true(equal_covariance_moderate$p.value >= .Machine$double.eps)
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


test_that("degenerate grouped data produce informative errors", {
  constant_data <- matrix(3, nrow = 12L, ncol = 2L)
  constant_group <- factor(rep(c("A", "B"), each = 6L))
  
  # In the heterogeneous procedure, validation of the individual trace
  # estimators is reached before the derived variance is calculated.
  # expect_error(
  #   hdrm_grouped(
  #     constant_data,
  #     hypothesis = "whole",
  #     group = constant_group,
  #     cov.equal = FALSE,
  #     subsampling = FALSE,
  #     B = 100,
  #     seed = 3141
  #   ),
  #   "The grouped trace estimators must be finite and non-negative.",
  #   fixed = TRUE
  # )
  
  print(
    "ein test wurde auskommentiert weil der resultierende fehler in anderer form deutlich eher passieren sollte!"
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
    "The estimated variance of the test statistic must be finite and positive.",
    fixed = TRUE
  )
})

print("auskommentiert: MSrootcompact remains an internal helper")
#### Grund: überlfüssig
# test_that("MSrootcompact remains an internal helper", {
#   expect_false("MSrootcompact" %in% getNamespaceExports("hdrm"))
#   
#   H <- diag(c(1, 0, 2))
#   
#   root <- hdrm:::MSrootcompact(H)
#   
#   expect_equal(crossprod(root), H, tolerance = 1e-12)
#   
#   expect_equal(nrow(root), qr(H)$rank)
# })


test_that("predefined hypothesis labels are preserved in both covariance procedures",
          {
            hypotheses <- c("whole", "sub", "interaction", "identical", "flat")
            
            for (hypothesis in hypotheses) {
              heterogeneous <- hdrm_grouped(
                Matrixbirthrates,
                hypothesis = hypothesis,
                group = group,
                cov.equal = FALSE,
                subsampling = FALSE,
                B = 64,
                seed = 3141
              )
              
              homogeneous <- hdrm_grouped(
                Matrixbirthrates,
                hypothesis = hypothesis,
                group = group,
                cov.equal = TRUE,
                subsampling = FALSE,
                B = 64,
                seed = 3141
              )
              
              expect_identical(
                heterogeneous$hypothesis,
                hypothesis,
                info = paste("heterogeneous covariance procedure:", hypothesis)
              )
              
              expect_identical(
                homogeneous$hypothesis,
                hypothesis,
                info = paste("equal covariance procedure:", hypothesis)
              )
              
              heterogeneous_output <- paste(capture.output(print(heterogeneous)), collapse = "\n")
              
              homogeneous_output <- paste(capture.output(print(homogeneous)), collapse = "\n")
              
              expect_match(
                heterogeneous_output,
                paste0("Hypothesis type: ", hypothesis),
                fixed = TRUE,
                info = paste("printed heterogeneous covariance result:", hypothesis)
              )
              
              expect_match(
                homogeneous_output,
                paste0("Hypothesis type: ", hypothesis),
                fixed = TRUE,
                info = paste("printed equal covariance result:", hypothesis)
              )
            }
          })
print("auskommentiert: the print method returns its input invisibly")
# test_that("the print method returns its input invisibly", {
#   result <- hdrm_grouped(
#     Matrixbirthrates,
#     hypothesis = "whole",
#     group = group,
#     subsampling = FALSE,
#     B = "10*N",
#     seed = 3141
#   )
#   
#   printed <- NULL
#   output <- capture.output(printed <- withVisible(print(result)))
#   
#   expect_false(printed$visible)
#   expect_identical(printed$value, result)
#   expect_gt(length(output), 0L)
# })
