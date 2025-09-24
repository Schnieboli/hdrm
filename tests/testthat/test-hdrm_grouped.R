data(birthrates)
M <- matrix(rnorm(1200), 40, 30)
L <- list(1:160, 1:40)
group <- factor(c(1,1,2,2,1,1,1,2,1,1,1,1,2,2,1,2), labels = c("west","east"))

test_that("perfect case works",{
  # hypothesis = whole
  expect_no_condition(
    hdrm_grouped(
      birthrates,
      hypothesis = "whole",
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # hypothesis = sub
  expect_no_condition(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # hypothesis = interaction
  expect_no_condition(
    hdrm_grouped(
      birthrates,
      hypothesis = "interaction",
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # hypothesis = identical
  expect_no_condition(
    hdrm_grouped(
      birthrates,
      hypothesis = "identical",
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # hypothesis = all_flat
  expect_no_condition(
    hdrm_grouped(
      birthrates,
      hypothesis = "flat",
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # subsampling = TRUE
  expect_no_condition(
    hdrm_grouped(
      birthrates,
      hypothesis = "whole",
      group = group,
      subsampling = TRUE,
      B = "10*N"
    )
  )

  # first element legal, second illegal
  expect_warning(
    expect_no_error(
      hdrm_grouped(
        birthrates,
        hypothesis = c("whole", "SUB"),
        group = group,
        subsampling = FALSE,
        B = "10*N"
      )
    )
  )

  # legal list to hypothesis
  expect_no_condition(
    hdrm_grouped(
      birthrates,
      hypothesis = list(TW = diag(2), TS = diag(34), TR = diag(120)),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )



  # argument to ...
  expect_no_condition(
    hdrm_grouped(
      birthrates,
      hypothesis = "whole",
      group = group,
      subsampling = FALSE,
      B = "100*N",
      a = 5
    )
  )
})


test_that("wrong input: data",{
  # data as list
  expect_error(
    hdrm_grouped(
      L,
      hypothesis = "whole",
      group = group,
      subsampling = TRUE,
      B = "10*N"
    )
  )

  # data as df
  expect_error(
    hdrm_grouped(
      Birthrates,
      hypothesis = "whole",
      group = group,
      subsampling = TRUE,
      B = "10*N"
    )
  )

})

test_that("missing values",{
  # NA in value
  df <- birthrates
  df[4,1] <- NA
  expect_warning(
    hdrm_grouped(
      df,
      hypothesis = "sub",
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )


})

test_that("wrong input: hypothesis",{

  # multiple arguments to hypothesis (default)
  expect_warning(
    hdrm_grouped(
      birthrates,
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # a number
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = 1,
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # illegal character
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = c("flart"),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )


  # list TS missing
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = list(TW = diag(4), ST = diag(40)),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  # TW missing
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = list(WT = diag(4), TS = diag(40)),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # wrong dimension of matrix TW
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = list(TW = diag(3), TS = diag(40)),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # wrong dimension of TS
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = list(TW = diag(2), TS = diag(51)),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # TW not symmetrical
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = list(TW = diag(2) + c(0,1,0,0), TS = diag(34)),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # TS not symmetrical
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = list(TW = diag(2), TS = diag(34) + c(1,0)),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # TW not idempotent
  expect_warning(
    hdrm_grouped(
      birthrates,
      hypothesis = list(TW = diag(1:2), TS = diag(34)),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # TS not idempotent
  expect_warning(
    hdrm_grouped(
      birthrates,
      hypothesis = list(TW = diag(2), TS = diag(1:34)),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # matrix
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = diag(34),
      group = group,
      subsampling = FALSE,
      B = "10*N"
    )
  )

})



test_that("false input: subsampling",{
  # works
  expect_no_condition(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = 1,
      B = "10*N"
    )
  )

  # works
  expect_no_condition(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = 0,
      B = "10*N"
    )
  )

  # logical as string
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = "TRUE",
      B = "10*N"
    )
  )

  # string
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = "notlogical",
      B = "10*N"
    )
  )

  # error because "maybe" makes it all a character;
  # warning because length 2
  expect_warning(
    expect_error(
      hdrm_grouped(
        birthrates,
        hypothesis = "sub",
        group = group,
        subsampling = c(TRUE, "maybe"),
        B = "10*N"
      )
    )
  )

  # works, but warning because length 2
  expect_warning(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = c(TRUE, FALSE),
      B = "10*N"
    )
  )

  # input logical matrix: works, because it can be coerced to vector
  # warning because length >1
  expect_warning(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = diag(c(TRUE, FALSE)),
      B = "10*N"
    )
  )

  # non logical matrix: error and warning, because "only first element used" comes first
  expect_warning(
    expect_error(
      hdrm_grouped(
        birthrates,
        hypothesis = "sub",
        group = group,
        subsampling = M,
        B = "10*N"
      )
    )
  )

  # input list: warning and error, because list is length 2
  expect_warning(
    expect_error(
      hdrm_grouped(
        birthrates,
        hypothesis = "sub",
        group = group,
        subsampling = L,
        B = "10*N"
      )
    )
  )


  # input data.frame: warning and error, because birthrates[1] is possible
  expect_warning(
    expect_error(
      hdrm_grouped(
        birthrates,
        hypothesis = "sub",
        group = group,
        subsampling = birthrates,
        B = "10*N"
      )
    )
  )
})


test_that("false input: B",{
  # character length 2
  expect_warning(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = FALSE,
      B = c("10*N", "20*N")
    )
  )

  # numeric length 2
  expect_warning(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = FALSE,
      B = c(100, 1000)
    )
  )

  # c(100, "10*N"): "100" works!
  expect_warning(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = FALSE,
      B = c(100, "10*N")
    )
  )

  # negative number
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = TRUE,
      B = -10
    )
  )

  # negative character
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = TRUE,
      B = "-10*N"
    )
  )

  # function with nonexistant argument
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = TRUE,
      B = "10 *asdfghjkl"
    )
  )

  # B wrong class
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = TRUE,
      B = M
    )
  )


})

test_that("false input: B",{
  # character length 2
  expect_warning(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = FALSE,
      B = c("10*N", "20*N")
    )
  )

  # numeric length 2
  expect_warning(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = FALSE,
      B = c(100, 1000)
    )
  )

  # c(100, "10*N"): "100" works!
  expect_warning(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = FALSE,
      B = c(100, "10*N")
    )
  )

  # negative number
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = TRUE,
      B = -10
    )
  )

  # negative character
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = TRUE,
      B = "-10*N"
    )
  )

  # function with nonexistant argument
  expect_error(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = TRUE,
      B = "10 *asdfghjkl"
    )
  )

  # B wrong class
  expect_error(
    hdrm_grouped(
      birthrates,
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
      birthrates,
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
      birthrates,
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
      birthrates,
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
      birthrates,
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
      birthrates,
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
      birthrates,
      hypothesis = "whole",
      group = group,
      subsampling = TRUE,
      B = "10*N",
      seed = 3141
    )$statistic,
    24.9872531
  )

  # AM = 0
  expect_equal(
    hdrm_grouped(
      birthrates,
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
      birthrates,
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
      birthrates,
      hypothesis = "whole",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$p.value,
    6.883383e-15
  )

  # hypothesis = sub
  expect_equal(
    hdrm_grouped(
      birthrates,
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
      birthrates,
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
      birthrates,
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
      birthrates,
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
      birthrates,
      hypothesis = "whole",
      group = group,
      subsampling = TRUE,
      B = "10*N",
      seed = 3141
    )$p.value,
    1.836742e-10
  )



  # AM = 0
  expect_equal(
    hdrm_grouped(
      birthrates,
      hypothesis = "whole",
      AM = 0,
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$p.value,
    6.883383e-15
  )

  # cov.equal = TRUE
  expect_equal(
    hdrm_grouped(
      birthrates,
      hypothesis = "whole",
      AM = 0,
      group = group,
      cov.equal = TRUE,
      B = "10*N",
      seed = 3141
    )$p.value,
    4.268119e-10
  )

})


test_that("hdrm_grouped f", {
  # hypothesis = whole
  expect_equal(
    hdrm_grouped(
      birthrates,
      hypothesis = "whole",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    3.78031367
  )

  # hypothesis = sub
  expect_equal(
    hdrm_grouped(
      birthrates,
      hypothesis = "sub",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    6.62627370
  )

  # hypothesis = interaction
  expect_equal(
    hdrm_grouped(
      birthrates,
      hypothesis = "interaction",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    6.09873285
  )

  # hypothesis = identical
  expect_equal(
    hdrm_grouped(
      birthrates,
      hypothesis = "identical",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    6.77215918
  )

  # hypothesis = all_flat
  expect_equal(
    hdrm_grouped(
      birthrates,
      hypothesis = "flat",
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    387.557345
  )

  # subsampling = TRUE
  expect_equal(
    hdrm_grouped(
      birthrates,
      hypothesis = "whole",
      group = group,
      subsampling = TRUE,
      B = "10*N",
      seed = 3141
    )$f,
    1.68699687
  )

  # AM = 0
  expect_equal(
    hdrm_grouped(
      birthrates,
      hypothesis = "whole",
      AM = 0,
      group = group,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    3.78031367
  )

  # cov.equal = TRUE
  expect_equal(
    hdrm_grouped(
      birthrates,
      hypothesis = "whole",
      AM = 0,
      group = group,
      cov.equal = TRUE,
      B = "10*N",
      seed = 3141
    )$f,
    16.8078464
  )

})





# Loading the dataset
data("EEG")

M <- matrix(rnorm(1200), 40, 30)
L <- list(1:160, 1:40)
test_that("perfect case works",{
  # hypothesis = whole
  expect_no_condition(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # hypothesis = sub
  expect_no_condition(
    hdrm_grouped(
      EEG$value,
      hypothesis = "sub",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # hypothesis = interaction
  expect_no_condition(
    hdrm_grouped(
      EEG$value,
      hypothesis = "interaction",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # hypothesis = identical
  expect_no_condition(
    hdrm_grouped(
      EEG$value,
      hypothesis = "identical",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # hypothesis = all_flat
  expect_no_condition(
    hdrm_grouped(
      EEG$value,
      hypothesis = "flat",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # subsampling = TRUE
  expect_no_condition(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = TRUE,
      B = "10*N"
    )
  )

  # first elemet legal, second illegal
  expect_warning(
    expect_no_error(
      hdrm_grouped(
        EEG$value,
        hypothesis = c("whole", "SUB"),
        group = EEG$group,
        subject = EEG$subject,
        subsampling = FALSE,
        B = "10*N"
      )
    )
  )

  # legal list to hypothesis
  expect_no_condition(
    hdrm_grouped(
      EEG$value,
      hypothesis = list(TW = diag(4), TS = diag(40), TR = diag(120)),
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )


  # argument to ...
  expect_no_condition(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "100*N",
      a = 5
    )
  )
})


test_that("wrong input: data",{
  # data as list
  expect_error(
    hdrm_grouped(
      L,
      hypothesis = "whole",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = TRUE,
      B = "10*N"
    )
  )

  # data as matrix
  expect_error(
    hdrm_grouped(
      M,
      hypothesis = "whole",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = TRUE,
      B = "10*N"
    )
  )

})

test_that("missing values",{
  # NA in value
  df <- EEG
  df$value[123] <- NA
  expect_warning(
    hdrm_grouped(
      df$value,
      hypothesis = "sub",
      group = df$group,
      subject = df$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # NA in group
  df <- EEG
  df$group[234] <- NA
  expect_error(
    hdrm_grouped(
      df$value,
      hypothesis = "sub",
      group = df$group,
      subject = df$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # NA in subject
  df <- EEG
  df$subject[145] <- NA
  expect_error(
    hdrm_grouped(
      df$value,
      hypothesis = "sub",
      group = df$group,
      subject = df$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )



  # NA in unrelated column
  df <- EEG
  df$variable[134] <- NA
  expect_no_condition(
    hdrm_grouped(
      df$value,
      hypothesis = "sub",
      group = df$group,
      subject = df$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # nonexistant column
  expect_error(
    hdrm_grouped(
      EEG$value,
      hypothesis = "sub",
      group = "nonexistant",
      subject = df$subject,
      subsampling = TRUE,
      B = "10*N"
    )
  )

})

test_that("wrong input: hypothesis",{

  # multiple arguments to hypothesis (default)
  expect_warning(
    hdrm_grouped(
      EEG$value,
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # a number
  expect_error(
    hdrm_grouped(
      EEG$value,
      hypothesis = 1,
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # illegal character
  expect_error(
    hdrm_grouped(
      EEG$value,
      hypothesis = c("flart"),
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )


  # list TS missing
  expect_error(
    hdrm_grouped(
      EEG$value,
      hypothesis = list(TW = diag(4), ST = diag(40)),
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )
  # TW missing
  expect_error(
    hdrm_grouped(
      EEG$value,
      hypothesis = list(WT = diag(4), TS = diag(40)),
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # wrong dimension of matrix TW
  expect_error(
    hdrm_grouped(
      EEG$value,
      hypothesis = list(TW = diag(3), TS = diag(40)),
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # wrong dimension of TS
  expect_error(
    hdrm_grouped(
      EEG$value,
      hypothesis = list(TW = diag(4), TS = diag(51)),
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # TW not symmetrical
  expect_error(
    hdrm_grouped(
      EEG$value,
      hypothesis = list(TW = diag(4) + c(0, 1), TS = diag(40)),
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # TS not symmetrical
  expect_error(
    hdrm_grouped(
      EEG$value,
      hypothesis = list(TW = diag(4), TS = diag(40) + c(1,0)),
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # TW not idempotent
  expect_warning(
    hdrm_grouped(
      EEG$value,
      hypothesis = list(TW = diag(1:4), TS = diag(40)),
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # TS not idempotent
  expect_warning(
    hdrm_grouped(
      EEG$value,
      hypothesis = list(TW = diag(4), TS = diag(1:40)),
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # matrix
  expect_error(
    hdrm_grouped(
      EEG$value,
      hypothesis = diag(160),
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

})



test_that("false input: subsampling",{
  # works
  expect_no_condition(
    hdrm_grouped(
      EEG$value,
      hypothesis = "sub",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = 1,
      B = "10*N"
    )
  )

  # works
  expect_no_condition(
    hdrm_grouped(
      EEG$value,
      hypothesis = "sub",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = 0,
      B = "10*N"
    )
  )

  # logical as string
  expect_error(
    hdrm_grouped(
      EEG$value,
      hypothesis = "sub",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = "TRUE",
      B = "10*N"
    )
  )

  # string
  expect_error(
    hdrm_grouped(
      EEG$value,
      hypothesis = "sub",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = "notlogical",
      B = "10*N"
    )
  )

  # error because "maybe" makes it all a character;
  # warning because length 2
  expect_warning(
    expect_error(
      hdrm_grouped(
        EEG$value,
        hypothesis = "sub",
        group = EEG$group,
        subject = EEG$subject,
        subsampling = c(TRUE, "maybe"),
        B = "10*N"
      )
    )
  )

  # works, but warning because length 2
  expect_warning(
    hdrm_grouped(
      EEG$value,
      hypothesis = "sub",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = c(TRUE, FALSE),
      B = "10*N"
    )
  )

  # input logical matrix: works, because it can be coerced to vector
  # warning because length >1
  expect_warning(
    hdrm_grouped(
      EEG$value,
      hypothesis = "sub",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = diag(c(TRUE, FALSE)),
      B = "10*N"
    )
  )

  # non logical matrix: error and warning, because "only first elemt used" comes first
  expect_warning(
    expect_error(
      hdrm_grouped(
        EEG$value,
        hypothesis = "sub",
        group = EEG$group,
        subject = EEG$subject,
        subsampling = M,
        B = "10*N"
      )
    )
  )

  # input list: warning and error, because list is length 2
  expect_warning(
    expect_error(
      hdrm_grouped(
        EEG$value,
        hypothesis = "sub",
        group = EEG$group,
        subject = EEG$subject,
        subsampling = L,
        B = "10*N"
      )
    )
  )


  # input data.frame: warning and error, because EEG[1] is possible
  expect_warning(
    expect_error(
      hdrm_grouped(
        EEG$value,
        hypothesis = "sub",
        group = EEG$group,
        subject = EEG$subject,
        subsampling = EEG$value,
        B = "10*N"
      )
    )
  )
})


test_that("false input: B",{
  # character length 2
  expect_warning(
    hdrm_grouped(
      EEG$value,
      hypothesis = "sub",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = c("10*N", "20*N")
    )
  )

  # numeric length 2
  expect_warning(
    hdrm_grouped(
      EEG$value,
      hypothesis = "sub",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = c(100, 1000)
    )
  )

  # c(100, "10*N"): "100" works!
  expect_warning(
    hdrm_grouped(
      EEG$value,
      hypothesis = "sub",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = c(100, "10*N")
    )
  )

  # negative number
  expect_error(
    hdrm_grouped(
      EEG$value,
      hypothesis = "sub",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = TRUE,
      B = -10
    )
  )

  # negative character
  expect_error(
    hdrm_grouped(
      EEG$value,
      hypothesis = "sub",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = TRUE,
      B = "-10*N"
    )
  )

  # function with nonexistant argument
  expect_error(
    hdrm_grouped(
      EEG$value,
      hypothesis = "sub",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = TRUE,
      B = "10 *asdfghjkl"
    )
  )

  # B wrong class
  expect_error(
    hdrm_grouped(
      EEG$value,
      hypothesis = "sub",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = TRUE,
      B = M
    )
  )


})


test_that("non continuous levels",{
  df <- EEG
  levels(df$subject) <- c(1, 3:161)

  expect_no_condition(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  df <- EEG
  levels(df$dimension) <- c(1, 3:50)

  expect_no_condition(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  df <- EEG
  levels(df$group) <- c(1,4,3,9)

  expect_no_condition(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N"
    )
  )
})




test_that("hdrm_grouped test statistics",{
  # hypothesis = whole
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$statistic,
    0.787384763
  )

  # hypothesis = sub
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "sub",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$statistic,
    3434.748
  )

  # hypothesis = interaction
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "interaction",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$statistic,
    2.620662063
  )

  # hypothesis = identical
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "identical",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$statistic,
    2.34961433
  )

  # hypothesis = all_flat
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "flat",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$statistic,
    1662.76162
  )

  # subsampling = TRUE
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = TRUE,
      B = "10*N",
      seed = 3141
    )$statistic,
    0.841567185
  )


  # AM=0
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      AM=0,
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$statistic,
    0.787384763
  )

  # cov.equal = TRUE,
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      AM = 0,
      group = EEG$group,
      cov.equal = TRUE,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$statistic,
    0.655452437
  )

})




test_that("hdrm_grouped p.value",{
  # hypothesis = whole
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$p.value,
    0.180241264
  )

  # hypothesis = sub
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "sub",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$p.value,
    2.220446e-16
  )

  # hypothesis = interaction
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "interaction",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$p.value,
    0.022748138
  )

  # hypothesis = identical
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "identical",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$p.value,
    0.03141055
  )

  # hypothesis = all_flat
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "flat",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$p.value,
    2.220446e-16
  )

  # subsampling = TRUE
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = TRUE,
      B = "10*N",
      seed = 3141
    )$p.value,
    0.160879089
  )

  # test.direction = one.sided
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      test.direction = "one.sided",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = TRUE,
      B = "10*N",
      seed = 3141
    )$p.value,
    0.160879089
  )

  # AM=0
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      AM=0,
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$p.value,
    0.180241264
  )

  # cov.equal = TRUE,
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      AM = 0,
      group = EEG$group,
      cov.equal = TRUE,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$p.value,
    0.184825068
  )

})


test_that("hdrm_grouped f",{
  # hypothesis = whole
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    3.494251930
  )

  # hypothesis = sub
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "sub",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    1.304193980
  )

  # hypothesis = interaction
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "interaction",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    3.73918263
  )

  # hypothesis = identical
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "identical",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    3.675343
  )

  # hypothesis = all_flat
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "flat",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    4.68698151
  )

  # subsampling = TRUE
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      group = EEG$group,
      subject = EEG$subject,
      subsampling = TRUE,
      B = "10*N",
      seed = 3141
    )$f,
    2.20681838
  )

  # AM=0
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      AM=0,
      group = EEG$group,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    3.494251930
  )

  # cov.equal = TRUE,
  expect_equal(
    hdrm_grouped(
      EEG$value,
      hypothesis = "whole",
      group = EEG$group,
      cov.equal = TRUE,
      subject = EEG$subject,
      subsampling = FALSE,
      B = "10*N",
      seed = 3141
    )$f,
    1.66609860
  )

})
