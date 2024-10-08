data(birthrates)
M <- matrix(rnorm(1200), 40, 30)
L <- list(1:160, 1:40)
group <- factor(c(1,1,2,2,1,1,1,2,1,1,1,1,2,2,1,2), labels = c("west","east"))

test_that("perfect case works",{
  # hypothesis = whole
  expect_no_condition(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "whole",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # hypothesis = sub
  expect_no_condition(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "sub",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # hypothesis = interaction
  expect_no_condition(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "interaction",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # hypothesis = identical
  expect_no_condition(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "identical",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # hypothesis = all_flat
  expect_no_condition(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "flat",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # subsampling = TRUE
  expect_no_condition(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "whole",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = TRUE,
      B = "10*N"
    )
  )

  # first elemet legal, second illegal
  expect_warning(
    expect_no_error(
      hdrm_grouped_widetable(
        birthrates,
        hypothesis = c("whole", "SUB"),
        group = group,
        value = "value",
        subject = "subject",
        dimension = "dimension",
        subsampling = FALSE,
        B = "10*N"
      )
    )
  )

  # legal list to hypothesis
  expect_no_condition(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = list(TW = diag(2), TS = diag(34), TR = diag(120)),
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # eine zahl
  expect_no_condition(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "whole",
      group = group,
      value = 2,
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # argument to ...
  expect_no_condition(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "whole",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "100*N",
      a = 5
    )
  )
})


test_that("wrong input: data",{
  # data as list
  expect_error(
    hdrm_grouped_widetable(
      L,
      hypothesis = "whole",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = TRUE,
      B = "10*N"
    )
  )

  # data as df
  expect_error(
    hdrm_grouped_widetable(
      EEG,
      hypothesis = "whole",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
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
    hdrm_grouped_widetable(
      df,
      hypothesis = "sub",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )


})

test_that("wrong input: hypothesis",{

  # multiple arguments to hypothesis (default)
  expect_warning(
    hdrm_grouped_widetable(
      birthrates,
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # a number
  expect_error(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = 1,
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # illegal character
  expect_error(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = c("flart"),
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )


  # list TS missing
  expect_error(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = list(TW = diag(4), ST = diag(40)),
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )
  # TW missing
  expect_error(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = list(WT = diag(4), TS = diag(40)),
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # wrong dimension of matrix TW
  expect_error(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = list(TW = diag(3), TS = diag(40)),
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # wrong dimension of TS
  expect_error(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = list(TW = diag(2), TS = diag(51)),
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # TW not symmetrical
  expect_error(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = list(TW = diag(2) + c(0,1,0,0), TS = diag(34)),
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # TS not symmetrical
  expect_error(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = list(TW = diag(2), TS = diag(34) + c(1,0)),
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # TW not idempotent
  expect_warning(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = list(TW = diag(1:2), TS = diag(34)),
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # TS not idempotent
  expect_warning(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = list(TW = diag(2), TS = diag(1:34)),
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # matrix
  expect_error(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = diag(34),
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

})


test_that("false input: colum indicators",{
  # column name doent exist
  expect_error(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "sub",
      group = "nonexistant",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # column name doent exist
  expect_error(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "sub",
      group = rep(1,16),
      value = "nonexistant",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

})

test_that("false input: subsampling",{
  # works
  expect_no_condition(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "sub",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = 1,
      B = "10*N"
    )
  )

  # works
  expect_no_condition(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "sub",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = 0,
      B = "10*N"
    )
  )

  # logical as string
  expect_error(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "sub",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = "TRUE",
      B = "10*N"
    )
  )

  # string
  expect_error(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "sub",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = "notlogical",
      B = "10*N"
    )
  )

  # error because "maybe" makes it all a character;
  # warning because length 2
  expect_warning(
    expect_error(
      hdrm_grouped_widetable(
        birthrates,
        hypothesis = "sub",
        group = group,
        value = "value",
        subject = "subject",
        dimension = "dimension",
        subsampling = c(TRUE, "maybe"),
        B = "10*N"
      )
    )
  )

  # works, but warning because length 2
  expect_warning(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "sub",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = c(TRUE, FALSE),
      B = "10*N"
    )
  )

  # input logical matrix: works, because it can be coerced to vector
  # warning because length >1
  expect_warning(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "sub",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = diag(c(TRUE, FALSE)),
      B = "10*N"
    )
  )

  # non logical matrix: error and warning, because "only first elemt used" comes first
  expect_warning(
    expect_error(
      hdrm_grouped_widetable(
        birthrates,
        hypothesis = "sub",
        group = group,
        value = "value",
        subject = "subject",
        dimension = "dimension",
        subsampling = M,
        B = "10*N"
      )
    )
  )

  # input list: warning and error, because list is length 2
  expect_warning(
    expect_error(
      hdrm_grouped_widetable(
        birthrates,
        hypothesis = "sub",
        group = group,
        value = "value",
        subject = "subject",
        dimension = "dimension",
        subsampling = L,
        B = "10*N"
      )
    )
  )


  # input data.frame: warning and error, because birthrates[1] is possible
  expect_warning(
    expect_error(
      hdrm_grouped_widetable(
        birthrates,
        hypothesis = "sub",
        group = group,
        value = "value",
        subject = "subject",
        dimension = "dimension",
        subsampling = birthrates,
        B = "10*N"
      )
    )
  )
})


test_that("false input: B",{
  # character length 2
  expect_warning(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "sub",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = c("10*N", "20*N")
    )
  )

  # numeric length 2
  expect_warning(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "sub",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = c(100, 1000)
    )
  )

  # c(100, "10*N"): "100" works!
  expect_warning(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "sub",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = c(100, "10*N")
    )
  )

  # negative number
  expect_error(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "sub",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = TRUE,
      B = -10
    )
  )

  # negative character
  expect_error(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "sub",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = TRUE,
      B = "-10*N"
    )
  )

  # function with nonexistant argument
  expect_error(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "sub",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = TRUE,
      B = "10 *asdfghjkl"
    )
  )

  # B wrong class
  expect_error(
    hdrm_grouped_widetable(
      birthrates,
      hypothesis = "sub",
      group = group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = TRUE,
      B = M
    )
  )


})
