data(EEG)
M <- matrix(rnorm(1200), 40, 30)
L <- list(1:160, 1:40)

test_that("perfect case works",{
  # hypothesis = whole
  expect_no_condition(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "whole",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # hypothesis = sub
  expect_no_condition(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )
  # hypothesis = interaction
  expect_no_condition(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "interaction",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )
  # subsampling = TRUE
  expect_no_condition(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "whole",
      group = "group",
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
      hdrm_grouped_longtable(
        EEG,
        hypothesis = c("whole", "SUB"),
        group = "group",
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
    hdrm_grouped_longtable(
      EEG,
      hypothesis = list(TW = diag(4), TS = diag(40), TR = diag(120)),
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # eine zahl
  expect_no_condition(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "whole",
      group = 1,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # argument to ...
  expect_no_condition(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "whole",
      group = "group",
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
    hdrm_grouped_longtable(
      L,
      hypothesis = "whole",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = TRUE,
      B = "10*N"
    )
  )

  # data as matrix
  expect_error(
    hdrm_grouped_longtable(
      M,
      hypothesis = "whole",
      group = "group",
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
  df <- EEG
  df$value[123] <- NA
  expect_warning(
    hdrm_grouped_longtable(
      df,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # NA in group
  df <- EEG
  df$group[234] <- NA
  expect_error(
    hdrm_grouped_longtable(
      df,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # NA in subject
  df <- EEG
  df$subject[145] <- NA
  expect_error(
    hdrm_grouped_longtable(
      df,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # NA in dimension
  df <- EEG
  df$dimension[165] <- NA
  expect_error(
    hdrm_grouped_longtable(
      df,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # NA in unrelated column
  df <- EEG
  df$variable[134] <- NA
  expect_no_condition(
    hdrm_grouped_longtable(
      df,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # nonexistant column
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "nonexistant",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = TRUE,
      B = "10*N"
    )
  )

})

test_that("wrong input: hypothesis",{

  # multiple arguments
  expect_warning(
    hdrm_grouped_longtable(
      EEG,
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

 # a number
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = 1,
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # illegal character
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = c("flart"),
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )


  # list TS missing
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = list(TW = diag(4), ST = diag(40)),
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )
  # TW missing
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = list(WT = diag(4), TS = diag(40)),
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # wrong dimension of matrix TW
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = list(TW = diag(3), TS = diag(40)),
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

 # wrong dimension of TS
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = list(TW = diag(4), TS = diag(51)),
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # TW not symmetrical
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = list(TW = diag(4) + c(0, 1), TS = diag(40)),
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # TS not symmetrical
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = list(TW = diag(4), TS = diag(40) + c(1,0)),
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # TW not idempotent
  expect_warning(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = list(TW = diag(1:4), TS = diag(40)),
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # TS not idempotent
  expect_warning(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = list(TW = diag(4), TS = diag(1:40)),
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # matrix
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = diag(160),
      group = "group",
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
    hdrm_grouped_longtable(
      EEG,
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
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "grouü",
      value = "nonexistant",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # column name doent exist
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "nonexistant",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # column name doent exist
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "nonexistant",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # number too high
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = 10,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # column name doent exist
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = "10",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # column name doent exist
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = 10,
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = 10,
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # vector
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = EEG$group,
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # vector
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = EEG$value,
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  # vector
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = EEG$subject,
      dimension = "dimension",
      subsampling = FALSE,
      B = "10*N"
    )
  )

  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = EEG$dimension,
      subsampling = FALSE,
      B = "10*N"
    )
  )

})

test_that("false input: subsampling",{
  # works
  expect_no_condition(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = 1,
      B = "10*N"
    )
  )

  # works
  expect_no_condition(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = 0,
      B = "10*N"
    )
  )

  # logical as string
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = "TRUE",
      B = "10*N"
    )
  )

  # string
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
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
      hdrm_grouped_longtable(
        EEG,
        hypothesis = "sub",
        group = "group",
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
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
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
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
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
      hdrm_grouped_longtable(
        EEG,
        hypothesis = "sub",
        group = "group",
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
      hdrm_grouped_longtable(
        EEG,
        hypothesis = "sub",
        group = "group",
        value = "value",
        subject = "subject",
        dimension = "dimension",
        subsampling = L,
        B = "10*N"
      )
    )
  )


  # input data.frame: warning and error, because EEG[1] is possible
  expect_warning(
    expect_error(
      hdrm_grouped_longtable(
        EEG,
        hypothesis = "sub",
        group = "group",
        value = "value",
        subject = "subject",
        dimension = "dimension",
        subsampling = EEG,
        B = "10*N"
      )
    )
  )
})


test_that("false input: B",{
  # character length 2
  expect_warning(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = c("10*N", "20*N")
    )
  )

  # numeric length 2
  expect_warning(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = c(100, 1000)
    )
  )

  # c(100, "10*N"): "100" works!
  expect_warning(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE,
      B = c(100, "10*N")
    )
  )

  # negative number
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = TRUE,
      B = -10
    )
  )

  # negative character
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = TRUE,
      B = "-10*N"
    )
  )

  # function with nonexistant argument
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = TRUE,
      B = "10 *asdfghjkl"
    )
  )

  # B wrong class
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = TRUE,
      B = M
    )
  )


})
