data(EEG)
M <- matrix(rnorm(1200), 40, 30)
L <- list(1:160, 1:40)


test_that("perfect case works",{
  # hypothesis = flat
  expect_no_condition(
    hdrm_single_longtable(
      EEG,
      hypothesis = "flat",
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )


  # first elemet legal, second illegal
  expect_warning(
    expect_no_error(
      hdrm_single_longtable(
        EEG,
        hypothesis = c("flat", "SUB"),
        value = "value",
        subject = "subject",
        dimension = "dimension"
      )
    )
  )

  # legal list to hypothesis
  expect_no_condition(
    hdrm_single_longtable(
      EEG,
      hypothesis = diag(40),
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )

  # eine zahl
  expect_no_condition(
    hdrm_single_longtable(
      EEG,
      hypothesis = "flat",
      value = 2,
      subject = "subject",
      dimension = "dimension"
    )
  )

  # argument to ...
  expect_no_condition(
    hdrm_single_longtable(
      EEG,
      hypothesis = "flat",
      value = "value",
      subject = "subject",
      dimension = "dimension",
      subsampling = FALSE
    )
  )
})


test_that("wrong input: data",{
  # data as list
  expect_error(
    hdrm_single_longtable(
      L,
      hypothesis = "flat",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )

  # data as matrix
  expect_error(
    hdrm_single_longtable(
      M,
      hypothesis = "flat",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )

})


test_that("missing values",{
  # NA in value
  df <- EEG
  df$value[123] <- NA
  expect_warning(
    hdrm_single_longtable(
      df,
      hypothesis = "flat",
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )


  # NA in subject
  df <- EEG
  df$subject[145] <- NA
  expect_error(
    hdrm_single_longtable(
      df,
      hypothesis = "flat",
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )

  # NA in dimension
  df <- EEG
  df$dimension[165] <- NA
  expect_error(
    hdrm_single_longtable(
      df,
      hypothesis = "flat",
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )

  # NA in unrelated column
  df <- EEG
  df$variable[134] <- NA
  expect_no_condition(
    hdrm_single_longtable(
      df,
      hypothesis = "flat",
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )

  # nonexistant column
  expect_error(
    hdrm_single_longtable(
      EEG,
      hypothesis = "flat",
      value = "value",
      subject = "nonexistant",
      dimension = "dimension"
    )
  )

})


test_that("wrong input: hypothesis",{


  # a number
  expect_error(
    hdrm_single_longtable(
      EEG,
      hypothesis = 1,
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )

  # illegal character
  expect_error(
    hdrm_single_longtable(
      EEG,
      hypothesis = c("flart"),
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )


  # list TS missing
  expect_error(
    hdrm_single_longtable(
      EEG,
      hypothesis = diag(1:41),
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )



  # wrong dimension of hypothesis
  expect_error(
    hdrm_single_longtable(
      EEG,
      hypothesis = matrix(0, 41, 39),
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )

  # TW not symmetrical
  expect_error(
    hdrm_single_longtable(
      EEG,
      hypothesis = list(TW = diag(4) + c(0, 1), TS = diag(40)),
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )

  # list
  expect_error(
    hdrm_single_longtable(
      EEG,
      hypothesis = list(diag(40)),
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )

})

test_that("false input: colum indicators",{


  # column name doent exist
  expect_error(
    hdrm_grouped_longtable(
      EEG,
      hypothesis = "sub",
      value = "nonexistant",
      subject = "subject",
      dimension = "dimension"
    )
  )

  # column name doent exist
  expect_error(
    hdrm_single_longtable(
      EEG,
      hypothesis = "flat",
      value = "value",
      subject = "nonexistant",
      dimension = "dimension"
    )
  )

  # column name doent exist
  expect_error(
    hdrm_single_longtable(
      EEG,
      hypothesis = "flat",
      value = "value",
      subject = "subject",
      dimension = "nonexistant"
    )
  )

  # number too high
  expect_error(
    hdrm_single_longtable(
      EEG,
      hypothesis = "flat",
      value = 10,
      subject = "subject",
      dimension = "dimension"
    )
  )

  # column name doent exist
  expect_error(
    hdrm_single_longtable(
      EEG,
      hypothesis = "flat",
      value = "10",
      subject = "subject",
      dimension = "dimension"
    )
  )

  # column index doent exist
  expect_error(
    hdrm_single_longtable(
      EEG,
      hypothesis = "flat",
      value = "value",
      subject = 10,
      dimension = "dimension"
    )
  )

  expect_error(
    hdrm_single_longtable(
      EEG,
      hypothesis = "flat",
      value = "value",
      subject = "subject",
      dimension = 10
    )
  )

  # vector
  expect_error(
    hdrm_single_longtable(
      EEG,
      hypothesis = "flat",
      value = EEG$value,
      subject = "subject",
      dimension = "dimension"
    )
  )



  # vector
  expect_error(
    hdrm_single_longtable(
      EEG,
      hypothesis = "flat",
      group = "group",
      value = "value",
      subject = EEG$subject,
      dimension = "dimension"
    )
  )

  expect_error(
    hdrm_single_longtable(
      EEG,
      hypothesis = "flat",
      value = "value",
      subject = "subject",
      dimension = EEG$dimension
    )
  )

})

test_that("non continuous levels",{
  df <- EEG
  levels(df$subject) <- c(1, 3:161)

  expect_no_condition(
    hdrm_single_longtable(
      EEG,
      hypothesis = "flat",
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )

  df <- EEG
  levels(df$dimension) <- c(1, 3:50)

  expect_no_condition(
    hdrm_single_longtable(
      EEG,
      hypothesis = "flat",
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )

  df <- EEG
  levels(df$group) <- c(1,4,3,9)

  expect_no_condition(
    hdrm_single_longtable(
      EEG,
      hypothesis = "flat",
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )
})
