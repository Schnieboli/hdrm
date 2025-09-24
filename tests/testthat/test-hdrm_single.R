#Tests for vector format

# Loading the dataset
data(EEG)

L <- list(1:160, 1:40)


test_that("perfect case works",{
  # hypothesis = flat
  expect_no_condition(
    hdrm_single(
      EEG$value,
      hypothesis = "flat",
      subject = EEG$subject
    )
  )


  # first elemet legal, second illegal
  expect_warning(
    expect_no_error(
      hdrm_single(
        EEG$value,
        hypothesis = c("flat", "SUB"),
        subject = EEG$subject
      )
    )
  )

  # legal list to hypothesis
  expect_no_condition(
    hdrm_single(
      EEG$value,
      hypothesis = diag(40),
      subject = EEG$subject
    )
  )


  # argument to ...
  expect_no_condition(
    hdrm_single(
      EEG$value,
      hypothesis = "flat",
      subject = EEG$subject,
      subsampling = FALSE
    )
  )
})


test_that("wrong input: data",{
  # data as list
  expect_error(
    hdrm_single(
      L1,
      hypothesis = "flat",
      subject = EEG$subject
    )
  )
})


test_that("missing values",{
  # NA in value
  vector <- EEG$value
  vector[123] <- NA
  expect_warning(
    hdrm_single(
      vector,
      hypothesis = "flat",
      subject = EEG$subject
    )
  )


  # NA in subject
  subjectvector <- EEG$subject
  subjectvector[145] <- NA
  expect_error(
    hdrm_single(
      EEG$value,
      hypothesis = "flat",
      subject = subjectvector,
    )
  )


  # different length of subject and data
  expect_error(
    hdrm_single(
      EEG$value[-1],
      hypothesis = "flat",
      subject = EEG$subject,
    )
  )

  # nonexistant column
  expect_error(
    hdrm_single(
      EEG$value,
      hypothesis = "flat",
      subject = "nonexistant",
    )
  )

})


test_that("wrong input: hypothesis",{


  # a number
  expect_error(
    hdrm_single(
      EEG$value,
      hypothesis = 1,
      subject = EEG$subject
    )
  )

  # illegal character
  expect_error(
    hdrm_single(
      EEG$value,
      hypothesis = c("flart"),
      subject = EEG$subject
    )
  )


  # list TS missing
  expect_error(
    hdrm_single(
      EEG$value,
      hypothesis = diag(1:41),
      subject = EEG$subject
    )
  )



  # wrong dimension of hypothesis
  expect_error(
    hdrm_single(
      EEG$value,
      hypothesis = matrix(0, 41, 39),
      subject = EEG$subject
    )
  )

  # TW not symmetrical
  expect_error(
    hdrm_single(
      EEG$value,
      hypothesis = list(TW = diag(4) + c(0, 1), TS = diag(40)),
      subject = EEG$subject
    )
  )

  # list
  expect_error(
    hdrm_single(
      EEG$value,
      hypothesis = list(diag(40)),
      subject = EEG$subject
    )
  )

})





test_that("hdrm_single test statistics",{
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


test_that("hdrm_single  p.values",{
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

test_that("hdrm_single f",{
          # hypothesis = flat
          expect_equal(
            hdrm_single(
              EEG$value,
              hypothesis = "flat",
              subject = EEG$subject
            )$f,
            1.0000
          )            })


#Tests for matrix format


data("birthrates")
Matrixbirthrates=as.matrix(birthrates)
L2 <- list(1:160, 1:40)


test_that("perfect case works",{
  # hypothesis = flat
  expect_no_condition(
    hdrm_single(
      Matrixbirthrates,
      hypothesis = "flat"
    )
  )


  # first element legal, second illegal
  expect_warning(
    expect_no_error(
      hdrm_single(
        Matrixbirthrates,
        hypothesis = c("flat", "SUB")
      )
    )
  )

  # legal list to hypothesis
  expect_no_condition(
    hdrm_single(
      Matrixbirthrates,
      hypothesis = diag(34)
    )
  )


  # argument to ...
  expect_no_condition(
    hdrm_single(
      Matrixbirthrates,
      hypothesis = "flat",
      subsampling = FALSE
    )
  )
})

test_that("wrong input: data",{
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
    subject= 1:34
  )
  )
  # data as df
  expect_error(
    hdrm_single(
      data = EEG,
      hypothesis = "flat"
    )
  )
})


test_that("missing values",{
  # NA in value
  M <- Matrixbirthrates
  M[4,7] <- NA
  expect_warning(
    hdrm_single(
      data = M,
      hypothesis = "flat"
    )
  )

  expect_warning(
    expect_true(
      all(dim(hdrm_single(
        data = M,
        hypothesis = "flat"
      )$data) == dim(M) - c(0,1))
    )
  )

})


test_that("wrong input: hypothesis",{


  # a number
  expect_error(
    hdrm_single(
      Matrixbirthrates,
      hypothesis = 1
    )
  )

  # illegal character
  expect_error(
    hdrm_single(
      Matrixbirthrates,
      hypothesis = c("flart")
    )
  )


  # wrong dimension of hypothesis
  expect_error(
    hdrm_single(
      Matrixbirthrates,
      hypothesis = matrix(0, 41, 39)
    )
  )

  # TW not symmetrical
  expect_warning(
    hdrm_single(
      Matrixbirthrates,
      hypothesis = diag(34) + c(0,1)
    )
  )

  # list
  expect_error(
    hdrm_single(
      Matrixbirthrates,
      hypothesis = list(diag(40))
    )
  )

})

test_that("hdrm_single test statistics",{
  # hypothesis = flat
  expect_equal(
    hdrm_single(
      Matrixbirthrates,
      hypothesis = "flat"
    )$statistic,
    7.481676
  )
})


test_that("hdrm_single p.values",{
  # hypothesis = flat
  expect_equal(
    hdrm_single(
      Matrixbirthrates,
      hypothesis = "flat"
    )$p.value,
    0.000666403,tolerance=0.0000002
  )
})


test_that("hdrm_single f",{
  # hypothesis = flat
  expect_equal(
    hdrm_single(
      Matrixbirthrates,
      hypothesis = "flat"
    )$f,
    1
  )
})




