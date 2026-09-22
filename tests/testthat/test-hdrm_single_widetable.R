data("birthrates")
data(EEG)
M <- matrix(rnorm(1200), 40, 30)
L <- list(1:160, 1:40)


test_that("perfect case works",{
  # hypothesis = flat
  expect_no_condition(
    hdrm_single_widetable(
      birthrates,
      hypothesis = "flat"
    )
  )


  # first elemet legal, second illegal
  expect_warning(
    expect_no_error(
      hdrm_single_widetable(
        birthrates,
        hypothesis = c("flat", "SUB")
      )
    )
  )

  # legal list to hypothesis
  expect_no_condition(
    hdrm_single_widetable(
      birthrates,
      hypothesis = diag(34)
    )
  )


  # argument to ...
  expect_no_condition(
    hdrm_single_widetable(
      birthrates,
      hypothesis = "flat",
      subsampling = FALSE
    )
  )
})

test_that("wrong input: data",{
  # data as list
  expect_error(
    hdrm_single_widetable(
      L,
      hypothesis = "flat",
      group = "group",
      value = "value",
      subject = "subject",
      dimension = "dimension"
    )
  )

  # data as df
  expect_error(
    hdrm_single_widetable(
      data = EEG,
      hypothesis = "flat",
    )
  )

})


test_that("missing values",{
  # NA in value
  M <- birthrates
  M[4,7] <- NA
  expect_warning(
    hdrm_single_widetable(
      data = M,
      hypothesis = "flat"
    )
  )

  expect_warning(
    expect_true(
      all(dim(hdrm_single_widetable(
        data = M,
        hypothesis = "flat"
      )$data) == dim(M) - c(0,1))
    )
  )

})


test_that("wrong input: hypothesis",{


  # a number
  expect_error(
    hdrm_single_widetable(
      birthrates,
      hypothesis = 1
    )
  )

  # illegal character
  expect_error(
    hdrm_single_widetable(
      birthrates,
      hypothesis = c("flart")
    )
  )


  # wrong dimension of hypothesis
  expect_error(
    hdrm_single_widetable(
      birthrates,
      hypothesis = matrix(0, 41, 39)
    )
  )

  # TW not symmetrical
  expect_warning(
    hdrm_single_widetable(
      birthrates,
      hypothesis = diag(34) + c(0,1)
    )
  )

  # list
  expect_error(
    hdrm_single_widetable(
      birthrates,
      hypothesis = list(diag(40))
    )
  )

})


