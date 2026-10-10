testthat::context("Testing 'filtering'")

testthat::test_that("filtering-se", {
  sacurine.se <- reading(system.file("extdata/sacurine", package = "phenomis"))

  assay.mn <- assay(sacurine.se)
  assay.mn[assay.mn < 1e5] <- NA
  assay(sacurine.se) <- assay.mn

  sacurine.se <- filtering(sacurine.se)

  expect_dims.vi <- c(85, 200)

  testthat::expect_equal(
    dim(sacurine.se),
    expect_dims.vi
  )
})

testthat::test_that("filtering-se-gender", {
  sacurine.se <- reading(system.file("extdata/sacurine", package = "phenomis"))

  assay.mn <- assay(sacurine.se)
  assay.mn[assay.mn < 1e5] <- NA
  assay(sacurine.se) <- assay.mn

  sacurine.se <- filtering(sacurine.se, class.c = "gender")

  expect_dims.vi <- c(90, 196)

  testthat::expect_equal(
    dim(sacurine.se),
    expect_dims.vi
  )
})

testthat::test_that("filtering-mae", {
  prometis.mae <- reading(system.file("extdata/prometis", package = "phenomis"))

  for (set.c in names(prometis.mae)) {
    set.se <- prometis.mae[[set.c]]
    set.mn <- assay(set.se)
    set.mn[set.mn < quantile(c(set.mn), 0.2)] <- NA
    assay(set.se) <- set.mn
    prometis.mae[[set.c]] <- set.se
  }

  prometis.mae <- filtering(prometis.mae)

  expect_dims.mn <- matrix(c(71, 28, 75, 28),
    nrow = 2,
    dimnames = list(NULL, c("metabo", "proteo"))
  )

  testthat::expect_equal(
    vapply(names(prometis.mae),
      function(set.c) dim(prometis.mae[[set.c]]),
      FUN.VALUE = integer(2)
    ),
    expect_dims.mn
  )
})
