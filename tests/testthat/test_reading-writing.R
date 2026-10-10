testthat::context("Testing 'read-write'")

#### .reading ####

testthat::test_that(".reading", {
  sacdir.c <- system.file("extdata/sacurine", package = "phenomis")

  sacurine.se1 <- .read_se(sacdir.c)

  testthat::expect_true(class(sacurine.se1) == "SummarizedExperiment")

  testthat::expect_equal(assay(sacurine.se1)[1, 1],
    477491,
    tolerance = 1e-3
  )
  # alternatively
  sacurine.se2 <- .read_se(
    NA,
    file.path(sacdir.c, "Galaxy1_dataMatrix.tabular"),
    file.path(
      sacdir.c,
      "Galaxy2_sampleMetadata.tabular"
    ),
    file.path(
      sacdir.c,
      "Galaxy3_variableMetadata.tabular"
    )
  )

  testthat::expect_true(class(sacurine.se2) == "SummarizedExperiment")

  testthat::expect_equal(assay(sacurine.se2)[1, 1],
    477491,
    tolerance = 1e-3
  )


  testthat::expect_error(.read_se(
    NA,
    file.path(sacdir.c, "Galaxy1_dataMatrix.tsv"),
    file.path(
      sacdir.c,
      "Galaxy2_sampleMetadata.tabular"
    ),
    file.path(
      sacdir.c,
      "Galaxy3_variableMetadata.tabular"
    )
  ))
})

#### reading_se ####

testthat::test_that("reading_SummarizedExperiment", {
  sacdir.c <- system.file("extdata/sacurine", package = "phenomis")

  sacurine.se2 <- reading(sacdir.c)

  testthat::expect_true(class(sacurine.se2) == "SummarizedExperiment")

  testthat::expect_equal(assay(sacurine.se2)[1, 1],
    477491,
    tolerance = 1e-3
  )
  # alternatively
  sacurine.se3 <- reading(NA,
    files.ls = list(
      dataMatrix = file.path(sacdir.c, "Galaxy1_dataMatrix.tabular"),
      sampleMetadata = file.path(sacdir.c, "Galaxy2_sampleMetadata.tabular"),
      variableMetadata = file.path(sacdir.c, "Galaxy3_variableMetadata.tabular")
    )
  )

  testthat::expect_true(class(sacurine.se3) == "SummarizedExperiment")

  testthat::expect_equal(assay(sacurine.se3)[1, 1],
    477491,
    tolerance = 1e-3
  )


  testthat::expect_error(reading(NA,
    files.ls = list(
      dataMatrix = file.path(sacdir.c, "Galaxy1_dataMatrix.tsv"),
      sampleMetadata = file.path(sacdir.c, "Galaxy2_sampleMetadata.tabular"),
      variableMetadata = file.path(sacdir.c, "Galaxy3_variableMetadata.tabular")
    )
  ))
})

#### reading_mae ####

testthat::test_that("reading_MultiAssayExperiment", {
  prometis_dir.c <- system.file("extdata/prometis", package = "phenomis")

  prometis.mae <- reading(prometis_dir.c)

  prometis.mae <- transforming(prometis.mae)

  testthat::expect_true(class(prometis.mae) == "MultiAssayExperiment")

  testthat::expect_identical(names(prometis.mae), c("metabo", "proteo"))

  testthat::expect_equal(MultiAssayExperiment::assays(prometis.mae)[["metabo"]]["isobutyric acid", "L818f"],
    24.64327,
    tolerance = 1e-3
  )


  met.mae <- reading(prometis_dir.c, subsets.vc = "metabo")
  met.mae <- transforming(met.mae)

  testthat::expect_true(class(met.mae) == "MultiAssayExperiment")

  testthat::expect_identical(names(met.mae), "metabo")

  testthat::expect_equal(MultiAssayExperiment::assays(met.mae)[["metabo"]]["isobutyric acid", "L818f"],
    24.64327,
    tolerance = 1e-3
  )
})
