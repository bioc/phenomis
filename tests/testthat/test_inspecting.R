testthat::context("Testing 'inspecting'")

testthat::test_that("inspecting-se", {
  sacurine.se <- reading(system.file("extdata/sacurine", package = "phenomis"))

  sacurine.se <- inspecting(sacurine.se,
    figure.c = "none",
    report.c = "none"
  )

  testthat::expect_equivalent(rowData(sacurine.se)["(2-methoxyethoxy)propanoic acid isomer", "pool_CV"],
    0.3160307,
    tolerance = 1e-6
  )

  proteo.se <- reading(system.file("extdata/prometis/proteo",
    package = "phenomis"
  ))

  proteo.se <- transforming(proteo.se)

  proteo.se <- inspecting(proteo.se,
    figure.c = "none",
    report.c = "none"
  )

  testthat::expect_equivalent(colData(proteo.se)["L818f", "deci_pval"],
    0.0525734,
    tolerance = 1e-6
  )
})

testthat::test_that("inspecting-mae", {
  prometis.mae <- reading(system.file("extdata/prometis", package = "phenomis"))
  prometis.mae <- transforming(prometis.mae)

  prometis.mae <- inspecting(prometis.mae,
    figure.c = "none",
    report.c = "none"
  )

  testthat::expect_equivalent(
    colData(prometis.mae[["metabo"]])[
      "W634m",
      "hotel_pval"
    ],
    0.818558,
    tolerance = 1e-6
  )
})
