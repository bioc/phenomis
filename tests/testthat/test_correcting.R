testthat::context("Testing 'correcting'")

testthat::test_that("correcting-se", {

  sacurine.se <- reading(system.file("extdata/sacurine",
                                     package = "phenomis"))
  sacurine.se <- correcting(sacurine.se, method.vc = "loess")

  testthat::expect_equal(assay(sacurine.se)["Testosterone glucuronide",
                                            "HU_neg_020"],
                         44136.83,
                         tolerance = 1e-6)

})

testthat::test_that("correcting-mae", {

  sacurine.se <- reading(system.file("extdata/sacurine", package = "phenomis"))
  sac_map.df <- data.frame(primary = colnames(sacurine.se),
                           colname = colnames(sacurine.se))
  map.ls <- list(sac1 = sac_map.df,
                 sac2 = sac_map.df)
  map.df <- MultiAssayExperiment::listToMap(map.ls)
  sac.mae <- MultiAssayExperiment::MultiAssayExperiment(experiments = list(sac1 = sacurine.se,
                                                                           sac2 = sacurine.se),
                                                        colData = colData(sacurine.se),
                                                        sampleMap = map.df)
  sac.mae <- correcting(sac.mae, method.vc = c("loess", "serrf"))

  testthat::expect_equal(assays(sac.mae)[["sac1"]]["Testosterone glucuronide",
                                                   "HU_neg_020"],
                         44136.83,
                         tolerance = 1e-6)

  if (Sys.info()["sysname"] != "Darwin") # Mac OS X
    testthat::expect_equal(assays(sac.mae)[["sac2"]]["Testosterone glucuronide",
                                                     "HU_neg_020"],
                           64426.17,
                           tolerance = 1e-6)


})

testthat::test_that("correcting-nazeros", {

  sacurine.se <- reading(system.file("extdata/sacurine", package = "phenomis"))

  ## introducing 0 in all QC of feature 1
  data.mn <- t(SummarizedExperiment::assay(sacurine.se))
  samp.df <- as.data.frame(SummarizedExperiment::colData(sacurine.se))
  pool.vi <- which(samp.df[, "sampleType"] == "pool")
  data.mn[pool.vi, 1] <- 0
  SummarizedExperiment::assay(sacurine.se) <- t(data.mn)

  ## running the correction
  normalized.se <- correcting(sacurine.se, method.vc = "loess")

  testthat::expect_equal(nrow(sacurine.se), 113)
  testthat::expect_equal(nrow(normalized.se), 112)

})
