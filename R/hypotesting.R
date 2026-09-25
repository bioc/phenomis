#### hypotesting (MultiAssayExperiment) ####

#' @rdname hypotesting
#' @export
setMethod(
  "hypotesting", signature(x = "MultiAssayExperiment"),
  function(x,
           test.c = c(
             "ttest", "limma", "wilcoxon",
             "anova", "kruskal",
             "pearson", "spearman",
             "limma2ways", "limma2waysInter",
             "anova2ways", "anova2waysInter"
           )[2],
           factor_names.vc,
           factor_levels.ls = list(
             factor1.vc = "default",
             factor2.vc = "default"
           ),
           adjust.c = c(
             "holm",
             "hochberg",
             "hommel",
             "bonferroni",
             "BH",
             "BY",
             "fdr",
             "none"
           )[5],
           adjust_thresh.n = 0.05,
           signif_maxprint.i = NA,
           title.c = NA,
           display_signif.l = FALSE,
           prefix.c = "",
           figure.c = c(
             "none",
             "interactive",
             "interactive_plotly",
             "myfile.pdf"
           )[2],
           report.c = c(
             "none",
             "interactive",
             "myfile.txt"
           )[2]) {
    if (!(report.c %in% c("none", "interactive"))) {
      sink(report.c, append = TRUE)
    }

    report_set.c <- report.c
    if (report_set.c != "none") {
      report_set.c <- "interactive"
    }

    if (!(figure.c %in% c(
      "none",
      "interactive",
      "interactive_plotly"
    ))) {
      grDevices::pdf(figure.c)
      figure_set.c <- "interactive"
    } else {
      figure_set.c <- figure.c
    }

    for (set.c in names(x)) {
      main.c <- paste0(
        ifelse(!is.na(title.c),
          paste0(title.c, ", "), ""
        ),
        set.c
      )

      if (report.c != "none") {
        message("Hypothesis testing for the '", set.c, "' dataset:")
      }

      set.se <- x[[set.c]]

      if (!(all(colnames(set.se) %in% rownames(colData(x))))) {
        stop("Identical naming of samples in the MultiAssayExperiment (e.g. in its colData slot) and in all SummarizedExperiment objects it contains is required in the current implementation of the hypotesting method")
      }

      for (factor.c in factor_names.vc) {
        if (!(factor.c %in% colnames(colData(x)))) {
          stop("Factor '", factor.c, "' was not found in the columns from the colData Data Frame of the MultiAssayExperiment object")
        }

        colData(set.se)[, factor.c] <- colData(x)[colnames(set.se), factor.c]
      }

      set.se <- hypotesting(
        x = set.se,
        test.c = test.c,
        factor_names.vc = factor_names.vc,
        factor_levels.ls = factor_levels.ls,
        adjust.c = adjust.c,
        adjust_thresh.n = adjust_thresh.n,
        signif_maxprint.i = signif_maxprint.i,
        title.c = main.c,
        display_signif.l = display_signif.l,
        prefix.c = prefix.c,
        set.c = set.c,
        figure.c = figure_set.c,
        report.c = report_set.c
      )

      for (factor.c in factor_names.vc) {
        colData(set.se)[, factor.c] <- NULL
      }

      x[[set.c]] <- set.se
    }

    if (!(figure.c %in% c(
      "none",
      "interactive",
      "interactive_plotly"
    ))) {
      grDevices::dev.off()
    }

    if (!(report.c %in% c("none", "interactive"))) {
      sink()
    }

    methods::validObject(x)

    return(invisible(x))
  }
)


#### hypotesting (SummarizedExperiment) ####

#' @rdname hypotesting
#' @export
setMethod(
  "hypotesting", signature(x = "SummarizedExperiment"),
  function(x,
           test.c = c(
             "ttest", "limma", "wilcoxon",
             "anova", "kruskal",
             "pearson", "spearman",
             "limma2ways", "limma2waysInter",
             "anova2ways", "anova2waysInter"
           )[2],
           factor_names.vc,
           factor_levels.ls = list(
             factor1.vc = "default",
             factor2.vc = "default"
           ),
           adjust.c = c(
             "holm",
             "hochberg",
             "hommel",
             "bonferroni",
             "BH",
             "BY",
             "fdr",
             "none"
           )[5],
           adjust_thresh.n = 0.05,
           signif_maxprint.i = NA,
           title.c = NA,
           display_signif.l = FALSE,
           prefix.c = "",
           set.c = NA_character_,
           figure.c = c(
             "none",
             "interactive",
             "interactive_plotly",
             "myfile.pdf"
           )[2],
           report.c = c(
             "none",
             "interactive",
             "myfile.txt"
           )[2]) {
    if (!(report.c %in% c("none", "interactive"))) {
      sink(report.c, append = TRUE)
    }

    # Checks

    factorVl <- factor_names.vc %in% colnames(colData(x))

    if (sum(!factorVl)) {
      unknown_factor.c <- paste(factor_names.vc[!factorVl],
        collapse = "', '"
      )
      coldata_column.c <- paste(colnames(colData(x)), collapse = "', '")
      stop(
        "The following factor(s) '",
        unknown_factor.c,
        "' was/were not found in the column names of colData(x): '",
        coldata_column.c, "'"
      )
    }

    testVc <- c(
      "ttest", "limma", "wilcoxon",
      "anova", "kruskal",
      "pearson", "spearman",
      "limma2ways", "limma2waysInter",
      "anova2ways", "anova2waysInter"
    )

    if (!test.c %in% testVc) {
      avail_test.c <- paste(testVc, collapse = "', '")
      stop(
        "'test.c = ", test.c, ": 'test.c' must be in '",
        avail_test.c, "'"
      )
    }

    if (!("transformed" %in% names(x@metadata))) {
      stop("Information about the transformation of the data is missing. Please use the 'transforming' method first.")
    }
    transformed.c <- x@metadata[["transformed"]]
    if (transformed.c == "none" & !(test.c %in% c(
      "ttest", "limma", "wilcoxon",
      "pearson", "spearman"
    ))) {
      stop("Anovas computation of fold changes are currently valid for log transformed data only. Please apply the 'transforming' method first with the 'log2' or 'log10' parameters.")
    }

    # Hypothesis testing

    if (report.c != "none") {
      message("Performing '", test.c, "'...")
    }

    if (test.c %in% c(
      "ttest", "limma", "wilcoxon",
      "pearson", "spearman"
    )) {
      testLs <- .twoSampCorTests(
        data.mn = t(assay(x)),
        samp.df = colData(x),
        feat.df = rowData(x),
        test.c = test.c,
        factorNameC = factor_names.vc,
        factorLevelsVc = factor_levels.ls[[1]],
        adjust.c = adjust.c,
        adjust_thresh.n = adjust_thresh.n,
        transformed.c = transformed.c,
        title.c = title.c,
        display_signif.l = display_signif.l,
        prefix.c = prefix.c,
        figure.c = figure.c
      )
    } else if (test.c %in% c("anova", "kruskal")) {
      testLs <- .anovas(
        data.mn = t(assay(x)),
        samp.df = colData(x),
        feat.df = rowData(x),
        test.c = test.c,
        factorNameC = factor_names.vc,
        factorLevelsVc = factor_levels.ls[[1]],
        adjust.c = adjust.c,
        adjust_thresh.n = adjust_thresh.n,
        transformed.c = transformed.c,
        title.c = title.c,
        prefix.c = prefix.c,
        figure.c = figure.c
      )
    } else if (test.c %in% c(
      "limma2ways", "limma2waysInter",
      "anova2ways", "anova2waysInter"
    )) {
      testLs <- .anovas2ways(
        data.mn = t(assay(x)),
        samp.df = colData(x),
        feat.df = rowData(x),
        test.c = test.c,
        factor_names.vc = factor_names.vc,
        factor_levels.ls = factor_levels.ls,
        adjust.c = adjust.c,
        adjust_thresh.n = adjust_thresh.n,
        transformed.c = transformed.c,
        title.c = title.c,
        prefix.c = prefix.c,
        figure.c = figure.c
      )
    }

    rowData(x) <- testLs[["feat.df"]]
    methods::validObject(x)
    metric.mn <- testLs[["metric.mn"]]

    ## Printing significant features

    signiColVl <- grepl("_signif", colnames(metric.mn))
    signiVl <- rowSums(metric.mn[, signiColVl, drop = FALSE],
      na.rm = TRUE
    ) > 0

    if (report.c != "none") {
      if (sum(signiVl) > 0) {
        metric.mn <- metric.mn[, !signiColVl, drop = FALSE]

        .signifPrint(
          metric.mn = metric.mn,
          signiVl = signiVl,
          signif_maxprint.i,
          adjust.c = adjust.c,
          adjust_thresh.n = adjust_thresh.n
        )
      } else {
        message(
          "\nNo significant variable found at the selected ",
          adjust_thresh.n, " level."
        )
      }
    }

    if (!(report.c %in% c("none", "interactive"))) {
      sink()
    }

    ## Plotting ("ttest", "limma", "wilcoxon", "pearson", "spearman")

    if (test.c %in% c(
      "ttest", "limma", "wilcoxon",
      "pearson", "spearman"
    ) && figure.c != "none") {
      plot_hypotesting(x,
        test.c = test.c,
        factor.c = factor_names.vc,
        adjust.c = adjust.c,
        adjust_thresh.n = adjust_thresh.n,
        display_signif.l = FALSE,
        title.c = ifelse(!is.na(set.c), set.c, NA_character_),
        figure.c = figure.c
      )
    }

    return(invisible(x))
  }
)


# ttest, wilcoxon, limma, pearson, spearman: test
.twoSampCorTests <- function(data.mn, ## data (matrix of numerics; samples x variables)
                             samp.df, ## sample metadata (dataframe; samples x metadata)
                             feat.df, ## feature metadata (dataframe; features x metadata)
                             test.c,
                             factorNameC,
                             factorLevelsVc,
                             adjust.c,
                             adjust_thresh.n,
                             transformed.c,
                             title.c,
                             display_signif.l,
                             prefix.c,
                             figure.c) {
  checkLs <- .oneFactorCheck(
    samp.df = samp.df,
    test.c = test.c,
    factorNameC = factorNameC,
    factorLevelsVc = factorLevelsVc
  )

  factorVc <- checkLs[["factorVc"]]
  factorFc <- checkLs[["factorFc"]]
  factorVn <- checkLs[["factorVn"]]

  statVn <- .twoSampCorStat(
    data.mn = data.mn,
    test.c = test.c,
    factorFc = factorFc,
    factorVn = factorVn,
    transformed.c = transformed.c
  )

  pvalVn <- .twoSampCorPval(
    data.mn = data.mn,
    test.c = test.c,
    factorFc = factorFc,
    factorVn = factorVn
  )

  pvalAdjustVn <- stats::p.adjust(pvalVn, adjust.c)

  adjustSigniVn <- as.numeric(pvalAdjustVn <= adjust_thresh.n)

  metric.mn <- cbind(
    statVn,
    pvalAdjustVn,
    adjustSigniVn
  )

  rownames(metric.mn) <- colnames(data.mn)
  colnames(metric.mn) <- .twoSampCorNames(
    test.c,
    factorNameC,
    factorFc,
    adjust.c,
    prefix.c,
    transformed.c
  )

  for (colC in colnames(metric.mn)) {
    feat.df[, colC] <- metric.mn[, colC]
  }

  metric.mn <- metric.mn[order(pvalVn), ,
    drop = FALSE
  ]


  return(list(
    feat.df = feat.df,
    metric.mn = metric.mn
  ))
}

# ttest, wilcoxon, limma, pearson, spearman, anova, kruskal: checking factor name and levels
.oneFactorCheck <- function(samp.df,
                            test.c,
                            factorNameC,
                            factorLevelsVc) {
  factorLevelsSelectedVc <- factorLevelsVc

  factorNameLengthN <- length(factorNameC)
  if (factorNameLengthN != 1) {
    stop(
      "'length(factorNameC) = ", factorNameLengthN,
      "': A single 'factorName' is required"
    )
  }

  factorVcn <- samp.df[, factorNameC]
  # either of 'character' or 'numeric' mode at this stage

  if (test.c %in% c(
    "ttest", "limma", "wilcoxon", "anova", "kruskal",
    "anova2ways",
    "anova2waysInter",
    "limma2ways",
    "limma2waysInter"
  )) {
    factorVc <- factorVcn
    factorVn <- NULL

    if (is.factor(factorVc)) {
      factorFc <- factorVc
      factorVc <- as.character(factorVc)
    } else {
      factorModeC <- mode(factorVc)
      if (factorModeC != "character") {
        stop(
          "'mode(factorVc) = ", factorModeC,
          "': Factor must be of 'character' mode."
        )
      }
      factorFc <- factor(factorVc)
    }

    factorLevelsDefaultVc <- levels(factorFc)
    levelDotDashVl <- grepl(".", factorLevelsDefaultVc, fixed = TRUE) |
      grepl("-", factorLevelsDefaultVc, fixed = TRUE)
    if (sum(levelDotDashVl)) {
      stop("'Factor levels must not contain '.' or '-' characters.")
    }

    factorLevelsDefaultN <- nlevels(factorFc)

    if (test.c %in% c("ttest", "limma", "wilcoxon") &&
      factorLevelsDefaultN != 2) {
      stop(
        "'nlevels(factorFc) = ", factorLevelsDefaultN,
        "': The factor must have 2 levels."
      )
    } else if (test.c %in% "anova" && factorLevelsDefaultN < 3) {
      stop(
        "'nlevels(factorFc) = ", factorLevelsDefaultN,
        "': The factor must have 3 levels or more."
      )
    }

    if (length(factorLevelsSelectedVc) > 1 ||
      factorLevelsSelectedVc != "default") {
      factorLevelsSelectedN <- length(factorLevelsSelectedVc)
      if (factorLevelsDefaultN != factorLevelsSelectedN) {
        stop(
          "'nlevels(factorFc) = ", factorLevelsDefaultN,
          "'and 'length(factorLevelsVc) = ", factorLevelsSelectedN,
          "': 'factorLevelsVc' must have the same number
             of levels as the factor"
        )
      }

      if (!identical(sort(factorLevelsVc), sort(factorLevelsDefaultVc))) {
        level_default.c <- paste(factorLevelsDefaultVc, collapse = ", ")
        level_select.c <- paste(factorLevelsSelectedVc, collapse = ", ")
        stop(
          "'levels(factorFc) = ", level_default.c,
          "' and 'factorLevelsVc = ", level_select.c,
          "': refLeveslVc' must match all level names of the factor."
        )
      }

      factorFc <- factor(factorVc, factorLevelsSelectedVc)
    }
  } else if (test.c %in% c("pearson", "spearman")) {
    factorVn <- factorVcn
    factorFc <- factorVc <- NULL

    factorModeC <- mode(factorVn)
    if (factorModeC != "numeric") {
      stop(
        "'mode(factorVn) = ", factorModeC,
        "': Factor must be of 'numeric' mode."
      )
    }
  }

  return(list(
    factorVc = factorVc,
    factorFc = factorFc,
    factorVn = factorVn
  ))
}

.twoSampDiffMean <- function(data.mn,
                             factorFc) {
  apply(data.mn, 2, function(y) {
    diff(tapply(y, factorFc, function(x) mean(x, na.rm = TRUE)))
  })
}

.twoSampDiffMedian <- function(data.mn,
                               factorFc) {
  apply(data.mn, 2, function(y) {
    diff(tapply(y, factorFc, function(x) stats::median(x, na.rm = TRUE)))
  })
}

.twoSampRatioMean <- function(data.mn,
                              factorFc) {
  apply(data.mn, 2, function(y) {
    means.vc <- tapply(y, factorFc, function(x) mean(x, na.rm = TRUE))
    means.vc[2] / means.vc[1]
  })
}

.twoSampRatioMedian <- function(data.mn,
                                factorFc) {
  apply(data.mn, 2, function(y) {
    medians.vc <- tapply(y, factorFc, function(x) stats::median(x, na.rm = TRUE))
    medians.vc[2] / medians.vc[1]
  })
}

.cor <- function(data.mn,
                 factorVn,
                 methodC) {
  apply(data.mn, 2, function(y) {
    stats::cor(factorVn, y, method = methodC, use = "pairwise.complete.obs")
  })
}

.twoSampCorStat <- function(data.mn,
                            test.c,
                            factorFc,
                            factorVn,
                            transformed.c) {
  if (test.c %in% c("ttest", "limma")) {
    if (transformed.c == "none") {
      return(.twoSampRatioMean(data.mn = data.mn, factorFc = factorFc))
    } else {
      return(.twoSampDiffMean(data.mn = data.mn, factorFc = factorFc))
    }
  } else if (test.c == "wilcoxon") {
    if (transformed.c == "none") {
      return(.twoSampRatioMedian(data.mn = data.mn, factorFc = factorFc))
    } else {
      return(.twoSampDiffMedian(data.mn = data.mn, factorFc = factorFc))
    }
  } else if (test.c %in% c("pearson", "spearman")) {
    return(.cor(data.mn = data.mn, factorVn = factorVn, methodC = test.c))
  }
}

.twoSampCorPval <- function(data.mn,
                            test.c,
                            factorFc,
                            factorVn) {
  if (test.c %in% c("ttest", "wilcoxon", "pearson", "spearman")) {
    if (test.c == "ttest") {
      hypotest <- function(y) stats::t.test(y ~ factorFc)[["p.value"]]
    } else if (test.c == "wilcoxon") {
      hypotest <- function(y) {
        stats::wilcox.test(y ~ factorFc,
          exact = FALSE
        )[["p.value"]]
      }
    } else if (test.c %in% c("pearson", "spearman")) {
      hypotest <- function(y) {
        stats::cor.test(factorVn, y,
          method = test.c,
          use = "pairwise.complete.obs",
          exact = FALSE
        )[["p.value"]]
      }
    }

    return(apply(data.mn, 2, hypotest))
  } else if (test.c == "limma") {
    designMN <- cbind(
      level2 = rep(1, nrow(data.mn)),
      level1.level2 = as.numeric(as.numeric(factorFc) == 1)
    )
    rownames(designMN) <- rownames(data.mn)
    limmaFitLs <- limma::lmFit(t(data.mn), designMN)
    limmaBayesLs <- limma::eBayes(limmaFitLs)

    return(limmaBayesLs[["p.value"]][, "level1.level2"])
  }
}


.twoSampCorNames <- function(test.c,
                             factorNameC,
                             factorFc,
                             adjust.c,
                             prefix.c,
                             transformed.c) {
  prefix.c <- paste0(prefix.c, test.c, "_", make.names(factorNameC), "_")

  if (test.c %in% c(
    "ttest", "limma", "wilcoxon",
    "limma2ways", "limma2waysInter",
    "anova2ways", "anova2waysInter"
  )) {
    prefix.c <- paste0(
      prefix.c,
      paste(levels(factorFc), collapse = "."), "_"
    )
    statC <- ifelse(transformed.c == "none", "ratio", "diff")
  } else if (test.c %in% c("pearson", "spearman")) {
    statC <- "cor"
  } else {
    stop()
  }

  paste0(prefix.c, c(statC, adjust.c, "signif"))
}


.anovas <- function(data.mn, ## data (matrix of numerics; samples x variables)
                    samp.df, ## sample metadata (dataframe; samples x metadata)
                    feat.df, ## feature metadata (dataframe; features x metadata)
                    test.c,
                    factorNameC,
                    factorLevelsVc,
                    adjust.c,
                    adjust_thresh.n,
                    transformed.c,
                    title.c,
                    prefix.c,
                    figure.c) {
  checkLs <- .oneFactorCheck(
    samp.df = samp.df,
    test.c = test.c,
    factorNameC = factorNameC,
    factorLevelsVc = factorLevelsVc
  )

  factorVc <- checkLs[["factorVc"]]
  factorFc <- checkLs[["factorFc"]]

  # metric and pairwise names
  namesLs <- .anovasNames(
    test.c = test.c,
    factorNameC = factorNameC,
    factorFc = factorFc,
    adjust.c = adjust.c,
    prefix.c = prefix.c
  )
  pairNamesVc <- namesLs[["pairNamesVc"]]
  metricNamesVc <- namesLs[["aovNamesVc"]]

  ## omnibus and post-hoc tests

  diffPvalLs <- .anovasDiffPval(
    data.mn = data.mn,
    test.c = test.c,
    factorFc = factorFc,
    pairNamesVc = pairNamesVc
  )

  pvalVn <- diffPvalLs[["pvalVn"]]
  pairDiffMN <- diffPvalLs[["pairDiffMN"]]
  pairPvalMN <- diffPvalLs[["pairPvalMN"]]

  # main ANOVA test
  pvalAdjustVn <- stats::p.adjust(pvalVn, method = adjust.c)
  pvalSigniVi <- as.integer(pvalAdjustVn <= adjust_thresh.n)

  # pairwise test
  pairAdjustMN <- apply(
    pairPvalMN, 2,
    function(pairPvalVn) {
      stats::p.adjust(pairPvalVn,
        method = adjust.c
      )
    }
  )
  pairSigniMI <- pairAdjustMN <= adjust_thresh.n
  mode(pairSigniMI) <- "integer"

  metric.mn <- cbind(
    pvalAdjustVn,
    pvalSigniVi
  )
  for (pairI in seq_len(ncol(pairAdjustMN))) {
    metric.mn <- cbind(
      metric.mn,
      pairDiffMN[, pairI],
      pairAdjustMN[, pairI],
      pairSigniMI[, pairI]
    )
  }
  rownames(metric.mn) <- colnames(data.mn)
  colnames(metric.mn) <- metricNamesVc

  for (colC in colnames(metric.mn)) {
    feat.df[, colC] <- metric.mn[, colC]
  }

  metric.mn <- metric.mn[order(pvalVn), , drop = FALSE]

  ## graphic

  if (figure.c != "none") {
    if (!grepl("interactive", figure.c)) {
      grDevices::pdf(figure.c)
    }

    .anovasPlot(
      data.mn = data.mn,
      factorNameC = factorNameC,
      factorFc = factorFc,
      metric.mn = metric.mn,
      pairNamesVc = pairNamesVc,
      adjust.c = adjust.c,
      adjust_thresh.n = adjust_thresh.n,
      title.c = title.c,
      test.c = test.c,
      figure.c = figure.c
    )

    if (!grepl("interactive", figure.c)) {
      grDevices::dev.off()
    }
  }

  return(list(
    feat.df = feat.df,
    metric.mn = metric.mn
  ))
}

.anovasDiffPval <- function(data.mn,
                            test.c,
                            factorFc,
                            pairNamesVc) {
  if (test.c == "anova") {
    aovMN <- t(apply(
      data.mn, 2,
      function(varVn) {
        aovModel <- stats::aov(varVn ~ factorFc)
        aovPvalN <- summary(aovModel)[[1]][1, "Pr(>F)"]
        tukeyHsdMN <- stats::TukeyHSD(aovModel)[["factorFc"]]

        padjVn <- diffVn <- rep(NA_real_, length(pairNamesVc))
        names(padjVn) <- names(diffVn) <- vapply(pairNamesVc,
          function(pairC) {
            paste(rev(unlist(strsplit(pairC, ".", fixed = TRUE))),
              collapse = "-"
            )
          },
          FUN.VALUE = character(1)
        )

        diffVn[rownames(tukeyHsdMN)] <- tukeyHsdMN[, "diff"]
        padjVn[rownames(tukeyHsdMN)] <- tukeyHsdMN[, "p adj"]
        c(
          aovPval = aovPvalN,
          diffVn,
          padjVn
        )
      }
    ))

    # main ANOVA test
    pvalVn <- aovMN[, "aovPval"]

    # pairwise differences
    pairI <- as.integer(nlevels(factorFc) * (nlevels(factorFc) - 1) / 2)
    pairDiffMN <- aovMN[, 1 + seq_len(pairI), drop = FALSE]

    # pairwise test
    pairPvalMN <- aovMN[, (2 + pairI):ncol(aovMN), drop = FALSE]
  } else if (test.c == "kruskal") {
    kruskalMN <- t(apply(data.mn, 2, function(varVn) {
      kruskalPvalN <- stats::kruskal.test(varVn ~ factorFc)[["p.value"]]
      nemenyiPvalMN <- PMCMRplus::kwAllPairsNemenyiTest(
        varVn,
        factorFc,
        "Tukey"
      )[["p.value"]]

      nemenyiPvalVn <- c(nemenyiPvalMN[lower.tri(nemenyiPvalMN, diag = TRUE)])
      names(nemenyiPvalVn) <- paste0(
        rep(
          rownames(nemenyiPvalMN),
          ncol(nemenyiPvalMN)
        ),
        "-",
        rep(colnames(nemenyiPvalMN),
          each = nrow(nemenyiPvalMN)
        )
      )[c(lower.tri(nemenyiPvalMN, diag = TRUE))]
      padjVn <- rep(NA_real_, length(pairNamesVc))
      names(padjVn) <- vapply(pairNamesVc, function(pairC) {
        paste(rev(unlist(strsplit(pairC, ".", fixed = TRUE))), collapse = "-")
      }, FUN.VALUE = character(1))

      padjVn[names(nemenyiPvalVn)] <- nemenyiPvalVn

      c(
        kruskalPval = kruskalPvalN,
        padjVn
      )
    }))

    # main test
    pvalVn <- kruskalMN[, "kruskalPval"]

    ## difference of the medians for each pairwise comparison

    pairDiffMN <- vapply(pairNamesVc, function(pairC) {
      pairVc <- unlist(strsplit(pairC, ".", fixed = TRUE))
      pairSamplesVi <- which(factorFc %in% pairVc)
      pairFactorFc <- factor(as.character(factorFc)[pairSamplesVi],
        levels = pairVc
      )
      apply(
        data.mn[pairSamplesVi, ], 2,
        function(varVn) {
          diff(as.numeric(tapply(
            varVn, pairFactorFc,
            function(x) {
              stats::median(x, na.rm = TRUE)
            }
          )))
        }
      )
    }, numeric(length(pvalVn)))

    # pairwise test
    pairPvalMN <- kruskalMN[, -1, drop = FALSE]
  }

  return(list(
    pvalVn = pvalVn,
    pairDiffMN = pairDiffMN,
    pairPvalMN = pairPvalMN
  ))
}

.anovasNames <- function(test.c,
                         factorNameC,
                         factorFc,
                         adjust.c,
                         prefix.c) {
  prefix.c <- paste0(prefix.c, test.c, "_", make.names(factorNameC), "_")

  # getting the names of the pairwise comparisons
  factLevVc <- levels(factorFc)
  pairMC <- matrix("", nrow = length(factLevVc), ncol = length(factLevVc))
  for (i in seq_along(factLevVc)) {
    for (j in i:length(factLevVc)) {
      pairMC[i, j] <- paste0(factLevVc[i], ".", factLevVc[j])
    }
  }
  pairMC <- t(pairMC)
  pairNamesVc <- pairMC[lower.tri(pairMC)]

  aovNamesVc <- paste0(
    prefix.c,
    c(adjust.c, "signif")
  )

  for (pairC in pairNamesVc) {
    aovNamesVc <- c(
      aovNamesVc,
      paste0(prefix.c, pairC, c(
        "_diff", paste0("_", adjust.c),
        "_signif"
      ))
    )
  }

  return(list(
    pairNamesVc = pairNamesVc,
    aovNamesVc = aovNamesVc
  ))
}

.anovasPlot <- function(data.mn,
                        factorNameC,
                        factorFc,
                        metric.mn,
                        pairNamesVc,
                        adjust.c,
                        adjust_thresh.n,
                        title.c,
                        test.c,
                        figure.c) {
  for (pairC in pairNamesVc) {
    pairmetric.mn <- metric.mn[, grep(pairC, colnames(metric.mn), fixed = TRUE)]

    pairFactorFc <- factorFc[
      factorFc %in% unlist(strsplit(pairC, ".",
        fixed = TRUE
      )),
      drop = TRUE
    ]

    pairAdjPvalVn <- pairmetric.mn[, grep(adjust.c, colnames(pairmetric.mn))]
    pairSigniVl <- pairAdjPvalVn <= adjust_thresh.n
    pairSigniVl[is.na(pairSigniVl)] <- FALSE

    mainC <- paste0(
      ifelse(!is.na(title.c), paste0(title.c, "\n"), ""),
      test.c, " '", factorNameC, "' (", pairC, ") ",
      adjust.c, " (signif.: ",
      format(sum(pairSigniVl),
        big.mark = ","
      ),
      ")"
    )

    if (!is.null(pairFactorFc)) {
      group.vc <- levels(pairFactorFc)
    } else {
      group.vc <- ""
    }

    label.vc <- rownames(pairmetric.mn)
    label.vc[!pairSigniVl] <- ""

    gg_volcanoplot(
      fold_change.vn = pairmetric.mn[, 1],
      adjusted_pvalue.vn = pairmetric.mn[, 2],
      adjust_method.c = adjust.c,
      adjust_thresh.n = adjust_thresh.n,
      label.vc = label.vc,
      title.c = mainC,
      xlab.c = "Fold Change",
      class_name.vc = group.vc,
      figure.c = figure.c
    )
  }

  pvalAdjVn <- metric.mn[, paste0(test.c, "_", factorNameC, "_", adjust.c)]
  pvalSigniVi <- which(metric.mn[
    ,
    paste0(
      test.c, "_",
      factorNameC, "_signif"
    )
  ] > 0)

  if (sum(pvalSigniVi)) {
    for (varI in pvalSigniVi) {
      varNameC <- rownames(metric.mn)[varI]

      mainC <- paste0(
        varNameC, " (", adjust.c, " = ",
        signif(pvalAdjVn[varI], 2), ")"
      )
      if (!is.na(title.c)) {
        mainC <- paste0(title.c, ", ", mainC)
      }

      graphics::boxplot(data.mn[, varNameC] ~ factorFc,
        main = mainC,
        xlab = factorNameC, ylab = ""
      )
    }
  }
}

# ANOVA 2 ways without interaction (2 levels for each factor)
.anovas2ways <- function(data.mn, ## data (matrix of numerics; samples x variables)
                         samp.df, ## sample metadata (dataframe; samples x metadata)
                         feat.df, ## feature metadata (dataframe; features x metadata)
                         test.c,
                         factor_names.vc,
                         factor_levels.ls,
                         adjust.c,
                         adjust_thresh.n,
                         transformed.c,
                         title.c,
                         prefix.c,
                         figure.c) {
  checkLs <- .twoFactorsCheck(
    samp.df,
    test.c,
    factor_names.vc,
    factor_levels.ls,
    adjust.c,
    adjust_thresh.n
  )

  factor1.vc <- checkLs[["factor1.vc"]]
  factor1.fc <- checkLs[["factor1.fc"]]
  factor2.vc <- checkLs[["factor2.vc"]]
  factor2.fc <- checkLs[["factor2.fc"]]

  diffPvalLs <- .anovas2waysDiffPval(
    data.mn = data.mn,
    test.c = test.c,
    factor1.fc = factor1.fc,
    factor2.fc = factor2.fc
  )

  diffMN <- diffPvalLs[["diffMN"]]
  pvalMN <- diffPvalLs[["pvalMN"]]

  pvalAdjustMN <- apply(
    pvalMN, 2,
    function(pvalVn) {
      stats::p.adjust(pvalVn, adjust.c)
    }
  )

  adjustSigniMN <- apply(
    pvalAdjustMN, 2,
    function(adjustVn) {
      as.numeric(adjustVn <= adjust_thresh.n)
    }
  )

  metric.mn <- NULL

  for (factorI in seq_len(2)) {
    metric.mn <- cbind(
      metric.mn,
      diffMN[, factorI],
      pvalAdjustMN[, factorI],
      adjustSigniMN[, factorI]
    )
  }

  if (grepl("Inter$", test.c)) {
    metric.mn <- cbind(metric.mn, pvalAdjustMN[, 3], adjustSigniMN[, 3])
  }

  rownames(metric.mn) <- colnames(data.mn)
  colnames(metric.mn) <- .anovas2waysNames(
    test.c = test.c,
    factor_names.vc = factor_names.vc,
    factor1.fc = factor1.fc,
    factor2.fc = factor2.fc,
    adjust.c = adjust.c,
    prefix.c = prefix.c,
    transformed.c = transformed.c
  )

  for (colC in colnames(metric.mn)) {
    feat.df[, colC] <- metric.mn[, colC]
  }

  metric.mn <- metric.mn[
    do.call(
      "order",
      as.list(as.data.frame(pvalAdjustMN))
    ), ,
    drop = FALSE
  ]

  ## graphic

  if (figure.c != "none") {
    .anovas2waysPlot(
      data.mn = data.mn,
      test.c = test.c,
      factor_names.vc = factor_names.vc,
      adjust.c = adjust.c,
      adjust_thresh.n = adjust_thresh.n,
      metric.mn = metric.mn,
      title.c = title.c,
      figure.c = figure.c
    )
  }

  return(list(
    feat.df = feat.df,
    metric.mn = metric.mn
  ))
}

# anova2ways: checking factor names and levels
.twoFactorsCheck <- function(samp.df,
                             test.c,
                             factor_names.vc,
                             factor_levels.ls,
                             adjust.c,
                             adjust_thresh.n) {
  if (length(factor_names.vc) != 2) {
    stop("Two factors are required.")
  }

  if (!is.list(factor_levels.ls) || length(factor_levels.ls) != 2) {
    stop("'factor_levels.ls' must be a list of length 2.", call. = FALSE)
  }

  factorLs <- vector(mode = "list", length = 2)

  for (i in seq_along(factor_names.vc)) {
    factorLs[[i]] <- .oneFactorCheck(
      samp.df = samp.df,
      test.c = test.c,
      factorNameC = factor_names.vc[i],
      factorLevelsVc = factor_levels.ls[[i]]
    )
  }

  factor1.vc <- factorLs[[1]][["factorVc"]]
  factor1.fc <- factorLs[[1]][["factorFc"]]
  factor2.vc <- factorLs[[2]][["factorVc"]]
  factor2.fc <- factorLs[[2]][["factorFc"]]

  if (nlevels(factor1.fc) != 2 || nlevels(factor2.fc) != 2) {
    stop("Only factors with 2 levels are handled by the current implementation
         of '2ways' tests")
  }

  return(list(
    factor1.vc = factor1.vc,
    factor1.fc = factor1.fc,
    factor2.vc = factor2.vc,
    factor2.fc = factor2.fc
  ))
}

.anovas2waysDiffPval <- function(data.mn,
                                 test.c,
                                 factor1.fc,
                                 factor2.fc) {
  # Difference of means

  anovas2waysDiffMN <- cbind(
    .twoSampDiffMean(
      data.mn = data.mn,
      factorFc = factor1.fc
    ),
    .twoSampDiffMean(
      data.mn = data.mn,
      factorFc = factor2.fc
    )
  )

  # Test statistic

  if (test.c == "anova2ways") {
    anova2waysPvalMN <- t(apply(
      data.mn, 2,
      function(varVn) {
        aovModel <- stats::lm(varVn ~ factor1.fc + factor2.fc)
        stats::anova(aovModel)[c(
          "factor1.fc",
          "factor2.fc"
        ), "Pr(>F)"]
      }
    ))
  } else if (test.c == "anova2waysInter") {
    anova2waysPvalMN <- t(apply(
      data.mn, 2,
      function(varVn) {
        aovModel <- stats::lm(varVn ~ factor1.fc * factor2.fc)
        stats::anova(aovModel)[c(
          "factor1.fc",
          "factor2.fc",
          "factor1.fc:factor2.fc"
        ), "Pr(>F)"]
      }
    ))
  } else if (test.c == "limma2ways") {
    # create design
    designMN <- stats::model.matrix(~ factor1.fc + factor2.fc)
    rownames(designMN) <- rownames(data.mn)

    # apply limma test
    limmaFitLs <- limma::lmFit(t(data.mn), designMN)
    limmaBayesLs <- limma::eBayes(limmaFitLs)

    anova2waysPvalMN <- limmaBayesLs[["p.value"]][, 2:3]
  } else if (test.c == "limma2waysInter") {
    # create the design
    designMN <- stats::model.matrix(~ 0 + factor1.fc:factor2.fc)
    colnames(designMN) <- make.names(gsub(
      "factor1.fc", "",
      gsub(
        "factor2.fc", "",
        colnames(designMN)
      )
    ))
    rownames(designMN) <- rownames(data.mn)

    # apply linear model
    limmaFit <- limma::lmFit(t(data.mn), designMN)

    # create contrast matrix
    factInterVc <- colnames(designMN)
    fact1AllC <- paste(paste(factInterVc[c(2, 4)], collapse = "+"),
      paste(factInterVc[c(1, 3)], collapse = "-"),
      sep = "-"
    )
    fact2AllC <- paste(paste(factInterVc[3:4], collapse = "+"),
      paste(factInterVc[seq_len(2)], collapse = "-"),
      sep = "-"
    )
    factInterC <- paste(paste(factInterVc[c(1, 4)], collapse = "+"),
      paste(factInterVc[c(2, 3)], collapse = "-"),
      sep = "-"
    )

    contrC <- c(fact1AllC, fact2AllC, factInterC)

    contrastMN <- limma::makeContrasts(
      contrasts = contrC,
      levels = factInterVc
    )

    # apply the contrast in linear model
    contrastFitLs <- limma::contrasts.fit(limmaFit, contrastMN)
    limmaBayesLs <- limma::eBayes(contrastFitLs)

    # recover signifcant variables
    anova2waysPvalMN <- limmaBayesLs[["p.value"]]
  }

  return(list(
    diffMN = anovas2waysDiffMN,
    pvalMN = anova2waysPvalMN
  ))
}

.anovas2waysNames <- function(test.c,
                              factor_names.vc,
                              factor1.fc,
                              factor2.fc,
                              adjust.c,
                              prefix.c,
                              transformed.c) {
  anova2waysNamesC <- c(
    .twoSampCorNames(
      test.c = test.c,
      factorNameC = factor_names.vc[1],
      factorFc = factor1.fc,
      adjust.c = adjust.c,
      prefix.c = prefix.c,
      transformed.c = transformed.c
    ),
    .twoSampCorNames(
      test.c = test.c,
      factorNameC = factor_names.vc[2],
      factorFc = factor2.fc,
      adjust.c = adjust.c,
      prefix.c = prefix.c,
      transformed.c = transformed.c
    )
  )

  if (grepl("Inter$", test.c)) {
    anova2waysNamesC <- c(
      anova2waysNamesC,
      paste0(
        paste0(
          test.c, "_",
          factor_names.vc[1], ":",
          factor_names.vc[2],
          "_"
        ),
        c(adjust.c, "signif")
      )
    )
  }

  anova2waysNamesC
}

.anovas2waysPlot <- function(data.mn,
                             test.c,
                             factor_names.vc,
                             adjust.c,
                             adjust_thresh.n,
                             metric.mn,
                             title.c,
                             figure.c) {
  metricSigniMN <- metric.mn[, grep("_signif$", colnames(metric.mn))]
  metricSigniVi <- which(rowSums(metricSigniMN) > 0)

  if (sum(metricSigniVi)) {
    metricSigniLs <- apply(
      metricSigniMN, 2,
      function(signiVn) which(signiVn > 0)
    )

    metricSigniNamesVc <- factor_names.vc
    if (grepl("Inter$", test.c)) {
      metricSigniNamesVc <- c(metricSigniNamesVc, paste(factor_names.vc,
        collapse = ":"
      ))
    }
    names(metricSigniLs) <- metricSigniNamesVc

    mainC <- test.c
    if (!is.na(title.c)) {
      mainC <- paste0(title.c, "\n", mainC)
    }

    if (!grepl("interactive", figure.c)) {
      file.tiffC <- gsub(".pdf", ".tiff", figure.c, fixed = TRUE)
    } else {
      file.tiffC <- "none"
    }

    vennplot(
      input.ls = metricSigniLs,
      title.c = mainC,
      figure.c = file.tiffC
    )
  }
}

# Print significant features
.signifPrint <- function(metric.mn,
                         signiVl,
                         signif_maxprint.i,
                         adjust.c,
                         adjust_thresh.n) {
  metricSigniMN <- metric.mn[signiVl, , drop = FALSE]

  if (is.na(signif_maxprint.i)) {
    signif_maxprint.i <- nrow(metricSigniMN)
  }

  message_thefirst.c <- paste0("The first ", signif_maxprint.i, " are")
  message(
    "\n", nrow(metricSigniMN), " variable",
    ifelse(nrow(metricSigniMN) > 1, "s", ""),
    " (", round(nrow(metricSigniMN) / nrow(metric.mn) * 100), "%) ",
    ifelse(nrow(metricSigniMN) > 1, "were", "was"),
    " found significant at the ", adjust_thresh.n,
    " level (after the '", adjust.c, "' correction).\n",
    ifelse(nrow(metricSigniMN) == 1,
      "It is",
      ifelse(nrow(metricSigniMN) <= signif_maxprint.i,
        "They are",
        message_thefirst.c
      )
    ),
    " displayed below",
    ifelse(nrow(metricSigniMN) > 1,
      " (sorted by increasing corrected p-values)",
      ""
    ),
    ":\n"
  )
  print(metricSigniMN[seq_len(min(signif_maxprint.i, nrow(metricSigniMN))), ,
    drop = FALSE
  ])
}


#### plot_hypotesting (SummarizedExperiment) ####

#' @rdname plot_hypotesting
#' @export
setMethod(
  "plot_hypotesting", signature(x = "SummarizedExperiment"),
  function(x,
           test.c = c(
             "ttest", "limma", "wilcoxon", "pearson", "spearman"
           )[1],
           factor.c,
           adjust.c = c(
             "holm",
             "hochberg",
             "hommel",
             "bonferroni",
             "BH",
             "BY",
             "fdr",
             "none"
           )[5],
           adjust_thresh.n = 0.05,
           title.c = NA_character_,
           class_color.vc = "",
           display_signif.l = TRUE,
           theme.c = c(
             "default",
             "bw",
             "classic",
             "dark",
             "gray",
             "linedraw",
             "light",
             "minimal",
             "void"
           )[2],
           size.ls = list(
             class.i = 5,
             lab.i = 16,
             point.i = 3,
             tick.i = 14,
             title.i = 20
           ),
           figure.c = c(
             "none",
             "interactive",
             "interactive_plotly",
             "myfile.pdf"
           )[2]) {
    if (!(test.c %in% c("ttest", "limma", "wilcoxon", "pearson", "spearman"))) {
      stop("'plot_hypotesting' method currently implemented only for the 'ttest', 'limma', 'wilcoxon', 'pearson', and 'spearman' tests.")
    }

    stopifnot(factor.c %in% colnames(colData(x)))
    factor.vcn <- colData(x)[, factor.c]

    test_colnames.vc <- grep(paste0(test.c, "_", factor.c), colnames(rowData(x)), value = TRUE)

    stopifnot(length(test_colnames.vc) == 3)

    stat.c <- tail(unlist(strsplit(test_colnames.vc[1], split = "_")), 1)

    stopifnot(stat.c %in% c("ratio", "diff", "cor"))

    fc_colname.c <- test_colnames.vc[1]

    adj_colname.c <- test_colnames.vc[2]

    stopifnot(tail(unlist(strsplit(adj_colname.c, split = "_")), 1) == adjust.c)

    fc.vn <- rowData(x)[, fc_colname.c]

    adj.vn <- rowData(x)[, adj_colname.c]
    signif.vl <- adj.vn <= adjust_thresh.n
    signif.vl[is.na(signif.vl)] <- FALSE

    main.c <- paste0(
      ifelse(!is.na(title.c), paste0(title.c, "\n"), ""),
      test.c, " '", factor.c, "', ",
      adjust.c, " (signif.: ",
      format(sum(signif.vl),
        big.mark = ","
      ),
      ")"
    )

    label.vc <- rownames(x)
    label.vc[!signif.vl] <- ""

    if (sum(signif.vl) > 30) {
      label.vc[adj.vn > max(head(sort(adj.vn), 30))] <- ""
    }

    figure_volcano.c <- figure.c
    if (grepl(".pdf$", figure_volcano.c)) {
      figure_volcano.c <- paste0(gsub(".pdf$", "", figure_volcano.c), "_", test.c, "_", factor.c, "_volcano.pdf")
    }

    gg_volcanoplot(
      fold_change.vn = fc.vn,
      adjusted_pvalue.vn = adj.vn,
      adjust_method.c = adjust.c,
      adjust_thresh.n = adjust_thresh.n,
      label.vc = label.vc,
      title.c = main.c,
      xlab.c = ifelse(stat.c == "cor", "Correlation", paste0("Fold Change (", ifelse(stat.c == "ratio", "ratio", "difference"), " between ", ifelse(test.c == "wilcoxon", "medians", "means"), ")")),
      class_name.vc = unlist(strsplit(unlist(strsplit(fc_colname.c, split = "_"))[3], split = ".", fixed = TRUE)),
      class_color.vc = class_color.vc,
      theme.c = theme.c,
      size.ls = size.ls,
      figure.c = figure_volcano.c
    )

    if (display_signif.l && sum(signif.vl) > 0 && (figure.c == "interactive" || grepl(".pdf$", figure.c))) {
      var_signif.vi <- which(signif.vl)

      if (test.c %in% c("pearson", "spearman")) {
        for (var.i in var_signif.vi) {
          var.c <- rownames(x)[var.i]

          main.c <- paste0(
            var.c, "\n", stat.c, " = ",
            signif(fc.vn[var.i], 2),
            ", ",
            adjust.c, " = ",
            signif(adj.vn[var.i], 2)
          )
          if (!is.na(title.c)) {
            main.c <- paste0(title.c, ", ", main.c)
          }

          figure_corplot.c <- figure.c
          if (grepl(".pdf$", figure.c)) {
            figure_corplot.c <- paste0(gsub(".pdf$", "", figure.c), "_", test.c, "_", factor.c, "_", make.names(var.c), "_corplot.pdf")
          }

          gg_corplot(
            data.frame(
              independent = factor.vcn,
              dependent = assay(x)[var.c, ]
            ),
            x.c = "independent",
            y.c = "dependent",
            title.c = main.c,
            xlab.c = factor.c,
            ylab.c = "intensity",
            label.vc = colnames(x),
            theme.c = theme.c,
            size.ls = size.ls,
            figure.c = figure_corplot.c
          )
        }
      } else { # "ttest", "limma", "wilcoxon"

        for (var.i in var_signif.vi) {
          var.c <- rownames(x)[var.i]

          main.c <- paste0(
            var.c, "\n(", stat.c, " ", ifelse(test.c == "wilcoxon", "medians", "means"), " = ",
            signif(fc.vn[var.i], 2),
            ", ",
            adjust.c, " = ",
            signif(adj.vn[var.i], 2), ")"
          )
          if (!is.na(title.c)) {
            main.c <- paste0(title.c, ", ", main.c)
          }

          figure_boxplot.c <- figure.c
          if (grepl(".pdf$", figure.c)) {
            figure_boxplot.c <- paste0(gsub(".pdf$", "", figure.c), "_", test.c, "_", factor.c, "_", make.names(var.c), "_boxplot.pdf")
          }

          gg_boxplot(
            data.frame(
              factor = factor.vcn,
              response = assay(x)[var.c, ]
            ),
            x.c = "factor",
            y.c = "response",
            color.c = "factor",
            title.c = main.c,
            xlab.c = factor.c,
            label.vc = colnames(x),
            theme.c = theme.c,
            size.ls = size.ls,
            figure.c = figure_boxplot.c
          )
        }
      }
    }
  }
)
