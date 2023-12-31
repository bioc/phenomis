#### clustering (MultiAssayExperiment) ####

#' @rdname clustering
#' @export
setMethod("clustering", signature(x = "MultiAssayExperiment"),
          function(x,
                   dissym.c = c("euclidean",
                                "maximum",
                                "manhattan",
                                "canberra",
                                "binary",
                                "minkowski",
                                "1-cor",
                                "1-abs(cor)")[7],
                   correl.c = c("pearson",
                                "kendall",
                                "spearman")[1],
                   agglo.c = c("ward.D",
                               "ward.D2",
                               "single",
                               "complete",
                               "average",
                               "mcquitty",
                               "median",
                               "centroid")[2],
                   clusters.vi = c(2, 2),
                   cex.vn = c(1, 1),
                   palette.c = c("blueOrangeRed",
                                 "redBlackGreen")[1],
                   scale_plot.l = TRUE,
                   title.c = NA,
                   figure.c = c("none", "interactive", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            report_set.c <- report.c
            if (report_set.c != "none")
              report_set.c <- "interactive"
            
            if (!(figure.c %in% c("none", "interactive")))
              grDevices::pdf(figure.c)
            
            figure_set.c <- figure.c
            if (figure_set.c != "none")
              figure_set.c <- "interactive"
            
            for (setC in names(x)) {
              
              if (report.c != "none")
                message("Performing hierarchical clustering on the '", setC,
                        "' dataset...")
              
              x[[setC]] <- clustering(x = x[[setC]],
                                      dissym.c = dissym.c,
                                      correl.c = correl.c,
                                      agglo.c = agglo.c,
                                      clusters.vi = clusters.vi,
                                      cex.vn = cex.vn,
                                      palette.c = palette.c,
                                      scale_plot.l = scale_plot.l,
                                      title.c = paste0("[", setC, "]"),
                                      figure.c = figure_set.c,
                                      report.c = report_set.c)
              
            }
            
            if (!(figure.c %in% c("none", "interactive")))
              grDevices::dev.off()
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })


#### clustering (SummarizedExperiment) ####

#' @rdname clustering
#' @export
setMethod("clustering", signature(x = "SummarizedExperiment"),
          function(x,
                   dissym.c = c("euclidean",
                                "maximum",
                                "manhattan",
                                "canberra",
                                "binary",
                                "minkowski",
                                "1-cor",
                                "1-abs(cor)")[7],
                   correl.c = c("pearson",
                                "kendall",
                                "spearman")[1],
                   agglo.c = c("ward.D",
                               "ward.D2",
                               "single",
                               "complete",
                               "average",
                               "mcquitty",
                               "median",
                               "centroid")[2],
                   clusters.vi = c(2, 2),
                   cex.vn = c(1, 1),
                   palette.c = c("blueOrangeRed",
                                 "redBlackGreen")[1],
                   scale_plot.l = TRUE,
                   title.c = NA,
                   figure.c = c("none", "interactive", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            metadata.ls <- .clustering(pro.mn = t(SummarizedExperiment::assay(x)),
                                       sam.Df = SummarizedExperiment::colData(x),
                                       var.Df = SummarizedExperiment::rowData(x),
                                       dissym.c = dissym.c,
                                       correl.c = correl.c,
                                       agglo.c = agglo.c,
                                       clusters.vi = clusters.vi,
                                       cex.vn = cex.vn,
                                       palette.c = palette.c,
                                       scale_plot.l = scale_plot.l,
                                       title.c = title.c,
                                       figure.c = figure.c)
            
            SummarizedExperiment::colData(x) <- metadata.ls[["sam.Df"]]
            SummarizedExperiment::rowData(x) <- metadata.ls[["var.Df"]] 
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })


#### clustering (MultiDataSet) ####

#' @rdname clustering
#' @export
setMethod("clustering", signature(x = "MultiDataSet"),
          function(x,
                   dissym.c = c("euclidean",
                                "maximum",
                                "manhattan",
                                "canberra",
                                "binary",
                                "minkowski",
                                "1-cor",
                                "1-abs(cor)")[7],
                   correl.c = c("pearson",
                                "kendall",
                                "spearman")[1],
                   agglo.c = c("ward.D",
                               "ward.D2",
                               "single",
                               "complete",
                               "average",
                               "mcquitty",
                               "median",
                               "centroid")[2],
                   clusters.vi = c(2, 2),
                   cex.vn = c(1, 1),
                   palette.c = c("blueOrangeRed",
                                 "redBlackGreen")[1],
                   scale_plot.l = TRUE,
                   title.c = NA,
                   figure.c = c("none", "interactive", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            report_set.c <- report.c
            if (report_set.c != "none")
              report_set.c <- "interactive"
            
            if (!(figure.c %in% c("none", "interactive")))
              grDevices::pdf(figure.c)
            
            figure_set.c <- figure.c
            if (figure_set.c != "none")
              figure_set.c <- "interactive"
            
            for (setC in names(x)) {
              
              if (report.c != "none")
                message("Performing hierarchical clustering on the '", setC,
                        "' dataset...")
              
              ese <- x[[setC]]
              
              ese <- clustering(x = ese,
                                dissym.c = dissym.c,
                                correl.c = correl.c,
                                agglo.c = agglo.c,
                                clusters.vi = clusters.vi,
                                cex.vn = cex.vn,
                                palette.c = palette.c,
                                scale_plot.l = scale_plot.l,
                                title.c = paste0("[", setC, "]"),
                                figure.c = figure_set.c,
                                report.c = report_set.c)
              
              x <- MultiDataSet::add_eset(x,
                                          ese,
                                          dataset.type = setC,
                                          GRanges = NA,
                                          overwrite = TRUE,
                                          warnings = FALSE)
              
            }
            
            if (!(figure.c %in% c("none", "interactive")))
              grDevices::dev.off()
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })


#### clustering (ExpressionSet) ####

#' @rdname clustering
#' @export
setMethod("clustering", signature(x = "ExpressionSet"),
          function(x,
                   dissym.c = c("euclidean",
                                "maximum",
                                "manhattan",
                                "canberra",
                                "binary",
                                "minkowski",
                                "1-cor",
                                "1-abs(cor)")[7],
                   correl.c = c("pearson",
                                "kendall",
                                "spearman")[1],
                   agglo.c = c("ward.D",
                               "ward.D2",
                               "single",
                               "complete",
                               "average",
                               "mcquitty",
                               "median",
                               "centroid")[2],
                   clusters.vi = c(2, 2),
                   cex.vn = c(1, 1),
                   palette.c = c("blueOrangeRed",
                                 "redBlackGreen")[1],
                   scale_plot.l = TRUE,
                   title.c = NA,
                   figure.c = c("none", "interactive", "myfile.pdf")[2],
                   report.c = c("none", "interactive", "myfile.txt")[2]) {
            
            if (!(report.c %in% c("none", "interactive")))
              sink(report.c, append = TRUE)
            
            metadata.ls <- .clustering(pro.mn = t(Biobase::exprs(x)),
                                       sam.Df = Biobase::pData(x),
                                       var.Df = Biobase::fData(x),
                                       dissym.c = dissym.c,
                                       correl.c = correl.c,
                                       agglo.c = agglo.c,
                                       clusters.vi = clusters.vi,
                                       cex.vn = cex.vn,
                                       palette.c = palette.c,
                                       scale_plot.l = scale_plot.l,
                                       title.c = title.c,
                                       figure.c = figure.c)
            
            Biobase::pData(x) <- metadata.ls[["sam.Df"]]
            Biobase::fData(x) <- metadata.ls[["var.Df"]]
            
            if (!(report.c %in% c("none", "interactive")))
              sink()
            
            methods::validObject(x)
            
            return(invisible(x))
            
          })

.clustering <- function(pro.mn, ## profiles (matrix; samples x variables)
                        sam.Df, ## sample metadata (data frame or DataFrame; samples x sample metadata)
                        var.Df, ## variable metadata (data frame or DataFrame; variables x variable metadata)
                        dissym.c,    ## dissimilarity
                        correl.c, ## correlation method
                        agglo.c, ## agglomeration method
                        clusters.vi,
                        cex.vn,
                        palette.c,    ## color scale
                        scale_plot.l,
                        title.c,
                        figure.c) {
  
  ncaN <- 14 ## Sample and variable name truncature for display
  
  hcl.ls <- .cluster_compute(pro.mn = pro.mn,
                             dissym.c = dissym.c,
                             correl.c = correl.c,
                             agglo.c = agglo.c)
  
  obsHcl <- hcl.ls[["obs"]]
  feaHcl <- hcl.ls[["fea"]]
  
  heaMN <- pro.mn[obsHcl[["order"]], feaHcl[["order"]]]
  
  ## figure
  
  if (figure.c != "none") {
    
    if (figure.c != "interactive")
      grDevices::pdf(figure.c)
    
    .cluster_plot(heaMN = heaMN,
                  scale_plot.l = scale_plot.l,
                  palette.c = palette.c,
                  obsHcl = obsHcl,
                  feaHcl = feaHcl,
                  clusters.vi = clusters.vi,
                  cex.vn = cex.vn,
                  ncaN = ncaN,
                  title.c = title.c)
    
    if (figure.c != "interactive")
      grDevices::dev.off()
    
  }
  
  ## Returning
  
  if (clusters.vi[1] > 1) ## number of sample clusters
    sam.Df[, "hclust"] <- stats::cutree(obsHcl, k = clusters.vi[1])
  
  if (clusters.vi[2] > 1) ## number of variable clusters
    var.Df[, "hclust"] <- stats::cutree(feaHcl, k = clusters.vi[2])
  
  return(list(sam.Df = sam.Df,
              var.Df = var.Df))
  
}


.cluster_compute <- function(pro.mn,
                             dissym.c,
                             correl.c,
                             agglo.c) {
  
  if (dissym.c %in% c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) {
    
    obsHcl <- stats::hclust(stats::dist(pro.mn, method = dissym.c),
                            method = agglo.c)
    
    feaHcl <- stats::hclust(stats::dist(t(pro.mn), method = dissym.c),
                            method = agglo.c)
    
  } else if (dissym.c == "1-cor") {
    
    obsHcl <- stats::hclust(stats::as.dist(1 - stats::cor(t(pro.mn),
                                                          method = correl.c,
                                                          use = "pairwise.complete.obs")),
                            method = agglo.c)
    
    feaHcl <- stats::hclust(stats::as.dist(1 - stats::cor(pro.mn,
                                                          method = correl.c,
                                                          use = "pairwise.complete.obs")),
                            method = agglo.c)
    
  } else if (dissym.c == "1-abs(cor)") {
    
    obsHcl <- stats::hclust(stats::as.dist(1 - abs(stats::cor(t(pro.mn),
                                                              method = correl.c,
                                                              use = "pairwise.complete.obs"))),
                            method = agglo.c)
    
    feaHcl <- stats::hclust(stats::as.dist(1 - abs(stats::cor(pro.mn,
                                                              method = correl.c,
                                                              use = "pairwise.complete.obs"))),
                            method = agglo.c)
    
  }
  
  return(list(obs = obsHcl,
              fea = feaHcl))
  
}


.cluster_plot <- function(heaMN,
                          scale_plot.l,
                          palette.c,
                          feaHcl,
                          obsHcl,
                          clusters.vi,
                          cex.vn,
                          ncaN,
                          title.c) {
  
  if (scale_plot.l)
    heaMN <- scale(heaMN)
  
  heaMN <- heaMN[, rev(seq_len(ncol(heaMN))), drop = FALSE]
  
  graphics::layout(matrix(seq_len(4), nrow = 2),
                   widths = c(1, 4), heights = c(1, 4))
  
  ## Color palette
  
  switch(palette.c,
         blueOrangeRed = {
           imaPalVn <- grDevices::colorRampPalette(c("blue", "orange", "red"),
                                                   space = "rgb")(5)[seq_len(5)]
         },
         redBlackGreen = {
           imaPalVn <- grDevices::colorRampPalette(c("red", "black", "green"),
                                                   space = "rgb")(5)[seq_len(5)]
         })
  
  ## Color scale
  
  scaN <- length(imaPalVn)
  
  graphics::par(mar = c(0.6, 0.6, 0.6, 4.1))
  
  ylimVn <- c(0, scaN)
  ybottomVn <- 0:(scaN - 1)
  ytopVn <- seq_len(scaN)
  
  graphics::plot(x = 0,
                 y = 0,
                 bty = "n",
                 font.axis = 2,
                 font.lab = 2,
                 type = "n",
                 xlim = c(0, 1),
                 ylim = ylimVn,
                 xlab = "",
                 ylab = "",
                 xaxs = "i",
                 yaxs = "i",
                 xaxt = "n",
                 yaxt = "n")
  
  graphics::rect(xleft = 0.8,
                 ybottom = ybottomVn,
                 xright = 1,
                 ytop = ytopVn,
                 col = imaPalVn,
                 border = NA)
  
  prtVn <- pretty(range(heaMN, na.rm = TRUE))
  graphics::axis(at = scaN / diff(range(prtVn)) * (prtVn - min(prtVn)),
                 font = 2,
                 font.axis = 2,
                 labels = prtVn,
                 las = 1,
                 lwd = 2,
                 lwd.ticks = 2,
                 side = 4,
                 xpd = TRUE)
  
  graphics::arrows(graphics::par("usr")[2],
                   graphics::par("usr")[4],
                   graphics::par("usr")[2],
                   graphics::par("usr")[3],
                   code = 0,
                   lwd = 2,
                   xpd = TRUE)
  
  ## Feature dendrogram
  
  graphics::par(mar = c(4.1, 0.6, 0, 0.1),
                lwd = 2)
  
  graphics::plot(rev(stats::as.dendrogram(feaHcl)), horiz = TRUE,
                 leaflab = "none",
                 main = "", xaxs = "i", yaxs = "i",
                 xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  
  revFeaHcl <- list(merge = cbind(feaHcl[["merge"]][, 2],
                                  feaHcl[["merge"]][, 1]),
                    height = feaHcl[["height"]],
                    order = rev(feaHcl[["order"]]),
                    labels = feaHcl[["labels"]])
  
  if (clusters.vi[2] > 1) {
    cluFeaVn <- stats::cutree(revFeaHcl,
                              k = clusters.vi[2])[revFeaHcl[["order"]]]
    cutFeaVn <- which(abs(diff(cluFeaVn)) > 0)
    cutFeaTxtVn <- c(cutFeaVn[1] / 2,
                     cutFeaVn + diff(c(cutFeaVn, length(cluFeaVn))) / 2) + 0.5
    cutFeaLinVn <- cutFeaVn + 0.5
    graphics::text(graphics::par("usr")[1] + 0.2 * diff(graphics::par("usr")[seq_len(2)]),
                   cutFeaTxtVn,
                   labels = unique(cluFeaVn),
                   cex = 2,
                   font = 2,
                   las = 2)
  }
  
  ## Observation dendrogram
  
  graphics::par(mar = c(0.1, 0, 0.6, 4.1),
                lwd = 2)
  
  graphics::plot(stats::as.dendrogram(obsHcl), leaflab = "none",
                 main = "", xaxs = "i", yaxs = "i",
                 yaxt = "n", xlab = "", ylab = "")
  
  if (clusters.vi[1] > 1) {
    cluObsVn <- stats::cutree(obsHcl, k = clusters.vi[1])[obsHcl[["order"]]]
    cutObsVn <- which(abs(diff(cluObsVn)) > 0)
    cutObsTxtVn <- c(cutObsVn[1] / 2,
                     cutObsVn + diff(c(cutObsVn, length(cluObsVn))) / 2) + 0.5
    cutObsLinVn <- cutObsVn + 0.5
    graphics::text(cutObsTxtVn,
                   0.8 * graphics::par("usr")[4],
                   labels =  unique(cluObsVn),
                   cex = 2,
                   font = 2)
  }
  
  ## Heatmap
  
  graphics::par(mar = c(4.1, 0, 0, 4.1))
  
  graphics::image(x = seq_len(nrow(heaMN)),
                  y = seq_len(ncol(heaMN)),
                  z = round(heaMN),
                  col = imaPalVn,
                  font.axis = 2,
                  font.lab = 2,
                  xaxt = "n",
                  yaxt = "n",
                  xlab = "",
                  ylab = "")
  
  obsOrdVc <- obsHcl[["labels"]][obsHcl[["order"]]]
  obsOrdLenVn <- vapply(obsOrdVc, nchar, FUN.VALUE = integer(1))
  obsOrdVc <- substr(obsOrdVc, 1, ncaN)
  obsOrdVc <- paste0(obsOrdVc, ifelse(obsOrdLenVn > ncaN, ".", ""), " ")
  
  graphics::mtext(obsOrdVc,
                  at = seq_len(nrow(heaMN)),
                  cex = cex.vn[1],
                  las = 2,
                  side = 1)
  
  feaOrdVc <- feaHcl[["labels"]][feaHcl[["order"]]]
  feaOrdLenVn <- vapply(feaOrdVc, nchar, FUN.VALUE = integer(1))
  feaOrdVc <- substr(feaOrdVc, 1, ncaN)
  feaOrdVc <- paste0(" ", feaOrdVc, ifelse(feaOrdLenVn > ncaN, ".", ""))
  
  graphics::mtext(feaOrdVc,
                  at = ncol(heaMN):1,
                  cex = cex.vn[2],
                  las = 2,
                  side = 4)
  
  if (clusters.vi[2] > 1)
    graphics::abline(h = cutFeaLinVn)
  if (clusters.vi[1] > 1)
    graphics::abline(v = cutObsLinVn)
  
  graphics::box()
  
  if (!is.na(title.c))
    graphics::title(paste0(title.c, " "), line = -1, adj = 1, outer = TRUE)
  
  
  
}
