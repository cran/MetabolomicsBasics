#' @title ClassificationHistogram.
#' @description \code{ClassificationHistogram} will plot the results of \link{ClassificationWrapper}.
#' @details No further details.
#' @param out_classific Output of \link{ClassificationWrapper}.
#' @param breaks Breaks for histogram.
#' @param ... Passed on to \code{par}. Useful to adjust \code{cex}.
#' @return Returns NULL invisibly.
#' @export
#' @importFrom graphics hist
ClassificationHistogram <- function(out_classific = NULL, breaks = seq(0, 1, 0.05), ...) {
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  par(mfrow = c(length(out_classific) * 2, 3))

  # if (!is.null(cex)) par("cex"=cex)
  par(...)

  x_title_set <- TRUE
  for (i in 1:length(out_classific)) {
    for (j in 1:2) {
      tmp.x <- sapply(out_classific[[i]][[j]], function(x) {
        x[["ConfusionMatrix"]][["overall"]]["Accuracy"]
      })
      if (any(is.na(tmp.x))) {
        for (k in 1:3) {
          graphics::plot(1, 1, ann = F, axes = F)
          graphics::text(1, 1, labels = "No valid models could be established.")
        }
      } else {
        # Accuracy histogram
        graphics::hist(tmp.x, breaks = breaks, las = 1, main = ifelse(x_title_set, "Accuracy", ""), ylab = paste(names(out_classific)[i], names(out_classific[[i]])[j], sep = "_"), xlab = "")
        mu <- round(mean(tmp.x), 2)
        mtext(text = paste("m =", mu), side = 3, line = -1.1, adj = ifelse(mu < 0.5, 0.975, 0.025))
        # Sensitivity histogram (with respect to first level
        tmp.x <- sapply(out_classific[[i]][[j]], function(x) {
          x[["ConfusionMatrix"]][["byClass"]]["Sensitivity"]
        })
        graphics::hist(tmp.x, breaks = breaks, las = 1, main = ifelse(x_title_set, "Sensitivity", ""), ylab = "", xlab = "")
        mu <- round(mean(tmp.x), 2)
        mtext(text = paste("m =", mu), side = 3, line = -1.1, adj = ifelse(mu < 0.5, 0.975, 0.025))
        # specificity histogram (with respect to first level
        tmp.x <- sapply(out_classific[[i]][[j]], function(x) {
          x[["ConfusionMatrix"]][["byClass"]]["Specificity"]
        })
        graphics::hist(tmp.x, breaks = breaks, las = 1, main = ifelse(x_title_set, "Specificity", ""), ylab = "", xlab = "")
        mu <- round(mean(tmp.x), 2)
        mtext(text = paste("m =", mu), side = 3, line = -1.1, adj = ifelse(mu < 0.5, 0.975, 0.025))
        x_title_set <- FALSE
      }
    }
  }

  invisible(NULL)
}
