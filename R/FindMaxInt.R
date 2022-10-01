#' @title FindMaxInt.
#' @description \code{FindMaxInt} will evaluate a list of xcmsRaw objects
#'   at a given time window (rt, rt_dev) for the maximum value of a single mz.
#'   Alternatively pre-processed narrow BPCs can be used if stored in an object
#'   of class xcmsEIC.
#' @details not yet
#' @param dat A list of xcmsRaws as imported by \code{\link{Load_mzML}} or \code{\link{xcmsRaw}}. Alternatively an object of class xcmsEIC containing pre-extracted data.
#' @param mz The target mass (m/z) to search maxima for and plot BPCs.
#' @param rt The target retention time.
#' @param rt_dev Allowed retention time window.
#' @param mz_dev Allowed lower and upper deviation mz deviation [Da].
#' @param test_plot Generate plot output if TRUE.
#' @param mfrow Specify mfrow explicitely.
#' @param ylim Specify ylim explicitely.
#' @param col Specify col explicitely.
#' @param ids Specify plot ids explicitely.
#' @param skip_plots Allows to block certain subplots in the mfrow matrix to bettern align replicates.
#' @param ann Names of result-vector may contain RT of maxInt-Scan OR Scan-ID OR MassDeviation.
#' @param output You can select to return a list of BPCs by selecting 'bpc' instead of 'maxint' here.
#' @return A named vector of intensity maximas with names being either Scan, RT or mass deviation.
#' @keywords internal
#' @noRd
#' @importFrom graphics lines axis abline text box
FindMaxInt <- function(dat = NULL, mz = NULL, rt = NULL, rt_dev = 5, mz_dev = 0.01, test_plot = FALSE, mfrow = NULL, skip_plots = NULL, ylim = c(100, 100), col = NULL, ids = NULL, ann = c("RT", "Scan", "dmz")[1], output = "maxint") {

  # check for xcms to be able to keep it in suggested packages
  if (!requireNamespace("xcms", quietly = TRUE)) {
    stop("The use of this function requires package 'xcms'. Please ",
         "install with 'Biobase::install(\"xcms\")'")
  }

  if (!inherits(dat, c("xcmsSet", "xcmsEIC", "list"))) {
    dat <- list(dat)
  }
  if (inherits(dat, "xcmsSet")) {
    # dat <- xcmsExtra::getBPC(dat, mzrange=matrix(mz+c(-1,1)*mz_dev, ncol=2), rtrange=matrix(rt+c(-1,1)*rt_dev, ncol=2))
    # dat <- getBPC(dat, mz=mz, mz_dev=mz_dev, rt=rt, rt_dev=rt_dev)
    cat("\nThe functionality to work on xcmsSets has been stripped on July 2017.\nModify the function to include it again if necessary.\n\n")
    dat <- lapply(dat@filepaths, xcms::xcmsRaw)
  }
  if (inherits(dat, "xcmsEIC")) {
    fm <- "eic"
    n_samples <- length(dat@eic)
  }
  if (inherits(dat, "list")) {
    fm <- "raw"
    n_samples <- length(dat)
  }
  if (test_plot) {
    opar <- par(no.readonly = TRUE)
    if (is.null(mfrow)) mfrow <- grDevices::n2mfrow(n_samples)
    if (is.null(col)) col <- rep(grey(0.8), n_samples)
    par(mfrow = mfrow, mar = c(2, 3, 3, 3), cex = max(c(0.5, 4 / max(mfrow))))
  }
  if (output == "maxint") res <- rep(NA, length(dat))
  if (output == "bpc") res <- vector("list", length(dat))
  cor_val <- 0
  for (i in 1:(n_samples + length(skip_plots))) {
    if (i %in% skip_plots) {
      cor_val <- cor_val + 1
      if (test_plot) plot(1, 1, ann = F, type = "n", axes = F)
    } else {
      j <- i - cor_val
      if (fm == "eic") {
        # nEIC <- which(abs(apply(dat@mzrange,1,mean)-mz)<mz_dev & abs(apply(dat@rtrange,1,mean)-rt)<rt_dev)
        nEIC <- which(dat@mzrange[, 1] < mz & dat@mzrange[, 2] > mz & dat@rtrange[, 1] < rt & dat@rtrange[, 2] > rt)
        if (length(nEIC) == 0) {
          warning("Couldn't find your specified mz/rt within provided xcmsEIC object.")
          # browser(TRUE)
          x <- NULL
        } else {
          if (length(nEIC) > 1) {
            # if several peaks with similar mass are co-located take closest to target rt
            nEIC <- nEIC[which.min(abs(apply(dat@rtrange[nEIC, ], 1, median) - rt))]
          }
          x <- dat@eic[[j]][[nEIC]]
          filt <- x[, "rt"] > rt - rt_dev & x[, "rt"] < rt + rt_dev
          x <- matrix(c(x[filt, "rt"], x[filt, "intensity"]), ncol = 2, dimnames = list(NULL, c("rt", "intensity")))
        }
      }
      if (fm == "raw") {
        x <- try(xcms::rawEIC(dat[[j]], mzrange = matrix(mz + c(-mz_dev, mz_dev), ncol = 2), rtrange = matrix(rt + c(-rt_dev, rt_dev), ncol = 2)), silent = TRUE)
        if (inherits(x, "try-error")) {
          warning(paste("Probably your defined scan time is out of range. mz =", mz, "rt =", rt), call. = FALSE)
          x <- NULL
        } else {
          x <- matrix(c(dat[[j]]@scantime[x[["scan"]]], x[["intensity"]], x[["scan"]]), ncol = 3, dimnames = list(NULL, c("rt", "intensity", "scan")))
        }
      }
      if (!is.null(x) && nrow(x) >= 1) {
        if (test_plot) {
          if (max(x[, 2]) > ylim[2]) tmp.ylim <- c(ylim[1], max(x[, 2])) else tmp.ylim <- ylim
          plot(x, main = ifelse(is.null(ids), j, ids[j]), type = "l", ylim = tmp.ylim, xlab = "", ylab = "", col = col[j], axes = F)
          axis(1)
          if (max(x[, 2]) > ylim[2]) {
            abline(h = ylim[2], col = grey(0.8), lty = 2)
            axis(4, at = ylim[2], col = grey(0.8))
          } else {
            lines(x = x[, 1], y = ylim[2] * (x[, 2] / max(x[, 2])), col = grey(0.8), lty = 2)
            axis(4, at = ylim[2], col = grey(0.8), labels = max(x[, 2]))
          }
          axis(2, at = max(x[, 2]), labels = formatC(max(x[, 2]), digits = 2, format = "e"))
          box()
        }
        f <- which.max(x[, 2])
        if (output == "maxint") res[j] <- x[f, 2]
        if (output == "bpc") res[[j]] <- x
        names(res)[j] <- round(x[f, ifelse(ann == "RT", 1, 3)], ifelse(ann == "RT", 1, 0))
      } else {
        if (test_plot) plot(1, 1, main = j, type = "n", ylim = ylim, xlab = "", ylab = "")
      }
    }
  }
  if (test_plot) par(opar)
  return(res)
}
