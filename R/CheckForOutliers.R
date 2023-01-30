#' @title CheckForOutliers.
#' @description \code{CheckForOutliers} will evaluate a numeric vector and check
#'   if outliers within groups based on group \eqn{mean \pm n \times sd}{mean +- n x sd}.
#' @details The numeric will be split by groups and each value will be evaluated
#'   with respect to its distance to the group mean (calculated out of the other
#'   values in the group). Distance here means the number of standard deviations
#'   the value is off the group mean. With different choices of \code{method}
#'   the output can be switched from the calculated fold-distances to a boolean
#'   of length(x) or and Index vector giving the outliers directly (see examples).
#' @param x Numeric vector.
#' @param group Factor vector of length(x).
#' @param n_sd Cutoff for outliers in E being mean(group)+-n_sd*sd(group) where
#'   group values are calculated without the outlier candidate.
#' @param method Different variants of the result value. See details.
#' @return Depending on the selected method. See details.
#' @examples
#' set.seed(0)
#' x <- runif(10)
#' x[1] <- 2
#' group <- gl(2, 5)
#' CheckForOutliers(x, group, method = "dist")
#' CheckForOutliers(x, group, method = "logical")
#' CheckForOutliers(x, group, method = "idx")
#' graphics::par(mfrow = c(1, 2))
#' bg <- c(3, 2)[1 + CheckForOutliers(x, group, method = "logical")]
#' graphics::plot(x = as.numeric(group), y = x, pch = 21, cex = 3,
#'   bg = bg, main = "n_sd=3", las = 1, xlim = c(0.5, 2.5))
#' bg <- c(3, 2)[1 + CheckForOutliers(x, group, n_sd = 4, method = "logical")]
#' graphics::plot(x = as.numeric(group), y = x, pch = 21, cex = 3,
#'   bg = bg, main = "n_sd=4", las = 1, xlim = c(0.5, 2.5))
#' graphics::par(mfrow = c(1, 1))
#'
#' # load raw data and sample description
#' raw <- MetabolomicsBasics::raw
#' sam <- MetabolomicsBasics::sam
#'
#' # no missing data in this matrix
#' all(is.finite(raw))
#'
#' # check for outliers (computing n-fold sd distance from group mean)
#' tmp <- apply(raw, 2, CheckForOutliers, group = sam$GT, method = "dist")
#' # plot a histogram of the observed distances
#' graphics::hist(tmp, breaks = seq(0, ceiling(max(tmp))), main = "n*SD from mean", xlab = "n")
#'
#' # Calculate the amount of values exceeding five-sigma and compare with a standard gaussian
#' table(tmp > 5)
#' round(100 * sum(tmp > 5) / length(tmp), 2)
#' \donttest{
#' gauss <- CheckForOutliers(x = rnorm(prod(dim(raw))), method = "dist")
#' sapply(1:5, function(i) {
#'   data.frame("obs" = sum(tmp > i), "gauss" = sum(gauss > i))
#' })
#'
#' # compare a PCA w/wo outliers
#' RestrictedPCA(
#'   dat = raw, sam = sam, use.sam = sam$GT %in% c("Mo17", "B73"), group.col = "GT",
#'   fmod = "GT+Batch+Order", P = 1, sign.col = "GT", legend.x = NULL, text.col = "Batch", medsd = TRUE
#' )
#' raw_filt <- raw
#' raw_filt[tmp > 3] <- NA
#' RestrictedPCA(
#'   dat = raw_filt, sam = sam, use.sam = sam$GT %in% c("Mo17", "B73"), group.col = "GT",
#'   fmod = "GT+Batch+Order", P = 1, sign.col = "GT", legend.x = NULL, text.col = "Batch", medsd = TRUE
#' )
#' }
#' @export
#' @importFrom stats na.omit
#' @importFrom stats sd
CheckForOutliers <- function(x = NULL, group = NULL, n_sd = 3, method = c("idx", "logical", "dist")) {
  method <- match.arg(method)
  if (is.null(group)) group <- gl(1, length(x))
  test <- sapply(1:length(x), function(i) {
    ig <- na.omit(x[setdiff(which(group == group[i]), i)])
    if (length(ig) >= 2) {
      return(abs(mean(ig) - x[i]) / sd(ig))
    } else {
      return(0)
    }
  })
  if (method == "idx") {
    return(which(test >= n_sd))
  }
  if (method == "logical") {
    return(test >= n_sd)
  }
  if (method == "dist") {
    return(test)
  }
}
