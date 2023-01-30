#' @title AdjustSymbols.
#' @description \code{AdjustSymbols} will generate plotting character and color vectors based
#'   on experimental factors.#'
#' @details Using a fixed color and symbol scheme indicating samples from different groups
#'   throughout all figures of a analysis workflow is a reasonable decision. This function
#'   allows to specify both and attach it to a sample table for further use.
#' @param cols Factor (color output) or numeric (grey-scale output) vector or NULL (omitted).
#' @param pchs Factor vector or NULL (omitted).
#' @param colorset Color definitions for the factor levels of `cols` (can be omitted to use default values).
#' @param symbolset Plotting character definitions for the factor levels of `pchs` (can be omitted to use default values).
#' @return Either a vector (if one parameter of `cols` and `pchs` remains NULL), a data frame
#'   with columns `cols` and `pchs` (if both are provided and of equal length) or a list of
#'   length 2 (if both are provided and of different length). Will be used by several plotting
#'   functions of the package internally.
#' @examples
#' # return color vector
#' x <- gl(6, 3)
#' y <- as.numeric(x)
#' plot(y, bg = AdjustSymbols(cols = x), pch = 21, cex = 2)
#' plot(y, bg = AdjustSymbols(cols = y), pch = 21, cex = 2)
#' plot(y, bg = AdjustSymbols(cols = x, colorset = 1:6), pch = 21, cex = 2)
#' plot(y, pch = AdjustSymbols(pchs = x), cex = 2)
#' plot(y, bg = 2, pch = AdjustSymbols(pchs = x, symbolset = 1:6), cex = 2)
#'
#' # load data and plot using provided color scheme
#' raw <- MetabolomicsBasics::raw
#' sam <- MetabolomicsBasics::sam
#' head(sam)
#' plot(y = raw[, 1], x = as.numeric(sam$GT), pch = sam$pchs, bg = sam$cols)
#'
#' # change colors to greyscale
#' head(AdjustSymbols(cols = sam$GT, pchs = sam$Origin))
#' tmp.set <- grDevices::rainbow(length(levels(sam$GT)))
#' head(AdjustSymbols(cols = sam$GT, pchs = sam$Batch, colorset = tmp.set))
#' plot(raw[, 1] ~ sam$GT, col = unique_labels(sam = sam, g = "GT")[, "cols"])
#' sam$cols <- AdjustSymbols(cols = as.numeric(sam$GT))
#' plot(raw[, 1] ~ sam$GT, col = unique_labels(sam = sam, g = "GT")[, "cols"]) #'
#' @export
#' @importFrom grDevices colors

AdjustSymbols <- function(cols = NULL, pchs = NULL, colorset = NULL, symbolset = NULL) {
  # specify colors if vector is provided
  if (!is.null(cols)) {
    stopifnot(is.vector(cols) | is.factor(cols))
    # if colorset is unspecified select favorite of 22 distinct depending on length(unique(cols))
    n <- length(unique(cols))
    if (is.null(colorset)) {
      colorset <- rep(c(rgb(222 / 255, 109 / 255, 80 / 255), rgb(40 / 255, 109 / 255, 222 / 255), rgb(64 / 255, 192 / 255, 0 / 255), rgb(255 / 255, 175 / 255, 0 / 255), grDevices::colors()[c(68, 134, 12, 38, 367, 31, 139, 128, 30, 100, 23, 142, 656, 53, 6, 26, 172, 11, 116, 47, 94, 58)]), length.out = n)
    } else {
      colorset <- rep(colorset, length.out = n)
    }
    if (is.numeric(cols)) {
      # specify greyscale
      cols <- (cols - min(cols, na.rm = T)) / (max(cols, na.rm = T) - min(cols, na.rm = T))
      ramp <- grDevices::colorRamp(colors = c(grDevices::grey(0.15), grDevices::grey(0.95)))
      cols <- sapply(cols, function(x) {
        grDevices::rgb(ramp(x) / 255)
      })
    } else {
      # specify colorscale
      cols <- colorset[as.numeric(factor(cols))]
    }
  }
  # specify symbols
  if (!is.null(pchs)) {
    stopifnot(is.vector(pchs) | is.factor(pchs))
    # ensure that pchs is a factor
    pchs <- factor(pchs)
    # if symbolset is unspecified select 5 favorite symbols (replicate if necessary)
    if (is.null(symbolset)) {
      symbolset <- rep(c(21:22, 24:25, 23), length.out = length(levels(pchs)))
    }
    pchs <- symbolset[as.numeric(pchs)]
  }
  if (!is.null(cols) && !is.null(pchs)) {
    if (length(cols) == length(pchs)) {
      return(data.frame(cols, pchs, stringsAsFactors = FALSE))
    } else {
      warning("'cols' and 'pchs' have not been of similar length. Returning list.")
      return(list(cols, pchs))
    }
  } else {
    if (!is.null(cols)) {
      return(cols)
    }
    if (!is.null(pchs)) {
      return(pchs)
    }
  }
  if (is.null(cols) && is.null(pchs)) {
    warning("Neither 'cols' nor 'pchs' have been provided. Returning NULL.")
    return(NULL)
  }
}
