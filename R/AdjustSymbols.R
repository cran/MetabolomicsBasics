#'@title AdjustSymbols.
#'
#'@description
#'\code{AdjustSymbols} will generate plotting character and color vectors based on experimental factors.
#'
#'@details
#'not yet
#'
#'@param cols Factor (color output) or numeric (greyscale output) vector or NULL (omitted).
#'@param pchs Factor vector or NULL (omitted).
#'@param colorset Can be selectively specified here. If NULL set automatically, else can be explicitly provided.
#'@param symbolset Can be selectively specified here. If NULL upt o 5 nice symbols are selected automatically where background can be colored.
#'
#'@return
#'data.frame with two columns (cols, pchs). Will be used by several plotting functions automatically.
#'
#'@examples
#'# load data and plot using provided color scheme
#'utils::data(raw, package = "MetabolomicsBasics")
#'utils::data(sam, package = "MetabolomicsBasics")
#'head(sam)
#'plot(y=raw[,1], x=as.numeric(sam$GT), pch=sam$pchs, bg=sam$cols)
#'
#'# change colors to greyscale
#'head(AdjustSymbols(cols=sam$GT, pchs=sam$Origin))
#'tmp.set <- grDevices::rainbow(length(levels(sam$GT)))
#'head(AdjustSymbols(cols=sam$GT, pchs=sam$Batch, colorset=tmp.set))
#'plot(raw[,1]~sam$GT, col=unique_labels(sam=sam, g="GT")[,"cols"])
#'sam$cols <- AdjustSymbols(cols=as.numeric(sam$GT))
#'plot(raw[,1]~sam$GT, col=unique_labels(sam=sam, g="GT")[,"cols"])
#'
#'@export
#'
#'@importFrom grDevices colors

AdjustSymbols <- function(cols=NULL, pchs=NULL, colorset=NULL, symbolset=NULL) {
  # specify colors
	if (!is.null(cols)) {
	  # if colorset is unspecified select 6 favorit of 22 distinct depending on length(unique(cols))
		if (is.null(colorset)) {
	    colorset <- rep(c(rgb(222/255,109/255,80/255),rgb(40/255,109/255,222/255),rgb(64/255,192/255,0/255),rgb(255/255,175/255,0/255), grDevices::colors()[c(68,134,12,38,367,31,139,128,30,100,23,142,656,53,6,26,172,11,116,47,94,58)]), length.out=length(unique(cols)))
		} else {
		  colorset <- rep(colorset, length.out=length(unique(cols)))
		}
		if (is.numeric(cols)) {
		  # specify greyscale
			cols <- (cols-min(cols,na.rm=T))/(max(cols,na.rm=T)-min(cols,na.rm=T))
			ramp <- grDevices::colorRamp(colors=c(grDevices::grey(0.15),grDevices::grey(0.95)))
			cols <- sapply(cols, function(x) { grDevices::rgb(ramp(x)/255) })
		} else {
		  # specify colorscale
			cols <- colorset[as.numeric(cols)]
		}
	}
  # specify symbols
  if (!is.null(pchs)) {
    # ensure that pchs is a factor
    pchs <- factor(pchs)
    # if symbolset is unspecified select 5 favorit symbols (replicate if necessary)
	  if (is.null(symbolset)) {
	    symbolset <- rep(c(21:22,24:25,23),length.out=length(levels(pchs)))
	  }
		pchs <- symbolset[as.numeric(pchs)]
	}
	if (!is.null(cols) && !is.null(pchs)) {
	  if (length(cols)==length(pchs)) {
	    return(data.frame(cols, pchs, stringsAsFactors = FALSE))
	  } else {
	    warning("'cols' and 'pchs' have not been of similar length. Returning list.")
	    return(list(cols, pchs))
	  }
	} else {
		if (!is.null(cols)) return(cols)
		if (!is.null(pchs)) return(pchs)
	}
  if (is.null(cols) && is.null(pchs)){
    warning("'cols' and 'pchs' have not been properly defined. Returning NULL.")
	  return(NULL)
  }
}
