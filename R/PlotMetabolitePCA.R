#'@title PlotMetabolitePCA
#'
#'@description
#'\code{PlotMetabolitePCA} will show PC1 and PC2 of a pcaMethods object and generate a flexible plot.
#'
#'@details
#'not yet
#'
#'@param pca_res A pcaRes object from the pcaMethods package.
#'@param sam Sample table including columns 'cols', 'pchs' (for data point color and shape) and 'ID' (to label data points) 'Group' (to split cols for legend) 'MP' (to adjust point size).
#'@param g Can be a factor vector of length=nrow(sam) and will influence legend and medsd.
#'@param medsd Calculate mean and sd for groups and overlay PCA plot with this information.
#'@param text.col Datapoints may be overlaid by textual information, e.g. sample ID and 'text.col' specifies the colummn name of sam to use for this purpose.
#'@param legend.x Position of a legend or NULL to omit it.
#'@param comm Will print commentary text to the bottom right of the plot (can be a character vector).
#'
#'@return
#'A vector fo similar length as input but with various name components removed.
#'
#'@examples
#'# load raw data and sample description
#'utils::data(raw, package = "MetabolomicsBasics")
#'utils::data(sam, package = "MetabolomicsBasics")
#'
#'# calculate pca Result using pcaMethods and plot
#'pca_res <- pcaMethods::pca(raw, method="rnipals", scale=c("none", "pareto", "uv")[2])
#'PlotMetabolitePCA(pca_res=pca_res, sam=sam, g=sam$GT)
#'# plot without legend and Group means instead
#'PlotMetabolitePCA(pca_res=pca_res, sam=sam, g=sam$GT, legend.x=NULL, text.col=NULL,
#'                  medsd=TRUE, comm=LETTERS[1:4])
#'
#'sam$Group <- interaction(sam$Origin, sam$Class, sep="_")
#'sam[,c("cols","pchs")] <- AdjustSymbols(cols=sam$Group, pchs=sam$Group)
#'PlotMetabolitePCA(pca_res=pca_res, sam=sam, g=sam$Group)
#'
#'@export
#'
#'@importFrom graphics title
#'@importFrom graphics mtext
#'

PlotMetabolitePCA <- function(pca_res=NULL, sam=NULL, g=NULL, medsd=FALSE, text.col="ID", legend.x="bottomleft", comm=NULL) {
# check if grouping variable exists
	stopifnot(!is.null(g) | ("Group" %in% colnames(sam)))
	if (!is.null(g)) {
		stopifnot(length(g)==nrow(sam))
		sam$Group <- factor(g)
	} else {
	  sam$Group <- factor(sam$Group)
	}
# # keep old plot parameters
#   opar <- par(no.readonly = TRUE)
#   on.exit(par(opar))
# check for existing col/pch values or provide respective cols
	if (!all(c("cols","pchs") %in% colnames(sam))) { sam[,c("cols","pchs")] <- AdjustSymbols(cols=sam$Group, pchs=sam$Group) }
	tmp.lab <- paste0("PC", 1:2, " (", round(100*pca_res@R2,2), "%)")
# plot pca result
	graphics::par(mar=c(4,4,0,0)+.5)
	# compute equal distance limits for x and y-axis
	tmp.lim <- range(pca_res@scores[,1:2])
	tmp.lim <- tmp.lim + c(-1,1) * 0.04*diff(tmp.lim)
	graphics::plot(pca_res@scores, las=1, xlim=tmp.lim, ylim=tmp.lim, ann=F, type="n")
	graphics::title(xlab=tmp.lab[1], ylab=tmp.lab[2], line=2.5)
# scale data points according to sam$MP if present or set to 2
	if ("MP" %in% colnames(sam)) {
		tmp.cex <- 1.5 + 1.5*(sam$MP-min(sam$MP))/(max(sam$MP)-min(sam$MP))
	} else {
		tmp.cex <- 2
	}
	graphics::points(pca_res@scores, bg=sam$cols, pch=sam$pchs, cex=tmp.cex, lend=1)
# annotate data points with values from column 'text.col' out of sam if present
	if (!is.null(text.col) && text.col %in% colnames(sam)) {
		graphics::text(pca_res@scores, labels=sam[,text.col], pos=3)
	}
# compute group wise symbols and colors
	l <- unique_labels(sam=sam, g=sam$Group)
# plot medians and sds if requested
	if (medsd) {
	  gt.med <- apply(pca_res@scores, 2, function(x) {sapply(split(x, sam$Group), stats::median, na.rm=TRUE)})
		gt.sdv <- apply(pca_res@scores, 2, function(x) {sapply(split(x, sam$Group), stats::sd, na.rm=TRUE)})
		graphics::segments(x0=gt.med[,1]-gt.sdv[,1], x1=gt.med[,1]+gt.sdv[,1], y0=gt.med[,2], col=grey(0.8))
		graphics::segments(y0=gt.med[,2]-gt.sdv[,2], y1=gt.med[,2]+gt.sdv[,2], x0=gt.med[,1], col=grey(0.8))
		graphics::points(gt.med, cex=stats::median(tmp.cex)*2, lwd=3, pch=l[,"pchs"], bg=l[,"cols"])
	}
# add a legend (!! the factor() function is necessary if not all groups are present within the plot)
	if (!is.null(legend.x)) {
	  graphics::legend(x=legend.x, pt.bg=l[,"cols"], pch=l[,"pchs"], legend=l[,"Level"], horiz=c(T,F)[2], inset=0.01)
	}
# add commentary text if provided
	for (i in 1:length(comm)) graphics::mtext(text=comm[i], side=1, line=-0.25-i, adj=0.99)
	invisible(NULL)
}
