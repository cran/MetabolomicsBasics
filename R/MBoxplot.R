#'@title MBoxplot
#'
#'@description
#'\code{MBoxplot} will generate an annotated boxplot.
#'A unifying function for MS-data Boxplots based on \'raw\' and \'sam\'.
#'
#'@details
#'not yet
#'
#'@param pk Colname of raw to plot if \code{pk} is character OR the colnum number if \code{pk} is numeric.
#'@param raw Plotting data as samples (rows) x metabolites (cols).
#'@param sam Sample table.
#'@param met Containing at minimum columns for annotation (see parameter \code{an}) and \code{nrow(met)} should be \code{ncol(raw)}.
#'@param g Grouping vector if \code{Group} not contained in \code{sam}.
#'@param flt Filter to exclude certain samples (T/F) vector.
#'@param an Switch to include annotation (from met) in the boxplot providing a character vector of colnames from \code{met}.
#'@param plot_sample_n Amend each box with the number of finite values which were a basis for plotting this group.
#'@param txt Character vector with information per sample to be plotted on top of the box as text.
#'@param cex.txt Specify size of annotation text.
#'@param plot_rel_axis Specify one level of \code{g} (or \code{sam$Group}) which to express the data relative against.
#'@param ... Further options parsed to \code{boxplot}.
#'
#'@return
#'Nothing. Will produce a plot (or file if specified).
#'
#'@examples
#'x <- data.frame("y"=runif(36), "GT"=gl(3,12), "TP"=factor(rep(rep(1:3,each=4),3)))
#'x <- cbind(x, AdjustSymbols(cols=x$GT, pchs=x$TP))
#'MBoxplot(pk="y", raw=x, sam=x, met=data.frame("Peak"="y", "Test"=I("info")),
#'         g=interaction(x$GT, x$TP), an="Test", plot_n_samples=TRUE, txt=rownames(x))
#'
#'@export
#'
#'@import grDevices

MBoxplot <- function(pk=pk, raw=NULL, sam=NULL, met=NULL, g=NULL, flt=NULL, an=NULL, plot_sample_n=FALSE, txt=NULL, cex.txt=0.5, plot_rel_axis=NULL, ...) {
	# check parameters
	if (is.character(pk)) {
	  stopifnot(pk %in% colnames(raw))
	  pk <- which(colnames(raw)==pk)[1]
	}
	if (is.numeric(pk)) {
		stopifnot(pk %in% 1:ncol(raw))
	}
	stopifnot(nrow(sam)==nrow(raw))
	stopifnot(any(is.null(an), all(an %in% colnames(met))))
	stopifnot(any(!is.null(g), "Group" %in% colnames(sam)))
	# prepare data
	if (is.null(flt)) flt <- rep(T,nrow(raw))
	stopifnot(length(flt)==nrow(raw))
	if (is.null(g)) { tmp.x <- factor(sam$Group[flt]) } else { tmp.x <- factor(g[flt]) }
	if (!("cols" %in% colnames(sam))) {
		tmp.col <- sapply(split(AdjustSymbols(cols=tmp.x),tmp.x),unique)
	} else {
		tmp.col <- sapply(split(sam$cols[flt],tmp.x),unique)
	}
	tmp.y <- raw[flt,pk]
  tmp.l <- unique(sapply(gregexpr("_",tmp.x),function(x){length(x[x!=(-1)])}))
  if (length(tmp.l)>1) tmp.l <- 0
	par(mar=c(tmp.l+2, 4, length(an), 0) + 0.5)
	if (!is.null(plot_rel_axis) && plot_rel_axis %in% levels(tmp.x)) {
	  if (!any(tmp.y<0,na.rm=T)) {
      tmp.y <- tmp.y/median(tmp.y[tmp.x==plot_rel_axis],na.rm=T)
    } else {
      warning("[MBoxplot] cant compute relative log2-axis as data contain negative values", .call=FALSE)
    }
	}
  if (!is.null(plot_rel_axis)) {
    boxplot(tmp.y ~ tmp.x, las=1, col=tmp.col, xaxt="n", log="y", ...)
  } else {
	  boxplot(tmp.y ~ tmp.x, las=1, col=tmp.col, xaxt="n", ...)
  }
  if (tmp.l==0) {
    axis(side=1, at=1:length(levels(tmp.x)), labels=levels(tmp.x))
  } else {
	  for (i in 1:(tmp.l+1)) {
			axis(side=1, at=1:length(levels(tmp.x)), labels=sapply(levels(tmp.x), function(x) {strsplit(x, "_")[[1]][i]}), line=i-1, tick=ifelse(i>1,F,T))
		}
  }
  if (plot_sample_n) {
    mtext(text=table(tmp.x[is.finite(tmp.y)]), side=1, at=1:length(levels(tmp.x)), line=-1, cex=cex.txt)
  }
	if (!is.null(txt)) {
	  text(x=jitter(as.numeric(tmp.x)), y=tmp.y, labels=txt, cex=cex.txt, col="lightblue")
	}
	if (length(an)>=1) {
		for (i in 1:length(an)) {
      txt <- met[pk,an[i]]
      txt <- ifelse((par("fin")[1]-0.4)<strwidth(txt,"inch"), paste(substr(txt,1,floor(nchar(txt)*((par("fin")[1]-0.4)/strwidth(txt,"inch")))),"..."), txt)
      mtext(side=3, line=length(an)-i, text=txt, adj=1)
		}
	}
	invisible(NULL)
}
