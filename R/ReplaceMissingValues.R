#'@title ReplaceMissingValues.
#'
#'@description
#'\code{ReplaceMissingValues} will replace missing values within a numeric matrix based on a principal component analysis.
#'
#'@details
#'The nipals algorithm is used to basically perform a PCA on the sparse matrix. Missing values are imputed based on the major components observed.
#'
#'@param x Numeric matrix.
#'@param ncomp Number of components to be used.
#'@param silent FALSE, suppress messages setting silent=TRUE.
#'
#'@return
#'Matrix without missing values.
#'
#'@examples
#'# load raw data and sample description
#'utils::data(raw, package = "MetabolomicsBasics")
#'utils::data(sam, package = "MetabolomicsBasics")
#'
# find outliers, store their values and substitute against NA within matrix
#'idx <- apply(raw, 2, CheckForOutliers, group=sam$GT, n_sd=5, method="logical")
#'sum(idx) # 215 values would be classified as outlier using a five-sigma band
#'old_vals <- raw[idx] # keep outlier values for comparison
#'raw_filt <- raw
#'raw_filt[idx] <- NA
#'raw_means <- apply(raw, 2, function(x) {
#'  sapply(split(x, sam$GT), mean, na.rm=TRUE)[as.numeric(sam$GT)]
#'})[idx]
#'raw_repl <- ReplaceMissingValues(x=raw_filt)
#'new_vals <- raw_repl[idx]
#'par(mfrow=c(2,1))
#'breaks <- seq(-0.7,1.3,0.05)
#'hist(raw_means-old_vals, breaks=breaks, main="", xlab="Outliers", las=1)
#'hist(raw_means-new_vals, breaks=breaks, main="", xlab="Replaced values", las=1)
#'
#'@export
#'
#'@importFrom utils flush.console
#'@importFrom mixOmics nipals
#'

ReplaceMissingValues <- function(x, ncomp=10, silent=FALSE) {
	if (!silent) cat(paste("\n...replacing missing values in a data matrix of m x n = ", nrow(x), " x ", ncol(x), "(=", prod(dim(x)), ")", sep="")); flush.console()
	nipals.x <- mixOmics::nipals(x, reconst = TRUE, ncomp = ncomp)
    id.na <- is.na(x) # only replace the imputation for the missing values
    n.na <- sum(id.na)
    i <- ncomp
	while (n.na > 0 && i>=3) {
		rec <- nipals.x$t[,1:i] %*% diag(nipals.x$eig[1:i],i,i) %*% t(nipals.x$p[,1:i])
	   x[id.na] <- rec[id.na]
		id.na <- x < 0 # check for missing values
		if (!silent) {
		  cat(paste("\n...replaced ", n.na-sum(id.na), " missing values using n=", i, " components.", sep=""))
		  utils::flush.console()
		}
    n.na <- sum(id.na)
    i <- i-1
	}
  if (!silent) cat("\n\n")
  return(x)
}
