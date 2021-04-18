#'@title find_boundaries.
#'
#'@description
#'\code{find_boundaries} will determine peak boundaries within a BPC or mass trace.
#'
#'@details
#'It is yet another peak finder or, more precisely, it is a function to identify two RT values which flank a intensity maximum which is required if one would like to integrate the peak area.
#'
#'@param int The measured intensity of the ion mass (obviously ordered according to consecutive RTs).
#'@param rt The respective retention times (can be omitted as currently not used).
#'@param p The anticipated peak position (as index of int) if several peaks are within the mass trace.
#'@param k The smoothing window parameter (provided to runmed).
#'@param bl The baseline value. Can be provided explicitly if automatic determination is insufficient.
#'@param local_min This is practically the upper end of the baseline. It can be set to avoid boundary detection at local minima (e.g. for peaks suffering ion suppression).
#'
#'@return
#'Numeric vector of length=2 specifying the start and end index of the peak.
#'
#'@examples
#'int <- sin(seq(-0.75*pi,1.75*pi,by=0.1))
#'plot(int)
#'abline(v=find_boundaries(int=int))
#'abline(v=find_boundaries(int=int, p=1))
#'
#'@importFrom stats runmed
#'
#'@export

find_boundaries <- function(int=NULL, rt=NULL, p=which.max(int), k=3, bl=min(int), local_min=int[p]) {
  # potential parameters
  min_inc <- 0.01 # at least 1% of peak height should be observed as difference between 2 consecutive values to differntiate noisebaseline from peak flan/increase
  # remove non-finite values and ensure positive values
  int[!is.finite(int)] <- 0
  int <- int-min(int)
  n <- length(int)
  # set up index
  idx <- 1:n
  # remove potential maxima from flanking peaks still present within int
  while (p %in% c(1,length(int)) && length(int)>5) {
    warning("Could not detect appropriate boundaries as intensity maximum was found on mass trace border.")
    int <- int[-p]
    idx <- idx[-p]
    n <- length(int)
    p <- which.max(int)
  }
  # default value for local_min
  #if (is.null(local_min)) local_min <- int[p]
  # default value for bl
  #if (is.null(bl)) bl <- min(int)
  # detect left/right boundaries
  if (p>=3) {
    dint <- runmed(diff(int[1:p]), k=k)
    lb <- dint<0 & int[1:(p-1)]<local_min
    lb <- ifelse(any(lb), (1+max(which(lb))), 1)
    test <- dint[lb:(p-1)]>min_inc*int[p]
    if (any(test)) lb <- (lb:(p-1))[min(which(test))]
  } else {
    lb <- 1
  }
  if (p<=n-3) {
    dint <- runmed(diff(int[p:n]), k=k)
    rb <- dint>0 & int[(p+1):n]<local_min
    rb <- ifelse(any(rb), (p+min(which(rb))-1), n)
    test <- abs(dint[1:(rb-p)])>min_inc*int[p]
    if (any(test)) rb <- ((p+1):rb)[max(which(test))]
  } else {
    rb <- n
  }
  return(idx[c(lb,rb)])
}

find_boundaries2 <- function(int=NULL, p=which.max(int), k=3, min_scans=3) {
  int[!is.finite(int)] <- 0
  idx <- 1:length(int)
  n <- length(idx)
  test_front <- diff(int[1:p])<=0
  test_front <- which(rev(cumsum(as.numeric(rev(test_front))))==k)
  test_front <- ifelse(length(test_front)>=1, max(test_front)+2, 1)
  lb <- idx[max(c(min(c(p-min_scans,test_front)),1))]
  test_tail <- diff(int[p:n])>=0
  test_tail <- which(cumsum(as.numeric(test_tail))==k)
  test_tail <- ifelse(length(test_tail)>=1, p+min(test_tail)-2, n)
  rb <- idx[min(c(max(c(p+min_scans,test_tail)),n))]
  return(idx[c(lb,rb)])
}
