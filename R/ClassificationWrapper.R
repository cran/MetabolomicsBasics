#'@title ClassificationWrapper.
#'
#'@description
#'\code{ClassificationWrapper} will do classification using SVM's and/or Decision Trees including cross validation.
#'
#'@details
#'not yet
#'
#'@param g Group-vector, factor.
#'@param d data, matrix or data.frame !! needs row/col-names.
#'@param n replicates of classifications.
#'@param k Fold cross validation.
#'@param method Currently \code{svm} and decison tree methods \code{C50} and \code{rpart} are supported.
#'@param train Either NULL (random permutations) or an index vector for a training subset out of \code{g}.
#'@param rpart.control Forwarded to rpart.
#'
#'@return
#'Classification results as list.
#'
#'@examples
#'utils::data(raw, package = "MetabolomicsBasics")
#'utils::data(sam, package = "MetabolomicsBasics")
#'gr <- sam$Origin
#'
#'# establish a basic rpart model and render a fancy plot including the accuracy
#'class_res <- ClassificationWrapper(d=raw, g=gr, method="rpart", n=10)
#'\donttest{
#'class_res <- ClassificationWrapper(d=raw, g=gr, method=c("rpart","svm"))
#'}
#'
#'@export
#'
#'@import C50
#'@import e1071
#'@import caret
#'@import graphics

ClassificationWrapper <- function(d=NULL, g=NULL, n=100, k=5, method=c("C50","svm","rpart"), train=NULL, rpart.control=list()) {

  implemented_methods <- c("C50","svm","rpart")

  # wrapper function only meaningful if replications are used
  stopifnot(n>=2, !is.null(d), !is.null(g))

  # limit requested methods to the ones available
  method <- method[method %in% implemented_methods]
  if (any(duplicated(method))) method <- unique(method)
  stopifnot(length(method)>=1)

  # classification methods ar often sensitive to missing data
  # replace missing values before classification (based on mixOmics package)
  if (!all(is.finite(d))) {
    warning("Missing values found. Replace using NIPALS algorithm.")
    d <- ReplaceMissingValues(x=d)
  }

  out_classific <- lapply(method, function(x) { 
    x <- vector("list",length=2)
    if (is.null(train)) {
      names(x) <- c("obs","perm")
    } else {
      names(x) <- c("train","test")
    }
    return(x) 
  })
  names(out_classific) <- method
  if (is.null(train)) {
    for (i in 1:length(method)) {
      out_classific[[i]][[1]] <- ClassificationCV(g=g, d=d, n=n, k=k, rand=F, method=method[i], rpart.control=rpart.control) # n replications of k-fold cross validated version of observed data
      out_classific[[i]][[2]] <- ClassificationCV(g=g, d=d, n=n, k=k, rand=T, method=method[i], rpart.control=rpart.control) # n replications of k-fold cross validated version of permuted data
    }    
  } else {
    for (i in 1:length(method)) {
      out_classific[[i]][[1]] <- ClassificationCV(g=g[train], d=d[train,], n=n, k=k, rand=F, method=method[i], rpart.control=rpart.control) # n replications of k-fold cross validated version of observed data
      out_classific[[i]][[2]] <- ClassificationCV(g=g[-train], d=d[-train,], n=n, k=k, rand=F, method=method[i], rpart.control=rpart.control) # n replications of k-fold cross validated version of permuted data
    }    
  }

  # return detailed results invisibly
  invisible(out_classific)
}
