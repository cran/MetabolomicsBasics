#' @title ClassificationWrapper.
#' @description \code{ClassificationWrapper} will do classification using SVM's and/or
#'   Decision Trees including cross validation.
#' @details Parameter `n_rand` will influence how permutation testing for robustness
#'   is conducted. If n_rand=1 than samples will be permuted exactly one time and
#'   subjected to n replications (with respect to fold splitting). If n_rand>1, samples
#'   will be permuted this many times but number of replications will be lowered to limit
#'   processing time. A good compromise is to balance both, using less replications than
#'   for observed data but on several randomizations.
#' @param g Group-vector, factor.
#' @param d data, matrix or data.frame !! needs row/col-names.
#' @param n replicates of classifications, i.e. number of different split into folds.
#' @param n_rand different number of randomizations, see Details.
#' @param k Fold cross validation.
#' @param method Currently \code{svm}, \code{ropls} and decision tree methods (\code{C50} and \code{rpart}) are supported.
#' @param train Either NULL (random permutations) or an index vector for a training subset out of \code{g}.
#' @param method.control A list of parameters, forwarded to the selected methods function.
#' @param silent Logical. Set TRUE to supress progress bar and warnings.
#' @return #' Classification results as list.
#' @examples
#' utils::data(raw, package = "MetabolomicsBasics")
#' utils::data(sam, package = "MetabolomicsBasics")
#' gr <- sam$Origin
#' # establish a basic rpart model and render a fancy plot including the accuracy
#' class_res <- ClassificationWrapper(d = raw, g = gr, method = c("rpart", "svm"), n = 3, k = 3)
#' ClassificationHistogram(class_res)
#' @export
ClassificationWrapper <- function(d = NULL, g = NULL, n = 100, n_rand = 1, k = 5, method = c("C50", "svm", "rpart", "ropls"), train = NULL, method.control = list(), silent = FALSE) {
  implemented_methods <- c("C50", "svm", "rpart", "ropls")

  # wrapper function only meaningful if replications are used
  stopifnot(n >= 2, !is.null(d), !is.null(g), n >= n_rand)

  # limit requested methods to the ones available
  method <- method[method %in% implemented_methods]
  if (any(duplicated(method))) method <- unique(method)
  stopifnot(length(method) >= 1)

  # classification methods ar often sensitive to missing data
  # replace missing values before classification (based on mixOmics package)
  if (!all(is.finite(d))) {
    warning("Missing values found. Replace using NIPALS algorithm.")
    d <- ReplaceMissingValues(x = d)
  }

  out_classific <- lapply(method, function(x) {
    x <- vector("list", length = 2)
    if (is.null(train)) {
      names(x) <- c("obs", "perm")
    } else {
      names(x) <- c("train", "test")
    }
    return(x)
  })
  names(out_classific) <- method

  if (is.null(train)) {
    for (i in 1:length(method)) {
      out_classific[[i]][["obs"]] <- ClassificationCV(g = g, d = d, n = n, k = k, rand = F, method = method[i], method.control = method.control) # n replications of k-fold cross validated version of observed data
      out_classific[[i]][["perm"]] <- sapply(1:n_rand, function(j) {
        ClassificationCV(g = g, d = d, n = ceiling(n / n_rand), k = k, rand = j, method = method[i], method.control = method.control)
      }) # n replications of k-fold cross validated version of permuted data
    }
  } else {
    for (i in 1:length(method)) {
      out_classific[[i]][["train"]] <- ClassificationCV(g = g[train], d = d[train, ], n = n, k = k, rand = F, method = method[i], method.control = method.control) # n replications of k-fold cross validated version of observed data
      out_classific[[i]][["test"]] <- sapply(1:n_rand, function(j) {
        ClassificationCV(g = g[-train], d = d[-train], n = ceiling(n / n_rand), k = k, rand = j, method = method[i], method.control = method.control)
      }) # n replications of k-fold cross validated version of permuted data
    }
  }

  # return detailed results invisibly
  invisible(out_classific)
}
