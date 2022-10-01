#' @title ClassificationCV.
#' @description \code{ClassificationCV} will perform a classification using SVM's
#'   and/or Decision Trees including cross validation on a data set according to
#'   a provided grouping vector.
#' @details This function allows to demonstrate the functionality of different
#'   classification tools with respect to building classifiers for metabolomics data.
#'   Check the examples in \code{\link{ClassificationWrapper}} for automatic
#'   multi-fold analysis.
#' @param d Data matrix or data.frame with named rows (samples) and columns (traits).
#' @param g Group-vector, factor.
#' @param n Replicates of classifications.
#' @param k Number of folds per replicate.
#' @param rand Randomize Group-vector (and apply according n and k to this randomization).
#' @param method Currently \code{svm}, \code{ropls} and decision tree methods \code{C50} and \code{rpart} are supported.
#' @param method.control A list of parameters, forwarded to the respective classification function.
#' @param silent Logical. Set TRUE to suppress progress bar and warnings.
#' @return A list of classification results which can be analyzed for accuracy,
#'   miss-classified samples and more.
#' @importFrom e1071 svm
#' @importFrom rpart rpart
#' @importFrom caret createMultiFolds
#' @importFrom caret confusionMatrix
#' @importFrom C50 C5.0
#' @importFrom C50 C5imp
#' @importFrom rlang call2
#' @importFrom rlang splice
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom stats runif
#' @importFrom stats predict
#' @importFrom stats as.formula
#' @export
ClassificationCV <- function(d = NULL, g = NULL, n = 1, k = 1, rand = F, method = c("svm", "C50", "rpart", "ropls")[1], method.control = list(), silent = FALSE) {

  # check for ropls to be able to keep it in suggested packages
  if (!requireNamespace("ropls", quietly = TRUE)) {
    stop("The use of this function requires package 'ropls'. Please ",
         "install with 'install.packages(\"opls\")'")
  }

  # warn and substitute in case of missing values
  if (method %in% c("C50", "svm") && !all(is.finite(d))) {
    warning("Missing values found. Replace using NIPALS algorithm.")
    d <- ReplaceMissingValues(x = d)
  }
  # without CV (returning only the classification model for the observed data (or a single random set))
  if (k == 1 && n > 1) {
    if (!silent) warning("...returning classification-model only as replications (n>1) are not meaningfull without crossvalidation (k=1).")
    n <- 1
  }
  # CV-fold must be <= number of replicates of smallest group
  stopifnot(1 <= min(table(g)) / k)
  stopifnot(!is.null(colnames(d)))
  stopifnot(!is.null(rownames(d)))
  # randomize group vector; !! if you want to do real permutations than use an outer loop and k=x, n=1 and rand=seed within this function
  if (is.numeric(rand)) {
    local.seed <- rand
    rand <- TRUE
  } else {
    local.seed <- n * k
  }
  if (rand) {
    set.seed(local.seed)
    g <- sample(factor(g))
  }
  # ensure that g is factor
  if (!is.factor(g)) {
    g <- as.factor(g)
    if (!silent) warning("Converting 'g' to factor.")
  }
  # ensure that d is data.frame
  d <- as.data.frame(d)
  # prepare function call
  # method.control <- list()
  if (method == "svm") {
    svm.control <- list("kernel" = "linear", "type" = "C-classification", "cost" = 1)
    for (arg in names(svm.control)) if (!arg %in% names(method.control)) method.control[[arg]] <- svm.control[[arg]]
    fnc <- rlang::call2(e1071::svm, rlang::splice(method.control))
  }
  if (method == "rpart") {
    rpart.control <- list("method" = "class")
    for (arg in names(rpart.control)) if (!arg %in% names(method.control)) method.control[[arg]] <- rpart.control[[arg]]
    fnc <- rlang::call2(rpart::rpart, rlang::splice(method.control))
  }
  if (method == "C50") {
    C50.control <- list("method" = "class")
    for (arg in names(C50.control)) if (!arg %in% names(method.control)) method.control[[arg]] <- C50.control[[arg]]
    fnc <- rlang::call2(C50::C5.0, rlang::splice(method.control))
  }
  if (method == "ropls") {
    ropls.control <- list("predI" = 1, "orthoI" = 1, printL = FALSE, plotL = FALSE)
    for (arg in names(ropls.control)) if (!arg %in% names(method.control)) method.control[[arg]] <- ropls.control[[arg]]
    fnc <- rlang::call2(ropls::opls, rlang::splice(method.control))
  }
  if (method == "rpart") {
    fnc$formula <- stats::as.formula("Group ~ .")
    fnc$data <- data.frame("Group" = g, d)
  } else {
    fnc$x <- d
    fnc$y <- g
  }
  if (k == 1) {
    # without CV
    tmp <- eval(fnc)
  } else {
    # with CV
    # sample balanced training sets in k-folds
    index <- caret::createMultiFolds(y = g, k = k, times = n)
    # set ProgressBar
    if (!silent) {
      print(paste("Doing", k, "fold classification with", method, "method using", n, "replications on", ifelse(rand, "randomized", "observed"), "data."))
      pb <- txtProgressBar(title = paste(n, "replicates"), label = "progress...", min = 0, max = n)
    }
    # do n replications of selected classification method and return interesting parameters
    tmp <- lapply(as.list(1:n), function(i) { # over all replications do...
      if (!silent) setTxtProgressBar(pb, value = i)
      pred <- vector("list", length = k)
      imp <- vector("list", length = k)
      misscl <- vector("list", length = k)
      numnodes <- vector("list", length = k)
      for (j in 1:k) { # for all k-folds of this replication 'n' do
        ind <- index[[(i - 1) * k + j]]
        if (method == "rpart") {
          fnc$formula <- stats::as.formula("Group ~ .")
          fnc$data <- data.frame("Group" = g[ind], d[ind, ])
        } else {
          fnc$x <- d[ind, ]
          fnc$y <- g[ind]
        }
        model <- try(eval(fnc))
        # try to catch some errors which may occur
        catch_err <- !inherits(model, c("opls", "svm", "rpart", "C5.0"))
        catch_err <- catch_err | (inherits(model, "C5.0") && model[["size"]] == 1)
        if (catch_err) {
          pred[[j]] <- NA
          misscl[[j]] <- NA
          imp[[j]] <- NA
        } else {
          if (method == "ropls") {
            pred[[j]] <- ropls::predict(model, newdata = d[-ind, ]) # prediction
          } else {
            pred[[j]] <- stats::predict(model, newdata = d[-ind, ]) # prediction
          }
          if (method == "rpart") {
            pred[[j]] <- apply(pred[[j]], 1, function(y) {
              colnames(pred[[j]])[which.max(y)]
            })
          }
          names(pred[[j]]) <- rownames(d[-ind, ])
          misscl[[j]] <- names(pred[[j]])[which(pred[[j]] != g[-ind])] # missclassification
          if (inherits(model, "svm")) imp[[j]] <- model[["variable.importance"]]
          if (inherits(model, "opls")) imp[[j]] <- model@"vipVn"
          if (inherits(model, "rpart")) numnodes[[j]] <- sum(model[["frame"]][, "var"] != "<leaf>")
          if (inherits(model, "C5.0")) imp[[j]] <- rownames(C50::C5imp(model)[C50::C5imp(model)$Overall != 0, , drop = F])
        }
      }
      pred <- factor(unlist(pred)[rownames(d)]) # combine the fold results
      if (any(is.na(pred))) {
        warning("Some decision tree models were of size 1 (no classification possible). Skip calculating confusion matrix and returning NA instead.")
        conf <- list("overall" = NA)
      } else {
        conf <- caret::confusionMatrix(data = pred, reference = g) # compute confusion matrix for this replicate
      }
      if (method == "ropls") {
        imp <- sapply(imp, function(x) {
          x
        })
      }
      return(list("Missclassified" = unlist(misscl), "ConfusionMatrix" = conf, "Importance" = unlist(imp), "NumberOfNodes" = unlist(numnodes)))
    })
    if (!silent) close(pb)
  }
  return(tmp)
}
