#'@title ClassificationCV.
#'
#'@description
#'\code{ClassificationCV} will perform a classification using SVM's and/or Decision Trees including cross validation on a data set according to a provided grouping vector.
#'
#'@details
#'not yet
#'
#'@param d Data matrix or data.frame with named rows (samples) and columns (traits).
#'@param g Group-vector, factor.
#'@param n Replicates of classifications.
#'@param k Number of folds per replicate.
#'@param rand Randomize Group-vector (and apply according n and k to this randomization).
#'@param method Currently \code{svm} and decison tree methods \code{C50} and \code{rpart} are supported.
#'@param svm.type Pass through options for svm method.
#'@param svm.kernel Pass through options for svm method.
#'@param tree.method Pass through options for tree method.
#'@param rpart.control Forwarded to rpart.
#'
#'@return
#'Nothing.
#'
#'@examples
#'# classify american and european parental lines
#'utils::data(raw, package = "MetabolomicsBasics")
#'utils::data(sam, package = "MetabolomicsBasics")
#'gr <- sam$Origin
#'
#'# establish a basic rpart model and render a fancy plot including the accuracy
#'class_res <- ClassificationCV(d=raw, g=gr, method="rpart")
#'rattle::fancyRpartPlot(class_res, main="Optimal rpart decision tree")
#'acc_opt <- caret::confusionMatrix(predict(class_res, newdata=data.frame(raw), type="class"), gr)
#'mtext(paste("Accuracy =", round(acc_opt[["overall"]]["Accuracy"],2)), adj = 0)
#'
#'# now repeat the analysis in a robust fashion (10-fold cross validation in 100 permutations)
#'\donttest{
#'class_res <- ClassificationCV(g=gr, d=raw, n=100, k=10, method="rpart")
#'acc_crossval <- sapply(class_res, function(x) {
#'  round(x$ConfusionMatrix[["overall"]]["Accuracy"],2)
#'})
#'hist(acc_crossval, sub=paste("Accuracy optimal =",round(acc_optacc_opt[["overall"]]["Accuracy"],2)))
#'}
#'
#'# check the examples in \code{\link{ClassificationWrapper}} for automatic multifold analysis
#'
#'@export
#'
#'@import C50
#'@import e1071
#'@import rpart
#'@import caret
#'
#'@importFrom utils setTxtProgressBar
#'@importFrom utils txtProgressBar
#'@importFrom stats runif
#'@importFrom stats predict

ClassificationCV <- function(d=NULL, g=NULL, n=1, k=1, rand=F, method=c("svm","C50","rpart")[1], svm.kernel="linear", svm.type="C-classification", tree.method="class", rpart.control=list()) {
  # warn and substitute in case of missing values
  if (method %in% c("C50","svm") && !all(is.finite(d))) {
    warning("Missing values found. Replace using NIPALS algorithm.")
    d <- ReplaceMissingValues(x=d)
  }
  if (!is.factor(g)) { g <- as.factor(g); warning("Converting 'g' to factor.")}
  # CV-fold must be <= number of replicates of smallest group
  stopifnot(1 <= min(table(g))/k)
  stopifnot(!is.null(colnames(d)))
  stopifnot(!is.null(rownames(d)))
  set.seed(n*k)
  # randomize group vector; !! if you want to do real permutations than use an outer loop and k=x n=1 within this function
  if (rand) g <- sample(factor(g))
  # sample balanced training sets in k-folds
  if (k==1) {
    # without CV (returning only the classification model for the observed data (or a single random set))
    if (n>1) {
      warning("...returning classification-model only as replications (n>1) are not meaningfull without crossvalidation (k=1).")
      n <- 1
    }
    if (method=="svm") tmp <- e1071::svm(x=d, y=g, kernel = svm.kernel, type = svm.type)
    if (method=="C50") tmp <- C50::C5.0(x=d, y=g, method=tree.method)
    if (method=="rpart") tmp <- rpart::rpart(Group ~ ., data=data.frame("Group"=g, d), method="class")
  } else {
    # with CV
    index <- caret::createMultiFolds(y=g, k=k, times=n)
    # set ProgressBar
    print(paste("Doing", k, "fold classification with", method, "method using", n, "replications on", ifelse(rand, "randomized", "observed"), "data."))
    pb <- txtProgressBar(title=paste(n, "replicates"), label="progress...", min=0, max=n)
    # do n replications of selected classification method and return interesting parameters
    tmp <- lapply(as.list(1:n), function(i) { # over all replications do...
      setTxtProgressBar(pb, value=i)
      pred <- vector("list",length=k)
      imp <- vector("list",length=k)
      misscl <- vector("list",length=k)
      numnodes <- vector("list",length=k)
      for (j in 1:k) { # for all k-folds of this replication 'n' do
        ind <- index[[(i-1)*k+j]]
        if (method=="svm") {
          model <- e1071::svm(x=d[ind,], y=g[ind], kernel = svm.kernel, type = svm.type) # svm model
          pred[[j]] <- predict(model, newdata=d[-ind,]) # prediction
          misscl[[j]] <- names(pred[[j]])[which(pred[[j]]!=g[-ind])] # missclassification
        }
        if (method=="C50") {
          ## decision tree model/prediction
          model <- C50::C5.0(x=d[ind,], y=g[ind], method=tree.method) # Tree model
          # if the model contains only one perfect node C5imp failes
          if (model$size==1) {
            pred[[j]] <- NA
            misscl[[j]] <- NA
            imp[[j]] <- matrix(NA,dimnames=list("NA","Overall"))
          } else {
            #browser()
            pred[[j]] <- predict(model, newdata=d[-ind,]); names(pred[[j]]) <- rownames(d[-ind,])
            misscl[[j]] <- names(pred[[j]])[which(pred[[j]]!=g[-ind])] # missclassification
            imp[[j]] <- C50::C5imp(model)[C50::C5imp(model)$Overall!=0,,drop=F] # importance; nodes from the tree
          }
        }
        if (method=="rpart") {
          #browser()
          model <- rpart::rpart(Group ~ ., data=data.frame("Group"=g[ind], d[ind,]), method="class", control=rpart.control)
          pred[[j]] <- predict(model, newdata=data.frame(d[-ind,]), type="class"); names(pred[[j]]) <- rownames(d[-ind,])
          misscl[[j]] <- names(pred[[j]])[which(pred[[j]]!=g[-ind])] # missclassification
          imp[[j]] <- model[["variable.importance"]]
          numnodes[[j]] <- sum(model[["frame"]][,"var"]!="<leaf>")
        }
      }
      pred <- unlist(pred)[rownames(d)] # combine the fold results
      if (any(is.na(pred))) {
        warning("Some decision tree models were of size 1 (no classification possible). Skip calculating confusion matrix and returning NA instead.")
        conf <- list("overall"=NA)
        #browser()
      } else {
        conf <- caret::confusionMatrix(data=pred, reference=g) # compute confusion matrix for this replicate
      }
      if (method=="svm") return(list("Missclassified"=unlist(misscl), "ConfusionMatrix"=conf))
      if (method=="C50") return(list("Missclassified"=unlist(misscl), "ConfusionMatrix"=conf, "Importance"=unlist(sapply(imp,rownames))))
      if (method=="rpart") return(list("Missclassified"=unlist(misscl), "ConfusionMatrix"=conf, "Importance"=unlist(imp), "NumberOfNodes"=unlist(numnodes)))
    })
    close(pb)
  }
  return(tmp)
}
