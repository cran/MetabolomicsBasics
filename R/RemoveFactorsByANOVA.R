#'@title RemoveFactorsByANOVA.
#'
#'@description
#'\code{RemoveFactorsByANOVA} will remove variance from data using an ANOVA model.
#'
#'@details
#'not yet
#'
#'@param y Data vector (or data matrix) to normalize (numeric + in same order as sam).
#'@param sam data.frame containing the factors or numerical vars for ANOVA model.
#'@param fmod Full model describing the experimental setting (provided as character string).
#'@param kmod Reduced model describing all the biological factors to keep (provided as character string).
#'@param output Should be \code{y_norm} in general but can be switched for testing.
#'@param remove_outliers Should be a numeric integer x (with $x=0$ : no effect; $x>=1$ remove all values which have error e with $e > abs(mean + x * sd)$ ).
#'
#'@return
#'Depends on \code{output}. Usually the normalized data vector (or matrix).
#'
#'@examples
#'# set up sample information
#'sam <- data.frame("GT"=gl(4,10),
#'                  "TR"=rep(gl(2,5),4),
#'                  "Batch"=sample(gl(2,20)),
#'                  "Order"=sample(seq(-1,1,length.out=40)))

#'# set up artificial measurement data
#'set.seed(1)
#'m1=c(5,6,2,9)[sam$GT]+c(-2,2)[sam$TR]+c(-3,3)[sam$Batch]+3*sam$Order+rnorm(nrow(sam), sd=0.5)
#'m2=c(5,-6,2,4)[sam$GT]+c(-2,2)[sam$TR]-5*sam$Order+rnorm(nrow(sam), sd=0.8)
#'dat <- data.frame(m1,m2)
#'
#'# apply function to remove variance
#'# full model incorporating all relevant factors defined in sample table
#'fmod="GT*TR+Batch+Order"
#'# reduced model: factors to be kept from full model; everything elso will be removed from the data
#'kmod="GT*TR"
#'RemoveFactorsByANOVA(y=dat[,"m1"], sam=sam, fmod=fmod, kmod=kmod, output="anova_y")
#'RemoveFactorsByANOVA(y=dat[,"m1"], sam=sam, fmod=fmod, kmod=kmod, output="anova_y_norm")
#'
#'@importFrom utils flush.console
#'@importFrom stats residuals
#'@importFrom stats coef
#'@importFrom stats model.matrix
#'@importFrom stats as.formula
#'@importFrom stats anova
#'@importFrom stats lm
#'@importFrom graphics par
#'
#'@export

RemoveFactorsByANOVA <- function(y=NULL, sam=NULL, fmod=NULL, kmod=NULL, output=c("y_norm","y_lm","anova_y","anova_y_norm","boxplot")[1], remove_outliers=0) {
  # Helper function
  print_info <- function(x, labove=1, lbelow=2, pdate=TRUE) {
    stopifnot(is.character(x) && length(x)==1)
    cat(paste(paste(rep("\n",labove),collapse=""), ifelse(pdate,date(),""), "\n", x, paste(rep("\n",lbelow),collapse=""), sep=""))
    flush.console()
  }
  facs <- strsplit(paste("y", fmod, sep=" ~ "), "[~+*: ]")[[1]]; facs <- facs[!(facs %in% c("","0","y"))]
  keep <- strsplit(paste("y", kmod, sep=" ~ "),"[~+*: ]")[[1]]; keep <- keep[!(keep %in% c("","0","y"))]
  fmod <- as.formula(paste("y", fmod, sep=" ~ "))
  if (length(keep)==0) {
    keep <- NULL
    kmod <- NULL
  } else {
    kmod <- as.formula(paste("y", kmod, sep=" ~ "))
  }
  stopifnot(length(facs)>=1)
  if (!is.null(keep) && !all(keep %in% facs)) {
    warning(paste("Factor(s) defined in 'kmod' need to be additionally specified in 'fmod' to give consistent ordering:", paste(keep[!(keep %in% facs)], collapse=", ")))
    stopifnot(all(keep %in% facs))
  }
  if (!all(facs %in% colnames(sam))) {
    warning(paste("Factor(s) undefined in 'sam':",paste(facs[!(facs %in% colnames(sam))], collapse=", ")))
    facs <- facs[facs %in% colnames(sam)]
    stopifnot(length(facs)>=1)
  }
  # initialize error counters
  count_unchanged <- 0
  count_outfiltfail <- 0
  outer_envir <- environment()
  # if y is vector, convert into data.frame before using apply and ensure rownames are present
  ydf <- as.data.frame(y); if (is.null(rownames(ydf))) rownames(ydf) <- 1:nrow(ydf)
  # apply normalization column wise
  out <- sapply(1:ncol(ydf), function(j) {
    tdf <- data.frame("y"=ydf[,j], sam[,facs,drop=F], row.names=rownames(ydf)) # set up dataframe for anova
    # check for incomplete sub groups (in case of more than 1 true factor) and return non-normalized values + warning
    true_fac_cols <- which(sapply(colnames(tdf),function(i){is.factor(tdf[,i])}))
    # [ToDo] check if test is correct for incomplete levels
    if (length(true_fac_cols)>=1 && any(sapply(split(tdf[,"y"], interaction(tdf[,true_fac_cols],drop=TRUE)),function(x){sum(is.finite(unlist(x)))==0}))) {
      # try to identify metabolite name for warning message
      incomplete_levels <- paste(names(which(sapply(split(tdf[,"y"], interaction(tdf[,true_fac_cols],drop=FALSE)),function(x){sum(is.finite(unlist(x)))==0}))), collapse="; ")
      print_info(paste("ANOVA Modell inapropriate. Return", colnames(ydf)[j], "without normalization. Incomplete levels:", incomplete_levels),0,0,F)
      yn <- ydf[,j]
      assign("count_unchanged", value=get("count_unchanged", envir=outer_envir) + 1, envir=outer_envir)
    } else {
      y.lm <- lm(fmod, data=tdf) # set up anova model
      re <- residuals(y.lm) # residuals
      # outlier removal
      if (is.numeric(remove_outliers) && length(remove_outliers)==1 && remove_outliers>0) {
        f <- abs(re-mean(re)) > remove_outliers*sd(re) # table(f)
        if (any(f)) {
          # if outliers are to be removed -- is model still valid? (i.e. no empty levels)
          f <- which(!is.na(tdf[,"y"]))[f]
          if (length(true_fac_cols)>=1 && any(sapply(split(tdf[-f,"y"], interaction(tdf[-f,true_fac_cols],drop=TRUE)),function(x){sum(is.finite(unlist(x)))==0}))) {
            print_info(paste("ANOVA Modell inapropriate after outlier removal. Skip this step for", colnames(ydf)[j],"."),0,0,F)
            assign("count_outfiltfail", value=get("count_outfiltfail", envir=outer_envir) + 1, envir=outer_envir)
          } else {
            tdf[f,"y"] <- NA
            y.lm <- lm(fmod, data=tdf) # set up anova model
            re <- residuals(y.lm) # residuals
          }
        }
      } else {
        f <- rep(T, nrow(tdf))
      }
      ce <- coef(y.lm) # these are the coefficients of the individual factors
      if (any(is.na(ce))) {
        ce[is.na(ce)] <- 0
        warning("Some coefficients were NA and had to be set to 0.\nYou probably have nested factors or to few degrees of freedom for your model. Please check.")
      }
      yn <- tdf[,"y"]; names(yn) <- rownames(tdf) # this maintaines any names y might have had
      fi <- is.finite(yn) # this preserves the NAs in y by restoring only the finite values
      if (is.null(kmod)) { # if no biological factors were specified
        tm <- rep(ce[1], nrow(tdf)) # this is the total mean
        yn[fi] <- tm[fi] + re
      } else {
        kmmat <- model.matrix(kmod, data=tdf)
        #if (length(yn)==nrow(kmmat)) yn[fi] <- (kmmat %*% ce[colnames(kmmat)])[fi] + re
        if (length(yn[fi])==nrow(kmmat)) yn[fi] <- (kmmat %*% ce[colnames(kmmat)]) + re else warning("check dimensions of 'kmmat'")
      }
      # correct for the offset
      yn[fi] <- yn[fi] + (mean(ydf[fi,j]) - mean(yn[fi]))
    }
    if (output=="y_lm") print(y.lm)
    if (output=="anova_y") print(stats::anova(y.lm))
    if (output=="anova_y_norm") {
      tdf <- data.frame(y=yn, sam[,facs])
      print(stats::anova(lm(fmod, data=tdf)))
    }
    if (output=="boxplot") {
      graphics::par(mfrow=c(1,2))
      ylm <- range(c(ydf[,j],yn),na.rm=T)
      ylb <- colnames(ydf)[j]
      if (is.null(keep)) {
        true_fac_cols <- which(colnames(sam) %in% facs[which(sapply(facs,function(i){is.factor(sam[,i])}))])
      } else {
        true_fac_cols <- which(colnames(sam) %in% keep[which(sapply(keep,function(i){is.factor(sam[,i])}))])
      }
      if (length(true_fac_cols)>=2) tmp.x <- interaction(sam[,true_fac_cols],sep="_",drop=T) else tmp.x <- sam[,true_fac_cols]
      graphics::plot(ydf[,j] ~ tmp.x, ylim=ylm, xlab=paste("Full:", paste(facs, collapse=", ")), ylab=ylb, las=1)
      if (remove_outliers & length(f>=1)) graphics::points(ydf[f,j] ~ tmp.x[f], col=2)
      graphics::plot(yn ~ tmp.x, ylim=ylm, xlab=paste("Keep:", paste(keep, collapse=", ")), ylab=paste(ylb,"norm",sep="_"), las=1)
      graphics::par(mfrow=c(1,1))
    }
    return(yn)
  })
  if (output=="y_norm") {
    if (is.vector(y)) {
      out <- as.vector(unlist(out))
      if (!is.null(names(y))) names(out) <- names(y)
    } else {
      colnames(out) <- colnames(ydf)
    }
  } else {
    out <- NULL
  }
  if(count_unchanged>=1) print_info(paste("In total", count_unchanged, "metabolites were returned unchanged."), 1, 1, F)
  if(count_outfiltfail>=1) print_info(paste("In total", count_outfiltfail, "metabolites were not filtered for outliers."), 1, 1, F)
  invisible(out)
}
