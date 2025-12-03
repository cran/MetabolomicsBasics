#' @title MetaboliteANOVA.
#' @description \code{MetaboliteANOVA} will perform an ANOVA on columns
#'   of a data matrix according to a specified model.
#' @details The function is a wrapper for \code{lm} including some sanity checks.
#'   It will accept a data matrix (traits in columns), sample information (data.frame)
#'   and a potential model as input, compute an ANOVA per column and return
#'   the respective P-values in a named matrix for further plotting or export.
#' @param dat Data  matrix (e.g. of metabolite).
#' @param sam Sample table (same number of row as 'dat' and containing all columns specified in 'model'.
#' @param model ANOVA model. May include +, * and : together with column names of sam (cf. Examples).
#' @param method The method to be used in column wise multiple testing adjustment, see \code{\link{p.adjust}}.
#' @param silent Logical. Shall the function print warnings to the console?
#' @return A named matrix of P-values (rows=metabolites/traits; cols=ANOVA factors).
#' @examples
#' # load raw data and sample description
#' raw <- MetabolomicsBasics::raw
#' sam <- MetabolomicsBasics::sam
#' # compute P-values according to specified ANOVA model (simple and complex)
#' head(m1 <- MetaboliteANOVA(dat = raw, sam = sam, model = "GT"))
#' head(m2 <- MetaboliteANOVA(dat = raw, sam = sam, model = "GT+Batch+Order+MP"))
#' # compare P-values for one factor determined in both models
#' hist(log10(m2[, "GT"]) - log10(m1[, "GT"]), main = "")
#' @importFrom plyr ldply
#' @importFrom stats p.adjust anova lm rnorm
#' @export
MetaboliteANOVA <- function(dat = NULL, sam = NULL, model = NULL, method = "none", silent = FALSE) {
  # extract factors from model
  facs <- strsplit(model, "[+*: ]")[[1]]
  facs <- facs[!(facs %in% "")]
  # convert model into standard R formula object
  fmod <- as.formula(paste("y", model, sep = " ~ "))
  # check if factors are defined in sam
  if (!all(facs %in% colnames(sam))) {
    warning(paste("Factor(s) undefined in 'sam':", paste(facs[!(facs %in% colnames(sam))], collapse = ", ")))
    facs <- facs[facs %in% colnames(sam)]
    stopifnot(length(facs) >= 1)
  }
  # if dat is vector, convert into data.frame before using apply
  ydf <- as.data.frame(dat)
  if (is.null(rownames(ydf))) {
    rownames(ydf) <- 1:nrow(ydf)
  }
  # set up model dataframe and remove empty levels if present
  mdf <- data.frame(sam[, facs, drop = F], row.names = rownames(ydf))
  # re-factor all true factors to remove empty levels if still present
  for (i in 1:ncol(mdf)) {
    if (is.factor(mdf[, i])) {
      mdf[, i] <- factor(mdf[, i])
    }
  }
  # evaluate which are the true factors in the model (not numeric) !! 1+ is necessary because
  true_fac_cols <- 1 + which(sapply(facs, function(i) {
    is.factor(sam[, i])
  }))
  # evaluate total number/names of terms from model (including interactions)
  fac_names <- rownames(stats::anova(stats::lm(fmod, data = data.frame("y" = stats::rnorm(nrow(mdf)), mdf))))
  fac_names <- fac_names[-length(fac_names)]
  # calculate ANOVA P-values consecutively for all columns of dat
  out <- plyr::ldply(1:ncol(ydf), function(j) {
    tdf <- data.frame("y" = ydf[, j], mdf) # set up dataframe for anova
    # check for incomplete sub groups (in case of more than 1 true factor) and return non-normalized values + warning
    if (length(true_fac_cols) >= 1 && any(sapply(split(tdf[, "y"], interaction(tdf[, true_fac_cols], drop = TRUE)), function(x) {
      sum(is.finite(unlist(x))) == 0
    }))) {
      # try to identify metabolite name for warning message
      incomplete_levels <- paste(names(which(sapply(split(tdf[, "y"], interaction(tdf[, true_fac_cols], drop = FALSE)), function(x) {
        sum(is.finite(unlist(x))) == 0
      }))), collapse = "; ")
      if (!silent) {
        message(paste("Return NAs for", colnames(ydf)[j], "(Incomplete levels:", incomplete_levels, ")"))
      }
      pval <- rep(NA, length(fac_names))
    } else {
      y.lm <- lm(fmod, data = tdf) # set up anova model
      pval <- anova(y.lm)[, "Pr(>F)"]
      names(pval) <- rownames(anova(y.lm))
      pval <- pval[-length(pval)]
    }
    return(pval)
  })
  rownames(out) <- colnames(ydf)
  colnames(out) <- fac_names
  out <- apply(out, 2, stats::p.adjust, method = method)
  attr(out, "p.adjust.method") <- method
  return(as.matrix(out))
}
