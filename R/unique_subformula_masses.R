#' @title unique_subformula_masses
#' @description \code{unique_subformula_masses} will generate a numeric vector of
#'   potential sub formula masses regarding a chemical formula as input.
#' @details tbd.
#' @param fml Chemical formula.
#' @param names Return respective sub formulas as vector names.
#' @param check_validity Filter for chemically valid formulas.
#' @return Numeric vector.
#' @examples
#' # specify a formula and calculate all potential combinatorial masses
#' fml <- c("C6H12O6", "C11H16NO4PS", "C24H51O4P")[1]
#' tmp <- unique_subformula_masses(fml = fml)
#' length(tmp); any(duplicated(tmp))
#' hist(tmp, breaks=seq(floor(min(tmp))-1, ceiling(max(tmp))), main=fml)
#' # do the same as above but check for chemical plausibility
#' tmp2 <- unique_subformula_masses(fml = fml, check_validity=TRUE)
#' length(tmp2)
#' hist(tmp2, breaks=seq(floor(min(tmp2))-1, ceiling(max(tmp2))), main=fml)
#' mz <- 147
#' tmp[abs(tmp-mz)<0.5]
#' tmp2[abs(tmp2-mz)<0.5]
#' @export
#' @importFrom InterpretMSSpectrum CountChemicalElements get_exactmass PlausibleFormula
unique_subformula_masses <- function(fml, names=TRUE, check_validity=FALSE) {
  ele <- InterpretMSSpectrum::CountChemicalElements(x=fml)
  n <- length(ele)
  idx <- matrix(NA, nrow=prod(ele+1), ncol=n)
  for (i in 1:n) {
    idx[,i] <- rep(rep(0:ele[i], each=prod(ele[-(1:i)]+1)), times=prod(ele[-(i:n)]+1))
  }
  # remove all zero row
  idx <- idx[-1,]
  m <- InterpretMSSpectrum::get_exactmass(names(ele))
  out <- apply(idx, 1, function(fac) {
    sum(m*fac)
  })
  if (names | check_validity) {
    sf <- sapply(1:nrow(idx), function(i) {
      flt <- idx[i,]>0; paste(names(ele)[flt], idx[i,flt], collapse="", sep="")
    })
  }
  if (names) names(out) <- sf
  if (check_validity) {
    test <- sapply(sf, InterpretMSSpectrum::PlausibleFormula)
    out <- out[test]
  }
  # if (!is.null(prec) && is.numeric(prec)) {
  #   out <- round(out, prec)
  # }
  return(out)
}
