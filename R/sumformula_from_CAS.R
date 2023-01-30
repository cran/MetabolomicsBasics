#' @title sumformula_from_CAS.
#' @description \code{sumformula_from_CAS} will try to get a chemical sum formula from a CAS ID.
#' @details tbd.
#' @param x Vector of CAS IDs.
#' @return A character vector of length input vector.
#' @examples
#' \dontrun{
#' x <- readLines("C:/Users/jlisec/Documents/Francesco Russo/RECTOX/RECTOX_GC-EI-MS_CASRN")
#' sf <- sumformula_from_CAS(x = x)
#' }
#' @export
#' @importFrom webchem pc_synonyms get_cid pc_prop
sumformula_from_CAS <- function(x = NULL) {
  flt <- webchem::is.cas(x)
  stopifnot(any(flt))
  out <- rep(NA, length(x))
  names(out) <- x
  # res <- vector(mode = "list", length = length(x))
  # for (i in 1:length(res)) {
  #   pc_syn_name_cas <- unlist(webchem::pc_synonyms(x[i], from = "name", match = "first"))
  #   ids <- webchem::get_cid(pc_syn_name_cas, from = "name", domain = "substance", match = "first")
  #   if (is.na(ids$cid)) {
  #     out <- as.data.frame(matrix(NA, nrow=1, ncol=4, dimnames=list(i, c("IUPACName", "MolecularFormula", "MonoisotopicMass", "InChIKey"))))
  #   } else {
  #     out <-  webchem::pc_prop(cid = ids$cid, properties = c("IUPACName", "MolecularFormula", "MonoisotopicMass", "InChIKey"))
  #   }
  # }
  out[flt] <- plyr::ldply(x[flt], function(cas) {
    #print(cas)
    pc_syn_name_cas <- unlist(webchem::pc_synonyms(cas, from = "name", match = "first"))
    if (!is.na(pc_syn_name_cas)) {
      ids <- webchem::get_cid(pc_syn_name_cas, from = "name", domain = "substance", match = "first")
      mf <- try(webchem::pc_prop(cid = ids$cid, properties = c("MolecularFormula")))
      ifelse(inherits(mf, "try-error"), NA, mf[,2])
    } else {
      NA
    }
  }, .progress = "text")[,1]
  return(out)
}
