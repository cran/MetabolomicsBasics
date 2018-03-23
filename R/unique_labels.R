#'@title unique_labels.
#'
#'@description
#'\code{unique_labels} will generate a dataframe with color and plotting character specification out of a sample table definition.
#'
#'@details
#'If a color/symbol specification exists for a sample set containing replicate groups this function will help in retrieving this information per group which is useful in boxplot or legend functions (cf. examples).
#'
#'@param sam Sample table.
#'@param g Either column name from sam containing factor column or factor of same length as sam.
#'
#'@return
#'Dataframe with group levels names and their color and plotting character specification.
#'
#'@examples
#'utils::data(raw, package = "MetabolomicsBasics")
#'utils::data(sam, package = "MetabolomicsBasics")
#'unique_labels(sam=sam, g="GT")
#'
#'@export
#'
unique_labels <- function(sam=NULL, g=NULL) {
  stopifnot(c("pchs","cols") %in% colnames(sam))
  if (length(g)==1 && is.character(g) && g%in%colnames(sam)) g <- factor(sam[,g])
  stopifnot(is.factor(g) & length(g)==nrow(sam))
  tmp <- data.frame("Level"=levels(g),
                    "pchs"=sapply(split(sam[,"pchs"], g), function(x) { sort(unique(x))[1] }),
                    "cols"=sapply(split(sam[,"cols"], g), function(x) { sort(unique(x))[1] }),
                    stringsAsFactors = FALSE)
  return(tmp)
}
