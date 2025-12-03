#' read_msdial
#'
#' @param file The alignment file as exported from MSDial.
#' @param sam An additional sample list (data.frame) where one column can be matched to the file names provided by MSDial.
#'
#' @return A list with at least 3 components (`sam`, `met` and `raw`) containing the relevant information on the experiment.
#' @export
#'
read_msdial <- function(file, sam = NULL) {
  # check if this is a GCMS or LCMS MSDial format
  # ToDo

  # extract sample information
  x <- plyr::ldply(readLines(file)[1:5], function(y) { strsplit(y, "\t")[[1]] })
  x <- t(as.matrix(x[,!(unlist(x[1,]) %in% c("", "NA"))]))
  colnames(x) <- x[1,]
  colnames(x)[3:5] <- c("Order","Batch","File")
  x <- as.data.frame(x[-1,])
  # match with sam if provided as parameter
  if (!is.null(sam)) {

  }
  sam <- cbind("ID" = paste0("S", formatC(1:nrow(x), width = ceiling(log10(nrow(x)+1)), flag = "0")), x)
  sam[,"Class"] <- factor(sam[,"Class"])
  sam[,"Order"] <- as.numeric(sam[,"Order"])
  sam[,"Batch"] <- factor(sam[,"Batch"])
  rownames(sam) <- sam[,"ID"]

  # extract metabolite and raw data information
  x <- utils::read.table(file, sep="\t", header=T, as.is=T, quote="", comment.char = "", skip=4, check.names = FALSE)
  met <- x[,1:(min(which(colnames(x) %in% sam[,"File"]))-1)]
  colnames(met)[2:4] <- c("rt", "mz", "name")
  met$mz <- round(met$mz, 4)
  met <- cbind("ID" = paste0("M", formatC(1:nrow(met), width = ceiling(log10(nrow(met)+1)), flag = "0")), met)
  rownames(met) <- met[,"ID"]


  raw <- t(as.matrix(x[,colnames(x) %in% sam[,"File"]]))
  rownames(raw) <- rownames(sam)
  colnames(raw) <- rownames(met)

  dat <- log10(raw)
  dat[!is.finite(dat)] <- NA

  return(list("sam" = sam, "met" = met, "raw" = raw, "dat" = dat))
}
