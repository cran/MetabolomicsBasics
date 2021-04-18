#'@title spectra_format_converter.
#'
#'@description
#'\code{spectra_format_converter} will generate a matrix with mz and int out of a text representation of a spectrum.
#'
#'@details
#'none.
#'
#'@param txt Sample table.
#'@param m_prec Mass precision of output spectrum.
#'@param i_prec Intensity precision of output spectrum.
#'
#'@return
#'Matrix with mz and int columns.
#'
#'@examples
#'spectra_format_converter(txt="57.1:100 58.0001:10")
#'spectra_format_converter(txt="58.0001:10 57.1:100", m_prec=4)
#'
#'@export
#'
spectra_format_converter <- function(txt=NULL, m_prec=3, i_prec=0) {
  s <- plyr::ldply(strsplit(txt, " ")[[1]], function(x) {
    x <- strsplit(x,":")[[1]]
    return(data.frame("mz"=as.numeric(x[1]), "int"=as.numeric(x[2])))
  })
  s <- as.matrix(s[order(s[,"mz"]),])
  s[,1] <- round(s[,1],m_prec)
  s[,2] <- round(s[,2],i_prec)
  return(s)
}
