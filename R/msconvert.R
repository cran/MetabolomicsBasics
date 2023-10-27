#' @title msconvert.
#' @description \code{msconvert} is calling ProteoWizards MSConvert as a command line tool on Windows.
#' @details It is a quick and dodgy function to show how to convert vendor MS data
#'   into an open format (mzML). You will have to download/install MSConvert prior
#'   to usage, and probably adjust the arguments according to your needs. Arguments
#'   are documented here \url{https://proteowizard.sourceforge.io/tools/msconvert.html}.
#'   If you don't know where the msconvert.exe is installed you can check for the correct
#'   path using \code{list.files(path="C:/", pattern="^msconvert.exe$", recursive = TRUE)}.
#' @param files A character vector of MS data files (wiff, raw, d, ...).
#' @param msc_exe The path to the installed `msconvert.exe`.
#' @param args The arguments passed to msconvert on the command line (see details for documentation).
#' @return Only some informative outputs printed to the console. The specified MS data
#'   files will be converted to mzML within the same folder.
#' @importFrom tools file_ext
#' @export
msconvert <- function(
    files = NULL,
    msc_exe = "C:\\Program Files\\ProteoWizard\\ProteoWizard 3.0.11856\\msconvert.exe",
    args = c(
      '--filter \"peakPicking cwt snr=0.01 peakSpace=0.1 msLevel=1\"',
      '--filter \"scanTime [0,3600]\"',
      '--filter \"metadataFixer\"',
      "--mzML",
      "--32",
      "--zlib"
    )) {
  if (!file.exists(msc_exe)) stop("Please provide the correct path to msconvert.exe.")
  if (!all(sapply(files, file.exists))) stop("Not all of your specified files exist.")
  in_file_ext <- unique(tools::file_ext(files))
  if (!length(in_file_ext) == 1) stop("You did provide different MS formats. That probably does not make sense regarding conversion parameters")
  outfiles <- gsub(paste0(in_file_ext, "$"), "mzML", files)
  for (i in 1:length(files)) {
    system2(command = msc_exe, args = c(files[i], paste0('--outfile \"', outfiles[i], '\"'), args))
  }
  invisible(NULL)
}
