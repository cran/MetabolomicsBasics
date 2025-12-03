#' @title msconvert.
#' @description \code{msconvert} is calling ProteoWizards MSConvert as a command line tool on Windows.
#' @details It is a quick and dodgy function to show how to convert vendor MS data
#'   into an open format (mzML). You will have to download/install MSConvert prior
#'   to usage, and probably adjust the arguments according to your needs. Arguments
#'   are documented here \url{https://proteowizard.sourceforge.io/tools/msconvert.html}.
#'   If you don't know where the msconvert.exe is installed you can check for the correct
#'   path using \code{list.files(path="C:/", pattern="^msconvert.exe$", recursive = TRUE)}.
#'   You may read on the potential arguments calling msconvert with the argument '--help',
#'   system2(command = 'path/to/msconvert.exe', args = "--help").
#' @param files A character vector of MS data files (wiff, raw, d, ...). If unspecified, function will return the args to msconvert.
#' @param msc_exe The path to the installed `msconvert.exe`.
#' @param out_path Specify valid path or keep NULL to write output in same folder as input.
#' @param args The arguments passed to msconvert on the command line (see details for documentation). Can be a keyword describing a predefined set (currently: 'default' and 'sciex_MS2' are available)
#' @return Only some informative outputs printed to the console. The specified MS data
#'   files will be converted to mzML within the input folder.
#' @examples
#' msconvert()
#'
#' @importFrom tools file_ext
#' @export
msconvert <- function(
  files = NULL,
  msc_exe = "C:\\Program Files\\ProteoWizard\\ProteoWizard 3.0.11856\\msconvert.exe",
  out_path = NULL,
  args = c(
    '--filter \"peakPicking cwt snr=0.01 peakSpace=0.1 msLevel=1\"',
    '--filter \"scanTime [0,3600]\"',
    '--filter \"metadataFixer\"',
    "--mzML",
    "--32",
    "--zlib"
  )
) {
  # predefined arg sets for msconvert; can be extended or could be substituted against a convenience function providing parameters more R friendly
  arg_sets <- list(
    "default" = c(
      '--filter \"peakPicking cwt snr=0.01 peakSpace=0.1 msLevel=1\"',
      '--filter \"scanTime [0,3600]\"',
      '--filter \"metadataFixer\"',
      "--mzML",
      "--32",
      "--zlib"
    ),
    "sciex_MS2" =  c(
      '--filter \"peakPicking vendor snr=0.01 peakSpace=0.1 msLevel=1-2\"',
      '--filter \"scanTime [0,3600]\"',
      '--filter \"metadataFixer\"',
      '--filter \"defaultArrayLength 2-\"',
      '--numpressLinearAbsTol 1e-4',
      "--mzML",
      "--32",
      "--zlib"
    ),
    "sciex_RTdrift" =  c(
      '--filter \"peakPicking vendor snr=0.01 peakSpace=0.1 msLevel=1-2\"',
      '--filter \"scanTime [0,3600]\"',
      '--filter \"metadataFixer\"',
      '--numpressLinearAbsTol 1e-4',
      "--mzML",
      "--32",
      "--zlib"
    )
  )
  if (length(args)==1 && as.character(args) %in% names(arg_sets)) {
    args <- arg_sets[[args]]
  }
  if (is.null(files)) {
    message("currently avialible arg sets are '", paste(names(arg_sets), collapse="' or '"), "'")
    message("Returning the current arg set parameters")
    return(args)
  }
  if (!file.exists(msc_exe)) {
    if (.Platform$OS.type == "windows") {
      message("You did not specify the path to 'msconvert.exe'. I try to find it for you...")
      msc_exe <- list.files(path = "C:/", pattern = "^msconvert.exe$", recursive = TRUE, full.names = TRUE)
      if (length(msc_exe)>=2) {
        message(paste(msc_exe, collapse = "\n"))
        message("Found several versions of 'msconvert.exe' and will use: ", msc_exe[1])
        msc_exe <- msc_exe[1]
      }
    } else {
      stop("Please provide a valid path to 'msconvert.exe'. You can use function 'list.files()' to search on your system.")
    }
  }
  if (!all(sapply(files, file.exists))) stop("Not all of your specified files exist.")
  in_file_ext <- unique(tools::file_ext(files))
  if (!length(in_file_ext) == 1) stop("You did provide different MS formats. That probably does not make sense regarding conversion parameters")
  if (is.null(out_path) || !file.exists(out_path)) {
    out_path <- sapply(files, function(x) { fs::path(gsub(basename(x), "", x)) })
  } else {
    out_path <- fs::path(out_path)
  }
  #outfiles <- fs::path(out_path, gsub(paste0(in_file_ext, "$"), "mzML", basename(files)))
  for (i in 1:length(files)) {
    message(files[i])
    system2(command = msc_exe, args = c(files[i], unname(ifelse(out_path[i]=="", "", paste0('--outdir \"', out_path[i], '\"'))), args))
  }
  invisible(NULL)
}
