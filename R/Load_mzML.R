#' @title Load_mzml.
#' @description \code{Load_mzml} will read and manipulate mz(X)ML data files
#'   using the xcmsRaw importer. The main purpose is to strip MS/MS signals
#'   as well as leading and/or trailing time ranges.
#' @details not yet
#' @param file_in A valid filename of a rawdata XML file (e.g. Bruker QTOF Export or anything that can be processed by xcmsRaw) or a vector of such files.
#' @param removePattern Removes all scans according to a T/F pattern. See examples.
#' @param removeStart Remove all Scans *before* this numeric time value.
#' @param removeEnd Remove all Scans *after* this numeric time value.
#' @param tc TimeConversion factor. Your specified time will be multiplied by tc before comparison with stored data. [ToDo:substitute against time scale evaluation based on XML]
#' @param export Boolean. Will write import result to mzData-File if TRUE.
#' @param correct_missing_scans T/F, will exclude the scan directly following a time gap as this seems to clear the collision cell and exagerate intensity data.
#' @param silent Suppress messages and warnings if TRUE.
#' @param profstep Passed to xcms to create ProfileMatrix (if >0) which is required to use some xcms based plotting functions like plotChrom.
#' @return A xcmsRaw object or a list of such if file_in is a vector.
#' @examples
#' # not yet
#' # data(test_mzXML)
#' # str(Load_mzML(file_in=test_mzXML))
#'
#' # removing avery secon scan (=stripping MS/MS)
#' # str(Load_mzML(file_in=test_mzXML, removePattern=c(T,F)))
#'
#' # removing 2 minutes from start and everything after 10 minutes (XML)
#' # str(Load_mzML(file_in=test_mzXML, removeStart=2, removeEnd=10))
#' @keywords internal
#' @noRd
Load_mzML <- function(file_in = NULL, removePattern = T, removeStart = NA, removeEnd = NA, correct_missing_scans = FALSE, tc = 60, export = FALSE, silent = FALSE, profstep = 0) {
  # check for xcms to be able to keep it in suggested packages
  if (!requireNamespace("xcms", quietly = TRUE)) {
    stop("The use of this function requires package 'xcms'. Please ",
         "install with 'Biobase::install(\"xcms\")'")
  }
  # Helper function
  removeScans <- function(x = NULL, s = NULL) {
    scannum <- as.integer(diff(c(x@scanindex, length(x@env$mz))))
    fs <- !((1:length(x@scantime)) %in% s) # scans to keep
    fm <- rep(fs, times = scannum)
    scannum <- scannum[fs]
    x@env$intensity <- x@env$intensity[fm]
    x@env$mz <- x@env$mz[fm]
    x@tic <- x@tic[fs]
    x@acquisitionNum <- 1:sum(fs)
    x@scantime <- x@scantime[fs]
    x@polarity <- x@polarity[fs]
    x@scanindex <- as.integer(c(0, cumsum(scannum)[-length(scannum)]))
    if ("profile" %in% names(x@env)) xcms::profStep(x) <- profstep
    x@scanrange <- c(1, length(x@scantime))
    return(x)
  }
  # main
  out <- vector("list", length(file_in))
  for (i in 1:length(out)) {
    if (file.exists(file_in[i])) {
      # import File and set scans to be removed to NULL
      if (!silent) cat(paste("\nProcessing File", file_in[i]))
      tmp <- xcms::xcmsRaw(file = file_in[i], profstep = profstep)
      s <- NULL

      # Establish filtering due to pattern (to e.g. remove BroadBand MS2 use: c(T,F))
      if (is.logical(removePattern) && any(removePattern) && any(removePattern == FALSE)) {
        # check for zero intensity scans
        aN <- tmp@acquisitionNum
        if (length(aN) != max(aN) && !all(removePattern)) {
          if (!silent) warning(paste0("Potentially ", max(tmp@acquisitionNum) - length(tmp@acquisitionNum), " zero-intensity scans have been removed at import of file ", file_in[i]))
        }
        if (min(aN) != 1) {
          if (!silent) warning(paste0("Acquisition does not start with scan 1 but with scan ", min(aN), " in file ", file_in[i], "which may lead to scan ordering problems."))
        }
        # remove certain scan-pattern manually
        s <- c(s, which(!rep(removePattern, length.out = max(aN))[aN]))
      }
      # Establish filtering for certain range
      if (is.numeric(removeStart) && length(removeStart) == 1 && removeStart > 0 && removeStart < (max(tmp@scantime) / tc)) {
        s <- c(s, which(tmp@scantime < removeStart * tc))
      }
      if (is.numeric(removeEnd) && length(removeEnd) == 1 && removeEnd < (max(tmp@scantime) / tc)) {
        s <- c(s, which(tmp@scantime > removeEnd * tc))
      }
      # Correct intensities after missing scan gaps by removing them
      if (correct_missing_scans) {
        dt <- diff(tmp@scantime)
        s <- c(s, which(c(FALSE, dt > 2 * median(dt))))
      }
      if (!is.null(s)) tmp <- removeScans(x = tmp, s = sort(unique(s)))
      # return including file name
      attr(tmp, "file_in") <- file_in[i]
      if (export) xcms::write.mzdata(tmp, gsub(paste0(tools::file_ext(file_in), "$"), "mzData", file_in[i]))
      out[[i]] <- tmp
    } else {
      if (!silent) cat(paste("\nFile", file_in[i], "not found..."))
      tmp <- NA
      attr(tmp, "file_in") <- file_in[i]
      out[[i]] <- tmp
    }
  }
  if (!silent) cat("\nFinished...\n\n")
  # return as list if many or unlisted if single
  if (length(out) == 1) {
    return(out[[1]])
  } else {
    return(out)
  }
}
