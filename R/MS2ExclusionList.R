#' @title MS2ExclusionList.
#' @description \code{MS2ExclusionList} will compute an exclusion list for MS2
#'     experiments based on an mzML file.
#' @param x A wiff file converted to mzML and imported by xcmsRaw.
#' @param ms2_previous A two column data.frame containing exclusion peaks from previous measurements or known contaminants.
#' @return A two column data.frame that can be written to txt and imported in Sciex-Analyst Software.
#' @noRd
#' @keywords internal
#' @examples
#' # example code
#' \dontrun{
#' fl <- paste0("C:/Users/jlisec/Documents/", c("B1-7-B34.wiff","Blank1_1.wiff")[2])
#' msc_exe <- paste0("C:/Users/jlisec/SoftwareLokal/",
#'   "pwiz-bin-windows-x86_64-vc143-release-3_0_23306_2250ca7.tar/",
#'   "pwiz-bin-windows-x86_64-vc143-release-3_0_23306_2250ca7/msconvert.exe")
#' msconvert(files = fl, msc_exe = msc_exe, out_path = "C:/Users/jlisec/Documents", args = "sciex_MS2")
#' x <- xcms::xcmsRaw(filename = gsub("wiff", "mzML", basename(fl)), profstep = 0, includeMSn = TRUE)
#' tmp.mzML <- "C:/Users/jlisec/Documents/B1-7-B34.mzML"
#' y <- xcms::xcmsRaw(filename = tmp.mzML, profstep = 0, includeMSn = TRUE)
#' el_blank <- MS2ExclusionList(x)
#' el_sample <- MS2ExclusionList(y)
#' }

MS2ExclusionList <- function(
  x = NULL, ms2_previous = NULL
) {
  # potential parameters
  mz_dev = 0.01
  ms2_int_limit <- 500

  # [1] extract previously obtained mz@rt
  if (length(x@msnPrecursorMz)>=1) {
    ms2_fl <- data.frame("mz_Da"=round(x@msnPrecursorMz, 4), "rt_min"=round(x@msnRt/60, 2))
  } else {
    ms2_fl <- NULL
  }
  mz_groups <- InterpretMSSpectrum::GetGroupFactor(x = ms2_fl[,1], gap = mz_dev)
  rt_groups <- InterpretMSSpectrum::GetGroupFactor(x = ms2_fl[,2], gap = 3/60)
  # remove redundancy of multiple MS2 scans
  ms2_fl <- plyr::ldply(split(ms2_fl, factor(paste(mz_groups, rt_groups, sep="_"))), function(x) {
    apply(x, 2, mean)
  }, .id = NULL)
  ms2_fl[,1] <- round(ms2_fl[,1], 4)
  ms2_fl[,2] <- round(ms2_fl[,2], 2)

  # [2] identify contaminants from solvents
  mz_test <- seq(x@mzrange[1], x@mzrange[2], 1.5*mz_dev)
  bpc <- HiResTEC::getMultipleBPC(x = x, mz = mz_test, mz_dev = mz_dev)
  md <- attr(bpc, "mass_defect")
  mz_test <- unique(round(na.omit(mz_test + md/1000),3))
  bpc <- HiResTEC::getMultipleBPC(x = x, mz = mz_test, mz_dev = mz_dev)
  flt <- apply(bpc, 2, function(x) { sum(is.finite(x))>=0.8*nrow(bpc) && max(x, na.rm=TRUE)>ms2_int_limit })
  md <- attr(bpc, "mass_defect")
  mz_test <- unique(round(na.omit(mz_test[flt] + md[flt]/1000),3))
  #bpc <- HiResTEC::getMultipleBPC(x = x, mz = mz_test, mz_dev = mz_dev/2)
  #HiResTEC::plotBPC(bpc = list(bpc))
  #table(flt)
  #bpc[1:4,flt][,1:10]
  ms2_cont <- data.frame("mz_Da"=mz_test, "rt_min"=0)

  # remove recorded from ms2_fl if present in ms2_cont
  flt <- sapply(ms2_fl$mz_Da, function(x) {
    any(abs(ms2_cont$mz_Da-x) < (mz_dev/2))
  })
  ms2_fl <- ms2_fl[flt,][order(ms2_fl[flt,1]),]

  # pdf("check_pos.pdf")
  # for (i in which(!flt)[order(ms2_fl[!flt,2])]){
  temp_file <- tempfile(pattern = paste(tools::file_path_sans_ext(as.character(x@filepath)), "cont_", sep="_"), fileext = ".pdf")
  if (any(flt)) {
    grDevices::pdf(temp_file)
      for (i in which(flt)[order(ms2_fl[flt,2])]) {
        mz <- ms2_fl[i,1]
        rt <- ms2_fl[i,2]*60
        HiResTEC::plotBPC(
          bpc = list(
            HiResTEC::getMultipleBPC(x = x, mz = mz+c(-1:3)*1.0034, mz_dev = mz_dev, rt = rt, rt_dev=20, smooth = 5),
            HiResTEC::getMultipleBPC(x = x, mz = mz+c(-1:3)*1.0034, mz_dev = mz_dev, smooth = 5)
          ), type = "both", col = c(grey(0.8),1:4), ids = paste0(i, " (", mz, "@", rt, ")")
        )
        scan_mz <- xcms::getScan(x, which.min(abs(x@scantime-rt)))
        #scan_mz[abs(scan_mz[,1]-mz)<0.01,]
        #scan_int_lim <- scan_mz[which.min(abs(scan_mz[,1]-mz)),2]
        scan_int_lim <- 0.1*max(scan_mz[,2])
        #scan_mz[scan_mz[,2]>=scan_int_lim,]
        check_mz <- scan_mz[scan_mz[,2]>=scan_int_lim,1]
        HiResTEC::plotBPC(bpc = list(HiResTEC::getMultipleBPC(x = x, mz = check_mz, mz_dev = mz_dev, rt = rt, rt_dev=20, smooth = 5)), ids = i, ann = "mz")
      }
    grDevices::dev.off()

    # Öffne die PDF-Datei mit dem Standardprogramm
    if (.Platform$OS.type == "windows") {
      shell.exec(temp_file)
    } else if (Sys.info()["sysname"] == "Darwin") {
      system2("open", temp_file)
    } else {
      system2("xdg-open", temp_file)
    }
  }

  temp_file2 <- tempfile(pattern = paste(tools::file_path_sans_ext(as.character(x@filepath)), "good_", sep="_"), fileext = ".pdf")
  if (any(!flt)) {
    grDevices::pdf(temp_file2)
    for (i in which(!flt)[order(ms2_fl[!flt,2])]) {
      mz <- ms2_fl[i,1]
      rt <- ms2_fl[i,2]*60
      HiResTEC::plotBPC(
        bpc = list(
          HiResTEC::getMultipleBPC(x = x, mz = mz+c(-1:3)*1.0034, mz_dev = mz_dev, rt = rt, rt_dev=20, smooth = 5),
          HiResTEC::getMultipleBPC(x = x, mz = mz+c(-1:3)*1.0034, mz_dev = mz_dev, smooth = 5)
        ), type = "both", col = c(grey(0.8),1:4), ids = paste0(i, " (", mz, "@", rt, ")")
      )
      scan_mz <- xcms::getScan(x, which.min(abs(x@scantime-rt)))
      #scan_mz[abs(scan_mz[,1]-mz)<0.01,]
      #scan_int_lim <- scan_mz[which.min(abs(scan_mz[,1]-mz)),2]
      scan_int_lim <- 0.1*max(scan_mz[,2])
      #scan_mz[scan_mz[,2]>=scan_int_lim,]
      check_mz <- scan_mz[scan_mz[,2]>=scan_int_lim,1]
      HiResTEC::plotBPC(bpc = list(HiResTEC::getMultipleBPC(x = x, mz = check_mz, mz_dev = mz_dev, rt = rt, rt_dev=20, smooth = 5)), ids = i, ann = "mz")
    }
    grDevices::dev.off()

    # Öffne die PDF-Datei mit dem Standardprogramm
    if (.Platform$OS.type == "windows") {
      shell.exec(temp_file2)
    } else if (Sys.info()["sysname"] == "Darwin") {
      system2("open", temp_file2)
    } else {
      system2("xdg-open", temp_file2)
    }
  }


  return(rbind(ms2_fl, ms2_previous, ms2_cont))
}
