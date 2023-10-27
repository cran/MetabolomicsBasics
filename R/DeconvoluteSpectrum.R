#' @title DeconvoluteSpectrum.
#' @description \code{DeconvoluteSpectrum} will evaluate a list of xcmsRaw objects
#'   at a given time (rt) and potential mass (mz1). The main purpose is
#'   to deconvolute the mass spectrum at rt including mz1.
#' @details Will test all mz at spectrum of base peak within range for co-apex,
#'   rt diff and ratio consistency/correlation over a set of samples.
#' @param dat A list of xcmsRaws or an xcmsSet object.
#' @param rt Retention time to search for maxima.
#' @param rt_dev Allowed retention time window.
#' @param mz1 If specified, ensure that this mass is included in the spectrum (assumed base peak). NULL otherwise.
#' @param mz_dev Allowed mz deviation [Da].
#' @param use.mz.adjust Will adjust mz on an experiment wide basis.
#' @param selected Numeric index vector. Use only this subset of files (to save time).
#' @param silent Supresses warnings and progress bar if TRUE.
#' @return A pseudo spectrum at rt (containing mz1 if specified). Effectively
#'   a 2-column matrix (mz, int) with rt as attribute.
#' @examples
#' # needs a list of xcmsRaws
#' @keywords internal
#' @noRd
#' @importFrom graphics par plot points legend
#' @importFrom stats median cor sd
#' @importFrom utils txtProgressBar setTxtProgressBar
DeconvoluteSpectrum <- function(dat = NULL, rt = NULL, rt_dev = 2, mz1 = NULL, mz_dev = 0.003, use.mz.adjust = FALSE, selected = NULL, silent = FALSE) {
  # Helper functions (internal)
  GetStableRatio <- function(dat = NULL, mz1 = NULL, mz2 = NULL, mz_dev = NULL, rt = NULL, rt_dev = NULL, cutoff = 200, plotting = FALSE) {
    # evaluate a stable int-ratio between two mz within a rt-range
    x <- xcms::rawEIC(dat, mzrange = matrix(mz1 + c(-mz_dev, mz_dev), ncol = 2), rtrange = matrix(rt + c(-rt_dev, rt_dev), ncol = 2))
    y <- xcms::rawEIC(dat, mzrange = matrix(mz2 + c(-mz_dev, mz_dev), ncol = 2), rtrange = matrix(rt + c(-rt_dev, rt_dev), ncol = 2))
    z <- matrix(c(dat@scantime[x[["scan"]]], x[["intensity"]], y[["intensity"]]), ncol = 3, dimnames = list(NULL, c("rt", "int_mz1", "int_mz2")))
    f <- z[, 3] > cutoff & z[, 2] > cutoff
    if (plotting) {
      opar <- par(no.readonly = TRUE)
      par(mfrow = c(1, 3))
      plot(x = z[f, 1], y = z[f, 2], ylim = c(cutoff, max(z[f, c(2, 3)])))
      points(x = z[f, 1], y = z[f, 3], col = 2)
      legend(x = "topleft", fill = c(1, 2), legend = round(c(mz1, mz2), 3))
      plot(x = z[f, 1], y = z[f, 2], ylim = c(cutoff, max(z[f, c(2, 3)])), log = "y")
      points(x = z[f, 1], y = z[f, 3], col = 2)
      # hist(z[f,3]/z[f,2])
      plot(x = z[f, 3], y = z[f, 2], log = "xy")
      par(opar)
    }
    if (sum(f) >= 6) {
      out <- median(z[f, 3] / z[f, 2])
      # names(out) <- round(sd(z[f,3]/z[f,2]),2)
      names(out) <- round(cor(log10(z[f, 3]), log10(z[f, 2]), use = "c", method = "p"), 2)
    } else {
      out <- NA
      names(out) <- NA
    }
    return(out)
  }
  GetMaxMzInRange <- function(dat = NULL, rt = NULL, rt_dev = NULL, mz_exclude = NULL, mz_exclude_dev = 0.003, ensure_apex = 3, plotting = FALSE) {
    sapply(dat, function(x) {
      s <- which(abs(x@scantime - rt) <= rt_dev)
      if (length(s) >= 1) {
        i <- x@env$intensity[(1 + x@scanindex[min(s)]):x@scanindex[max(s)]]
        m <- x@env$mz[(1 + x@scanindex[min(s)]):x@scanindex[max(s)]]
        # determine if mz is within exclude masses
        if (!is.null(mz_exclude)) {
          fe <- !sapply(m, function(m_local) {
            any(abs(m_local - mz_exclude) < mz_exclude_dev)
          })
        } else {
          fe <- rep(T, length(m))
        }
        # determine if mz is likely of coeluting peak
        if (!is.null(ensure_apex)) {
          count <- 0
          w <- ensure_apex # number of scans from search window boarder which shall not contain the max int
          border_scans <- c((1 + x@scanindex[min(s)]):x@scanindex[min(s) + w], (1 + x@scanindex[max(s) - w]):x@scanindex[max(s)])
          border_scans <- border_scans - min(border_scans) + 1
          border_scans <- 1:length(m) %in% border_scans
          if (max(i[fe]) >= 1000) {
            fb <- which.max(i[fe])
            while (border_scans[fe][fb]) {
              test <- !sapply(m[fe], function(m_local) {
                any(abs(m_local - m[fe][fb]) < mz_exclude_dev)
              })
              fe[fe] <- fe[fe] & test
              count <- count + 1
              fb <- which.max(i[fe])
            }
          }
          # print(count)
        }
        if (max(i[fe]) >= 1000) {
          # determine if mz fulfills certain S/N ratio
          if (plotting) {
            m_test <- m[fe][which.max(i[fe])][1]
            test <- sapply(s, function(s_local) {
              ind <- ((1 + x@scanindex[s_local]):x@scanindex[s_local + 1]) - x@scanindex[min(s)]
              wind <- abs(m[ind] - m_test) < mz_exclude_dev
              ifelse(any(wind), i[ind][wind], NA)
            })
            plot(y = test, x = x@scantime[s], main = round(m_test, 4), ylab = "int")
          }
          # i[abs(m-m[which.max(i)][1])<mz_exclude_dev]
          return(m[fe][which.max(i[fe])][1])
        } else {
          return(NA)
        }
      } else {
        return(NA)
      }
    })
  }
  GetSpectrum <- function(x = NULL, rt = NULL, cutoff = 200) {
    s <- which.min(abs(x@scantime - rt))
    i <- x@env$intensity[(1 + x@scanindex[s]):x@scanindex[s + 1]]
    m <- x@env$mz[(1 + x@scanindex[s]):x@scanindex[s + 1]]
    return(m[i > cutoff])
  }
  GetMeanMz <- function(dat = NULL, mz = NULL, mz_dev = NULL, rt = NULL, rt_dev = NULL, method = NULL, cutoff = 200, plotting = FALSE) {
    median(sapply(dat, function(x) {
      s <- which(abs(x@scantime - rt) <= rt_dev)
      if (length(s) >= 1) {
        a <- sapply(s, function(k) {
          i <- x@env$intensity[(1 + x@scanindex[k]):x@scanindex[k + 1]]
          m <- x@env$mz[(1 + x@scanindex[k]):x@scanindex[k + 1]]
          f <- which.min(abs(m - mz))
          c(m[f], i[f])
        })
        # f <- abs(a[1,]-mz)<mz_dev & a[2,]>(cutoff+max(a[2,])/1000)
        f <- abs(a[1, ] - mz) < mz_dev & a[2, ] > cutoff & a[2, ] < 10^5
        # f <- abs(a[1,]-mz)<mz_dev & a[2,]>200 & a[2,]<10^5
        return(median(a[1, f]))
      } else {
        return(NA)
      }
    }), na.rm = T)
  }

  # check for xcms to be able to keep it in suggested packages
  verify_suggested("xcms")

  # convert xcmsSet to list of XCMSraw's
  if (inherits(dat, "xcmsSet")) {
    file_in <- dat@filepaths
    if (!all(grep("mzXML$", file_in) == (1:length(file_in)))) {
      stop("The filenames stored in dat@filepaths do not all end on mzXML. Stop processing.")
    }
    if (is.null(selected)) selected <- 1:length(file_in)
    if (!silent) pb <- txtProgressBar(min = 0, max = length(selected), title = "Loading XCMSraw's", style = 3)
    dat <- vector("list", length = length(selected))
    for (i in 1:length(selected)) {
      if (!silent) setTxtProgressBar(pb, i)
      dat[[i]] <- Load_mzML(file_in = file_in[selected[i]], removeStart = rt - 15, removeEnd = rt + 15, tc = 1)
    }
    if (!silent) close(pb)
  }

  if (is.null(selected)) selected <- 1:length(dat)

  # prepare output result
  out <- list()

  # check if rt/rt_dev definition is valid
  # sapply(dat, function(x) {all(((range(x@scantime)-(rt+c(-1,1)*rt_dev)) >= c(0,0)) == c(F,T))})
  # check_rt <- sapply(dat, function(x) {all(min(x@scantime) <= rt-rt_dev & max(x@scantime) >= rt+rt_dev)})
  # if (!all(check_rt)) warning(paste(sum(!check_rt), "datasets do not cover desired rt window."))

  # find max mz
  if (is.null(mz1)) mz1 <- median(GetMaxMzInRange(dat = dat, rt = rt, rt_dev = rt_dev))

  # determine rt and int of max mz over exp
  out[[1]] <- FindMaxInt(dat = dat, "mz" = mz1, "rt" = rt, rt_dev = rt_dev, mz_dev = mz_dev, test_plot = FALSE)
  # remove files which are empty in this region (e.g. no data recorded within time window or for defined mz-window)
  filt_files <- which(out[[1]] == 0 | is.na(out[[1]]))
  if (length(filt_files) >= 1) {
    if (length(filt_files) == length(dat)) {
      warning("[DeconvoluteSpectrum] No valid intensities for this mz/rt combination.")
      return(NULL)
    } else {
      if (!silent) warning(paste("[DeconvoluteSpectrum] removed", length(filt_files), "of", length(dat), "data files being empty at rt =", rt, "for mz =", mz1), call. = FALSE)
      out[[1]] <- out[[1]][-filt_files]
      dat <- dat[-filt_files]
    }
  }
  # readjust rt and rt_dev based on data
  rt <- median(as.numeric(names(out[[1]]))[out[[1]] > 0], na.rm = T)
  # ==========================================================================================
  # [ToDo] think about readjusting rt and rt_dev based on this result
  # ==========================================================================================
  rt_dev <- ceiling(max(abs(as.numeric(names(out[[1]]))[out[[1]] > 0] - rt)))
  # browser()
  # if rt_dev is still high (shouldn't necessary be larger than 2 for well aligned samples) try to fix
  if (rt_dev > 2) {
    if (!silent) warning("Time window of peak apexes was larger than 2sec.")
    rt_dev <- 2
    out[[1]] <- FindMaxInt(dat = dat, "mz" = mz1, "rt" = rt, rt_dev = rt_dev, mz_dev = mz_dev, test_plot = FALSE)
    # rt <- median(as.numeric(names(out[[1]]))[out[[1]]>0],na.rm=T)
  }
  if (rt_dev == 0) {
    # happens if single samples are processed
    rt_dev <- 0.5
    out[[1]] <- FindMaxInt(dat = dat, "mz" = mz1, "rt" = rt, rt_dev = rt_dev, mz_dev = mz_dev, test_plot = FALSE)
    # rt <- median(as.numeric(names(out[[1]]))[out[[1]]>0],na.rm=T)
  }

  # remove files which are empty in this region (e.g. no data recorded within time window)
  filt_files <- which(out[[1]] == 0 | is.na(out[[1]]))
  if (length(filt_files) >= 1) {
    n <- length(filt_files)
    if (!silent) warning(paste("[DeconvoluteSpectrum] empty data file", ifelse(n > 1, "s", ""), paste(filt_files, collapse = ", "), "at rt =", rt, collapse = " "), call. = FALSE)
    out[[1]] <- out[[1]][-filt_files]
    dat <- dat[-filt_files]
  }

  # determine masses from spectrum
  bpm <- which.max(out[[1]])
  mz2 <- GetSpectrum(x = dat[[bpm]], rt = as.numeric(names(out[[1]]))[bpm], cutoff = min(c(500, 0.1 * max(out[[1]]))))

  if (use.mz.adjust) {
    mz2 <- sapply(mz2, function(mz2local) {
      GetMeanMz(dat = dat, mz = mz2local, mz_dev = mz_dev, rt = rt, rt_dev = rt_dev, cutoff = min(c(200, 0.1 * max(out[[1]]))))
    })
  }

  # determine maxima for these masses from all files
  if (any(is.na(mz2))) warning("func_DeconvoluteSpectrum :: NAs in mz2 object")
  for (j in 1:length(mz2)) {
    out[[1 + j]] <- FindMaxInt(dat = dat, "mz" = mz2[j], "rt" = rt, rt_dev = rt_dev, mz_dev = mz_dev, test_plot = FALSE)
  }

  # get CV(intensity)
  # cv <- function(x) {sd(x,na.rm=T)/mean(x,na.rm=T)}
  # dcv_int <- round(sapply(out[-1], function(x) {cv(x)})-cv(out[[1]]),2) # get diffs of CV for max int

  # get diffs for rt at max int between candidates (mz2 and base peak mz1)
  d_rt <- round(sapply(out[-1], function(x) {
    median(as.numeric(names(x)) - as.numeric(names(out[[1]])), na.rm = T)
  }), 2)

  # get stable ratios and correlation between mz1 and all mz2's
  if (length(names(dat)) >= 1) names(dat) <- NULL # [JL: names have to be removed not to be concatenated with the desired output]
  test <- lapply(mz2, function(mz2local) {
    sapply(dat, function(x) {
      GetStableRatio(dat = x, mz1 = mz1, mz2 = mz2local, mz_dev = mz_dev, rt = rt, rt_dev = rt_dev, cutoff = 200)
    })
  })
  # rng_rat <- round(apply(test, 2, function(x) { diff(range(x, na.rm=T)) }), 4) # get diff(range) for int-ratios to main mz
  # rng_rat <- round(sapply(test, function(x) { median(as.numeric(names(x)),na.rm=T) }), 4) # get cv for int-ratios to main mz
  cor_rat <- round(sapply(test, function(x) {
    median(as.numeric(names(x)), na.rm = T)
  }), 4)

  # get mz intensity levels by scaling mz2 with median ratios
  mn_int <- round(out[[1]][bpm] * sapply(test, function(x) {
    median(x, na.rm = T)
  }))

  # combine these infos into dataframe
  tmp <- cbind(mz2, mn_int, d_rt, cor_rat)
  rownames(tmp) <- 1:nrow(tmp)
  # print(tmp)

  # filter for mz2 belonging to base peak (=mz1) and return pseudo spectrum
  f <- abs(d_rt) < 0.15 & cor_rat > 0.7
  spec <- matrix(c(mz2[which(f)], mn_int[which(f)]), ncol = 2, dimnames = list(NULL, c("mz", "int")))

  attr(spec, "rt") <- rt
  return(spec)
}
