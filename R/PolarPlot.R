#' @title PolarCoordHeterPlot.
#' @description \code{PolarCoordHeterPlot} will draw a plot in polar coordinates
#'    visualizing heterosis effects according to a layout by Swanson-Wagner,
#'    where plot radius represents log2 of fold change between lowest and highest
#'    genotype and plot angle represents the ratio between lowest, intermediate
#'    and highest genotype.
#' @details See examples.
#' @param x Data matrix with measurement values (traits in rows and genotypes in columns).
#' @param gt Character vector of length=3 indicating P1, F and P2. These are used to  filter by column name from x.
#' @param rev_log If you've log transformed your data, you might want to revert the log transformation.
#' @param exp_fac Expansion factor to increase figure size.
#' @param thr Alpha level used in ANOVA to filter insignificant rows. Keep thr=1 to include all matrix rows.
#' @param plot_lab Show 'text' or 'graph' style labels of the polar sections (or keep 'none' to omit).
#' @param col Provide a color vector of length nrow(x).
#' @return Will generate a plot in polar coordinates and return the x/y coordinates of the data points invisibly.
#' @examples
#' # using the provided experimental data
#' raw <- MetabolomicsBasics::raw
#' sam <- MetabolomicsBasics::sam
#' x <- t(raw)
#' colnames(x) <- sam$GT
#' gt <- c("B73","B73xMo17","Mo17")
#' PolarCoordHeterPlot(x=x, gt=gt, plot_lab="graph", thr=0.01, rev_log=exp(1))
#'
#' coord <- PolarCoordHeterPlot(x=x, gt=gt, thr=0.01, rev_log=exp(1))
#' points(x=coord$x[3], coord$y[3], pch=22, cex=4, col=2)
#' # using random data
#' gt <- c("P1","P1xP2","P2")
#' set.seed(0)
#' x <- matrix(rnorm(150), nrow = 10, dimnames = list(paste0("M",1:10), sample(rep(gt, 5))))
#' x[1:4,1:6]
#' PolarCoordHeterPlot(x=x, gt=gt)
#' # using text style labels for the sections
#' PolarCoordHeterPlot(x=x, gt=gt, plot_lab="text", exp_fac=0.75)
#' # reverting the order of parental genotypes
#' PolarCoordHeterPlot(x=x, gt=c("P2","P1xP2","P1"), plot_lab="text", exp_fac=0.75)
#' # using graph style labels for the sections
#' PolarCoordHeterPlot(x=x, gt=c("P2","P1xP2","P1"), plot_lab="graph")
#' # coloring data points
#' PolarCoordHeterPlot(x=x, gt=gt, col=1:10)
#' # applying ANOVA P value threshold to input rows
#' PolarCoordHeterPlot(x=x, gt=gt, col=1:10, thr=0.5)
#' PolarCoordHeterPlot(x=x, gt=gt, plot_lab="graph", col=1:10, thr=0.5)
#' @export
#' @importFrom graphics points
#' @importFrom stats anova lm
PolarCoordHeterPlot <- function(
    x,
    gt = c("P1", "P1xP2", "P2"),
    rev_log = NULL,
    exp_fac = 1,
    thr = 1,
    plot_lab = c("none", "text", "graph"),
    col = NULL)
{

  plot_lab <- match.arg(plot_lab)

  if (is.null(col)) {
    col <- rep(0, nrow(x))
  } else {
    col <- as.vector(col)
    if (length(col)!=nrow(x)) {
      message("Parameter 'col' should be a color vector of length nrow(x)")
      col <- rep(col, length.out=nrow(x))
    }
  }

  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  # Helper Functions
  PlotCircle <- function(r=1, n=100, x0=0, y0=0, phi1=0, phi2=360, ...) {
    tmp.x <- x0 + r * cos((2 * pi / 360) * seq(phi1, phi2, length.out=n))
    tmp.y <- y0 + r * sin((2 * pi / 360) * seq(phi1, phi2, length.out=n))
    graphics::points(tmp.x, tmp.y, pch=".", ...)
  }
  PlotRadLine <- function(r=1, phi=0, x0=0, y0=0, ...) {
    # 0 <= phi <= 180
    tmp.x <- x0 + r * cos((2 * pi / 360) * (c(0, 180) + phi))
    tmp.y <- y0 + r * sin((2 * pi / 360) * (c(0, 180) + phi))
    graphics::lines(tmp.x, tmp.y, pch=".", ...)
  }
  PlotRadLegend <- function(r=1, phi=0, x0=0, y0=0, legend="test    ", ...) {
    tmp.x <- x0 + r * cos((2 * pi / 360) * phi)
    tmp.y <- y0 + r * sin((2 * pi / 360) * phi)
    graphics::legend(tmp.x, tmp.y, legend=legend, ...)
  }
  CombineGTsForLabel <- function(gt = c("B73", "B73xMo17", "Mo17")) {
    tmp.lab <- paste(rep(gt, times=c(5,4,3)),
                     rep(c("=", "<", "<", "<"), 3),
                     gt[c(3,3,2,2,2,1,1,3,3,2,1,1)],
                     c("<", "=", ">")[c(1,1,2,1,1,1,2,1,1,1,2,1)],
                     gt[c(2,2,3,3,3,3,3,1,1,1,2,2)], sep=" ")
    return(tmp.lab)
  }
  PlotSideLegend <- function(tmp.text = "B73 < Mo17 = B73xMo17", tmp.ord = c(3,2,4), tmp.cex = 2, plot.f1.shortname = TRUE) {
    par("mar"=rep(3,4))
    plot(1:5, 1:5, type="n", axes=F, ann=F)
    par("usr"=c(1,5,1,5))
    abline(v=1:5); abline(h=1:5);
    axis(1, at=2, labels=gt[1], tick=FALSE, cex.axis=tmp.cex, hadj=1)
    axis(1, at=3, labels=ifelse(plot.f1.shortname, "F1", gt[2]), tick=FALSE, cex.axis=tmp.cex, hadj=0.5)
    axis(1, at=4, labels=gt[3], tick=FALSE, cex.axis=tmp.cex, hadj=0)
    axis(3, at=3, labels=tmp.text, tick=FALSE, cex.axis=tmp.cex)
    points(x=2:4, y=tmp.ord, pch=16, cex=2)
    lines(x=2:4, y=tmp.ord, lwd=2)
  }
  main_fig <- function() {
    # compute trait mean values per genotype
    p1 <- apply(x[,grep(gt[1], colnames(x))], 1, mean, na.rm=T)
    f1 <- apply(x[,grep(gt[2], colnames(x))], 1, mean, na.rm=T)
    p2 <- apply(x[,grep(gt[3], colnames(x))], 1, mean, na.rm=T)
    dat <- data.frame(p1, f1, p2)

    # [option] revert log conversion of data
    if (!is.null(rev_log) && is.numeric(rev_log) && is.finite(rev_log)) {
      dat <- rev_log^dat
    }

    # [option] compute ANOVA P-values per row to filter rows
    if (thr<1) {
      filt <- apply(x[, colnames(x) %in% gt], 1, function(y) {
        anova(lm(unlist(y) ~ factor(colnames(x)[colnames(x) %in% gt])))$Pr[1]
      })
      if (sum(filt<=thr)==0) {
        stop("No trait below this significance level. Can't plot anything")
      }
      # only show metabolites with ANOVA P below 'thr'
      dat <- dat[filt<=thr,,drop=FALSE]
      col <- col[filt<=thr]
    }

    # compute radius for polar plot (see top)
    rad <- apply(dat, 1, function(x) {max(x)/min(x)})

    # compute angles as given in Swanson-Wagner Example, check for uniqueness first
    if ( any(apply(dat, 1, function(x) {duplicated(x)})) ) {
      warning("Identical trait values for related genotypes found.")
    }
    posseb <- matrix(c(1,1,2,2,3,3,3,2,1,3,2,1,2,3,3,1,1,2), nrow=6)
    ang <- apply(dat, 1, function(x) {
      fac <- (max(x)-median(x))/(max(x)-min(x))*60
      tmp <- which(apply(posseb, 1, function(y) {all(order(x)==y)}))
      abs(c(0, -120, -180, 180, 240, -360)[tmp] + fac)
    })

    # main plot
    # transform angles (Swanson-Wagner: clockwise from 12; R: counter-clockwise from 3)
    tmp.ang <- ang; tmp.ang[tmp.ang <= 90] <- 360 + tmp.ang[tmp.ang <= 90]; ang_norm <- 360 - abs(-90 + tmp.ang)
    x <- rad * cos((2*pi/360) * ang_norm)
    y <- rad * sin((2*pi/360) * ang_norm)

    tmp.lim <- 1/exp_fac * c(-1,1) * max(abs((c(x,y))))

    par(mar=c(0,0,0,0))
    if (plot_lab != "graph") par(pty="s")
    plot(x, y, ylim=tmp.lim, xlim=tmp.lim, axes=F, xlab="", ylab="", main="", type="n", asp=ifelse(plot_lab != "graph", 1, 0))

    tmp.r <- ceiling(max(abs(rad)))

    #plot dotted circle
    for (i in 1:tmp.r) PlotCircle(r=i, n=i*200, col=grey(0.8))

    # plot separating lines
    for (i in c(0, 60, 120)) { PlotRadLine(r=tmp.r, phi=i, x0=0, y0=0, lty=2) }
    for (i in c(30, 90, 150)){ PlotRadLine(r=tmp.r, phi=i, x0=0, y0=0) }

    # plot legends
    if (plot_lab == "text") {
      tmp.lab <- CombineGTsForLabel(gt)
      tmp.xjust <- c(0.5, 0, 0, 0, 0, 0, 0.5, 1, 1, 1, 1, 1)
      tmp.yjust <- c(0, 0, 0, 0.5, 1, 1, 1, 1, 1, 0.5, 0, 0)
      tmp.seq <- c(seq(90, 0, -30), seq(330, 120, -30))
      for (i in 1:12) {
        PlotRadLegend(r=tmp.r + .1, phi=tmp.seq[i], legend=tmp.lab[i], bg="white", bty="o", xjust=tmp.xjust[i], yjust=tmp.yjust[i], cex=0.6)
      }
    }

    # plot data points
    points(x, y, pch=21, bg=col, cex=ifelse(plot_lab != "graph", 1, 5/3))

    invisible(data.frame(x, y))
  }

  # use fancy outer legend
  if (plot_lab == "graph") {
    layout(mat=matrix(c(1,16:14,2,17,5,5,5,13,6,5,5,5,12,7,5,5,5,11,3,8,9,10,4), nrow=5))
    for (j in 1:4) { plot(1, type="n", axes=F, xlab="", ylab="") }
    out <- main_fig()
    for (j in 1:12) {
      PlotSideLegend(
        tmp.text = CombineGTsForLabel(gt = c(gt[1], "F1", gt[3]))[j],
        tmp.ord=matrix(c(2,4,2,2,4,3,2,4,4,2,3,4,2,2,4,3,2,4,4,2,4,4,2,3,4,2,2,4,3,2,4,4,2,3,4,2), nrow=3)[,j],
        tmp.cex=1, plot.f1.shortname=TRUE
      )
    }
  } else {
    out <- main_fig()
  }

  invisible(out)

}
