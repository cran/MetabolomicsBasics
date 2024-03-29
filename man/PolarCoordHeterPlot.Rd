% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PolarPlot.R
\name{PolarCoordHeterPlot}
\alias{PolarCoordHeterPlot}
\title{PolarCoordHeterPlot.}
\usage{
PolarCoordHeterPlot(
  x,
  gt = c("P1", "P1xP2", "P2"),
  rev_log = NULL,
  exp_fac = 1,
  thr = 1,
  plot_lab = c("none", "text", "graph"),
  col = NULL
)
}
\arguments{
\item{x}{Data matrix with measurement values (traits in rows and genotypes in columns).}

\item{gt}{Character vector of length=3 indicating P1, F and P2. These are used to  filter by column name from x.}

\item{rev_log}{If you've log transformed your data, you might want to revert the log transformation.}

\item{exp_fac}{Expansion factor to increase figure size.}

\item{thr}{Alpha level used in ANOVA to filter insignificant rows. Keep thr=1 to include all matrix rows.}

\item{plot_lab}{Show 'text' or 'graph' style labels of the polar sections (or keep 'none' to omit).}

\item{col}{Provide a color vector of length nrow(x).}
}
\value{
Will generate a plot in polar coordinates and return the x/y coordinates of the data points invisibly.
}
\description{
\code{PolarCoordHeterPlot} will draw a plot in polar coordinates
   visualizing heterosis effects according to a layout by Swanson-Wagner,
   where plot radius represents log2 of fold change between lowest and highest
   genotype and plot angle represents the ratio between lowest, intermediate
   and highest genotype.
}
\details{
See examples.
}
\examples{
# using the provided experimental data
raw <- MetabolomicsBasics::raw
sam <- MetabolomicsBasics::sam
x <- t(raw)
colnames(x) <- sam$GT
gt <- c("B73","B73xMo17","Mo17")
PolarCoordHeterPlot(x=x, gt=gt, plot_lab="graph", thr=0.01, rev_log=exp(1))

coord <- PolarCoordHeterPlot(x=x, gt=gt, thr=0.01, rev_log=exp(1))
points(x=coord$x[3], coord$y[3], pch=22, cex=4, col=2)
# using random data
gt <- c("P1","P1xP2","P2")
set.seed(0)
x <- matrix(rnorm(150), nrow = 10, dimnames = list(paste0("M",1:10), sample(rep(gt, 5))))
x[1:4,1:6]
PolarCoordHeterPlot(x=x, gt=gt)
# using text style labels for the sections
PolarCoordHeterPlot(x=x, gt=gt, plot_lab="text", exp_fac=0.75)
# reverting the order of parental genotypes
PolarCoordHeterPlot(x=x, gt=c("P2","P1xP2","P1"), plot_lab="text", exp_fac=0.75)
# using graph style labels for the sections
PolarCoordHeterPlot(x=x, gt=c("P2","P1xP2","P1"), plot_lab="graph")
# coloring data points
PolarCoordHeterPlot(x=x, gt=gt, col=1:10)
# applying ANOVA P value threshold to input rows
PolarCoordHeterPlot(x=x, gt=gt, col=1:10, thr=0.5)
PolarCoordHeterPlot(x=x, gt=gt, plot_lab="graph", col=1:10, thr=0.5)
}
