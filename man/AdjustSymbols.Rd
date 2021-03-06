% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AdjustSymbols.R
\name{AdjustSymbols}
\alias{AdjustSymbols}
\title{AdjustSymbols.}
\usage{
AdjustSymbols(cols = NULL, pchs = NULL, colorset = NULL,
  symbolset = NULL)
}
\arguments{
\item{cols}{Factor (color output) or numeric (greyscale output) vector or NULL (omitted).}

\item{pchs}{Factor vector or NULL (omitted).}

\item{colorset}{Can be selectively specified here. If NULL set automatically, else can be explicitly provided.}

\item{symbolset}{Can be selectively specified here. If NULL upt o 5 nice symbols are selected automatically where background can be colored.}
}
\value{
data.frame with two columns (cols, pchs). Will be used by several plotting functions automatically.
}
\description{
\code{AdjustSymbols} will generate plotting character and color vectors based on experimental factors.
}
\details{
not yet
}
\examples{
# load data and plot using provided color scheme
utils::data(raw, package = "MetabolomicsBasics")
utils::data(sam, package = "MetabolomicsBasics")
head(sam)
plot(y=raw[,1], x=as.numeric(sam$GT), pch=sam$pchs, bg=sam$cols)

# change colors to greyscale
head(AdjustSymbols(cols=sam$GT, pchs=sam$Origin))
tmp.set <- grDevices::rainbow(length(levels(sam$GT)))
head(AdjustSymbols(cols=sam$GT, pchs=sam$Batch, colorset=tmp.set))
plot(raw[,1]~sam$GT, col=unique_labels(sam=sam, g="GT")[,"cols"])
sam$cols <- AdjustSymbols(cols=as.numeric(sam$GT))
plot(raw[,1]~sam$GT, col=unique_labels(sam=sam, g="GT")[,"cols"])

}
