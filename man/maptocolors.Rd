\name{maptocolors}
\alias{maptocolors}
\title{maps numeric vector to color palette}
\usage{
maptocolors(vec, col, num = 100, range = NULL)
}
\arguments{
  \item{vec}{numeric vector to map to color}

  \item{col}{color palette to which to be mapped}

  \item{num}{number of bins of colors}

  \item{range}{range of values to map}
}
\description{
maps numeric vector to color palette
}
\examples{
plot((1:10),col=maptocolors(vec=(1:10),colorRampPalette(c("blue","red"))),pch=19,cex=4)
}

