\name{plotBedgraph}
\alias{plotBedgraph}
\title{plots data stored in bed file format}
\usage{
plotBedgraph(signal, chrom, chromstart, chromend, range = NULL,
  color = SushiColors(2)(2)[1], lwd = 1, linecolor = NA,
  addscale = FALSE, overlay = FALSE, rescaleoverlay = FALSE,
  transparency = 1, flip = FALSE, xaxt = "none", yaxt = "none",
  xlab = "", ylab = "", xaxs = "i", yaxs = "i", bty = "n",
  ymax = 1.04, colorbycol = NULL, ...)
}
\arguments{
  \item{signal}{signal track data to be plotted (in
  bedgraph format)}

  \item{chrom}{chromosome of region to be plotted}

  \item{chromstart}{start position}

  \item{chromend}{end position}

  \item{range}{y-range to plpt ( c(min,max) )}

  \item{color}{color of signal track}

  \item{lwd}{color of line outlining signal track.  (only
  valid if linecol is not NA)}

  \item{linecolor}{color of line outlining signal track.
  use NA for no outline}

  \item{addscale}{TRUE/FALSE whether to add a y-axis}

  \item{overlay}{TRUE / FALSE whether this data should be
  plotted on top of an existing plot}

  \item{rescaleoverlay}{TRUE/FALSE whether the new plot
  shold be rescaled based on the maximum value to match the
  existing plot (only valid when overlay is set to 'TRUE')}

  \item{transparency}{Value between 0 and 1 indication the
  degree of transparency of the plot}

  \item{flip}{TRUE/FALSE whether the plot should be flipped
  over the x-axis}

  \item{xaxt}{A character which specifies the x axis type.
  See \code{\link{par}}}

  \item{yaxt}{A character which specifies the y axis type.
  See \code{\link{par}}}

  \item{xlab}{Label for the x-axis}

  \item{ylab}{Label for the y-axis}

  \item{xaxs}{Must be set to 'i' for appropriate
  integration into Sushi plots.  See \code{\link{par}}}

  \item{yaxs}{Must be set to 'i' for appropriate
  integration into Sushi plots.  See \code{\link{par}}
  plottype}

  \item{bty}{A character string which determined the type
  of box which is drawn about plots.  See
  \code{\link{par}}}

  \item{ymax}{fraction of max y value to set as height of
  plot.}

  \item{colorbycol}{palette to use to shade the signal
  track plot.  Only applicable when overlay is set to
  FALSE.}

  \item{...}{values to be passed to \code{\link{plot}}}
}
\description{
plots data stored in bed file format
}
\examples{
data(Sushi_ChIPSeq_CTCF.bedgraph)
data(Sushi_DNaseI.bedgraph)

chrom            = "chr11"
chromstart       = 1955000
chromend         = 1965000

plotBedgraph(Sushi_ChIPSeq_CTCF.bedgraph,chrom,chromstart,chromend,transparency=.50,flip=FALSE,color="blue",linecol="blue")
plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart,chromend,transparency=.50,flip=FALSE,color="#E5001B",linecol="#E5001B",overlay=TRUE,rescaleoverlay=TRUE)
labelgenome(chrom,chromstart,chromend,side=1,scipen=20,n=3,line=.18,chromline=.5,scaleline=0.5,scale="Mb")

transparency = 0.5
col1 = col2rgb("blue")
finalcolor1 = rgb(col1[1],col1[2],col1[3],alpha=transparency * 255,maxColorValue = 255)
col2 = col2rgb("#E5001B")
finalcolor2 = rgb(col2[1],col2[2],col2[3],alpha=transparency * 255,maxColorValue = 255)

legend("topright",inset=0.025,legend=c("DnaseI","ChIP-seq (CTCF)"),fill=c(finalcolor1,finalcolor2),border=c("blue","#E5001B"),text.font=2,cex=0.75)
}

