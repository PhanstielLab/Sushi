\name{addlegend}
\alias{addlegend}
\title{adds a legend to a Sushi plot}
\usage{
addlegend(range, title = "", labels.digits = 1, palette = topo.colors,
  side = "right", labelside = "left", xoffset = 0.1, width = 0.05,
  bottominset = 0.025, topinset = 0.025, tick.num = 5,
  tick.length = 0.01, txt.font = 1, txt.cex = 0.75, title.offset = 0.05,
  title.font = 2, title.cex = 1)
}
\arguments{
  \item{range}{the rang of values to be plotted.  ie
  c(min,max)}

  \item{title}{title of values to be mapped}

  \item{labels.digits}{Number of digits after the decimal
  point to include in labels}

  \item{palette}{color palette to use}

  \item{side}{side of plot to place legend
  ('right','left')}

  \item{labelside}{side of legend to place legend title}

  \item{xoffset}{fraction of plot to offset the legend}

  \item{width}{width as a fraction of the plot width}

  \item{bottominset}{inset from the bottom of the blot as a
  fraction of the plot width}

  \item{topinset}{inset from the top of the blot as a
  fraction of the plot width}

  \item{tick.num}{desired number of tickmarks}

  \item{tick.length}{length of tick marks}

  \item{txt.font}{font type of legend text}

  \item{txt.cex}{font size of legned text}

  \item{title.offset}{offset of title from the key}

  \item{title.font}{font type of legend title}

  \item{title.cex}{font size of legned text}
}
\description{
This function adds a legend to Sushi plots that have a
colorby function (e.g. plotHic, plotGenes, and plotBedpe)
}
\examples{
data(Sushi_HiC.matrix)

chrom            = "chr11"
chromstart       = 500000
chromend         = 5050000

phic = plotHic(Sushi_HiC.matrix,chrom,chromstart,chromend,max_y = 20,zrange=c(0,28),palette = topo.colors,flip=FALSE)

labelgenome(chrom,chromstart,chromend,side=1,scipen=20,n=4,scale="Mb",edgeblankfraction=0.20,line=.18,chromline=.5,scaleline=0.5)

addlegend(phic[[1]],palette=phic[[2]],title="score",side="right",bottominset=0.4,topinset=0,xoffset=-.035,labelside="left",width=0.025,title.offset=0.035)
}

