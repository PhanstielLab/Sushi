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

