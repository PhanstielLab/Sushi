\name{plotManhattan}
\alias{plotManhattan}
\title{plots a Manhattan plot}
\usage{
plotManhattan(bedfile, chrom = NULL, chromstart = NULL, chromend = NULL,
  pvalues, genome = NULL, col = SushiColors(5), space = 0.01,
  ymax = 1.04, ...)
}
\arguments{
  \item{bedfile}{bedfile for Manhattan plot}

  \item{chrom}{chromosome of region to be plotted}

  \item{chromstart}{start position}

  \item{chromend}{end position}

  \item{pvalues}{pvalues to be used for plotting (will be
  converted to -log(10) space)}

  \item{genome}{A genome object (2 columns: column 1 =
  chromosome name, column 2 = length of chromosome).
  Required if plotting multiple chromosomes at once.}

  \item{col}{single colors, vector of colors, or color
  palette for coloring points}

  \item{space}{the space in between each chromosome as a
  fraction of the width of the plot}

  \item{ymax}{fraction of max y value to set as height of
  plot.}

  \item{...}{Arguments to be passed to methods such as
  \code{\link{plot}}}
}
\description{
plots a Manhattan plot
}
\examples{
data(Sushi_GWAS.bed)
data(Sushi_hg18_genome)

chrom1            = "chr11"
chromstart1       = 500000
chromend1         = 5050000

plotManhattan(bedfile=Sushi_GWAS.bed,pvalues=Sushi_GWAS.bed[,5],genome=Sushi_hg18_genome,col=topo.colors,cex=0.75)
labelgenome(genome=Sushi_hg18_genome,side=1,scipen=20,n=4,scale="Mb",edgeblankfraction=0.20,line=.18,chromline=.5,scaleline=0.5)
axis(side=2,las=2,tcl=.2)
mtext("log10(P)",side=2,line=1.75,cex=.75,font=2)
}

