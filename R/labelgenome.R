#' Adds genome coordinates to the x-axis of a Sushi plot
#'
#' @param chrom chromosome to plot
#' @param chromstart start position
#' @param chromend end position
#' @param genome a genome object (2 columns: column 1 = chromosome name, column 2 = length of chromosome). Only for multi chromosomal plots
#' @param space  the space in between each chromosome as a fraction of the width of the plot.  Only for multi chromosomal plots
#' @param scale Scale of the plot ('bp','Kb','Mb') 
#' @param side Side of the scale to add the plot to.  Only tested for sides 1 and 3.
#' @param scipen higher values decrease the likelihood of using scientific for the position labels.
#' @param n Desired number of ticks 
#' @param chromfont font type of chromosome label
#' @param chromadjust position, as a fraction of the width of the plot, of the chomosome label
#' @param chromcex font size of the chomosome label
#' @param chromline vertical offset of the chomosome label
#' @param scalefont font type of scale label
#' @param scaleadjust  position, as a fraction of the width of the plot, of the scale label
#' @param scalecex font size of the scale label
#' @param scaleline vertical offset of the scale label
#' @param line vertical offset of position labels
#' @param edgeblankfraction percent of the edges to leave black for chromosome and scale labels
#' @param ... values to be passed to \code{\link{axis}}
#' @export 
#' @examples
#' data(Sushi_DNaseI.bedgraph)
#' # set the genomic regions
#' 
#' plotBedgraph(Sushi_DNaseI.bedgraph,chrom="chr11",chromstart=1650000,chromend=2350000,colorbycol=SushiColors(7))
#' labelgenome(chrom="chr11",chromstart=1650000,chromend=2350000,side=1,n=4,scale="Mb")
#' axis(side=2,las=2,tcl=.2)
#' mtext("Read Depth",side=2,line=1.75,cex=.75,font=2)
#' 
labelgenome <-
function(chrom, chromstart,chromend,genome=NULL,space=0.01,scale="bp",side=1,scipen=20,n=5,
                        chromfont=2,chromadjust=0.015,chromcex=1,chromline=0.5,
                        scalefont=2,scaleadjust=0.985,scalecex=1,scaleline=0.5,line=0.18,edgeblankfraction=0.10,...)
{
  
  minxpos = par('usr')[1]
  maxxpos = par('usr')[2]
  range   = abs(maxxpos - minxpos)
  curxpd =  par()$xpd
  par(xpd = FALSE)
  curscipen =  getOption("scipen")
  options("scipen" = 20)
  
  if(is.null(genome) == FALSE)
  {
    chromoffsets = chromOffsets(genome,space)
    
    # label the chromsomes
    chromcenters = (chromoffsets[,3] + chromoffsets[,4]) / 2
    labels = gsub("chr", "", chromoffsets[,1]) 
    axis(side=side,at=chromcenters,labels=labels,...)  
}

  if(is.null(genome) ==TRUE)
  {
    chromstart = chromstart + edgeblankfraction * range
    chromend   = chromend   - edgeblankfraction * range
    
    if(scale == "bp")
    {
      chromstartlabel = chromstart
      chromendlabel   = chromend
      fact = 1
    }
    
    if(scale == "Mb")
    {
      chromstartlabel = chromstart /1000000
      chromendlabel   = chromend   /1000000
      fact = 1000000
    }
    if(scale == "Kb")
    {
      chromstartlabel = chromstart /1000
      chromendlabel   = chromend   /1000
      fact = 1000
    }
    
    # determine the label names and positions
    labels = pretty(c(chromstartlabel,chromendlabel),n=n)
    labels = labels[2:(length(labels)-1)]
    
    #plot the axis
    axis(side=side,at = c(par('usr')[1]-range,labels*fact,par('usr')[2]+range), labels=c("",labels,""),line=line,...)
    mtext(chrom,side=side,font=chromfont,adj=chromadjust,line=chromline + line,cex=chromcex)
    mtext(scale,side=side,font=chromfont,adj=scaleadjust,line=scaleline + line,cex=scalecex)
  }
  
  # reset the scipen value ot its original settings
  options("scipen" = curscipen)
  par(xpd = curxpd)
}
