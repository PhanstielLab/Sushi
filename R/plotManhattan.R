#' plots a Manhattan plot
#'
#' @param bedfile bedfile for Manhattan plot
#' @param chrom chromosome of region to be plotted
#' @param chromstart start position
#' @param chromend end position
#' @param pvalues pvalues to be used for plotting (will be converted to -log(10) space)
#' @param genome A genome object (2 columns: column 1 = chromosome name, column 2 = length of chromosome). Required if plotting multiple chromosomes at once.
#' @param col single colors, vector of colors, or color palette for coloring points
#' @param space the space in between each chromosome as a fraction of the width of the plot
#' @param ymax fraction of max y value to set as height of plot.
#' @param ... Arguments to be passed to methods such as \code{\link{plot}}
#' @export 
#' @examples
#' 
#' data(Sushi_GWAS.bed)
#' data(Sushi_hg18_genome)
#' 
#' chrom1            = "chr11"
#' chromstart1       = 500000
#' chromend1         = 5050000
#' 
#' plotManhattan(bedfile=Sushi_GWAS.bed,pvalues=Sushi_GWAS.bed[,5],genome=Sushi_hg18_genome,col=topo.colors,cex=0.75)
#' labelgenome(genome=Sushi_hg18_genome,side=1,scipen=20,n=4,scale="Mb",edgeblankfraction=0.20,line=.18,chromline=.5,scaleline=0.5)
#' axis(side=2,las=2,tcl=.2)
#' mtext("log10(P)",side=2,line=1.75,cex=.75,font=2)
#' 
#' 
#' 
plotManhattan <-
function(bedfile, chrom=NULL,chromstart=NULL,chromend=NULL,pvalues,genome=NULL,col=SushiColors(5),space=0.01,ymax=1.04,...)
{
  if (is.null(genome) == FALSE)
  {
    chromoffsets = chromOffsets(genome,space)
    
    if (class(col) == "function")
    {
      col = col(nrow(chromoffsets))
    }
    else
    {
      col = rep(col,ceiling(nrow(chromoffsets) / length(col)   ))
    }
    
    
    
    # remove data from chroms not in genome
    bedfile = bedfile[bedfile[,1]%in%chromoffsets[,1],]
    
    # add offsets for chromosome and colors
    columber = 0
    columbers = rep(0,nrow(bedfile))
    for (i in (1:nrow(chromoffsets)))
    {
      columber = columber + 1
      if (columber > length(col))
      {
        columber = 1
      }

      rowsofinterest = which(bedfile[,1] == chromoffsets[i,1])
      bedfile[rowsofinterest,2] =  bedfile[rowsofinterest,2] + chromoffsets[i,3]
      columbers[rowsofinterest] = columber
    }
    
    cumsums = cumsum(as.numeric(genome[,2]))
    spacer  = cumsums[length(cumsum(as.numeric(genome[,2])))] * space
    
    # make the plot
    yrange = c(min(-log10(bedfile[,5])),max(-log10(bedfile[,5])) * ymax)
    plot(bedfile[,2],-log10(bedfile[,5]),col=col[columbers],xlim=c(0,max(chromoffsets[,4])+spacer),ylim=yrange,pch=19,xaxt='n',xlab='',yaxt='n',ylab='',xaxs = 'i',yaxs='i',...) 
  }
  
  if (is.null(genome) == TRUE)
  {
    bedfile = bedfile[which( bedfile[,1] == chrom),]
    
    # make the plot
    plot(bedfile[,2],-log10(bedfile[,5]),col=col,xlim=c(chromstart,chromend),pch=19,xaxt='n',xlab='',yaxt='n',ylab='',xaxs = 'i',yaxs='i',...) 
  }
  
}
