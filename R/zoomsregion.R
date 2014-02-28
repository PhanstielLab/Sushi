#' Adds a zoom region to a plot
#' 
#' This function is used on the first plot of a zoom in
#'
#' @param region chromosome start and stop to zoom in on
#' @param chrom chromosome of region to be plotted
#' @param genome A genome object (2 columns: column 1 = chromosome name, column 2 = length of chromosome). Set to NULL if adding zoom to a plot with only a singe chromosome.
#' @param space the space in between each chromosome as a fraction of the width of the plot.  Only used when adding a zoomsregion to a plot with multiple chromosomes (e.g. a Manhattan plot)
#' @param padding The minimum size of a zoom region  (as a fraction of the plot width).  If the specified zoom region is too small it will zoom on a region twice this wide cerntered on te specified zoom region.
#' @param col Color of the zoom region
#' @param zoomborder Color of the border of the zoom region
#' @param lty line type of zoom region border.  See \code{\link{plot}}
#' @param lwd line type of zoom region border.  See \code{\link{plot}}
#' @param extend single value or vector of 2 values specifying how far the zoom region extend above and below the plot region (as a fraction of the plot height).  Note this valu only applies to the narrow portion of the zoom region.
#' @param wideextend  Value specifying how below the plot region (as a fraction of the plot height) the wide portion of the zoom window starts. Only applicable if highlight is set to FALSE.
#' @param offsets vector of 2 values specifying offsets to the left and right side of the wide portion of the zoom window.  It may be neccesary to adjust these by trial and error for more complicated layouts. Only applicable if highlight is set to FALSE.
#' @param highlight TRUE/FALSE indicating if you are adding a highlight region as opposed to a zoom in.  Highlight regions simply draw a box around the region of interest
#' @export
#' @examples
#'
#' data(Sushi_DNaseI.bedgraph)
#' data(Sushi_ChIPSeq_CTCF.bedgraph)
#' 
#' # make a layout for all of the plots
#' layout(matrix(c(1,1,
#'                 2,2)
#'               ,2, 2, byrow = TRUE))
#' par(mgp=c(3, .3, 0))
#' 
#' par(mar=c(3,4,2,1))
#' chrom            = "chr11"
#' chromstart       = 1650000
#' chromend         = 2350000
#' zoomregion1      = c(1955000,1965000)
#' 
#' plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart,chromend,transparency=1.0,color="#5900E5",lwd=1,linecol="#5900E5")
#'
#' zoomsregion(zoomregion1,col=NA,zoomborder="black",lty=2,lwd=1,extend=c(0.01,0.09),wideextend=0.10,offsets=c(0,0))
#' 
#' labelgenome(chrom,chromstart,chromend,side=1,scipen=20,n=4,line=.18,chromline=.5,scaleline=0.5,scale="Mb")
#' 
#' axis(side=2,las=2,tcl=.2)
#' mtext("Read Depth",side=2,line=1.75,cex=.75,font=2)
#' 
#' # plot dnaseI data
#' plotBedgraph(Sushi_DNaseI.bedgraph,chrom,zoomregion1[1],zoomregion1[2],transparency=.50,flip=FALSE,color="#E5001B",linecol="#E5001B")
#' 
#' # plot chip-seq data
#' plotBedgraph(Sushi_ChIPSeq_CTCF.bedgraph,chrom,zoomregion1[1],zoomregion1[2],transparency=.30,flip=FALSE,color="blue",linecol="blue",overlay=TRUE,rescaleoverlay=TRUE)
#' 
#' # add zoombox
#' zoombox(zoomregion = NULL,lwd = 1,col="black")
#' 
#' axis(side=2,las=2,tcl=.2)
#' mtext("Read Depth",side=2,line=1.75,cex=.75,font=2)
#' 
#' # add the genome labels
#' labelgenome(chrom,zoomregion1[1],zoomregion1[2],side=1,scipen=20,n=3,line=.18,chromline=.5,scaleline=0.5,scale="Mb")
#' 
#' # set the legend colors
#' transparency = 0.5
#' col1 = col2rgb("blue")
#' finalcolor1 = rgb(col1[1],col1[2],col1[3],alpha=transparency * 255,max = 255)
#' col2 = col2rgb("#E5001B")
#' finalcolor2 = rgb(col2[1],col2[2],col2[3],alpha=transparency * 255,max = 255)
#' 
#' # add legend
#' legend("topright",inset=0.025,legend=c("DnaseI","ChIP-seq (CTCF)"),fill=c(finalcolor1,finalcolor2),border=c("blue","#E5001B"),text.font=2,cex=0.75)
#' 
zoomsregion <-
function(region,chrom=NULL,genome=NULL,space=0.01,padding=0.005,col=NA,zoomborder="black",lty=2,lwd=1,extend=0,wideextend=0.1,offsets=c(0,0),highlight=FALSE)
{
  if (is.null(genome) == FALSE)
  {
    chromoffsets = chromOffsets(genome,space)
    region[1] = region[1] + chromoffsets[which(chromoffsets[,1]==chrom),3]
    region[2] = region[2] + chromoffsets[which(chromoffsets[,1]==chrom),3]

    if (abs(region[2]-region[1]) < 2*padding*(max(chromoffsets[,4])))
    {
      center = mean(region)
      region[1] = center - padding*max(chromoffsets[,4])
      region[2] = center + padding*max(chromoffsets[,4])
    }
  }
  
  orgimal_xpd = par()$xpd
  par(xpd=TRUE)
  
  # get plot limits
  plotlef = par('usr')[1]
  plotrig = par('usr')[2]
  plotbot = par('usr')[3]
  plottop = par('usr')[4]
  xrange  = abs(plotrig - plotlef)
  yrange  = abs(plotbot - plottop)
  
  # determine upper and lower extension
  if (length(extend) == 1)
  {
    extend.upper  = extend * yrange
    extend.lower  = extend * yrange
  }
  if (length(extend) > 1)
  {
    extend.upper  = extend[1] * yrange
    extend.lower  = extend[2] * yrange
  }
  
  # determine bottom xoffsets
  plotlef = plotlef + xrange * offsets[1]
  plotrig = plotrig - xrange * offsets[2]
  
  # Set the lower values
  miny    = plotbot - extend.lower - (wideextend * yrange)
  lowy    = plotbot - extend.lower - (5*yrange)
  
  if (is.null(region)==FALSE)
  {
    # find the original xpd value
    currentxpd = par()$xpd
    
    # set xpd to allwo for plotting outside of borders
    par(xpd=TRUE)
    
    # add padding
    if ((abs(region[1] - region[2] ) / xrange) < (2* padding))
    {
      center = mean(region)
      region[1] = center - padding*xrange
      region[2] = center + padding*xrange
    }
    
    # highlight
    if (highlight == TRUE)
    {
      plotrig = region[2]
      plotlef = region[1]
      miny = plotbot-extend.lower
      lowy = plotbot-extend.lower
    }
    
    # draw a polygon
    polygon(x = c(region[1],region[1],region[2],region[2],plotrig,plotrig,plotlef,plotlef,region[1]),
            y = c(plotbot-extend.lower,plottop+extend.upper,plottop+extend.upper,plotbot-extend.lower,miny,lowy,lowy,miny,plotbot-extend.lower),
            col=col,
            border=zoomborder,
            lty = lty,
            lwd = lwd
    )
    
    # reset the original xpd value
    par(xpd=currentxpd)
  }

  par(xpd = orgimal_xpd)
}