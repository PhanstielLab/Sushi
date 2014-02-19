#' Adds a zoom box to a plot
#' 
#' This function is used on the second plot of a zoom in
#'
#' @param zoomregion  Region of another zoom on this plot.  Only required if this plot has another zoomregion on it.
#' @param lty line type for box.  See \code{\link{par}}
#' @param lwd line width.  See See \code{\link{par}}
#' @param col  Color for zoombox line
#' @param topextend How far to exted the lines above the current plot (as a fraction of the plot height)
#' @param passthrough TRUE / FALSE whether or not to pass the zoom though this plot.  If set to FALSE no horizontal line is drawn on the botoom of the plot
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
zoombox <- function(zoomregion = NULL,lty=2,lwd=1,col="black",topextend=2,passthrough=FALSE)
{

    #determine the y range
    yrange = abs(par('usr')[4]-par('usr')[3])
    topextend = topextend * yrange
    
    originalxpd = par()$xpd
    par(xpd = TRUE)
    
    if (passthrough == TRUE)
    {
      segments(x0=par('usr')[1], y0=par('usr')[3]-topextend, x1 = par('usr')[1], y1 = topextend,lty=lty,lwd=lwd,col=col)
      segments(x0=par('usr')[2], y0=par('usr')[3]-topextend, x1 = par('usr')[2], y1 = topextend,lty=lty,lwd=lwd,col=col)
    }
    
    if (passthrough == FALSE)
    {
      segments(x0=par('usr')[1], y0=par('usr')[3], x1 = par('usr')[1], y1 = topextend,lty=lty,lwd=lwd,col=col)
      segments(x0=par('usr')[2], y0=par('usr')[3], x1 = par('usr')[2], y1 = topextend,lty=lty,lwd=lwd,col=col)
      if (is.null(zoomregion) == TRUE)
      {
        segments(x0=par('usr')[1], y0=par('usr')[3], x1 = par('usr')[2], y1 = par('usr')[3],lty=lty,lwd=lwd,col=col)
      }
      if (is.null(zoomregion) == FALSE)
      {
        segments(x0=par('usr')[1], y0=par('usr')[3], x1 = min(zoomregion), y1 = par('usr')[3],lty=lty,lwd=lwd,col=col)
        segments(x0=max(zoomregion), y0=par('usr')[3], x1 = par('usr')[2], y1 = par('usr')[3],lty=lty,lwd=lwd,col=col)
      }
    }
    par(xpd = originalxpd)
}
