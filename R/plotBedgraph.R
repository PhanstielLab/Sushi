#' plots data stored in bed file format
#' 
#'
#' @param signal signal track data to be plotted (in bedgraph format)
#' @param chrom chromosome of region to be plotted
#' @param chromstart start position
#' @param chromend end position
#' @param range y-range to plpt ( c(min,max) )
#' @param color color of signal track  
#' @param lwd color of line outlining signal track.  (only valid if linecol is not NA)
#' @param linecolor color of line outlining signal track.  use NA for no outline
#' @param addscale TRUE/FALSE whether to add a y-axis
#' @param overlay TRUE / FALSE whether this data should be plotted on top of an existing plot
#' @param rescaleoverlay TRUE/FALSE whether the new plot shold be rescaled based on the maximum value to match the existing plot (only valid when overlay is set to 'TRUE')   
#' @param transparency Value between 0 and 1 indication the degree of transparency of the plot
#' @param flip TRUE/FALSE whether the plot should be flipped over the x-axis
#' @param xaxt A character which specifies the x axis type.  See \code{\link{par}}
#' @param yaxt A character which specifies the y axis type.  See \code{\link{par}}
#' @param xlab Label for the x-axis
#' @param ylab Label for the y-axis
#' @param xaxs Must be set to 'i' for appropriate integration into Sushi plots.  See \code{\link{par}}
#' @param yaxs Must be set to 'i' for appropriate integration into Sushi plots.  See \code{\link{par}}
#' plottype
#' @param bty A character string which determined the type of box which is drawn about plots.  See \code{\link{par}}
#' @param ymax fraction of max y value to set as height of plot.
#' @param colorbycol palette to use to shade the signal track plot.  Only applicable when overlay is set to FALSE.
#' @param ... values to be passed to \code{\link{plot}}
#' @export
#' @examples
#'
#' data(Sushi_ChIPSeq_CTCF.bedgraph)
#' data(Sushi_DNaseI.bedgraph)
#'
#' chrom            = "chr11"
#' chromstart       = 1955000
#' chromend         = 1965000
#' 
#' plotBedgraph(Sushi_ChIPSeq_CTCF.bedgraph,chrom,chromstart,chromend,transparency=.50,flip=FALSE,color="blue",linecol="blue")
#' plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart,chromend,transparency=.50,flip=FALSE,color="#E5001B",linecol="#E5001B",overlay=TRUE,rescaleoverlay=TRUE)
#' labelgenome(chrom,chromstart,chromend,side=1,scipen=20,n=3,line=.18,chromline=.5,scaleline=0.5,scale="Mb")
#' 
#' transparency = 0.5
#' col1 = col2rgb("blue")
#' finalcolor1 = rgb(col1[1],col1[2],col1[3],alpha=transparency * 255,maxColorValue = 255)
#' col2 = col2rgb("#E5001B")
#' finalcolor2 = rgb(col2[1],col2[2],col2[3],alpha=transparency * 255,maxColorValue = 255)
#' 
#' legend("topright",inset=0.025,legend=c("DnaseI","ChIP-seq (CTCF)"),fill=c(finalcolor1,finalcolor2),border=c("blue","#E5001B"),text.font=2,cex=0.75)
plotBedgraph <-
function(signal,chrom,chromstart,chromend,range=NULL,color=SushiColors(2)(2)[1],
                         lwd=1,linecolor=NA,addscale=FALSE,overlay=FALSE,rescaleoverlay=FALSE,transparency=1.0,
                         flip=FALSE, xaxt='none',yaxt='none',xlab="",ylab="",xaxs="i",yaxs="i",bty='n',ymax=1.04,
                         colorbycol=NULL,...)
{
  if (overlay == TRUE)
  {
    colorbycol = NULL
  }
  
  if(is.na(linecolor ) == TRUE)
  {
    linecolor = color
  }
  
  # ensure that the chromosome is a character
  signal[,1] = as.character(signal[,1])
  
  # filter for desired region
  signaltrack = signal[which(signal[,1] == chrom & ((signal[,2] > chromstart & signal[,2] < chromend) |  (signal[,3] > chromstart & signal[,3] < chromend))),(2:4)]
  
  # exit if overlay is TRUE and there isn't enough data
  if (overlay ==TRUE && nrow(signaltrack) < 2)
  {
    return ("not enough data within range to plot")
  }
  
  # exit if overlay is FALSE and there isn't enough data
  if (nrow(signaltrack) < 2)
  {
    if (is.null(range) == TRUE)
    {
      range = c(0,1)
    }
    
    # make blank plot
    plot(0,0,xlim=c(chromstart,chromend),type='n',ylim=range,xaxt=xaxt,yaxt=yaxt,ylab=ylab,xaxs=xaxs,yaxs=yaxs,bty=bty,xlab=xlab,...) 
    return ("not enough data within range to plot")
  }
  
  # downsample for plotting
  while (nrow(signaltrack) > 8000)
  {
    # downsample for plotting if neccesary
    if (nrow(signaltrack) %% 2 != 0)
    {
      signaltrack = signaltrack[1:(nrow(signaltrack)-1),]
    }
    
    starts = signaltrack[seq(1, nrow(signaltrack), 2),1]
    stops  = signaltrack[seq(2, nrow(signaltrack), 2),2]
    meanval = apply(cbind(signaltrack[seq(1, nrow(signaltrack), 2),3],signaltrack[seq(2, nrow(signaltrack), 2),3]), 1, mean)
    signaltrack = cbind(starts,stops,meanval)
  }

  # add col names
  names(signaltrack)[(1:3)]    = c("V1","V2","V3")
  
  # make linking regions if neccesary
  linkingregions = cbind(signaltrack[1:(nrow(signaltrack)-1),2], signaltrack[2:nrow(signaltrack),1])
  linkingregions = matrix(linkingregions[which(linkingregions[,1] != linkingregions[,2]),],ncol=2)
  
  if (nrow(linkingregions) > 0)
  {
    linkingregions = cbind(linkingregions,0)
    
    # make col names the same
    names(linkingregions)[(1:3)] = c("V1","V2","V3")
    
    # add linking regions to signaltrack
    signaltrack = rbind(signaltrack,linkingregions)
  }
  
  # sort data
  signaltrack = signaltrack[order(signaltrack[,1]),]
  
  # convert two columns to one
  signaltrack = cbind(as.vector(t(signaltrack[,c(1,2)])),as.vector(t(signaltrack[,c(3,3)])))
  
  # add slighltly negative vaue to both ends to ensure proper polygon plotting
  signaltrack = rbind(c(min(signaltrack[,1]),-.00001),signaltrack)
  signaltrack = rbind(signaltrack, c(max(signaltrack[,1]),-.00001))
  
  if (flip == TRUE)
  {
    signaltrack[,2] = signaltrack[,2]*-1
  }
  
  # determine the y-limits
  if (is.null(range) == TRUE)
  {
    range = c(0,ymax*max(signaltrack[,2]))
    if (flip == TRUE)
    {
      range = c(ymax*min(signaltrack[,2]),0)
    }
  }

  if (overlay == FALSE)
  {
    # make blank plot
    plot(signaltrack,xlim=c(chromstart,chromend),type='n',ylim=range,xaxt=xaxt,yaxt=yaxt,ylab=ylab,xaxs=xaxs,yaxs=yaxs,bty=bty,xlab=xlab,...) 
  }
    
  # rescale the overlay plot for comparative purposes
  if (rescaleoverlay == TRUE)
  {
    if (flip == FALSE)
    {
      signaltrack[,2] =   par('usr')[4] * signaltrack[,2] / max(abs(signaltrack[,2]) )
    }
    if (flip == TRUE)
    {
      
      signaltrack[,2] =   abs(par('usr')[3]) * signaltrack[,2] / max(abs(signaltrack[,2]) )
    }
  }

  # set the transparency
  rgbcol = col2rgb(color)
  finalcolor = rgb(rgbcol[1],rgbcol[2],rgbcol[3],alpha=transparency * 255,maxColorValue = 255)
  
  if (is.null(colorbycol) == FALSE)
  {
    # add the gradient to the background
    if (is.null(colorbycol) == FALSE)
    {
      plotlef = par('usr')[1]
      plotrig = par('usr')[2]
      plotbot = par('usr')[3]
      plottop = par('usr')[4]
      
      bgcol = maptocolors((1:100),colorbycol)
      
      if (flip == FALSE)
      {
        tops = seq(plotbot,plottop,length.out=101)[2:101]  
        bots =  seq(plotbot,plottop,length.out=101)[1:100]
      }
      if (flip == TRUE)
      {
        bots =  rev(seq(plotbot,plottop,length.out=101)[2:101])
        tops =  rev(seq(plotbot,plottop,length.out=101)[1:100] )
      }
      
      for (i in (3:(nrow(signaltrack)-1)))
      {
        if (flip == FALSE)
        {
          ybots = bots[which(bots <= signaltrack[i,2])]
          ytops = tops[which(tops <  signaltrack[i,2])]
        }
        if (flip == TRUE)
        {
          ybots = bots[which(bots >= signaltrack[i,2])]
          ytops = tops[which(tops >  signaltrack[i,2])]
        }
        ytops = c(ytops,signaltrack[i,2] )
        xleft = rep(signaltrack[i-1,1],length(ytops))
        xrigh = rep(signaltrack[i,1],length(ytops))
        xybgcol = bgcol[1:length(ytops)]
        
        rect(xleft=xleft, ybottom=ybots, xright=xrigh, ytop=ytops,col= xybgcol,border=bgcol,lwd=lwd) 
      } 
    }
  }

  if (is.null(colorbycol) == TRUE)
  {
    # plot the signal track
    polygon(signaltrack,col=finalcolor,border=linecolor,lwd=lwd)
  }
  
  # add scale to upper right corner
  if (addscale == TRUE)
  {
    mtext(paste(range[1],range[2],sep="-"),side=3,font=1,adj=1.00,line=-1)
  }
}
