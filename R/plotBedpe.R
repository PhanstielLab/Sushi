#' plots data stored in bed file format
#' 
#'
#' @param bedpedata bed paired end data to be plotted
#' @param chrom chromosome of region to be plotted
#' @param chromstart start position
#' @param chromend end position
#' @param heights single value or vector specifying the height of the arches to be plotted (only valid when plottype is set to "loops" )
#' @param color single value or vector specifying colors of bedpe elements
#' @param colorby vector to scale colors by
#' @param colorbycol palette to apply color scale to (only valid when colorby is not NULL)
#' @param colorbyrange the range of values to apply the color scale to.  Values outside that range will be set to the limits of the range.
#' @param lwdby vector to scale line widths by
#' @param lwdrange the range of values to apply the line width scale to.  Values outside that range will be set to the limits of the range.
#' @param offset offset of bedpe elements from the x-axis
#' @param flip TRUE/FALSE whether the plot should be flipped over the x-axis
#' @param lwd linewidth for bedpe elements (only valid when colorby is not NULL)
#' @param xaxt A character which specifies the x axis type.  See \code{\link{par}}
#' @param yaxt A character which specifies the y axis type.  See \code{\link{par}}
#' @param plottype type of plot (acceptable values are 'loops' and 'lines')
#' @param maxrows The maximum number of rows to plot on the y-axis
#' @param ymax fraction of max y value to set as height of plot. Only applies when plottype is set to 'loops'
#' @param height the height of the boxes at either end of a bedpe element if plottype is set to 'lines'. Typical vaues range form 0 to 1.  (only valid when plottype is set to 'lines')
#' @param bty A character string which determined the type of box which is drawn about plots.  See \code{\link{par}}
#' @param ... values to be passed to \code{\link{plot}}
#' @examples
#' 
#' data(Sushi_5C.bedpe)
#' 
#' chrom            = "chr11"
#' chromstart       = 1650000
#' chromend         = 2350000
#' pbpe = plotBedpe(Sushi_5C.bedpe,chrom,chromstart,chromend,heights = Sushi_5C.bedpe$score,offset=0,flip=FALSE,bty='n',
#' lwd=1,plottype="loops",colorby=Sushi_5C.bedpe$samplenumber,colorbycol=topo.colors)
#' labelgenome(chrom, chromstart,chromend,side=1,scipen=20,n=3,scale="Mb",line=.18,chromline=.5,scaleline=0.5)
#' legend("topright",inset =0.01,legend=c("K562","HeLa","GM12878"),col=c(topo.colors(3)),pch=19,bty='n',text.font=2)
#' axis(side=2,las=2,tcl=.2)
#' mtext("Z-score",side=2,line=1.75,cex=.75,font=2)
#' 


plotBedpe <-
function(bedpedata,chrom,chromstart,chromend,heights,
                      color="black",colorby=NULL,colorbycol=NULL,colorbyrange=NULL,
                      lwdby    =  NULL,lwdrange =  c(1,5),
                      offset=0,flip=FALSE,lwd=1,xaxt='n',yaxt='n',bty='n',
                      plottype="loops",maxrows=10000,height=.3,ymax=1.04,...)
{
  
  if (nrow(bedpedata) == 0)
  {
    plot(c(1,1),type='n',xlab="",ylab="",xaxs = 'i',yaxs='i',xlim=c(chromstart,chromend), ylim=c(0,1),xaxt=xaxt,yaxt=yaxt,bty=bty,...)
    return(list(colorbyrange,colorbycol))
  }
  
  # convert strand info 
  if (ncol(bedpedata) >= 10)
  {
    bedpedata[,9] = convertstrandinfo(bedpedata[,9])
    bedpedata[,10] = convertstrandinfo(bedpedata[,10])
  }
  
  # Define a function that determines which row to plot a gene on
  checkrow <- function(data,alldata,maxrows,wiggle=0)
  {
    for (row in (1:maxrows))
    {
      thestart = min(as.numeric(data[c(2,3,5,6)])) - wiggle
      thestop  = max(as.numeric(data[c(2,3,5,6)])) + wiggle
      
      if (nrow(alldata[which(alldata$plotrow == row &
                               ((thestart >= alldata[,3] & thestart <= alldata[,6]) |
                                  (thestop >= alldata[,3] & thestop <= alldata[,6]) |
                                  (thestart <= alldata[,3] & thestop >= alldata[,6]))),]) == 0)
      {
        return (row)
      }
    }
    return (NA)
  }
  
  # Define a function that plots a looping interaction on a graph
  plotpair <- function(start,end,height,offset=0,flip=FALSE,...)
  {
    x = NULL
    
    posneg = 1
    if (flip == TRUE)
    {
      posneg = -1
      offset = offset * -1
    }
    
    # input variables
    x1 = min(start,end)
    x2 = max(start,end)
    vert_y  = height
    
    # get x coordinate of vertex
    vert_x= (x1 + x2) / 2
    
    # solve for 'a'
    a = vert_y / ((vert_x - x1) * (vert_x - x2))
    
    # plot curve
    curve(offset+posneg * a * (x - x1) * (x - x2),from=x1, to=x2, add=TRUE,...)
  }
  
  # convert data to data frame
  bedpedata = data.frame(bedpedata[,1:6])
  names(bedpedata) = c("chrom1","start1","stop1","chrom2","start2","stop2")
  
  # make sure data is organized correctly
  start1 = apply(bedpedata[,c(2,3)],1,min)
  stop1  = apply(bedpedata[,c(2,3)],1,max)
  start2 = apply(bedpedata[,c(5,6)],1,min)
  stop2  = apply(bedpedata[,c(5,6)],1,max)
  bedpedata$start1 = start1
  bedpedata$stop1  = stop1
  bedpedata$start2 = start2
  bedpedata$stop2  = stop2
  
  if (plottype == "loops")
  {
    # add height column
    bedpedata$heights = heights
  }
  
  # add color column
  bedpedata$color = rep("black",nrow(bedpedata))
  
  # add colorby column
  if (is.null(colorby) == FALSE)
  {
    bedpedata$colorbyvalue = colorby
  }
  
  # add a lwdby column
  if (is.null(lwdby) == FALSE )
  {
    bedpedata$lwdby = lwdby
  }
  
  # change colors
  if (length(color) == 1)
  {
    bedpedata$color = rep(color,nrow(bedpedata))
  }
  if (length(color) > 1)
  {
    bedpedata$color = color
  }
  
  # add lwd column
  bedpedata$lwd = rep(1,nrow(bedpedata))
  
  # change lwd
  if (length(lwd) == 1)
  {
    bedpedata$lwd = rep(lwd,nrow(bedpedata))
  }
  if (length(color) > 1)
  {
    bedpedata$lwd = lwd
  }
  
  # add position columns
  bedpedata$pos1 = apply(bedpedata[,c(2,3)],1,mean)
  bedpedata$pos2 = apply(bedpedata[,c(5,6)],1,mean)
  
  # filter for region of interest 
  bedpedata = bedpedata[which(
    (bedpedata$chrom1 == chrom & (bedpedata$pos1 > chromstart & bedpedata$pos1 < chromend))
    &
      (bedpedata$chrom2 == chrom & (bedpedata$pos2 > chromstart & bedpedata$pos2 < chromend))
  ),]  
  
  if (nrow(bedpedata) == 0)
  {
    plot(c(1,1),type='n',xlab="",ylab="",xaxs = 'i',yaxs='i',xlim=c(chromstart,chromend), ylim=c(0,1),xaxt=xaxt,yaxt=yaxt,bty=bty,...)
    return(list(colorbyrange,colorbycol))
  }
  
  
  # color by
  if (is.null(colorby) == FALSE)
  {
    bedpedata$color = maptocolors(bedpedata$colorbyvalue,colorbycol,range=colorbyrange)
    if (is.null(colorbyrange))
    {
      colorbyrange = c(min(bedpedata$colorbyvalue),max(bedpedata$colorbyvalue))
    }
    
  }
  
  # lwd by
  if (is.null(lwdby) == FALSE )
  {
    bedpedata$lwd = maptolwd(bedpedata$lwdby,range=lwdrange)
  }
  
  if (plottype == "lines")
  {
    
    # sort by distance for prettier plotting
    bedpedata$distance = abs(bedpedata$pos2 - bedpedata$pos1)
    bedpedata = bedpedata[order(bedpedata$distance,decreasing=TRUE),]
    
    # set wiggle
    wiggle = abs( chromend - chromstart) * .04
    
    bedpedata$plotrow = 0
    
    # get row info for bedpes
    for (l in (1:nrow(bedpedata)))
    {
      # figure out which row to plot it in
      bedpedata[l,]$plotrow = checkrow(bedpedata[l,],bedpedata,maxrows=maxrows,wiggle=wiggle)
    }
    
    # count total lines
    totallines = min(maxrows,max(bedpedata$plotrow))
    
    # first make empty plot
    if (flip == FALSE)
    {
      plot(c(1,1),type='n',xlab="",ylab="",xaxs = 'i',yaxs='i',xlim=c(chromstart,chromend), ylim=c(0,totallines+1),xaxt=xaxt,yaxt=yaxt,bty=bty,...)
    }
    if (flip == TRUE)
    {
      plot(c(1,1),type='n',xlab="",ylab="",xaxs = 'i',yaxs='i',xlim=c(chromstart,chromend), ylim=c(-totallines-1,0),xaxt=xaxt,yaxt=yaxt,bty=bty,...)
    } 
    
    # now plot all of the pairs
    for (i in (1:nrow(bedpedata)))
    {
      line = bedpedata$plotrow[i]
      if (flip == TRUE)
      {
        line = -bedpedata$plotrow[i]
        height = -1 * height
      }
      
      rect(xleft=min(bedpedata[i,2] , bedpedata[i,3]),
           ybottom=line-height,
           xright=max(bedpedata[i,2] , bedpedata[i,3]), 
           ytop=line+height,
           col=bedpedata$color[i],
           lwd=bedpedata$lwd[i],
           border=bedpedata$color[i])
      
      rect(xleft=min(bedpedata[i,5] , bedpedata[i,6]),
           ybottom=line-height,
           xright=max(bedpedata[i,5] , bedpedata[i,6]), 
           ytop=line+height,
           col=bedpedata$color[i],
           lwd=bedpedata$lwd[i],
           border=bedpedata$color[i])
      
      segments(x0 =  min(bedpedata[i,2] , bedpedata[i,3],bedpedata[i,5] , bedpedata[i,6]),
               x1 =  max(bedpedata[i,2] , bedpedata[i,3],bedpedata[i,5] , bedpedata[i,6]),
               y0 = line,
               y1 = line,
               lwd = bedpedata$lwd[i],
               col=bedpedata$color[i])
    }
  }
  
  if (plottype == "loops")
  {
    # first make empty plot
    if (flip == FALSE)
    {
      plot(c(1,1),type='n',xlab="",ylab="",xaxs = 'i',yaxs='i',xlim=c(chromstart,chromend), ylim=c(0,(max(bedpedata$heights) + offset )*ymax),xaxt=xaxt,yaxt=yaxt,bty=bty,...)
    }
    if (flip == TRUE)
    {
      plot(c(1,1),type='n',xlab="",ylab="",xaxs = 'i',yaxs='i',xlim=c(chromstart,chromend), ylim=c(-max((bedpedata$heights)-offset)*ymax,0  ),xaxt=xaxt,yaxt=yaxt,bty=bty,...)
    }

    
    # plot the data
    for (row in (1:nrow(bedpedata)))
    {
      x1     = bedpedata$pos1[row]
      x2     = bedpedata$pos2[row]
      height = bedpedata$heights[row]
      color  = bedpedata$color[row]
      lwd    = bedpedata$lwd[row]
      plotpair(x1,x2,height,offset=offset,flip=flip,col=as.character(color),lwd=lwd)
    }
  }
  
  # return color by range and palette
  if (is.null(colorby) == FALSE)
  {
    list(colorbyrange,colorbycol)
  }
}
