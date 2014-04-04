#' adds a legend to a Sushi plot
#'
#' This function adds a legend to Sushi plots that have a colorby function (e.g. plotHic, plotGenes, and plotBedpe)
#'
#' @param range the rang of values to be plotted.  ie c(min,max)
#' @param title title of values to be mapped
#' @param labels.digits Number of digits after the decimal point to include in labels 
#' @param palette color palette to use
#' @param side side of plot to place legend ('right','left')
#' @param labelside side of legend to place legend title
#' @param xoffset fraction of plot to offset the legend
#' @param width width as a fraction of the plot width 
#' @param bottominset inset from the bottom of the blot as a fraction of the plot width
#' @param topinset inset from the top of the blot as a fraction of the plot width
#' @param tick.num desired number of tickmarks
#' @param tick.length length of tick marks
#' @param txt.font font type of legend text
#' @param txt.cex font size of legned text
#' @param title.offset offset of title from the key
#' @param title.font font type of legend title
#' @param title.cex font size of legned text
#' @export
#' @examples
#' 
#' data(Sushi_HiC.matrix)
#' 
#' chrom            = "chr11"
#' chromstart       = 500000
#' chromend         = 5050000
#' 
#' phic = plotHic(Sushi_HiC.matrix,chrom,chromstart,chromend,max_y = 20,zrange=c(0,28),palette = topo.colors,flip=FALSE)
#' 
#' labelgenome(chrom,chromstart,chromend,side=1,scipen=20,n=4,scale="Mb",edgeblankfraction=0.20,line=.18,chromline=.5,scaleline=0.5)
#' 
#' addlegend(phic[[1]],palette=phic[[2]],title="score",side="right",bottominset=0.4,topinset=0,xoffset=-.035,labelside="left",width=0.025,title.offset=0.035)
#'
addlegend <-
  function(range,title="",labels.digits=1,palette=topo.colors,side="right",labelside="left",
           xoffset=0.1,width=0.05,bottominset=0.025,topinset=0.025,tick.num = 5,tick.length=0.01,
           txt.font=1,txt.cex = .75,
           title.offset=0.05,title.font=2,title.cex = 1)
  {
    # allow plotting outside the plotting area
    originalxpd = par()$xpd
    par(xpd=TRUE)
    
    # get colors
    scalecol = maptocolors(seq(range[1],range[2],length.out=99),palette,num=100)
    
    #plot boxes
    minxpos = par('usr')[1]
    maxxpos = par('usr')[2]
    minypos = par('usr')[3]
    maxypos = par('usr')[4]
    xrange = abs(maxxpos - minxpos)
    yrange = abs(maxypos - minypos)
    
    if (side == "right")
    {
      xleft  = maxxpos + xrange * xoffset
      xright = xleft   + xrange * width
    }
    
    if (side == "left")
    {
      xright = minxpos - xrange * xoffset
      xleft  = xright  - xrange * width
    }
    
    # determine the y-coordinates of box edges
    rectedges = seq(minypos+yrange*bottominset,maxypos-yrange*topinset,length.out=length(scalecol) +1)
    
    # get the top and bottom y-values for each legend box
    recttops  = rectedges[2:length(rectedges)]
    rectbots  = rectedges[1:length(rectedges)- 1]
    
    # draw legend
    for (i in (1:length(scalecol)))
    {
      rect(xleft =  xleft, ybottom = rectbots[i], xright = xright, ytop = recttops[i],col=scalecol[i],border=NA)
    }
    
    # draw legend border
    rect(xleft= xleft , ybottom = rectbots[1], xright = xright, ytop = recttops[length(recttops)],col=NA,border="black")
    
    
    if (labelside == "left")
    {
      tickleft  = xleft-tick.length*xrange
      tickright = xleft
      txt.x     = xleft-tick.length*xrange
      txt.adj=1.2
      # add legend label
      text(x= xleft-tick.length*xrange - title.offset*xrange,y=mean(rectedges),labels = title,srt=90,cex=title.cex,font=title.font)
    }
    if (labelside == "right")
    {
      tickleft  =  xright
      tickright = xright + tick.length*xrange
      txt.x     = xright + tick.length*xrange
      txt.adj=-.2
      # add legend label
      text(x= xright+tick.length*xrange + title.offset*xrange,y=mean(rectedges),labels = title,srt=270,cex=title.cex,font=title.font)
    }
    
    # add ticks
    labels = seq(range[1],range[2],length.out=tick.num)
    labels = round(labels,digits=labels.digits)
    j = 0 
    for (i in seq(rectedges[1],rectedges[length(rectedges)],length.out = tick.num))
    {
      j = j + 1
      segments(x0 = tickleft, y0 = i, x1 = tickright, y1 = i,col="black")
      text(x= txt.x,y=i,adj=txt.adj,labels =  labels[j],cex=txt.cex,font=txt.font)
    }
    
    # go back to original xpd setting
    par(xpd=originalxpd)
  }
