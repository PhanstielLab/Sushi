#' plots data stored in bed file format
#' 
#'
#' @param beddata genomic data to be plotted (in bed format)
#' @param chrom chromosome of region to be plotted
#' @param chromstart start position
#' @param chromend end position
#' @param type type of plot ('region','circles','density')
#' @param colorby vector to scale colors by
#' @param colorbycol palette to apply color scale to (only valid when colorby is not NULL)
#' @param colorbyrange the range of values to apply the color scale to.  Values outside that range will be set to the limits of the range.
#' @param rownumber vector giving the row numbers of each bed element to be plotted.
#' @param row How row number should be determined.  Appropriate values are 'auto' or 'supplied'
#' @param height Value, typically between 0 and 1, that sets the height of each bed element
#' @param plotbg The background color of the plot
#' @param wiggle the fraction of the plot to leave blank on either side of each element to avoid overcrowding.
#' @param splitstrand TRUE/FALSE indicating whether reverse strnad bed elements shold be plotted below the x axis.  (only valid when row is set to 'auto')
#' @param numbins  The number of bins to divide the region into when type is set to density (only valid when type is set to 'density')
#' @param binsmoothing  umber of bins to sum together when type is set to density  (only valid when type is set to 'density')
#' @param palettes list of color palettes used for density plots.  Each row can have a unique palette.  number of palettes is less than the number of rows then only the first palette is used (only valid when type is set to 'density')
#' @param rowlabels labels for the y-axis
#' @param rowlabelcol color of the y-axis labels
#' @param rowlabelfont font of the y-axis labels
#' @param rowlabelcex font size of the y-axis labels
#' @param maxrows The maximum number of rows to plot on the y-axis
#' @param color single color or vector of colors to use to plot the points or regions (not valid when type is set to 'density')
#' @param xaxt A character which specifies the x axis type.  See \code{\link{par}}
#' @param yaxt A character which specifies the y axis type.  See \code{\link{par}}
#' @param xlab Label for the x-axis
#' @param ylab Label for the y-axis
#' @param xaxs Must be set to 'i' for appropriate integration into Sushi plots.  See \code{\link{par}}
#' @param yaxs Must be set to 'i' for appropriate integration into Sushi plots.  See \code{\link{par}}
#' @param bty A character string which determined the type of box which is drawn about plots.  See \code{\link{par}}
#' @param border border color drawn around each bed element or density bin.  Set to 'n' for none.
#' @param ... values to be passed to other functions
#' @examples
#' data(Sushi_ChIPSeq_severalfactors.bed)
#' chrom            = "chr15"
#' chromstart      = 72800000
#' chromend         = 73100000
#' Sushi_ChIPSeq_severalfactors.bed$color = heat.colors(max(Sushi_ChIPSeq_severalfactors.bed$row))[Sushi_ChIPSeq_severalfactors.bed$row]
#' plotBed(beddata    = Sushi_ChIPSeq_severalfactors.bed,chrom = chrom,chromstart = chromstart,chromend =chromend,
#'         rownumber  = Sushi_ChIPSeq_severalfactors.bed$row, type = "circles",color=Sushi_ChIPSeq_severalfactors.bed$color,row="given",plotbg="grey95",
#'         rowlabels=unique(Sushi_ChIPSeq_severalfactors.bed$name),rowlabelcol=unique(Sushi_ChIPSeq_severalfactors.bed$color),rowlabelcex=0.75)
#' 
#' Sushi_ChIPSeq_severalfactors.bed$color = heat.colors(max(Sushi_ChIPSeq_severalfactors.bed$row))[Sushi_ChIPSeq_severalfactors.bed$row]
#' 
#' plotBed(beddata    = Sushi_ChIPSeq_severalfactors.bed,chrom = chrom,chromstart = chromstart,chromend =chromend,
#'         rownumber  = Sushi_ChIPSeq_severalfactors.bed$row, type = "region",color=Sushi_ChIPSeq_severalfactors.bed$color,row="given",plotbg="grey95",
#'         rowlabels=unique(Sushi_ChIPSeq_severalfactors.bed$name),rowlabelcol=unique(Sushi_ChIPSeq_severalfactors.bed$color),rowlabelcex=0.75)
#' 
#' colors = c("dodgerblue1","firebrick2","violet","yellow",
#'         "dodgerblue1","firebrick2","violet","yellow",
#'         "dodgerblue1","firebrick2","violet")
#' plotBed(beddata    = Sushi_ChIPSeq_severalfactors.bed,chrom = chrom,chromstart = chromstart,chromend =chromend,
#'          rownumber  = Sushi_ChIPSeq_severalfactors.bed$row, type = "density",row="supplied",
#'          rowlabels=unique(Sushi_ChIPSeq_severalfactors.bed$name),rowlabelcol=colors,rowlabelcex=0.75,
#'          palettes=list(
#'          colorRampPalette(c("black",colors[1])),
#'          colorRampPalette(c("black",colors[2])),
#'          colorRampPalette(c("black",colors[3])),
#'          colorRampPalette(c("black",colors[4])),
#'          colorRampPalette(c("black",colors[5])),
#'          colorRampPalette(c("black",colors[6])),
#'          colorRampPalette(c("black",colors[7])),
#'          colorRampPalette(c("black",colors[8])),
#'          colorRampPalette(c("black",colors[9])),
#'          colorRampPalette(c("black",colors[10])),
#'          colorRampPalette(c("black",colors[11]))))

plotBed <-
function(beddata,chrom,chromstart,chromend,type="region",
                    colorby = NULL,colorbycol = NULL,colorbyrange=NULL,
                    rownumber = NULL,
                    row="auto",height=0.4,
                    plotbg     = "white",
                    wiggle = 0.02,
                    splitstrand = FALSE,
                    numbins = 200, binsmoothing=10,palettes = topo.colors,
                    rowlabels = NULL,
                    rowlabelcol= "dodgerblue2",
                    rowlabelfont = 2,
                    rowlabelcex = 1,
                    maxrows=1000000,color="dodgerblue4",xaxt='none',yaxt='none',xlab="",ylab="",xaxs="i",yaxs="i",bty='n',border=NA,...)
{
  
  # convert strand info 
  if (ncol(beddata) >= 6)
  {
    beddata[,6] = convertstrandinfo(beddata[,6])
  }
  
  # Define a function that determines which row to plot a gene on
  checkrow <- function(data,alldata,maxrows,wiggle=0)
  {
    for (row in (1:maxrows))
    {
      thestart = min(as.numeric(data[c(2,3)])) - wiggle
      thestop  = max(as.numeric(data[c(2,3)])) + wiggle
      
      if (nrow(alldata[which(alldata$plotrow == row &
                               ((thestart >= alldata[,2] & thestart <= alldata[,3]) |
                                  (thestop >= alldata[,2] & thestop <= alldata[,3]) |
                                  (thestart <= alldata[,2] & thestop >= alldata[,3]))),]) == 0)
      {
        return (row)
      }
    }
    return (NA)
  }
  
  if (is.null(color) == FALSE)
  {
    beddata$plotcolor = color
  }
  
  # color by
  if (is.null(colorby) == FALSE)
  {
    beddata$plotcolor = colorby
  }

  
  if (is.null(rownumber) == FALSE)
  {
    beddata$plotrow = rownumber
  }
  if (is.null(rownumber) == TRUE)
  {
    beddata$plotrow = rep(1,nrow(beddata))
  }
  
  # determine the number of rows
  numberofrows = length(unique(beddata$plotrow))
  
  # filter for region of interest 
  beddata = beddata[which(beddata[,1] == chrom & ((beddata[,2] > chromstart & beddata[,2] < chromend) |  (beddata[,3] > chromstart & beddata[,3] < chromend))),]

  
  # color by
  if (is.null(colorby) == FALSE)
  {
    beddata$plotcolor = maptocolors(beddata$plotcolor,colorbycol,range=colorbyrange)
  }
  
  wiggle = abs( chromend - chromstart) * wiggle
  
  if (row == "auto")
  {
    # randomize the order
    beddata=beddata[sample(nrow(beddata)),]
    
    beddata$plotrow = rep(0,nrow(beddata))
    
    if (splitstrand == TRUE)
    {
      plusstrand  = beddata[which(beddata[,6] > 0),]
      minusstrand = beddata[which(beddata[,6] < 0),]
      for (i in (1:nrow(plusstrand)))
      {
        plusstrand[i,]$plotrow = checkrow(plusstrand[i,],plusstrand,maxrows,wiggle=wiggle)
      }
      for (i in (1:nrow(minusstrand)))
      {
        minusstrand[i,]$plotrow = checkrow(minusstrand[i,],minusstrand,maxrows,wiggle=wiggle)
      }
      minusstrand$plotrow = -1 * minusstrand$plotrow
      beddata = rbind(plusstrand,minusstrand)
    }
    if (splitstrand == FALSE)
    {
      for (i in (1:nrow(beddata)))
      {
        beddata[i,]$plotrow = checkrow(beddata[i,],beddata,maxrows,wiggle=wiggle)
      }
    }
    numberofrows = max(beddata$plotrow)
  }
  
  if (nrow(beddata) == 0)
  {
    print ("Warning: No bed elements in selected range")
  }
  
  # plot the bed
  if (type == "region")
  {
    # make the empty plot
    plot(1,1,xlim=c(chromstart,chromend),ylim=c((min(beddata$plotrow)-(height*1.2)),min(maxrows, max(numberofrows))+(height*1.2)),type ='n',bty='n',xaxt='n',yaxt='n',ylab="",xlab="",xaxs="i")
    
    # set background color
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = plotbg,border=NA)
    
    if (splitstrand == TRUE)
    {
      abline(h=0)
    }
      
    for (bed in (1:nrow(beddata)))
    {
      line = beddata$plotrow[bed]
      rect(xleft=beddata[bed,2] , ybottom=line-height, xright=beddata[bed,3], ytop=line+height,col=beddata[bed,]$plotcolor,border=beddata[bed,]$plotcolor,...)
    }
  }
  
  # plot the bed
  if (type == "circles")
  {
    # make the empty plot
    plot(1,1,xlim=c(chromstart,chromend),ylim=c((min(beddata$plotrow)-(height*1.2)),min(maxrows, max(beddata$plotrow))+(height*1.2)),type ='n',bty='n',xaxt='n',yaxt='n',ylab="",xlab="",xaxs="i")
    
    # set background color
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = plotbg,border=NA)
    
    if (splitstrand == TRUE)
    {
      abline(h=0)
    }
    
    beddata$xval = apply(beddata[,c(2,3)],1,mean)
    points(x=beddata$xval,y=beddata$plotrow,pch=19,col=beddata$plotcolor,new=FALSE,bty='n',xaxt='n',yaxt='n',ylab="",xlab="",xaxs="i")
  }
  
  
  # type == density
  if (type == "density")
  {
    # Define a function that draws a rectangle for heatmap plotting
    drawRect <- function(vector,row)
    {
      xleft    = vector[1]
      ybottom  = row - 0.5
      xright   = vector[2]
      ytop     = row + 0.5
      color    = vector[3]
      rect(xleft, ybottom, xright, ytop,col=color,border=NA)
    }
    
    # check for enough palettes
    if (length(palettes) < numberofrows)
    {
      print ("less palettes than samples")
      print ("Using first palette for all samples")
      for (i in (1:numberofrows))
      {
        palettes[[i]] = palettes[[1]]
      }
    }
    
    # make the empty plot
    plot(c(1,1),xlim=c(chromstart,chromend),ylim=c(0.5,numberofrows+0.5),type ='n',bty='n',xaxt='n',yaxt='n',ylab="",xlab="",xaxs="i",yaxs="i")
    
    # add each row
    for (row in (1:numberofrows))
    {
      # filter for sample
      samplebeddata = beddata[which(beddata$plotrow == row),]
      
      # establish break points
      bins = seq(chromstart,chromend,length.out=numbins)
      
      # get center of peaks
      samplebeddata$xval = apply(samplebeddata[,c(2,3)],1,mean)
      
      # adjust numbers to fit in range
      samplebeddata$xval[which(samplebeddata$xval < min(chromstart,chromend))] = min(chromstart,chromend)
      samplebeddata$xval[which(samplebeddata$xval > max(chromstart,chromend))] = max(chromstart,chromend)
      
      # get values in bins
      bindata = hist(samplebeddata$xval,breaks=bins,plot = FALSE)$counts
      # smooth bins
      binsdatasmooth = rollapply(bindata, binsmoothing,sum,partial=TRUE)
      
      # build plotting matrix 
      xyvalues = data.frame(x = bins[1:(length(bins)-1)],
                            y = bins[2:length(bins)],
                            c = maptocolors(binsdatasmooth,col=palettes[[row]],num=100)
      )
      
      # draw heatmap
      apply(xyvalues,1, drawRect,row=row)
      
    }
  }
  
  # Define a function that labels rows of a plotBed plot
  labelRowsBed <- function(rows,labels,font=2,adj=1.2,col="dodgerblue4",cex=1)
  {
    orginalxpd = par()$xpd
    par(xpd=TRUE)
    text(x=rep(par('usr')[1],length(rows)),y=rows, labels=labels,adj=adj,col=col,font=font,cex=cex)
    par(xpd=orginalxpd)
  }
  
  # add lables
  labelRowsBed(rows=seq(1:numberofrows),labels=rowlabels,col=rowlabelcol,font=rowlabelfont,cex=rowlabelcex)  
}
