#' plots HiC interactio matrix
#' 
#'
#' @param hicdata interaction matrix representing HiC data.  Row and column names should be postions along a chromosome
#' @param chrom chromosome of region to be plotted
#' @param chromstart start position
#' @param chromend end position
#' @param max_y The maximum bin distance to plot
#' @param zrange The range of interaction scores to plot (more extreme value will be set to the max or min)
#' @param palette color palette to use for representing interaction scores
#' @param flip TRUE/FALSE whether plot should be flipped over the x-axis
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
plotHic <-
function(hicdata,chrom,chromstart,chromend,max_y = 30,zrange=NULL,palette = SushiColors(7),flip=FALSE)
{
  # grab only regions of interest
  rows = as.numeric(rownames(hicdata))
  cols = as.numeric(colnames(hicdata))
  
  hicregion        = hicdata[which(rows >= chromstart & rows <= chromend) ,
                         which(cols >= chromstart & cols <= chromend)]
  
  # determine number of bins
  nbins = nrow(hicregion)
  stepsize = abs(chromstart-chromend) / (2*nbins)
  
  # convert to matrix if it isnt't already
  hicm = as.matrix(hicregion)
  
  if (is.null(zrange) == TRUE)
  {
    max_z = max(hicm)
    min_z = min(hicm)
  }
  if (is.null(zrange) == FALSE)
  {
    min_z = zrange[1]
    max_z = zrange[2]
  }
  
  # map to colors
  hicmcol = matrix(maptocolors(hicm,palette,num=100,range=zrange),nrow=nrow(hicm))
  
  # make an empty plot
  if (flip == FALSE)
  {
    plot(1,1,xlim=c(chromstart,chromend),ylim=c(0,max_y),type='n',xaxs='i',yaxs='i',bty='n',xaxt='n',yaxt='n',xlab="",ylab="")
  }
  
  if (flip == TRUE)
  {
    plot(1,1,xlim=c(chromstart,chromend),ylim=c(-max_y,0),type='n',xaxs='i',yaxs='i',bty='n',xaxt='n',yaxt='n',xlab="",ylab="")
  }
  
  # fill plot
  for (rownum in (1:nrow(hicm)))
  {
    y = -.5
    if (flip == TRUE)
    {
      y = .5
    }
    
    x = chromstart + (rownum * 2 * stepsize) - (stepsize * 2)
    for (colnum in (rownum:ncol(hicm)))
    {
      x = x + stepsize
      
      if (flip == FALSE)
      {
        y = y + .5
      }
      if (flip == TRUE)
      {
        y = y - .5
      }
      xs = c(x-stepsize,x,x+stepsize,x,x-stepsize)
      ys = c(y,y+.5,y,y-.5,y)
      polygon(xs,ys,border=NA,col=hicmcol[colnum,rownum])
    }
  }
  
  # return color by range and palette
  return(list(c(min_z,max_z),palette))
}
