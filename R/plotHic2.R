#' plots HiC interaction matrix
#' 
#'
#' @param hicdata interaction matrix representing HiC data.  Row and column names should be postions along a chromosome
#' @param chrom chromosome of region to be plotted
#' @param chromstart start position
#' @param chromend end position
#' @param max_dist The maximum distance of interaction to plot
#' @param zrange The range of interaction scores to plot (more extreme value will be set to the max or min)
#' @param palette color palette to use for representing interaction scores
#' @param flip TRUE/FALSE whether plot should be flipped over the x-axis
#' @param format The input format of the HiC data.  Options are "sparse" (all data in first 3 colums) or "full" (large matrix)
#' @param labeltype options are "bin" or "bp"
#' @param plottype options are "triangle" or "square"
#' @param resolution the width in bp of each pixel
#' @param add TRUE/FALSE whether plot should be added to an existing plot
#' @param half Should the plot be on the top, bottom, or both (onlly applies when plottype == "square)"
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
plotHic2 <-
  function(hicdata,chrom,chromstart,chromend,max_dist = 100000,zrange=NULL,
           palette = SushiColors(7),flip=FALSE,
           format = "auto",labeltype = "bin",resolution=NULL,plottype="triangle",
           half="both",add=FALSE
  )
  {
    
    # define a function that takes midpoints, resolution, and color and draws a ploygon
    drawpoly <- function(df,res,flip=FALSE,plottype="triangle",half="both")
    {
      col = rgb(df[4],df[5],df[6],maxColorValue=255)
      if (plottype=="triangle")
      {
        x = (df[1] + df[2])/2
        if (flip == FALSE)
        {
          y = abs(df[2]-df[1])/2
        }
        if (flip == TRUE)
        {
          y = -abs(df[2]-df[1])/2
        }
        
        xs = unlist(c(x,x-.5*res,x,x+.5*res,x))
        ys = unlist(c(y-.5*res,y,y+.5*res,y,y-.5*res))
        polygon (xs,ys,col=col,border=NA)
      }
      
      if (plottype=="square")
      {
        if (half == "both")
        {
          x = df[1]
          y = df[2]
          rect(xleft=x-.5*res, ybottom=y-.5*res, xright=x+.5*res, ytop=y+.5*res,col=col,border=NA)
        }
        if (half == "top")
        {
          x = df[1]
          y = df[2]
          if (y >x)
          {
            rect(xleft=x-.5*res, ybottom=y-.5*res, xright=x+.5*res, ytop=y+.5*res,col=col,border=NA)
          }
          if (y == x)
          {
            xs = c(x-.5*res,x-.5*res,x+.5*res,x-.5*res)
            ys = c(y-.5*res,y+.5*res,y+.5*res,y-.5*res)
            polygon (xs,ys,col=col,border=NA)
          }
        }
        if (half == "bottom")
        {
          x = df[1]
          y = df[2]
          if (y < x)
          {
            rect(xleft=x-.5*res, ybottom=y-.5*res, xright=x+.5*res, ytop=y+.5*res,col=col,border=NA)
          }
          if (y == x)
          {
            xs = c(x-.5*res,x+.5*res,x+.5*res,x-.5*res)
            ys = c(y-.5*res,y+.5*res,y-.5*res,y-.5*res)
            polygon (xs,ys,col=col,border=NA)
          }
        }
        
      }
    }
    
    # convert bins to bp
    if (is.null(resolution) == TRUE)
    {
      return ("Error: resolution is a required parameter")
    }
    
    # 
    if (plottype=="square")
    {
      max_dist = Inf
    }
    
    # establish format type
    if (format == "auto")
    {
      if (ncol(head(hicdata)) <= 6)
      {
        format = "sparse"
      }
      if (ncol(head(hicdata)) > 6)
      {
        format = "full"
      }
      
      print (paste("HiC data format auto detected and set to",format))
    }
    
    # Deal with full format
    if (format == "full")
    {
      if (labeltype == "bin")
      {
        
        # convert bins to bp
        rownames(hicdata) = ceiling(as.numeric(rownames(hicdata)) * resolution - (resolution/2))
        colnames(hicdata) = ceiling(as.numeric(colnames(hicdata)) * resolution - (resolution/2))
      }
      
      # filter for region of interest
      rows = as.numeric(rownames(hicdata))
      cols = as.numeric(colnames(hicdata))
      paddedchromstart = chromstart - resolution
      paddedchromend   = chromend + resolution
      hicregion        = hicdata[which(rows >= chromstart & rows <= chromend) ,
                                 which(cols >= chromstart & cols <= chromend)]
      
      # convert to 3 column format
      hicregion = data.frame(col = rep(colnames(hicregion), each = nrow(hicregion)), 
                             row = rep(rownames(hicregion), ncol(hicregion)), 
                             value = as.vector(as.matrix(hicregion)))
    }
    
    
    # Deal with sparse format
    if (format == "sparse")
    {
      if (labeltype == "bin")
      {
        # convert bins to bp
        hicdata[,1] = ceiling(as.numeric(hicdata[,1]) * resolution - (resolution/2))
        hicdata[,2] = ceiling(as.numeric(hicdata[,2]) * resolution - (resolution/2))
      }

      # filter for region
      hicregion        = hicdata[which(hicdata[,1] >= chromstart & hicdata[,1] <= chromend &
                                 hicdata[,2] >= chromstart & hicdata[,2] <= chromend & abs(hicdata[,1]-hicdata[,2]) <= max_dist ),]
    }
    
    # make sure columns are numbers
    hicregion[,1] = as.numeric(as.character(hicregion[,1]))
    hicregion[,2] = as.numeric(as.character(hicregion[,2]))
    hicregion[,3] = as.numeric(as.character(hicregion[,3]))
 
    # set the zrange
    if (is.null(zrange) == TRUE)
    {
      max_z = max(hicregion[,3])
      min_z = min(hicregion[,3])
    }
    if (is.null(zrange) == FALSE)
    {
      min_z = zrange[1]
      max_z = zrange[2]
    }
    
    # map to colors and covnert to rgb
    hicregioncolors = maptocolors(hicregion[,3],palette,num=100,range=zrange)   
    hicregion       = cbind(hicregion,t(col2rgb(hicregioncolors)))
    
    # plot polygons
    if (plottype =="triangle")
    {
      # make an empty plot
      if (flip == FALSE)
      {
        plot(1,1,xlim=c(chromstart,chromend),ylim=c(0,max_dist/2),type='n',xaxs='i',yaxs='i',bty='n',xaxt='n',yaxt='n',xlab="",ylab="")
      }
      
      if (flip == TRUE)
      {
        plot(1,1,xlim=c(chromstart,chromend),ylim=c(-max_dist/2,0),type='n',xaxs='i',yaxs='i',bty='n',xaxt='n',yaxt='n',xlab="",ylab="")
      }
      
      # plot the polygons
      apply(hicregion,1,drawpoly,res=resolution,flip=flip,plottype)
    }
    
    # plot polygons
    if (plottype =="square")
    {
      # make symmetric
      hicregionrev = hicregion[,c(2,1,3,4,5,6)]
      hicregion = rbind(hicregion,setNames(hicregionrev,names(hicregion)))
      hicregion = hicregion[!duplicated(hicregion),]
      
      if (add == FALSE)
      {
        # make an empty plot
        plot(1,1,xlim=c(chromstart,chromend),ylim=c(chromstart,chromend),type='n',xaxs='i',yaxs='i',bty='n',xaxt='n',yaxt='n',xlab="",ylab="")
      }
        
      # plot the polygons
      apply(hicregion,1,drawpoly,res=resolution,flip=flip,plottype,half=half)
    }
    
    # return color by range and palette
    return(list(c(min_z,max_z),palette))
  }
