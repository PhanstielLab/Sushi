#' highlight loci on a HiC plot
#' 
#'
#' @param x bin or midpoint of one interacting loci
#' @param y bin or midpoint of other interacting loci
#' @param labeltype the units for x and y options are 'bin' or 'bp'
#' @param flip TRUE/FALSE whether plot should be flipped over the x-axis
#' @param diameter of circle as a fraction of the x-axis
#' @param lwd The line width, a positive number, defaulting to 1.
#' @param resolution the width in bp of each pixel
#' @param plottype options are "triangle" or "square"
#' @param half if plottype is set to "square" should the highlight be "above" or "below" the diagonal
#' @param ... values to be passed to \code{\link{plot}}
#' @export
#' @examples
#'

highlighthic <-
function(x,y,labeltype = "bin",flip=FALSE,diameter=.025,lwd=1,
         resolution=NULL,plottype="triangle",half="above",...)
{
  
  # establish coordinates
  if (labeltype == "bin")
  {
    newx = x * resolution - (resolution/2)
    newy = y * resolution - (resolution/2)
  }
  if (labeltype == "bp")
  {
    newx = x
    newy = y
  }
  if (plottype=="triangle")
  {
    finalx = (newx + newy) / 2
    if (flip == FALSE)
    {
      finaly = abs(newx-newy)/2
    }
    if (flip == TRUE)
    {
      finaly = -abs(newx-newy)/2
    }
  }
  
  if (plottype=="square")
  {
    finalx = x
    if (flip == FALSE)
    {
      finaly = y
    }
    if (flip == TRUE)
    {
      finaly = y
    }
    
    if(half == "below")
    {
      finaly = x
      finalx = y
    }
  }
  
  # establish diameter
  if (length(diameter) == 1)
  {
    diameter = rep(diameter,length(finalx))
  }
  size = (par("usr")[2] - par("usr")[1])*diameter
  
  # plot the data
  symbols(x=finalx,y=finaly,circles=size,inches=FALSE,add=TRUE,lwd=lwd,...)
}
