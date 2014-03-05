#' makes colors transparent (or opaque)
#' 
#'
#' @param color color or colors to make opaque
#' @param transparency value between 0 and 1 indicating desired opaqueness
#' @export
#' @examples
#' plot((1:10),col="red",pch=19)
#' points((10:1),col=opaque("red",transparency=0.3),pch=19)
opaque <- function(color=SushiColors(7)(7),transparency=0.5)
{
  if (length(color) == 1)
  {
    colrgb = col2rgb(color)
    newcol = rgb(colrgb[1],colrgb[2],colrgb[3],alpha=transparency * 255,maxColorValue = 255)
    return(newcol)
  }
  else
  {
    newcols = c()
    for (col in color)
    {
      colrgb = col2rgb(col)
      newcols = c(newcols,rgb(colrgb[1],colrgb[2],colrgb[3],alpha=transparency * 255,maxColorValue = 255))
    }
    return (newcols)
  }
}
  



