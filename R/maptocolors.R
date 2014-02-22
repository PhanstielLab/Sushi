#' maps numeric vector to color palette
#'
#'
#' @param vec numeric vector to map to color
#' @param col color palette to which to be mapped
#' @param num number of bins of colors
#' @param range range of values to map
#' @export
#' @examples
#' plot((1:10),col=maptocolors(vec=(1:10),colorRampPalette(c("blue","red"))),pch=19,cex=4)
maptocolors <-
function(vec,col,num=100,range=NULL)
{
  if (is.null(range) == TRUE)
  {
    breaks <- seq(min(vec), max(vec),length.out=num)
  }
  if (is.null(range) == FALSE)
  {
    vec[which(vec < range[1])] = range[1]
    vec[which(vec > range[2])] = range[2]
    breaks <- seq(range[1], range[2],length.out=num)
  }
  
  cols <- col(length(breaks) + 1)
  colvec = as.character(cut(vec, c(-Inf, breaks, Inf), labels=cols))
  return(colvec)
}
