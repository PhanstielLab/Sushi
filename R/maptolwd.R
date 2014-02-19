#' maps numeric vector to line widths
#'
#'
#' @param lwdby numeric vector to map to line widths
#' @param range range of values to map
#' @examples
#' plot((1:10),lwd=maptolwd(lwdby=(1:10)))
maptolwd <-
function(lwdby,range=c(1,5))
{
  lwdby =  lwdby - min(lwdby)
  fractionofmax = lwdby / max(lwdby)
  lwds = round(fractionofmax * abs(range[2] - range[1]),digits = 1) + range[1]
  return(lwds) 
}
