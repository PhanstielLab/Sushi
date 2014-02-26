#' Generates a Sushi color palette
#'
#' @param palette The name of the Sushi palette to return.  For list of available palettes try (SushiColors(list))
#' @examples
#' plot(1,xlab='',xaxt='n',ylab='',yaxt='n',xlim=c(0,8),ylim=c(2,8),type='n',bg="grey")
#' for (i in (2:7))
#' {
#'   points(x=(1:i),y=rep(i,i),bg=SushiColors(i)(i),cex=3,pch=21)
#'   
#' }
#' 
#' axis(side=2,at=(2:7),labels=(2:7),las=2)
#' axis(side=1,at=(1:7),labels=(1:7))
#' mtext("SushiColors",side=3,font=2, line=1, cex=1.5)
#' mtext("colors",side=1,font=2, line=2)
#' mtext("palette",side=2,font=2, line=2)
SushiColors <- function(palette="fire")
{
  if (palette == 7 )
  {
    return ( colorRampPalette(c("black","blue","#5900E5","#E5001B","orange","yellow","white")))
  }
  if (palette == 6 )
  {
    return ( colorRampPalette(c("black","blue","#5900E5","#E5001B","orange","yellow")))
  }
  if (palette == 5 )
  {
    return ( colorRampPalette(c("black","blue","#5900E5","#E5001B","orange")))
  }
  if (palette == 4 )
  {
    return ( colorRampPalette(c("blue","#5900E5","#E5001B","orange")))
  }
  if (palette == 3 )
  {
    return ( colorRampPalette(c("blue","#E5001B","orange")))
  }
  
  if (palette == 2 )
  {
    return ( colorRampPalette(c("blue","#E5001B")))
  }

  if (palette == "list")
  {
    print ((2:7))
  }
  
}


