#' Generates a Sushi color palette
#'
#' @param palette The name of the Sushi palette to return.  For list of available palettes try (SushiColors(list))
#' @examples
#' plot((1:10),col=SushiColors("fire")(10),pch=19,cex=3)
SushiColors <- function(palette="fire")
{
  if (palette == "fire" )
  {
    return ( colorRampPalette(c("black","blue","#5900E5","#E5001B","orange","yellow","white")))
  }
  if (palette == "firenowhite" )
  {
    return ( colorRampPalette(c("black","blue","#5900E5","#E5001B","orange","yellow")))
  }
  if (palette == "firedark" )
  {
    return ( colorRampPalette(c("black","blue","#5900E5","#E5001B","orange")))
  }
  if (palette == "two" )
  {
    return ( colorRampPalette(c("blue","#E5001B")))
  }
  if (palette == "three" )
  {
    return ( colorRampPalette(c("blue","#E5001B","orange")))
  }

  if (palette == "list")
  {
    print (c("fire","firedark","two","three"))
  }
  
}

