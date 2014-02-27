#' Converts strand info to 1 / -1
#'
#' @param strandvector vector of strand information to convert from +/- to 1/-1 if neccesary
convertstrandinfo <- function(strandvector)
{
  if (class((strandvector)) == "character")
  {
    print ("yes")
    b = rep(1,length((strandvector)))
    b[which((strandvector) == "-")] = -1
    strandvector = b
  }
  return (strandvector)
}