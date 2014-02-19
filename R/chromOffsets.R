#' defines chromosome offsets for plotting multi chromosomal plot (eg plotManhattan)
#'
#' @param genome A genome object to be used (2 columns: column 1 = chromosome name, column 2 = length of chromosome)
#' @param space the space in between each chromosome as a fraction of the width of the plot
chromOffsets <-
function(genome,space=0.01)
{
  #### determine the offset for each chromosome ####
  genome.sorted = sortChrom(genome)
  
  # add cum sum column
  cumsums =  cumsum(as.numeric(genome.sorted[,2]))
  spacer  = cumsums[length(cumsums)] * space
  aditionalspace = (1:length(cumsums)-0) * spacer
  
  # start pos
  startpos = c(0,cumsums[1:length(cumsums)-1])
  #startpos = cumsums[1:length(cumsums)]
  startpos = startpos + aditionalspace
  genome.sorted[,3] = startpos  
  
  # stop pos
  stoppos = cumsums + (1:(length(cumsums)))*spacer
  
  return(data.frame(chrom=genome.sorted[,1],size=genome.sorted[,3],start=startpos,stop=stoppos))
}
