#' sort chromosome files by chom name
#'
#' @param genome A genome object to be used (2 columns: column 1 = chromosome name, column 2 = length of chromosome)
sortChrom <-
function(genome)
{
  # first remove the 'chr'
  data_nochr = genome
  data_nochr[,1] = gsub("chr", "", data_nochr[,1])
  
  # sort the data
  numericvalues = suppressWarnings(as.numeric(data_nochr$V1))
  
  # sort the numeric ones
  data_nochr_num = data_nochr[which(is.na(numericvalues) == FALSE),]
  data_nochr_num_sorted = data_nochr_num[order(numericvalues[which(is.na(numericvalues) == FALSE)]),]
  
  # sort the non-numeric ones
  data_nochr_nonnum = data_nochr[which(is.na(numericvalues) == TRUE),]
  data_nochr_nonnum_sorted = data_nochr_nonnum[order(data_nochr_nonnum[,1]),]
  
  result = rbind(data_nochr_num_sorted,data_nochr_nonnum_sorted)
  result[,1] =  paste("chr",result[,1],sep="")
  
  return (result)
}
