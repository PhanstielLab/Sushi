#' defines chromosome offsets for plotting multi chromosomal plot (eg plotManhattan)
#'
#' @param biomart biomart database to connect to
#' @param host biomart host to connect to if different then www.biomart.org 
#' @param dataset biomart dataset you want to use
#' @param attributes biomart attributes you want to retrieve.
#' @param filters biomart filters (one or more) that should be used in the query.
#' @param values values for the filters
getgeneinfo <-
function(biomart="ensembl",host=NULL,
                        dataset="hsapiens_gene_ensembl",
                        attributes = c("chromosome_name","start_position","end_position","exon_chrom_start","exon_chrom_end","external_gene_id","strand"),
                        filters=c("chromosome_name","start","end","with_ox_refseq_mrna"),
                        values)
{
  
  # set the mart
  if (is.null(host)==TRUE)
  {
    mart=useMart(biomart=biomart, dataset=dataset)
  }
  if (is.null(host)==FALSE)
  {
    mart=useMart(host=host, biomart=biomart, dataset=dataset)
  }
  
  
  # retrieve the vaues
  geneinfo = getBM(attributes = attributes,
                   filters=filters,
                   values=values,
                   mart=mart)
  
  # add a score column
  geneinfo$score = rep(".",nrow(geneinfo))
  
  # reorder the columns into bedpe format
  geneinfo = geneinfo[,c(1,2,3,1,4,5,6,8,7,7)]
  
  # add generic column names
  names(geneinfo) = c("chrom_gene","start_gene","stop_gene","chrom_exon",
                      "start_exon","stop_exon","gene_name","score","strand_gene","strand_transcript")
  
  # return data frame
  return (geneinfo)
}
