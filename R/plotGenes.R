#' plots gene structure or transcript structures
#' 
#'
#' @param geneinfo gene info stored in a bed-like format.  If NULL it will look up genes in the region using biomart (with biomart="ensembl" and dataset="hsapiens_gene_ensembl"). See also \code{\link{useMart}}
#' @param chrom chromosome of region to be plotted
#' @param chromstart start position
#' @param chromend end position
#' @param col single value or vector specifying colors of gene structures
#' @param colorby vector to scale colors by
#' @param colorbycol palette to apply color scale to (only valid when colorby is not NULL)
#' @param colorbyrange the range of values to apply the color scale to.  Values outside that range will be set to the limits of the range.
#' @param bheight the height of the boxes drawn for exons
#' @param lheight the height of the bent line is bent is set to TRUE
#' @param bentline TRUE/FALSE indicating whether lines between exons should be bent
#' @param packrow TRUE / FALSE indicating whether genes should be packed or whether each gene should be plotted on its own row
#' @param types single value or vector specifying types of elements (acceptable values are 'exon','utr')
#' @param plotgenetype String specifying whether the genes should resemble a 'box' or a 'arrow'
#' @param arrowlength value (between 0 and 1) specifying the length of the tail of each arrow as a fraction of the total plot width (only valid when plotgenetype is set to "arrow")
#' @param wigglefactor the fraction of the plot to leave blank on either side of each element to avoid overcrowding.
#' @param maxrows The maximum number of rows to plot on the y-axis
#' @param labeltext TRUE/FALSE indicating whether genes should be labeled
#' @param labeloffset value (between 0 and 1) specifying the vertical offset of gene labels
#' @param fontsize font size of gene labels
#' @param fonttype font type of gene labels
#' @param labelat postion along gene to place labels (acceptable values are "middle","start",and "end")
#' @param ... values to be passed to \code{\link{plot}}
#' @export
#' @examples
#' 
#' data(Sushi_genes.bed)
#' 
#' chrom            = "chr15"
#' chromstart       = 72998000
#' chromend         = 73020000
#' chrom_biomart    = 15
#' 
#' plotGenes(Sushi_genes.bed,chrom_biomart,chromstart,chromend ,types=Sushi_genes.bed$type,
#'      maxrows=1,height=0.5,plotgenetype="arrow",bentline=FALSE,col="blue",
#'      labeloffset=1,fontsize=1.2)
#' 
#' labelgenome( chrom, chromstart,chromend,side=1,scipen=20,n=3,scale="Mb",line=.18,chromline=.5,scaleline=0.5)
#'
plotGenes <-
function(geneinfo=NULL, chrom=NULL, chromstart=NULL,chromend=NULL,
                      col=SushiColors(2)(2)[1],bheight=0.3,lheight=0.3,bentline=TRUE,
                      packrow=TRUE,maxrows=10000,
                      colorby=NULL,colorbyrange=NULL, colorbycol=colorRampPalette(c("blue","red")),
                      types="exon",plotgenetype = "box",
                      arrowlength=.005,wigglefactor=0.05,
                      labeltext=TRUE,labeloffset=0.4,fontsize=.7,fonttype=2,labelat="middle",...)
{
  
  # if the gene info is nullusing current human annotations
  if (is.null(geneinfo)==TRUE)
  {
    # grab info
    mart=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    chrom_biomart    = gsub("chr","",chrom)
    geneinfo = getBM(attributes = c("chromosome_name","exon_chrom_start","exon_chrom_end","external_gene_id","strand"),
                           filters= c("chromosome_name","start","end"),
                           values=list(chrom_biomart,chromstart,chromend),
                           mart=mart)
    
    # make names the same
    names(geneinfo) = c("chrom","start","stop","gene","strand")
    
    # reorder and make proper bed format
    geneinfo$score = "."
    geneinfo = geneinfo[,c(1,2,3,4,6,5)]    
  }

  # convert strand info 
  if (ncol(geneinfo) >= 6)
  {
    geneinfo[,6] = convertstrandinfo(geneinfo[,6])
  }
  
  # Define a function that merges overlapping regions
  mergetypes <- function(exons)
  {
    newexons = c()
    for (subtype in names(table(exons[,7])))
    {
      # get subtype
      subexons = exons[which(exons[,7] == subtype),]
      subexons = subexons[order(subexons[,2]),]
      
      # skip if none
      if (nrow(subexons) == 0)
      {
        next
      }
      
      # check for overlap
      i = 0
      for (j in (1:nrow(subexons)))
      {
        i = i + 1
        if (i > 1)
        {
          
          if (subexons[j,2] <= curstop)
          {
            subexons[j,2] = curstart
            subexons[j,3] = max(subexons[j,3],curstop)
          }
        }
        curstart = subexons[j,2]
        curstop = subexons[j,3]
      }
      
      # make sure all starts and stops are the same
      for (startpos in names(table(subexons[,2])))
      {
        subexons[which(subexons[2] == startpos),3] = max(subexons[which(subexons[2] == startpos),3])
      }
      
      # remove duplicate rows
      subexons = subexons[!duplicated(subexons),]
      newexons = rbind(newexons,subexons)
    }
    return (newexons)
  }
  

  # Define a function that plots the gene structure given a y-value and all exons
  plottranscript <- function(exons,col,yvalue,bheight,lheight,bentline=TRUE,border="black",
                             arrowlength,bprange,strandlength=0.04,strandarrowlength=0.10,plotgenetype="box",
                             labeltext=TRUE,labeloffset=0.4,fontsize=.7,fonttype=2,labelat="middle",...)
  {
    strand = exons[1,6]
    
    # if label is true add the label
    
   labellocation = mean(c(max(exons[,c(2,3)]),min(exons[,c(2,3)])))
   if (labeltext == TRUE)
   {
     if ((labelat == "start" & strand == 1) || (labelat == "end" & strand == -1)  )
     {
       labellocation = min(exons[,c(2,3)])
     }
     if ((labelat == "start" & strand == -1) || (labelat == "end" & strand == 1)  )
     {
       labellocation =  max(exons[,c(2,3)])
     }
      
      adj = 0.0
      if (strand == 1)
      {
        adj = 1.0
      }
      text(labellocation,yvalue+labeloffset,labels=exons[1,4],adj=adj,cex=fontsize,font=fonttype)
      arrows(labellocation+strand*bprange*strandlength/4,yvalue+labeloffset,
             labellocation+strand*bprange*strandlength,yvalue+labeloffset,
             length=strandarrowlength)
   }
    
    # make sure coordinates are in the correct order
    min = apply(exons[,c(2,3)],1,min)
    max = apply(exons[,c(2,3)],1,max)
    exons[,2] = min
    exons[,3] = max
    
    # merge types
    exons = mergetypes(exons) 
    
    # sort by order
    exons = exons[order(exons[,2],exons[,3]),]
    
    # combine all types of exons for proper line plotting
    allexons =  exons
    allexons$types = "exon"
    allexons = mergetypes(allexons)
    
    if (nrow(allexons) > 1)
    {
      # organize line info
      linecoords = cbind(allexons[(1:nrow(allexons)-1),3],allexons[(2:nrow(allexons)),2])
      linecoords = cbind(linecoords,apply(linecoords,1,mean))
      
      # plot lines
      top = yvalue 
      if (bentline == TRUE)
      {
        top = yvalue + lheight
      }
      
      for (i in (1:nrow(linecoords)))
      {
        if ((linecoords[i,2] - linecoords[i,1]) > 0 )
        {
          segments(x0 = linecoords[i,1],
                   x1 = linecoords[i,3],
                   y0 = yvalue,
                   y1 = top,
                   col     = col,
                   ...)
          segments(x0 = linecoords[i,3],
                   x1 = linecoords[i,2],
                   y0 = top,
                   y1 = yvalue,
                   col     = col,
                   ...)
        }
      }
    }
    
    # add boxes
    for (i in (1:nrow(exons)))
    {
      height = bheight
      if (exons[i,7] == "utr")
      {
        height = bheight/2
      }
      
      # boxes
      if (plotgenetype == "box")
      {
        rect(xleft    = exons[i,2],
             ybottom  = yvalue - height,
             xright   = exons[i,3],
             ytop    = yvalue + height,
             col     = col,
             border  = col,
             ...)
      }
      
      # arrows
      if (plotgenetype == "arrow")
      {
        # strand
        strand = exons[i,6]
        
        offset = bprange * arrowlength * strand * -1
        
        x= c(exons[i,2],
             exons[i,2]+offset,
             exons[i,3]+offset,
             exons[i,3],
             exons[i,3]+offset,
             exons[i,2]+offset,
             exons[i,2]
        )
        
        y= c(yvalue,
             yvalue+height,
             yvalue+height,
             yvalue,
             yvalue-height,
             yvalue-height,
             yvalue)
        
        polygon(x=x , y=y ,col=col,border=col, ...)
      }     
    }
  }
  
  
  
  # Define a function th determines which row to plot a gene on
  checkrow <- function(data,alldata,maxrows,strand,wiggle=0,plotgenetype="box",arrowlength=0.005,bprange=0)
  {

    startcol = 2
    stopcol  = 3

    
    strand = data[1,5]
    
    for (row in (1:maxrows))
    {
      
      if (plotgenetype == "arrow")
      {
        if (strand == -1)
        {
          thestart = as.numeric(data[startcol]) - wiggle - (bprange * arrowlength * 3)
          thestop  = as.numeric(data[stopcol]) + wiggle
        }
        if (strand == 1)
        {
          thestart = as.numeric(data[startcol]) - wiggle
          thestop  = as.numeric(data[stopcol]) + wiggle + (bprange * arrowlength * 3)
        }
      }
      
      if (plotgenetype == "box")
      {
        thestart = as.numeric(data[startcol]) - wiggle
        thestop  = as.numeric(data[stopcol])  + wiggle
      }
      
      if (nrow(alldata[which(alldata$plotrow == row &
                               ((thestart >= alldata[,startcol] & thestart <= alldata[,stopcol]) |
                                  (thestop >= alldata[,startcol] & thestop <= alldata[,stopcol]) |
                                  (thestart <= alldata[,startcol] & thestop >= alldata[,stopcol]))),]) == 0)
      {
        return (row)
      }
    }
    return (NA)
  }
  
  # remove unwanted columns
  geneinfo = geneinfo[,seq(1,6)]
  
  # establish start and stop columns
  startcol = 2
  stopcol  = 3
  
  # add types
  geneinfo$types = types
  
  # remove lines with NA
  geneinfo = geneinfo[which(is.na(geneinfo[,2]) ==FALSE),]
  
  # color by 
  if (is.null(colorby) == FALSE)
  {
    if (is.null(colorbyrange) == TRUE)
    {
      colorbyrange = c(min(colorby),max(colorby))
    }
    geneinfo$colors = maptocolors(colorby,col=colorbycol,num=100,range=colorbyrange)
  }
  
  if (is.null(colorby) == TRUE)
  {
    geneinfo$colors = col
  }
  
  # set xlim
  if (is.null(chrom) == TRUE | is.null(chromstart) == TRUE | is.null(chromend) == TRUE)
  {
    chrom      = geneinfo[1,1]
    chromstart = min(geneinfo[,c(startcol,stopcol)])
    chromend   = max(geneinfo[,c(startcol,stopcol)])
    extend     = abs(chromend - chromend) * 0.04
    chromstart = chromstart - extend
    chromend   = chromend   + extend
  }
  
  # get bprange
  bprange = abs(chromend - chromstart)
  
  # set wiggle
  wiggle = bprange * wigglefactor
  
  # get number of geneinfo
  numberofgeneinfo = length(names(table(geneinfo[,4])))
  namesgeneinfo    = names(table(geneinfo[,4]))
  
  # sort by length
  starts = c()
  stops  = c()
  sizes  = c()
  strands = c()
  
  
  # collect the info for each transcript
  for (i in (1:numberofgeneinfo))
  {
    subgeneinfo  = geneinfo[which(geneinfo[,4] == namesgeneinfo[i]),]
    starts = c(starts,min(subgeneinfo[,2:3]))
    stops  = c(stops, max(subgeneinfo[,2:3]))
    sizes  = c(sizes, stops[i] - starts[i])
    strands = c(strands,subgeneinfo[1,6])
  }
  
  transcriptinfo = data.frame(names=namesgeneinfo,starts=starts,stops=stops,sizes=sizes,strand=strands)
  transcriptinfo = transcriptinfo[order(sizes,decreasing=TRUE),]
  

  
  # get row information
  if (packrow == TRUE)
  {
    transcriptinfo$plotrow = 0
    for (i in (1:nrow(transcriptinfo)))
    {
      transcriptinfo$plotrow[i] = checkrow(transcriptinfo[i,],transcriptinfo,maxrows=maxrows,wiggle=wiggle,plotgenetype=plotgenetype,arrowlength=arrowlength)
    }
  }
  if (packrow == FALSE)
  {
    transcriptinfo$plotrow = seq(1:nrow(transcriptinfo))
  }
  

  
  # make the empty plot
  offsettop = 0.5
 
  # filter out rows above max row 
  transcriptinfo = transcriptinfo[which(is.na(transcriptinfo$plotrow)==FALSE),]
  
  # filter out transcrits that don't overlap region
  transcriptinfo = transcriptinfo[which((transcriptinfo[,2] > chromstart & transcriptinfo[,2] < chromend)
                       | (transcriptinfo[,3] > chromstart & transcriptinfo[,3] < chromend)),]

  if (nrow(transcriptinfo) == 0)
  {
    toprow = 1
  }
  if (nrow(transcriptinfo) > 0)
  {
    toprow = max(transcriptinfo$plotrow)
  }
  
  plot(c(1,1),xlim=c(chromstart,chromend),ylim=c(0.5,(toprow  + offsettop)),type ='n',bty='n',xaxt='n',yaxt='n',ylab="",xlab="",xaxs="i")
  
  
  if (nrow(transcriptinfo) > 0)
  {
    for (i in (1:nrow(transcriptinfo)))
    {
      subgeneinfo  = geneinfo[which(geneinfo[,4] == transcriptinfo[i,1]),]
      plottranscript(subgeneinfo,col=subgeneinfo$colors[1],yvalue=transcriptinfo$plotrow[i],bheight=bheight,lheight=lheight,bentline=bentline,border=col,
                     bprange=bprange,arrowlength=arrowlength,plotgenetype=plotgenetype,
                     labeltext=labeltext,labeloffset=labeloffset,fontsize=fontsize,fonttype=fonttype,labelat=labelat)
    }
  }
  
  # return color by range and palette
  return(list(colorbyrange,colorbycol))

}
