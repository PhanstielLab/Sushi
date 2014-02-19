library('Sushi')
pdfname = "~/Dropbox/Sushi/Figure_dhp_new.pdf"
Sushi_data = data(package = 'Sushi')
data(list = Sushi_data$results[,3]) 
makepdf = FALSE

############################################################################
#                                                                          #
#                             DEFAULT PALETTES                             #
#                                                                          #
############################################################################

# Default plotHic palette
palette_fireNice = colorRampPalette(c("#3362A5", "dodgerblue1", "deepskyblue", "gray99", "gold", "firebrick2","#A31D1D"))
palette_fire = colorRampPalette(c("black","blue","#5900E5","#E5001B","orange","yellow","white"))
palette_fire_dark = colorRampPalette(c("black","blue","#5900E5","#E5001B","orange"))

# Default PlotBedpe palettes
palette_plotBed = colorRampPalette(c("firebrick3","orange","deepskyblue2","dodgerblue2","dodgerblue4"))
palette_bedpe = colorRampPalette(rev(c("#3362A5", "dodgerblue1", "deepskyblue", "gold", "firebrick2","#A31D1D")), space="Lab")

# Default PlotBed palettes
palette_plotBedDensity =list(
  colorRampPalette(c("black",palette_plotBed(5)[1])),
  colorRampPalette(c("black",palette_plotBed(5)[2])),
  colorRampPalette(c("black",palette_plotBed(5)[3])),
  colorRampPalette(c("black",palette_plotBed(5)[4])),
  colorRampPalette(c("black",palette_plotBed(5)[5])))


############################################################################
#                                                                          #
#                                   CODE                                   #
#                                                                          #
############################################################################


if (makepdf ==TRUE)
{
  pdf(pdfname,height=10,width=12)
}


# make a layout for all of the plots
layout(matrix(c(1,1,1,1,
                1,1,1,1,
                2,2,8,8,
                2,2,9,9,
                3,3,10,10,
                3,3,10,10,
                4,4,11,11,
                4,4,11,11,
                5,5,12,12,
                5,5,12,12,
                6,7,13,13,
                6,7,14,14
),12, 4, byrow = TRUE))
par(mgp=c(3, .3, 0))


#################################################################################################
#                                                                                               #
#                                 (A) manhattan plot                                            #
#                                                                                               #
#################################################################################################

# set the margins
par(mar=c(3,4,3,2))

# set the genomic regions
chrom1            = "chr11"
chromstart1       = 500000
chromend1         = 5050000

chrom2            = "chr15"
chromstart2       = 73000000
chromend2         = 89500000


# make the manhattan plot
plotManhattan(bedfile=Sushi_GWAS.bed,pvalues=Sushi_GWAS.bed[,5],genome=Sushi_hg18_genome,col=palette_fire_dark(nrow(Sushi_hg18_genome)),cex=0.75)

# add zoom 1
zoomsregion(region=c(chromstart1,chromend1),chrom=chrom1,genome=Sushi_hg18_genome, zoomborder = "black", lty=2,lwd = 1,extend=c(0.07,0.2),wideextend=0.2,offsets=c(0,.535))

# add zoom 2
zoomsregion(region=c(chromstart2,chromend2),chrom=chrom2,genome=Sushi_hg18_genome, zoomborder = "black", lty=2,lwd = 1,extend=c(0.07,0.2),wideextend=0.2,offsets=c(.535,0))

# add labels
labelgenome(genome=Sushi_hg18_genome,side=1,scipen=20,n=4,scale="Mb",edgeblankfraction=0.20,line=.18,chromline=.5,scaleline=0.5)

# add y-axis
axis(side=2,las=2,tcl=.2)
mtext("log10(P)",side=2,line=1.75,cex=.75,font=2)

# Add plot label
mtext("A)   GWAS",side=3, adj=-.025,line=1,font=2)


#################################################################################################
#                                                                                               #
#                                 (B) Hi-C                                                      #
#                                                                                               #
#################################################################################################

# set the margins
par(mar=c(3,4,2,2))

# set the genomic regions
chrom            = "chr11"
chromstart       = 500000
chromend         = 5050000
zoomregion       = c(1700000,2350000)

# plot the HiC data
phic = plotHic(Sushi_HiC.matrix,chrom,chromstart,chromend,max_y = 20,zrange=c(0,28),palette = palette_fire,flip=FALSE)

# add labels
labelgenome(chrom,chromstart,chromend,side=1,scipen=20,n=4,scale="Mb",edgeblankfraction=0.20,line=.18,chromline=.5,scaleline=0.5)

# add the legend
addlegend(phic[[1]],palette=phic[[2]],title="score",side="right",bottominset=0.4,topinset=0,xoffset=-.035,labelside="left",width=0.025,title.offset=0.035)

# add zoom
zoomsregion(region=zoomregion, zoomborder = "black", lty=2,lwd = 1,extend=c(0.05,0.25))

# add zoombox
zoombox(zoomregion=zoomregion,lty=2,lwd = 1,col = "black")

# Add plot label
mtext("B)    HiC",side=3, adj=-0.065,line=0.5,font=2)


#################################################################################################
#                                                                                               #
#                                 (C) 5C                                                        #
#                                                                                               #
#################################################################################################

# set the margins
par(mar=c(3,4,2,2))

# set the genomic regions
chrom            = "chr11"
chromstart       = 1650000
chromend         = 2350000


# plot the loops
pbpe = plotbedpe(Sushi_5C.bedpe,chrom,chromstart,chromend,heights = Sushi_5C.bedpe$score,offset=0,flip=FALSE,bty='n',
                 lwd=1,plottype="loops",colorby=Sushi_5C.bedpe$samplenumber,colorbycol=colorRampPalette(c("#5900E5","#E5001B","orange")))

# add zoombox
zoombox(passthrough=TRUE,lty=2,lwd = 1,col="black")

# add the genome labels
labelgenome(chrom, chromstart,chromend,side=1,scipen=20,n=3,scale="Mb",line=.18,chromline=.5,scaleline=0.5)

# add the legend
legend("topright",inset =0.01,legend=c("K562","HeLa","GM12878"),col=c("#5900E5","#E5001B","orange"),pch=19,bty='n',text.font=2)

# add y-axis
axis(side=2,las=2,tcl=.2)
mtext("Z-score",side=2,line=1.75,cex=.75,font=2)

# Add plot label
mtext("C)    5C",side=3, adj=-0.065,line=0.5,font=2)

#################################################################################################
#                                                                                               #
#                                 (D) ChIA PET (PolII)                                          #
#                                                                                               #
#################################################################################################

# set the margins
par(mar=c(3,4,2,2))

# set the genomic regions
chrom            = "chr11"
chromstart       = 1650000
chromend         = 2350000

# plot the loops
pbpe = plotbedpe(Sushi_ChIAPET_pol2.bedpe,chrom,chromstart,chromend,flip=TRUE,bty='n',lwd=1,plottype="lines",
                 colorby=abs(Sushi_ChIAPET_pol2.bedpe$start1-Sushi_ChIAPET_pol2.bedpe$start2),colorbycol=palette_fire_dark)

# add the genome labels
labelgenome(chrom, chromstart,chromend,side=1,scipen=20,n=4,scale="Mb",line=.18,chromline=.5,scaleline=0.5)

# add the legend
addlegend(pbpe[[1]],palette=pbpe[[2]],title="distance (bp)",side="right",bottominset=0.05,topinset=0.35,xoffset=-.035,labelside="left",width=0.025,
          title.offset=0.08,labels.digits=0,)

# add zoombox
zoombox(passthrough=TRUE,lty=2,lwd = 1,col="black")

# Add plot label
mtext("D)    ChIA-PET (Pol2)",side=3, adj=-0.065,line=0.5,font=2)

#################################################################################################
#                                                                                               #
#                                 (E) DnaseI                                                    #
#                                                                                               #
#################################################################################################

# set the margins
par(mar=c(3,4,2,2))

# set the genomic regions
chrom            = "chr11"
chromstart       = 1650000
chromend         = 2350000
zoomregion1      = c(1955000,1965000)
zoomregion2      = c(2281200,2282200)

# overlapping, transparent, and rescaled
plotBedgraph(Sushi_DNaseI.bedgraph,chrom,chromstart,chromend,transparency=1.0,color="#5900E5",lwd=1,linecol="#5900E5")

# add zoom 1
zoomsregion(zoomregion1,col=NA,zoomborder="black",lty=2,lwd=1,extend=c(0.01,0.18),wideextend=0.10,offsets=c(0,.577))

# add the genome labels
labelgenome(chrom,chromstart,chromend,side=1,scipen=20,n=4,line=.18,chromline=.5,scaleline=0.5,scale="Mb")

# add zoombox
zoombox(zoomregion=zoomregion,lty=2,lwd = 1,col="black")

# add zoom 2
zoomsregion(zoomregion2,col=NA,zoomborder="black",lty=2,lwd=1,extend=c(0.01,0.18),wideextend=0.10,offsets=c(.577,0))

# add y-axis
axis(side=2,las=2,tcl=.2)
mtext("Read Depth",side=2,line=1.75,cex=.75,font=2)

# Add plot label
mtext("E)    DnaseI",side=3, adj=-0.065,line=0.5,font=2)


#################################################################################################
#                                                                                               #
#                                 (F) ChIP-Seq ChIP Exo                                         #
#                                                                                               #
#################################################################################################

# set the genomic regions
chrom            = "chr11"
chromstart       = 1650000
chromend         = 2350000
zoomregion1      = c(1955000,1965000)
zoomregion2      = c(2281200,2282200)

# plot chip-seq data
plotBedgraph(Sushi_ChIPSeq_CTCF.bedgraph,chrom,zoomregion1[1],zoomregion1[2],transparency=.50,flip=FALSE,color="blue",linecol="blue")

# Add plot label
mtext("F)    DnaseI / ChIP-Seq",side=3, adj=-0.6,line=0.5,font=2)

# plot dnaseI data
plotBedgraph(Sushi_DNaseI.bedgraph,chrom,zoomregion1[1],zoomregion1[2],transparency=.50,flip=FALSE,color="#E5001B",linecol="#E5001B",overlay=TRUE,rescaleoverlay=TRUE)

# add the genome labels
labelgenome(chrom,zoomregion[1],zoomregion[2],side=1,scipen=20,n=3,line=.18,chromline=.5,scaleline=0.5,scale="Mb")

# add zoombox
zoombox(zoomregion = NULL,lwd = 1,col="black")

# set the legend colors
transparency = 0.5
col1 = col2rgb("blue")
finalcolor1 = rgb(col1[1],col1[2],col1[3],alpha=transparency * 255,max = 255)
col2 = col2rgb("#E5001B")
finalcolor2 = rgb(col2[1],col2[2],col2[3],alpha=transparency * 255,max = 255)

# add legend
legend("topright",inset=0.025,legend=c("DnaseI","ChIP-seq (CTCF)"),fill=c(finalcolor1,finalcolor2),border=c("blue","#E5001B"),text.font=2,cex=0.75)



#################################################################################################
#                                                                                               #
#                                 (G) Bed Pile up                                               #
#                                                                                               #
#################################################################################################

# set the genomic regions
chrom            = "chr11"
chromstart       = 1650000
chromend         = 2350000
zoomregion1      = c(1955000,1965000)
zoomregion2      = c(2281200,2282200)

# plt the chip-seq data as a pile-up
plotBed(beddata    = Sushi_ChIPSeq_pol2.bed,chrom = chrom,chromstart = zoomregion2[1],chromend =zoomregion2[2],
        colorby    = Sushi_ChIPSeq_pol2.bed$strand,colorbycol = colorRampPalette(c("blue","#E5001B")),row  = "auto",wiggle=0.001,splitstrand=FALSE, type = "region")

# add the genome labels
labelgenome(chrom,zoomregion2[1],zoomregion2[2],side=1,scipen=20,n=2,line=.18,chromline=.5,scaleline=0.5,scale="Mb")

# add zoombox
zoombox(zoomregion = NULL,lwd = 1,col="black")

# add legend
legend("topright",inset=0.025,legend=c("forward","reverse"),fill=c("#E5001B","blue"),border=c("#E5001B","blue"),text.font=2,cex=0.75)

# Add plot label
mtext("G)    ChIP-Seq",side=3, adj=-.275,line=0.5,font=2)

#################################################################################################
#                                                                                               #
#                                 (H) manhattan plot zoomed                                     #
#                                                                                               #
#################################################################################################

# set the margins
par(mar=c(0.1,4,2,2))

# set the genomic regions
chrom            = "chr15"
chromstart       = 60000000
chromend         = 80000000
chromstart2       = 72000000
chromend2         = 74000000

# make the manhattan plot
plotManhattan(bedfile=Sushi_GWAS.bed,chrom=chrom2,chromstart= chromstart, chromend = chromend,    pvalues=Sushi_GWAS.bed$pval.GC.DBP,col=palette_fire_dark(nrow(Sushi_hg18_genome))[15],cex=0.75)

# add zoom in
zoomsregion(region=c(chromstart2,chromend2),chrom=chrom2,genome=NULL, zoomborder = "black", lty=2,lwd = 1,extend=c(0.05,1),offsets=c(0.0,0))

# add zoom box
zoombox(passthrough=TRUE,lty=2,lwd = 1,col="black",topextend=5)

# add y-axis
axis(side=2,las=2,tcl=.2)
mtext("Z-score",side=2,line=1.75,cex=.75,font=2)

# Add plot label
mtext("H)    GWAS",side=3, adj=-0.065,line=0.5,font=2)


#################################################################################################
#                                                                                               #
#                                 (I) Gene density                                               #
#                                                                                               #
#################################################################################################

# set the margins
par(mar=c(3,4,1.8,2))
par(mar=c(3,4,1.8,2))

# set the genomic regions
chrom            = "chr15"
chromstart       = 60000000
chromend         = 80000000
chrom_biomart    = gsub("chr","",chrom)

# set the mart (since we want hg18 coordinates)
mart=useMart(host='may2009.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')

# get just gene info
geneinfobed = getBM(attributes = c("chromosome_name","start_position","end_position"),
                    filters= c("chromosome_name","start","end"),
                    values=list(chrom_biomart,chromstart,chromend),
                    mart=mart)

# add "chr" to the chrom column
geneinfobed[,1] = paste("chr",geneinfobed[,1],sep="")

# plot gene density
plotBed(beddata = geneinfobed[!duplicated(geneinfobed),],chrom = chrom,chromstart = chromstart,row='supplied',
        chromend =chromend,palettes = palette_plotBedDensity, type = "density")

#label genome
labelgenome(chrom=chrom, chromstart, chromend  ,side=1,scipen=20,n=4,scale="Mb",edgeblankfraction=0.10,line=.18,chromline=.5,scaleline=0.5)

# add zoom in
zoomsregion(region=c(chromstart2,chromend2),chrom=chrom2,genome=NULL, zoomborder = "black", lty=2,lwd = 1,extend=c(2,1.0),wideextend=.75,offsets=c(0.0,0))

# add zoombox
zoombox(zoomregion=c(chromstart2,chromend2),lty=2,lwd = 1,col="black",topextend=5)

# Add plot label
mtext("I)    Gene Density",side=3, adj=-0.065,line=0.20,font=2)


#################################################################################################
#                                                                                               #
#                                 (J) RNA seq                                                    #
#                                                                                               #
#################################################################################################

# set the margins
par(mar=c(3,4,2,2))

# set the genomic regions
chrom2            = "chr15"
chromstart2       = 72800000
chromend2         = 73100000
zoomregion       = c(72998000,73020000)
chrom2_biomart   = 15

# plot transcripts
pg = plotgenes(Sushi_transcripts.bed,chrom2,chromstart2,chromend2 ,types = Sushi_transcripts.bed$type,
          colorby=Sushi_transcripts.bed$score,colorbycol= palette_fire_dark ,labeltext=FALSE,
          maxrows=50,height=0.4,plotgenetype="box",
          zoomregion = zoomregion1, zoomborder = "grey", lty=2,lwd = 3,zoomextend=c(0.01,0.15))

# label genome
labelgenome( chrom2, chromstart2,chromend2,side=1,scipen=20,n=3,scale="Mb",line=.18,chromline=.5,scaleline=0.5)

# add the legend
addlegend(pg[[1]],palette=pg[[2]],title="FPKM",side="right",bottominset=0.4,topinset=0,xoffset=-.035,labelside="left",width=0.025,title.offset=0.055)

# add zoombox
zoombox(passthrough=TRUE,lty=2,lwd = 1,col="black")

# add zoom
zoomsregion(region=zoomregion,lty=2,lwd = 1,extend=c(-.025,1),zoomborder="black")

# Add plot label
mtext("J)    RNA-seq",side=3, adj=-0.065,line=0.5,font=2)

#################################################################################################
#                                                                                               #
#                                 (K) ChIP Seq peaks                                            #
#                                                                                               #
#################################################################################################

# set the margins
par(mar=c(3,4,2,2))

# set the genomic regions
chrom            = "chr15"
chromstart      = 72800000
chromend         = 73100000
zoomregion       = c(72998000,73020000)

Sushi_ChIPSeq_severalfactors.bed$color = palette_fire_dark(max(Sushi_ChIPSeq_severalfactors.bed$row))[Sushi_ChIPSeq_severalfactors.bed$row]

# plot it
plotBed(beddata    = Sushi_ChIPSeq_severalfactors.bed,chrom = chrom,chromstart = chromstart,chromend =chromend,
        rownumber  = Sushi_ChIPSeq_severalfactors.bed$row, type = "circles",color=Sushi_ChIPSeq_severalfactors.bed$color,row="given",plotbg="grey95",
        rowlabels=unique(Sushi_ChIPSeq_severalfactors.bed$name),rowlabelcol=unique(Sushi_ChIPSeq_severalfactors.bed$color),rowlabelcex=0.75)

# label genome
labelgenome(chrom,chromstart,chromend,side=1,scipen=20,n=3,line=.18,chromline=.5,scaleline=0.5,scale="Mb")

# add zoom
zoombox(zoomregion = zoomregion,lwd = 1,col="black")

# add zoom in
zoomsregion(region=zoomregion,chrom=chrom, zoomborder = "black", lty=2,lwd = 1,extend=c(0.5,.22),wideextend=0.15,offsets=c(0.0,0))

# Add plot label
mtext("K)    ChIP-seq",side=3, adj=-0.065,line=0.5,font=2)



#################################################################################################
#                                                                                               #
#                                 (L) Pol2 bedgrpah                             #
#                                                                                               #
#################################################################################################

# set the margins
par(mar=c(3,4,2,2))

# set the genomic regions
chrom            = "chr15"
chromstart       = 72998000
chromend         = 73020000

# plot the Pol2 bedgraph data
plotBedgraph(Sushi_ChIPSeq_pol2.bedgraph, chrom,chromstart,chromend,transparency=1.0,flip=FALSE,color="blue",linecol="blue")

# label genome
labelgenome(chrom,chromstart,chromend,side=1,scipen=20,n=3,line=.18,chromline=.5,scaleline=0.5,scale="Mb")

# add zoombox
zoombox(passthrough=TRUE,lty=2,lwd = 1,col="black")

# add y-axis
axis(side=2,las=2,tcl=.2)
mtext("Read Depth",side=2,line=1.75,cex=.75,font=2)

# Add plot label
mtext("L)    Chip-Seq (Pol2)",side=3, adj=-.065,line=0.5,font=2)

#################################################################################################
#                                                                                               #
#                                 (M) RNA-seq bedgrpah                             #
#                                                                                               #
#################################################################################################

# set the margins
par(mar=c(2,4,.5,2))

# set the genomic regions
chrom            = "chr15"
chromstart       = 72998000
chromend         = 73020000

# plot the K562 RNAseq bedgraph data
plotBedgraph(Sushi_RNASeq_K562.bedgraph, chrom,chromstart,chromend,transparency=1.0,flip=FALSE,color="red",linecol="red")

# label genome
labelgenome(chrom,chromstart,chromend,side=1,scipen=20,n=3,line=.18,chromline=.5,scaleline=0.5,scale="Mb")

# add zoombox
zoombox(passthrough=TRUE,lty=2,lwd = 1,col="black")

# Add plot label
mtext("M)    RNA-seq",side=3, adj=-0.065,line=-1,font=2)

#################################################################################################
#                                                                                               #
#                                 (N) Gene Structures                                           #
#                                                                                               #
#################################################################################################

# set the margins
par(mar=c(3,4,.5,2))

# set the genomic region
chrom            = "chr15"
chromstart       = 72998000
chromend         = 73020000

# plot gene structures
plotgenes(Sushi_genes.bed,chrom_biomart,chromstart,chromend ,types=Sushi_genes.bed$type,
          maxrows=1,height=0.5,plotgenetype="arrow",bentline=FALSE,col="blue",
          labeloffset=1,fontsize=1.2)

# label genome
labelgenome( chrom, chromstart,chromend,side=1,scipen=20,n=3,scale="Mb",line=.18,chromline=.5,scaleline=0.5)

# add zoombox
zoombox(lty=2,lwd = 1,col="black")

# Add plot label
mtext("N)    Gene Structures",side=3, adj=-0.065,line=-.5,font=2)


if (makepdf ==TRUE)
{
  dev.off()
}



