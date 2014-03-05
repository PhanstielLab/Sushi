#' adds a letter and a title to a plot
#'
#' This function adds a letter and a title (both are optional) to the top of a plot.  Udeful for generating paper figures.
#'  
#' @param letter A string, typically a letter or number (eg 'A', 'A)', '1', etc) to lable the plot with
#' @param title A string for a plot title
#' @param letteradj adj of letter. See \code{\link{par}}
#' @param titleadj adj of title. See \code{\link{par}}
#' @param letterfont font of letter. See \code{\link{par}}
#' @param titlefont font of title See \code{\link{par}}
#' @param lettercex cex of letter. See \code{\link{par}}
#' @param titlecex cex of title See \code{\link{par}}
#' @param letterline line of letter. See \code{\link{par}}
#' @param titleline line of title See \code{\link{par}}
#' @param lettercol color of letter. See \code{\link{par}}
#' @param titlecol color of title See \code{\link{par}}
#' @export
#' @examples
#' 
#' par(mar=c(3,3,3,3))
#' plot((1:10),col=maptocolors(vec=(1:10),colorRampPalette(c("blue","red"))),pch=19,cex=4)
#' labelplot("A)"," sample plot",lettercex=2,titlecex=2,titlecol="blue")
#' 
labelplot <- function(letter=NULL,title=NULL,letteradj=-0.05,titleadj=0.0,letterfont=2,titlefont=2,
                      lettercex=1.2,titlecex=1,letterline=0.5,titleline=0.5,lettercol="black",titlecol="black")
{
  if (is.null(letter) == FALSE)
  {
    mtext(letter,side=3, adj=letteradj,line=letterline,font=letterfont,cex=lettercex,col=lettercol)
  }
  if (is.null(title) == FALSE)
  {
    mtext(title,side=3, adj=titleadj,line=titleline,font=titlefont,cex=titlecex,col=titlecol)
  }
}


