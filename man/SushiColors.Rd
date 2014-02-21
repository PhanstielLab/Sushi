\name{SushiColors}
\alias{SushiColors}
\title{Generates a Sushi color palette}
\usage{
SushiColors(palette = "fire")
}
\arguments{
  \item{palette}{The name of the Sushi palette to return.
  For list of available palettes try (SushiColors(list))}
}
\description{
Generates a Sushi color palette
}
\examples{
plot((1:10),col=SushiColors("fire")(10),pch=19,cex=3)
}

