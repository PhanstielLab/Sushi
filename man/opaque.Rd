\name{opaque}
\alias{opaque}
\title{makes colors transparent (or opaque)}
\usage{
opaque(color = SushiColors(7)(7), transparency = 0.5)
}
\arguments{
  \item{color}{color or colors to make opaque}

  \item{transparency}{value between 0 and 1 indicating
  desired opaqueness}
}
\description{
makes colors transparent (or opaque)
}
\examples{
plot((1:10),col="red",pch=19)
points((10:1),col=opaque("red",transparency=0.3),pch=19)
}

