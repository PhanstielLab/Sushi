\name{maptolwd}
\alias{maptolwd}
\title{maps numeric vector to line widths}
\usage{
maptolwd(lwdby, range = c(1, 5))
}
\arguments{
  \item{lwdby}{numeric vector to map to line widths}

  \item{range}{range of values to map}
}
\description{
maps numeric vector to line widths
}
\examples{
plot((1:10),lwd=maptolwd(lwdby=(1:10)))
}

