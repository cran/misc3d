\name{slices3d}
\alias{slices3d}
\title{Interactive Image Slices of 3D or 4D Volume Data}
\description{
  Uses \pkg{tkrplot} to create an interactive slice view of three or four
  dimensional volume data.
}
\usage{
slices3d(vol, main, scale = 0.8, col=gray.colors(512), cross = TRUE)
}
\arguments{
  \item{vol}{a three or four dimensional real array.}
  \item{main}{a character vector; main title for the plot.} 
  \item{scale}{real value for scaling embedded plot size.}
  \item{col}{a list of colors.}
  \item{cross}{logical; if \code{TRUE}, show cross hairs of current slices.}
}
\details{
  Shows slices of 3D array along the axes as produced by \code{image},
  along with sliders for controlling which slices are shown.  For 4D
  data an additional slider selects the value of the fourth index.
  This is useful for viewing medical imaging data (e.g. PET scans and
  fMRI data).
}
\examples{
  #Example 1: View of a mixture of three tri-variate normal densities
  nmix3 <- function(x, y, z, m, s) {
      0.4 * dnorm(x, m, s) * dnorm(y, m, s) * dnorm(z, m, s) +
      0.3 * dnorm(x, -m, s) * dnorm(y, -m, s) * dnorm(z, -m, s) +
      0.3 * dnorm(x, m, s) * dnorm(y, -1.5 * m, s) * dnorm(z, m, s)
  }
  x<-seq(-2, 2, len=40)
  g<-expand.grid(x = x, y = x, z = x)
  v<-array(nmix3(g$x,g$y,g$z, .5,.5), c(40,40,40))
  slices3d(v,main="View of a mixture of three tri-variate normals", col=heat.colors(256))

\dontrun{
  #Example 2: View of FMRI data.
  library(AnalyzeFMRI)
  a<-f.read.analyze.volume(system.file("example.img", package="AnalyzeFMRI"))
  slices3d(aperm(a), main="View of FMRI")
  }
}
\keyword{hplot}
