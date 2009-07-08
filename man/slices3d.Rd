\name{slices3d}
\alias{slices3d}
\title{Interactive Image Slices of 3D or 4D Volume Data}
\description{
  Uses \pkg{tkrplot} to create an interactive slice view of three or four
  dimensional volume data.
}
\usage{
slices3d(vol1, vol2=NULL, rlim1, rlim2, col1, col2, main,
         scale = 0.8, alpha=1, cross = TRUE,
         layout=c("counterclockwise", "clockwise"))
}
\arguments{
  \item{vol1}{a three or four dimensional real array. If two images are
    overlaid, then this is the one at bottom.}
  \item{vol2}{a three or four dimensional real array. If two images are
    overlaid, then this is the one on top. The default value is
    \code{NULL}, when only \code{vol1} is drawn.}
  \item{rlim1}{the minimum and maximum \code{vol1} values for which
    colors should be plotted, defaulting to the range of the values of
    \code{vol1}.}
  \item{rlim2}{the minimum and maximum \code{vol2} values for which
    colors should be plotted, defaulting to the range of the values of
    \code{vol2}, if two images are overlaid.}
  \item{col1}{a list of colors for \code{vol1}.}
  \item{col2}{a list of colors for \code{vol2}.}
  \item{main}{a character vector; main title for the plot.} 
  \item{scale}{real value for scaling embedded plot size.}
  \item{alpha}{real value for transparency level, if two images are
    overlaid. The default value is 1.}
  \item{cross}{logical; if \code{TRUE}, show cross hairs of current
	slices.}
  \item{layout}{a character string specifying the layout. It must be
	either "counterclockwise" or "clockwise", and may be abbreviated.
	The default is "counterclockwise". Images corresponding to
	the x-y planes are always displayed in the third quadrant.
	If \code{layout} is counterclockwise, then the first quadrant
	shows images from the y-z planes and the second quadrant the x-z
	planes.	Otherwise, the images in the first and the second quadrant
	are switched. The fourth quadrant is left for the slider used to
	select the value of the fourth index (if any) of input array(s).}
}
\details{
  Shows slices of 3D array along the axes as produced by \code{image},
  along with sliders for controlling which slices are shown.  For 4D
  data an additional slider selects the value of the fourth index.
  Two images can be overlaid. This is useful for viewing medical imaging
  data (e.g. PET scans and fMRI data).
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
  slices3d(vol1=v, main="View of a mixture of three tri-variate normals", col1=heat.colors(256))

\dontrun{
  #Example 2: Put a z-map from fMRI data on top of a structure
  #           image. The threshold value of the z-map is 2.
  library(AnalyzeFMRI)
  temp<-f.read.analyze.volume("standard.img")
  z<-f.read.analyze.volume("z-map.img")
  slices3d(vol1=temp, vol2=z[,,,1], rlim2=c(2,Inf),col2=heat.colors(20),
          main="Regions above threshold values.")
  }
}
\keyword{hplot}
