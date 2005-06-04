\name{lines3d}
\alias{lines3d}
\title{Draw 3D Lines}
\description{
  Draw connected line segments in three dimensions. 
}
\usage{
  lines3d(x,y,z,add=FALSE,...)
}
\arguments{
  \item{x,y,z}{coordinates of points to connect.}
  \item{add}{logical; if \code{TRUE}, add to current \code{rgl} graph.}
  \item{...}{material and texture properties. See \code{rgl.material}
    for details.}
}
\references{ 
   Daniel Adler, Oleg Nenadic and Walter Zucchini (2003)
   RGL: A R-library for 3D visualization with OpenGL
}
\seealso{
  \code{\link[rgl]{rgl.lines}}.
}
\examples{
  #Plot a three-dimensional helix.
  t = seq(0,10*pi,by=pi/50)
  lines3d(4*sin(t),4*cos(t),t,color="red")
}
\keyword{hplot}
