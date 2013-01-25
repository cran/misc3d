\name{computeContour3d}
\alias{computeContour3d}
\title{Compute Isosurface, a Three Dimension Contour}
\description{
  Computes a 3D contours or isosurface by the marching cubes algorithm.
}
\usage{
computeContour3d(vol, maxvol = max(vol), level,
                 x = 1:dim(vol)[1],
                 y = 1:dim(vol)[2],
                 z = 1:dim(vol)[3], mask)
}
\arguments{
  \item{vol}{a three dimensional array.}
  \item{maxvol}{maximum of the \code{vol} array.}
  \item{level}{The level at which to construct the contour surface.}
  \item{x,y,z}{locations of grid planes at which values in \code{vol} are
    measured.}
  \item{mask}{a function of 3 arguments returning a logical array, a
    three dimensional logical array, or \code{NULL}.  If not
    \code{NULL}, only cells for which \code{mask} is true at all eight
    vertices are used in forming the contour.}
}
\value{
  A matrix of three columns representing the triangles making up the
  contour surface. Each row represents a vertex and goups of three
  rows represent a triangle.
  }
\details{
  Uses the marching-cubes algorithm, with adjustments for dealing with
  face and internal ambiguities, to compute an  isosurface.
  See references for the details. The function
  \code{\link[misc3d]{contour3d}} provides a higher-level interface.
}
\references{
 Chernyaev E. (1995)  
 Marching Cubes 33: Construction of Topologically Correct Isosurfaces  
 \emph{Technical Report CN/95-17, CERN}

 Lorensen W. and Cline H. (1987)  
 Marching Cubes: A High Resolution 3D Surface Reconstruction Algorithm  
 \emph{Computer Graphics} \bold{vol. 21, no. 4}, 163-169
 
 Nielson G. and Hamann B. (1992)
 The Asymptotic Decider: Resolving the Ambiguity in Marching Cubes  
 \emph{Proc. IEEE Visualization}  \bold{92}, 83-91
} 
\seealso{
  \code{\link[misc3d]{contour3d}}
}
\examples{
  x <- seq(-2,2,len=50)
  g <- expand.grid(x = x, y = x, z = x)
  v <- array(g$x^4 + g$y^4 + g$z^4, rep(length(x),3))
  con <- computeContour3d(v, max(v), 1)
  drawScene(makeTriangles(con))
}
\keyword{hplot}



