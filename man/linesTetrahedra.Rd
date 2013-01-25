\name{linesTetrahedra}
\alias{linesTetrahedra}
\title{Create a Set of Lines with Tetrahetra Centered at Points along
  the Lines}
\description{
  Creates a scene consisting of lines made up of small tetrahedra
  centered at points along them.
}
\usage{
linesTetrahedra(x, y, z, delta=c(min(x[,2]-x[,1])/10,
                                    min(y[,2]-y[,1])/10,
                                    min(z[,2]-z[,1])/10),
                lwd = 0.01, color = "black", ...)
}
\arguments{
  \item{x, y, z}{numeric vectors of length two or matrices with two
    columns representing coordinates of starting and ending points 
    of line(s).}
  \item{delta}{numeric; increase in each dimension used to locate points
    along the lines; recycled to length 3.}
  \item{lwd}{numeric; used for the size of the tetrahedron in each
    dimension; recycled to length 3.}
  \item{color}{color to use for the tetrahedra.}
  \item{...}{additional arguments to be passed on to
    \code{makeTriangles}.}
}
\value{
  Returns a triangle mesh scene representing the lines.
}
\details{
  The function uses the Bresenham's line algorithm to locate points along
  lines and then creates a triangle mesh scene representing tetrahedra
  centered at those points.
}
\seealso{
  \code{\link[rgl]{lines3d}}.
}
\examples{
p <- pointsTetrahedra(x=c(100,100, 257, 257), 
                      y=c(100,100, 257, 257),
                      z=c(100,257, 257, 100), size=1)
l <- linesTetrahedra(x=matrix(c(100,257,
                                100,257), nrow=2, byrow=TRUE),
                     y=matrix(c(100,257,
                                100,257), nrow=2, byrow=TRUE),
                     z=matrix(c(100,257,
                                257,100), nrow=2, byrow=TRUE),
                     lwd=0.4,
                     col="red")
drawScene.rgl(list(p, l))
}
\keyword{hplot}
