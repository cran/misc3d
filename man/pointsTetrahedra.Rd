\name{pointsTetrahedra}
\alias{pointsTetrahedra}
\title{Create a Set of Tetrahetra Centered at Data Points}
\description{
  Creates a scene consisting of small tetrahedra centered at
  specified data points in three dimensions.
}
\usage{
pointsTetrahedra(x, y, z, size = 0.01, color = "black", ...)
}
\arguments{
  \item{x, y, z}{numeric vectors representing point coordinates.}
  \item{size}{numeric; multiple of data range to use for the size of the
    tetrahedron in each dimension; recycled to length 3.}
  \item{color}{color to use for the tetrahedra.}
  \item{...}{additional arguments to be passed on to
    \code{makeTriangles}.}
}
\value{
  Returns a triangle mesh scene representing the tetrahedra.
}
\details{
  This function is useful, for example, for incorporating raw data along
  with a density estimate surface in a scene rendered using standard or
  grid graphics.  For \pkg{rgl} rendering \code{\link[rgl]{points3d}} is
  an alternative.
}
\seealso{
  \code{\link[rgl]{points3d}}.
}
\examples{
  with(quakes, {
      d <- kde3d(long, lat, -depth, n = 40)
      v <- contour3d(d$d, exp(-12),d$x/22, d$y/28, d$z/640,
                     color="green", color2="gray", draw=FALSE)
      p <- pointsTetrahedra(long/22, lat/28, -depth/640,
                            size = 0.005)
      drawScene(list(v, p))
  })
}
\keyword{hplot}
