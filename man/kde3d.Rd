\name{kde3d}
\alias{kde3d}
\title{Compute a Three Dimension Kernel Density Estimate}
\description{
  Evaluates a three dimensional kernel density estimate using a Gaussian
  kernel with diagonal covariance matrix on a regular grid.
}
\usage{
   kde3d(x, y, z, h, n = 20, lims = c(range(x), range(y), range(z)))
}
\arguments{
  \item{x,y,z}{\code{x}, \code{y}, and \code{z} coordinates of the data.}
  \item{h}{vector of three bandwidths for the density estimate;
    recycled if length is less than three; default is based on the
    normal reference bandwidth (see \code{\link[MASS]{bandwidth.nrd}}).}
  \item{n}{numbers of grid points to use for each dimension; recycled if
    length is less than three.}
  \item{lims}{lower and upper limits on the region for which the density
    estimate is to be computed, provides as a vector of length 6,
    corresponding to low and high values of \code{x}, \code{y}, and
    \code{z}; recycled if only two values are supplied.}
}
\value{A list of four components, \code{x}, \code{y}, \code{z}, and
  \code{d}. \code{x}, \code{y}, and \code{z} are the coordinates of the
  grid points at which the density estimate has been evaluated, and
  \code{d} is a three dimensional array of the estimated density values.
}
\references{Based on the function \code{\link[MASS]{kde2d}} in package
  \pkg{MASS}.
} 
\seealso{
  \code{\link[MASS]{kde2d}}.
}
\examples{
  with(quakes, {
      d <- kde3d(long, lat, -depth, n = 40)
      contour3d(d$d, exp(-12), d$x/22, d$y/28, d$z/640,
                color = "green", color2 = "gray", scale=FALSE,
                engine = "standard")
  })
}
\keyword{dplot}



