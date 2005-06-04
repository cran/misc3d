\name{parametric3d}
\alias{parametric3d}
\title{Draw a 3D Parametric Plot}
\description{
  Plot a two-parameter surface in three dimensions.
}
\usage{
  parametric3d(fx, fy, fz, umin, umax, vmin, vmax, n = 100, add = FALSE, ...)
}
\arguments{
  \item{fx,fy,fz}{vectorized functions of u and v to compute the
    \code{x}, \code{y}, and \code{z} coordinates.}
  \item{umin}{numeric; the minimum value of u.}
  \item{umax}{numeric; the maximum value of u.}
  \item{vmin}{numeric; the minimum value of v.}
  \item{vmax}{numeric; the maximum value of v.}
  \item{n}{the number of equally spaced \code{u} and \code{v} values to use.}
  \item{add}{logical; if \code{TRUE}, add to current \code{rgl} graph.}
  \item{...}{material and texture properties. See \code{rgl.material}
    for details.}
}
\details{
  Analogous to Mathematica's \code{Param3D}.  Evaluates the
  functions \code{fx}, \code{fy}, and \code{fz} specifying the
  coordinates of the surface at an equally spaced grid of values for the
  parameters \code{u} and \code{v}.
}
\references{ 
   Daniel Adler, Oleg Nenadic and Walter Zucchini (2003)
   RGL: A R-library for 3D visualization with OpenGL
}
\seealso{
  \code{\link[rgl]{rgl.quads}},\code{\link[rgl]{rgl.surface}},
  \code{\link[rgl]{rgl.material}},\code{\link[scatterplot3d]{scatterplot3d}}.
}
\examples{
  #Example 1: Ratio-of-Uniform sampling region of bivariate normal
  parametric3d(fx = function(u, v) u * exp(-0.5 * (u^2 + v^2 -
                      2 * 0.75 * u * v)/sqrt(1-.75^2))^(1/3),
               fy = function(u, v) v * exp(-0.5 * (u^2 + v^2 -
                      2 * 0.75 * u * v)/sqrt(1-.75^2))^(1/3),
               fz = function(u, v) exp(-0.5 * (u^2 + v^2 - 2 * 0.75 * u * 
                      v)/sqrt(1-.75^2))^(1/3),
               umin = -20, umax = 20, vmin = -20, vmax = 20, 
               n = 100)	

  #Example 2: Ratio-of-Uniform sampling region of Bivariate t      
  parametric3d(fx = function(u,v) u*(dt(u,2) * dt(v,2))^(1/3), 
               fy = function(u,v) v*(dt(u,2) * dt(v,2))^(1/3),
               fz = function(u,v) (dt(u,2) * dt(v,2))^(1/3), 
               umin = -20, umax = 20, vmin = -20, vmax = 20, 
               n = 100, color = "green", alpha=1.0)

  #Example 3: Surface of revolution
  parametric3d(fx = function(u,v) u,
               fy = function(u,v) sin(v)*(u^3+2*u^2-2*u+2)/5,
               fz = function(u,v) cos(v)*(u^3+2*u^2-2*u+2)/5,
               umin = -2.3, umax = 1.3, vmin = 0, vmax = 2*pi)

}
}
\keyword{hplot}
