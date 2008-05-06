\name{parametric3d}
\alias{parametric3d}
\title{Draw a 3D Parametric Plot}
\description{
  Plot a two-parameter surface in three dimensions.
}
\usage{
  parametric3d(fx, fy, fz, u, v, umin, umax, vmin, vmax, n = 100,
               color = "white", color2 = NA, alpha = 1,
               fill = TRUE, col.mesh = if (fill) NA else color,
               smooth = 0, material = "default", 
               add = FALSE, draw = TRUE, engine = "rgl", ...)
}
\arguments{
  \item{fx,fy,fz}{vectorized functions of u and v to compute the
    \code{x}, \code{y}, and \code{z} coordinates.}
  \item{u}{numeric vector of u values.}
  \item{v}{numeric vector of v values.}
  \item{umin}{numeric; the minimum value of u. Ignored in \code{u} is
    supplied.}
  \item{umax}{numeric; the maximum value of u. Ignored in \code{u} is
    supplied.}
  \item{vmin}{numeric; the minimum value of v. Ignored in \code{v} is
    supplied.}
  \item{vmax}{numeric; the maximum value of v. Ignored in \code{v} is
    supplied.}
  \item{n}{the number of equally spaced \code{u} and \code{v} values to
    use.  Ignored in \code{u} and \code{v} are supplied.}
  \item{color}{color to use for the surface. Can also be a function of
    three arguments.  This is called with three arguments, the
    coordinates of the midpoints of the triangles making up the surface.
    The function should return a vector of colors to use for the
    triangles.}
  \item{color2}{opposite face color.}
  \item{alpha}{alpha channel level, a number between 0 and 1..}
  \item{fill}{logical; if \code{TRUE}, drawing should use filled
    surfaces; otherwise a wire frame should be drawn.}
  \item{col.mesh}{color to use for the wire frame.}
  \item{smooth}{integer specifying Phong shading level; currently only
    used by "standard" and "grid" engines.}
  \item{material}{material specification; currently only used by
    "standard" and "grid" engines.  Currently possible values are the
    character strings "dull", "shiny", "metal", and "default".}
  \item{add}{logical; if \code{TRUE}, add to current graph.}
  \item{draw}{logical; if \code{TRUE}, draw the results; otherwise,
    return contour triangles.}
  \item{engine}{character; currently "rgl", "standard", "grid" or "none";
    for "none" the computed triangles are returned.}
  \item{...}{additional rendering arguments, e.g. material and texture
    properties for the "rgl" engine.   See documentation for
    \code{\link[misc3d]{drawScene}} and \code{\link[misc3d]{drawScene.rgl}}}
}
\value{
  For the "rgl" engine the returned value is \code{NULL}. For the
  "standard" and "grid" engines the returned value is the viewing
  transformation as returned by \code{persp}.  For the engine "none", or
  when \code{draw} is not true, the returned value is a structure
  representing the triangles making up the surface.
  }
\details{
  Analogous to Mathematica's \code{Param3D}.  Evaluates the
  functions \code{fx}, \code{fy}, and \code{fz} specifying the
  coordinates of the surface at an equally spaced grid of values for the
  parameters \code{u} and \code{v}.
}
\note{
  The "rgl" engine now uses the standard rgl coordinates instead of
  negating \code{y} and swapping \code{y} and \code{z}. If you need to
  reproduce the previous behavior you can use
  \code{options(old.misc3d.orientation=TRUE)}.

  Transparency only works properly in the "rgl" engine.  For standard or
  grid graphics on pdf or quartz devices using alpha levels less than 1
  does work but the triangle borders show as a less transparent mesh.
  }
\references{ 
   Daniel Adler, Oleg Nenadic and Walter Zucchini (2003)
   RGL: A R-library for 3D visualization with OpenGL
}
\seealso{
  \code{\link[rgl]{rgl.surface}},
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
  parametric3d(fx = function(u, v) u * exp(-0.5 * (u^2 + v^2 -
                      2 * 0.75 * u * v)/sqrt(1-.75^2))^(1/3),
               fy = function(u, v) v * exp(-0.5 * (u^2 + v^2 -
                      2 * 0.75 * u * v)/sqrt(1-.75^2))^(1/3),
               fz = function(u, v) exp(-0.5 * (u^2 + v^2 - 2 * 0.75 * u * 
                      v)/sqrt(1-.75^2))^(1/3),
               u = qcauchy((1:100)/101), v = qcauchy((1:100)/101))	
  parametric3d(fx = function(u, v) u * exp(-0.5 * (u^2 + v^2 -
                      2 * 0.75 * u * v)/sqrt(1-.75^2))^(1/3),
               fy = function(u, v) v * exp(-0.5 * (u^2 + v^2 -
                      2 * 0.75 * u * v)/sqrt(1-.75^2))^(1/3),
               fz = function(u, v) exp(-0.5 * (u^2 + v^2 - 2 * 0.75 * u * 
                      v)/sqrt(1-.75^2))^(1/3),
               u = qcauchy((1:100)/101), v = qcauchy((1:100)/101),
               engine = "standard", scale = FALSE, screen = list(x=-90, y=20))

  #Example 2: Ratio-of-Uniform sampling region of Bivariate t      
  parametric3d(fx = function(u,v) u*(dt(u,2) * dt(v,2))^(1/3), 
               fy = function(u,v) v*(dt(u,2) * dt(v,2))^(1/3),
               fz = function(u,v) (dt(u,2) * dt(v,2))^(1/3), 
               umin = -20, umax = 20, vmin = -20, vmax = 20, 
               n = 100, color = "green")
  parametric3d(fx = function(u,v) u*(dt(u,2) * dt(v,2))^(1/3),
               fy = function(u,v) v*(dt(u,2) * dt(v,2))^(1/3),
               fz = function(u,v) (dt(u,2) * dt(v,2))^(1/3),
               u = qcauchy((1:100)/101), v = qcauchy((1:100)/101),
               color = "green")
  parametric3d(fx = function(u,v) u*(dt(u,2) * dt(v,2))^(1/3),
               fy = function(u,v) v*(dt(u,2) * dt(v,2))^(1/3),
               fz = function(u,v) (dt(u,2) * dt(v,2))^(1/3),
               u = qcauchy((1:100)/101), v = qcauchy((1:100)/101),
               color = "green", engine = "standard", scale = FALSE)


  #Example 3: Surface of revolution
  parametric3d(fx = function(u,v) u,
               fy = function(u,v) sin(v)*(u^3+2*u^2-2*u+2)/5,
               fz = function(u,v) cos(v)*(u^3+2*u^2-2*u+2)/5,
               umin = -2.3, umax = 1.3, vmin = 0, vmax = 2*pi)
  parametric3d(fx = function(u,v) u,
               fy = function(u,v) sin(v)*(u^3+2*u^2-2*u+2)/5,
               fz = function(u,v) cos(v)*(u^3+2*u^2-2*u+2)/5,
               umin = -2.3, umax = 1.3, vmin = 0, vmax = 2*pi,
               engine = "standard", scale = FALSE,
               color = "red", color2 = "blue", material = "shiny")

}
\keyword{hplot}
