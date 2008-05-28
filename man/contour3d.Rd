\name{contour3d}
\alias{contour3d}
\title{Draw an Isosurface, a Three Dimension Contour Plot}
\description{
  Computes and renders 3D contours or isosurfaces computed by the
  marching cubes algorithm.
}
\usage{
   contour3d(f, level, x, y, z, mask = NULL, color = "white", color2 = NA,
             alpha = 1, fill = TRUE, col.mesh = if (fill) NA else color,
             material = "default", smooth = 0,
             add = FALSE, draw = TRUE, engine = "rgl", ...)
}
\arguments{
  \item{f}{a function of 3 arguments or a three dimensional array.}
  \item{level}{The level or levels at which to construct contour surfaces.}
  \item{x,y,z}{locations of grid planes at which values in \code{f} are
    measured or \code{f} is to be evaluated. Can be omitted if \code{f}
    is an array.}
  \item{mask}{a function of 3 arguments returning a logical array, a
    three dimensional logical array, or \code{NULL}.  If not
    \code{NULL}, only cells for which \code{mask} is true at all eight
    vertices are used in forming the contour.  Can also be a list of
    functions the same length as \code{level}.}
  \item{color}{color to use for the contour surface. Recycled to the
    length of \code{'levels'}. Can also be a function, or list of
    functions, of three arguments.  These are called for each level with
    three arguments, the coordinates of the midpoints of the triangles
    making up the surface.  They should return a vector of colors to use
    for the triangles.}
  \item{color2}{opposite face color. Recycled to the length of
    \code{'levels'}.}
  \item{alpha}{alpha channel level, a number between 0 and 1. Recycled to the
    length of \code{'levels'}.}
  \item{fill}{logical; if \code{TRUE}, drawing should use filled
    surfaces; otherwise a wire frame should be drawn. Recycled to the
    length of \code{'levels'}.}
  \item{col.mesh}{color to use for the wire frame. Recycled to the
    length of \code{'levels'}.}
  \item{smooth}{integer or logical specifying Phong shading level for
    "standard" and "grid" engines or whether or not to use shading for
    the "rgl" engine. Recycled to the length of \code{'levels'}.}
  \item{material}{material specification; currently only used by
    "standard" and "grid" engines. Currently possible values are the
    character strings "dull", "shiny", "metal", and "default".  Recycled
    to the length of \code{'levels'}.}
  \item{add}{logical; if \code{TRUE}, add to current \code{rgl} graph.}
  \item{draw}{logical; if \code{TRUE}, draw the results; otherwise,
    return contour triangles.}
  \item{engine}{character; currently "rgl", "standard", "grid" or "none";
    for "none" the computed triangles are returned.}
  \item{...}{additional rendering arguments, e.g. material and texture
    properties for the "rgl" engine.  See documentation for
    \code{\link[misc3d]{drawScene}} and \code{\link[misc3d]{drawScene.rgl}}}
}
\value{
  For the "rgl" engine the returned value is \code{NULL}. For the
  "standard" and "grid" engines the returned value is the viewing
  transformation as returned by \code{persp}.  For the engine "none", or
  when \code{draw} is not true, the returned value is a structure
  representing the triangles making up the contour, or a list of such
  structures for multiple contours.
  }
\details{
  Uses the marching-cubes algorithm, with adjustments for dealing with
  face and internal ambiguities, to draw isosurfaces.
  See references for the details.
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
 Chernyaev E. (1995)  
 Marching Cubes 33: Construction of Topologically Correct Isosurfaces  
 \emph{Technical Report CN/95-17, CERN}

 Daniel Adler, Oleg Nenadic and Walter Zucchini (2003)
 RGL: A R-library for 3D visualization with OpenGL

 Lorensen W. and Cline H. (1987)  
 Marching Cubes: A High Resolution 3D Surface Reconstruction Algorithm  
 \emph{Computer Graphics} \bold{vol. 21, no. 4}, 163-169
 
 Nielson G. and Hamann B. (1992)
 The Asymptotic Decider: Resolving the Ambiguity in Marching Cubes  
 \emph{Proc. IEEE Visualization}  \bold{92}, 83-91
} 
\seealso{
  \code{\link[rgl]{triangles3d}}, \code{\link[rgl]{material3d}},
  \code{\link[rgl]{surface3d}}.
}
\examples{
  #Example 1: Draw a ball
  f <- function(x, y, z)x^2+y^2+z^2
  x <- seq(-2,2,len=20)
  contour3d(f,4,x,x,x)
  contour3d(f,4,x,x,x, engine = "standard")

  # ball with one corner removed.
  contour3d(f,4,x,x,x, mask = function(x,y,z) x > 0 | y > 0 | z > 0)
  contour3d(f,4,x,x,x, mask = function(x,y,z) x > 0 | y > 0 | z > 0,
            engine="standard", screen = list(x = 290, y = -20),
            color = "red", color2 = "white")

  # ball with computed colors
  w <- function(x,y,z) {
      v <- sin(x) + cos(2 * y) * sin(5 * z)
      r <- range(v)
      n <- 100
      i <- pmax(pmin(ceiling(n * (v - r[1]) / (r[2] - r[1])), n), 1)
      terrain.colors(n)[i]
  }
  contour3d(f,4,x,x,x, color = w)

  #Example 2: Nested contours of mixture of three tri-variate normal densities
  nmix3 <- function(x, y, z, m, s) {
      0.4 * dnorm(x, m, s) * dnorm(y, m, s) * dnorm(z, m, s) +
      0.3 * dnorm(x, -m, s) * dnorm(y, -m, s) * dnorm(z, -m, s) +
      0.3 * dnorm(x, m, s) * dnorm(y, -1.5 * m, s) * dnorm(z, m, s)
  }
  f <- function(x,y,z) nmix3(x,y,z,.5,.5)
  g <- function(n = 40, k = 5, alo = 0.1, ahi = 0.5, cmap = heat.colors) {
      th <- seq(0.05, 0.2, len = k)
      col <- rev(cmap(length(th)))
      al <- seq(alo, ahi, len = length(th))
      x <- seq(-2, 2, len=n)
      contour3d(f,th,x,x,x,color=col,alpha=al)
      bg3d(col="white")
  }
  g(40,5)
  gs <- function(n = 40, k = 5, cmap = heat.colors, ...) {
      th <- seq(0.05, 0.2, len = k)
      col <- rev(cmap(length(th)))
      x <- seq(-2, 2, len=n)
      m <- function(x,y,z) x > .25 | y < -.3
      contour3d(f,th,x,x,x,color=col, mask = m, engine = "standard",
                scale = FALSE, ...)
      bg3d(col="white")
  }
  gs(40, 5, screen=list(z = 130, x = -80), color2 = "lightgray", cmap=rainbow)

\dontrun{
  #Example 3: Nested contours for FMRI data.
  library(AnalyzeFMRI)
  a <- f.read.analyze.volume(system.file("example.img", package="AnalyzeFMRI"))
  a <- a[,,,1]
  contour3d(a, 1:64, 1:64, 1.5*(1:21), lev=c(3000, 8000, 10000),
            alpha = c(0.2, 0.5, 1), color = c("white", "red", "green"))

  # alternative masking out a corner
  m <- array(TRUE, dim(a))
  m[1:30,1:30,1:10] <- FALSE
  contour3d(a, 1:64, 1:64, 1.5*(1:21), lev=c(3000, 8000, 10000),
            mask = m, color = c("white", "red", "green"))
  contour3d(a, 1:64, 1:64, 1.5*(1:21), lev=c(3000, 8000, 10000),
            color = c("white", "red", "green"),
            color2 = c("gray", "red", "green"),
            mask = m, engine="standard",
            scale = FALSE, screen=list(z = 60, x = -120))
  }
}
\keyword{hplot}



