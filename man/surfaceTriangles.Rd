\name{surfaceTriangles}
\alias{surfaceTriangles}
\title{Create a Triangle Mesh Representing a Surface}
\description{
  Creates a triangle mesh object representing a surface over a
  rectangular grid.
}
\usage{
surfaceTriangles(x, y, f, color = "red", color2 = NA, alpha = 1,
                 fill = TRUE, col.mesh = if (fill) NA else color,
                 smooth = 0, material = "default")
}
\arguments{
  \item{x, y}{numeric vectors.}
  \item{f}{numeric matrix of dimension \code{length(x)} by
    \code{length(y)} or vectorized function of two arguments.}
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
}
\value{
  Returns a triangle mesh object representing the surface.
  }
\seealso{
  \code{\link{persp}}, \code{\link[rgl]{rgl.surface}},
  \code{\link[rgl]{surface3d}}.
}
\examples{
  drawScene(surfaceTriangles(seq(-1,1,len=30), seq(-1,1,len=30),
                             function(x, y) (x^2 + y^2), color2 = "green"))
  drawScene.rgl(surfaceTriangles(seq(-1,1,len=30), seq(-1,1,len=30),
                             function(x, y) (x^2 + y^2), color2 = "green"))
}
}
\keyword{hplot}
