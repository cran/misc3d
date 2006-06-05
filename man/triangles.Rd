\name{triangles}
\alias{makeTriangles}
\alias{updateTriangles}
\alias{scaleTriangles}
\alias{translateTriangles}
\title{Triangle Mesh Functions}
\description{
  Functions to create and modify triangle mesh objects representing 3D
  surfaces..
}
\usage{
makeTriangles(v1, v2, v3, color = "red", color2 = NA, alpha = 1,
              fill = TRUE, col.mesh = if (fill) NA else color,
              smooth = 0, material = "default")
updateTriangles(triangles, color, color2, alpha, fill, col.mesh,
                material, smooth)
translateTriangles(triangles, x = 0, y = 0, z = 0)
scaleTriangles(triangles, x = 1, y = x, z = x)
}
\arguments{
  \item{v1,v2,v3}{specification of triangle coordinates. If all three
    are provided then they should be matrices with three columns
    representing coordinates of the first, second, and third vertices of
    the triangles.  If only \code{v1} and \code{v2} are provided then
    \code{v1} should be a numeric matrix with three rows specifying
    coordinates of vertices, and \code{v2} should be an integer matrix
    with three rows specifying the indexes if the vertices in the
    triangles.  If only \code{v1} is provided then it should be a matrix
    with three columns and number of rows divisible by three specifying
    the vertices of the triangles in groups of three.}
  \item{triangles}{triangle mesh object.}
  \item{x,y,z}{numeric of length one. Amounts by which to translate or
    scale corresponding coordinates.}
  \item{color}{color to use for the surface. Can also be a function of
    three arguments.  This is called with three arguments, the
    coordinates of the midpoints of the triangles making up the surface.
    The function should return a vector of colors to use for the
    triangles.}
  \item{color2}{opposite face color; only used for "standard" and "grid"
    engines.}
  \item{alpha}{alpha channel level, a number between 0 and 1.}
  \item{fill}{logical; if \code{TRUE}, drawing should use filled
    surfaces; otherwise a wire frame should be drawn.}
  \item{col.mesh}{color to use for the wire frame.}
  \item{smooth}{integer specifying Phong shading level; currently only
    used by "standard" and "grid" engines.}
  \item{material}{material specification; currently only used by
    "standard" and "grid" engines. Currently possible values are the
    character strings "dull", "shiny", "metal", and "default".}
}
\value{
  A triangle mesh object of class \code{Triangles3D}.
}
\details{
  \code{makeTriangles} creates a triangle mesh object.
  \code{updateTriangles} modifies fields of such an object.  Both may
  perform some consistency checks.

  \code{translateTriangles} and \code{scaleTriangles} translate or scale
  the vertex locations of triangle mesh objects by specified amounts.
}
}
\keyword{hplot}



