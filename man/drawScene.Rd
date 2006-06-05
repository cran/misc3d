\name{drawScene}
\alias{drawScene}
\alias{drawScene.rgl}
\title{Rendering of Triangular Mesh Surface Data}
\description{
  Draw scenes consisting of one or more surfaces described by triangular
  mesh data structures.
}
\usage{
drawScene(scene, light = c(0, 0, 1),
          screen = list(z = 40, x = -60), scale = TRUE, R.mat = diag(4),
          perspective = FALSE, distance = if (perspective) 0.2 else 0, 
          fill = TRUE, xlim = NULL, ylim = NULL, zlim = NULL,
          aspect = c(1, 1), col.mesh = if (fill) NA else "black",
          polynum = 100, lighting = phongLighting, add = FALSE,
          engine = "standard", col.bg = "transparent", depth = 0)
drawScene.rgl(scene, add = FALSE, ...)
}
\arguments{
  \item{scene}{a triangle mesh object of class \code{Triangles3D} or a
    list of such objects representing the scene to be rendered.}
  \item{light}{numeric vector of length 3 or 4.  The first three
    elements represent the direction to the light in viewer coordinates;
    the viewer is at \code{(0, 0, 1 / distance)} looking down along the
    positive z-axis.  The fourth element, if present, represents light
    intensity; the default is 1.}
  \item{screen}{as for \code{panel.3dwire}, a list giving sequence of
    rotations to be applied to the scene before being rendered. The
    initial position starts with the viewing point along the positive
    z-axis, and the x and y axes in the usual position. Each component
    of the list should be named one of "x", "y" or "z"; repetitions are
    allowed.  The values values indicate the amount of rotation about
    that axis in degrees.}
  \item{scale}{logical. Before viewing the x, y and z coordinates of the
    scene defining the surface are transformed to the interval
    [-0.5,0.5].  If \code{scale} is true the x, y and z coordinates are
    transformed separately.  Otherwise, the coordinates are scaled so
    that aspect ratios are retained.}
  \item{R.mat}{initial rotation matrix in homogeneous coordinates, to be
    applied to the data before \code{screen} rotates the view further.}
  \item{perspective}{logical, whether to render a perspective
    view. Setting this to \code{FALSE} is equivalent to setting
    \code{distance} to 0}
  \item{distance}{numeric, between 0 and 1, controls amount of
    perspective. The distance of the viewing point from the origin (in
    the transformed coordinate system) is \code{1 / distance}. This is
    described in a little more detail in the documentation for
    \code{cloud}.}
  \item{fill}{logical; if \code{TRUE}, drawing should use filled
    surfaces or wire frames as indicated by the object properties.
    Otherwise all objects in the scene should be rendered as wire
    frames.}
  \item{xlim,ylim,zlim}{x-, y-  and z-limits.  The scene is rendered so that
    the rectangular volume defined by these limits is visible.}
  \item{aspect}{vector of length 2. Gives the relative aspects of the
    y-size/x-size and z-size/x-size of the enclosing cube.}
  \item{col.mesh}{color to use for the wire frame if \code{frames} is
    true.}
  \item{polynum}{integer.  Number of triangles to pass in batches to
    grid primitives for the "grid" engine.  The default should be
    adequate.}
  \item{lighting}{a lighting function.  Current options are
    \code{phongLighting} and \code{perspLighting}.}
  \item{add}{logical; if \code{TRUE}, add to current \code{rgl} graph.}
  \item{engine}{character; currently "standard" or "grid".}
  \item{col.bg}{background dolor to use in color depth cuing.}
  \item{depth}{numeric, between 0 and 1.  Controls the amount of color
    blending to \code{col.bg} for objects farther from the
    viewer. \code{depth} equal to zero means no depth cuing.}
  \item{...}{rgl material and texture properties; see documentation for
    \code{\link[rgl]{rgl.material}}}
}
\value{
  \code{drawScene.rgl} returns \code{NULL}. The return value of
  \code{drawScene} is the viewing transformation as returned by
  \code{persp}.
}
\details{
  \code{drawScene} renders a scene consisting of one or more triangle
  mesh objects using standard or grid graphics.  Object-specific
  rendering features such as smoothing and material are controlled by
  setting in the objects.  Arguments to \code{drawScene} control global
  factors such as viewer and light position.

  \code{drawScene.rgl} renders the scene in an rgl window.
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
\seealso{\code{\link[rgl]{rgl.material}}}
\examples{
  vtri <- local({
      z <- 2 * volcano
      x <- 10 * (1:nrow(z))
      y <- 10 * (1:ncol(z))
      surfaceTriangles(x, y, z, color="green3")
  })
  drawScene(vtri, scale = FALSE)
  drawScene(vtri, screen=list(x=40, y=-40, z=-135), scale = FALSE)
  drawScene(vtri, screen=list(x=40, y=-40, z=-135), scale = FALSE,
            perspective = TRUE)
  drawScene(vtri, screen=list(x=40, y=-40, z=-135), scale = FALSE,
            perspective = TRUE, depth = 0.4)
}
\keyword{hplot}



