\name{exportScene}
\alias{exportScene}
\title{Writing Out Triangular Mesh Scenes}
\description{
  Writing out scenes consisting of one or more surfaces represented by triangular
  mesh data structures to textual files.
}
\usage{
exportScene(scene, filename, format=c("OFF", "IDTF", "ASY"))
}
\arguments{
  \item{scene}{a triangle mesh object of class \code{Triangles3D} or a
    list of such objects representing the scene to be exported.}
  \item{filename}{the name of the exported textual file.}
  \item{format}{the format of the exported textual file.
    It must be one of "OFF", "IDTF", or "ASY" and can be abbreviated.
    The default is "OFF".}
}
\value{
  Textual files representing triangular mesh scenes.
}
\details{
  \code{exportScene} writes out scenes to textual files,
  which can be used for other purposes, for example the
  generation of U3d and PRC files for interactive 3D visualization in
  a PDF.
}
\examples{
  nmix3 <- function(x, y, z, m, s) {
      0.4 * dnorm(x, m, s) * dnorm(y, m, s) * dnorm(z, m, s) +
      0.3 * dnorm(x, -m, s) * dnorm(y, -m, s) * dnorm(z, -m, s) +
      0.3 * dnorm(x, m, s) * dnorm(y, -1.5 * m, s) * dnorm(z, m, s)
  }
 
  f <- function(x,y,z) nmix3(x,y,z,.5,.5)
 
  gs1 <- function(n = 40, k = 5, cmap = heat.colors, ...) {
      th <- seq(0.05, 0.2, len = k)
      col <- rev(cmap(length(th)))
      x <- seq(-2, 2, len=n)
      m <- function(x,y,z) x > .25 | y < -.3
      contour3d(f,th,x,x,x,color=col, mask = m, engine = "none",
               scale = FALSE, ...)
  }
 
  conts <- gs1(40, 5, screen=list(z = 130, x = -80),
               color2 = "lightgray", cmap=rainbow)
  filename <- file.path(tempdir(), "nmix")
  exportScene(conts, filename, "OFF")
}
\keyword{hplot}



