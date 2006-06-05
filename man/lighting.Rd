\name{lighting}
\alias{perspLighting}
\alias{phongLighting}
\title{Lighting Functions}
\description{
  Functions to compute colors modified for lighting effects.
}
\usage{
phongLighting(normals, view, light, color, color2, alpha, material = "default")
perspLighting(normals, view, light, color, color2, alpha, material = "default")
}
\arguments{
  \item{normals}{numeric matrix with three columns representing surface
    normal vectors.}
  \item{view}{numeric vector of length 3 representing the direction to
    the viewer.}
  \item{light}{numeric vector of length 3 or 4.  The first three
    elements represent the direction to the light.  The fourth element,
    if present, represents light intensity; the default is 1.}
  \item{color}{colors to use for faces in the direction of the normal
    vectors.}
  \item{color2}{opposite face color.}
  \item{alpha}{alpha channel level, a number between 0 and 1.}
  \item{material}{material specification. Currently possible values are
    the character strings "dull", "shiny", "metal", and "default".}  }
\value{
  Vector of color specifications.
}
\details{
  \code{phongLighting} uses the Phong lighting model to compute colors
  modified for view direction, light direction, and material properties.
  \code{perspLighting} implements approximately the same lighting model
  as the \code{persp} function.
}
}
\keyword{hplot}



