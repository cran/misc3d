\name{contour3d}
\alias{contour3d}
\title{Draw an Isosurface, a Three Dimension Contour Plot}
\description{
  Uses rgl to render an isosurface computed by the marching cubes algorithm.
}
\usage{
   contour3d(f, level, x, y, z, add=FALSE,...)
}
\arguments{
  \item{f}{a function of 3 arguments or a three dimensional array.}
  \item{level}{numeric scalar. The value according to which to draw the
    isosurface.}
  \item{x,y,z}{locations of grid planes at which values in \code{f} are
    measured or \code{f} is evaluated. Can be omitted if \code{f} is an
    array.}
  \item{add}{logical; if \code{TRUE}, add to current \code{rgl} graph.}
  \item{...}{material and texture properties. See \code{rgl.material}
    for details.}
}
\details{
  Uses the marching-cubes algorithm, with adjustments for dealing with
  face and internal ambiguities, to draw the isosurface.
  See references for the details.
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
  \code{\link[rgl]{rgl.triangles}}, \code{\link[rgl]{rgl.material}},
  \code{\link[rgl]{rgl.surface}}.
}
\examples{
  #Example 1: Draw a ball
  f <- function(x, y, z)x^2+y^2+z^2
  x <- seq(-2,2,len=20)
  contour3d(f,4,x,x,x)

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
      for (i in seq(along = th)) {
           contour3d(f,th[i],x,x,x,color=col[i],alpha=al[i],
                     add = (i != 1))
      }
      rgl.bg(col="white")
  }
  g(40,5)

\dontrun{
  #Example 3: Nested contours for FMRI data.
  library(AnalyzeFMRI)
  a<-f.read.analyze.volume(system.file("example.img", package="AnalyzeFMRI"))
  contour3d(aperm(a[,,,1],c(1,3,2)),1:64,1.5*(1:21),1:64,lev=3000,alpha=0.5)
  contour3d(aperm(a[,,,1],c(1,3,2)),1:64,1.5*(1:21),1:64,lev=8000,
            col="red",add=TRUE, alpha=0.5)
  contour3d(aperm(a[,,,1],c(1,3,2)),1:64,1.5*(1:21),1:64,lev=10000,
            col="green",add=TRUE)
  }
}
\keyword{hplot}



