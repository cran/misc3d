findNonFinite3d <- function(x, y, z, group) {
  n <- length(x) / group
  nf <- ! (is.finite(x) & is.finite(y) & is.finite(z))
  if (any(nf)) {
    idxg <- group * (1 : n)
    if (group == 2)
      rep(nf[idxg] | nf[idxg - 1], rep(2, n))
    else if (group == 3)
      rep(nf[idxg] | nf[idxg - 1] | nf[idxg - 2], rep(3, n))
    else if (group == 4)
      rep(nf[idxg] | nf[idxg - 1] | nf[idxg - 2] | nf[idxg - 3], rep(4, n))
    else stop("unknown group size for NA grouping")
  }
  else nf
}

lines3d <- function(x,y,z,add=FALSE,...){
    n <- length(x)
    if(n != length(y) || n != length(z))
      stop("Lengths of x,  y, z do not match")

    if(!add) rgl.clear()

    off <- c(1,1)
    x <- kronecker(x, off)[-c(1, 2 * n)]
    y <- kronecker(y, off)[-c(1, 2 * n)]
    z <- kronecker(z, off)[-c(1, 2 * n)]

    nf <- findNonFinite3d(x, y, z, 2)
    if (any(nf)) {
      x <- x[! nf]
      y <- y[! nf]
      z <- z[! nf]
    }
    rgl.lines(x, z, -y,...)
}
