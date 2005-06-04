parametric3d <- function(fx, fy, fz, umin, umax, vmin, vmax, n=100,
                         add = FALSE,...){
    origion.u <- seq(umin, umax, len=n)
    origion.v <- seq(vmin, vmax, len=n)

    if (add==FALSE) rgl.clear()
        
    off1 <- c(0, 0, 1, 1)
    off2 <- c(0, 1, 1, 0)
    ppm1 <- n - 1
    gridu <- rep(1 : ppm1, rep(4 * (n - 1), n - 1)) + off1
    gridv <- rep(rep(1 : (n - 1), n - 1), rep(4, (n - 1) ^ 2)) + off2
    u <- origion.u[gridu]
    v <- origion.v[gridv]
      
    x <- fx(u,v)
    y <- fy(u,v)
    z <- fz(u,v)

    nf <- findNonFinite3d(x, y, z, 4)
    if (any(nf)) {
      x <- x[! nf]
      y <- y[! nf]
      z <- z[! nf]
    }
    rgl.quads(x, z, -y, ...)
}
