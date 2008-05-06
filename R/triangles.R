## These functions work with collections of n triangles.  A collection of
## triangles is a list with components v1, v2, v3 representing the
## coordinates of the three vertices; each of these components is an n by
## 3 matrix.

makeTriangles <- function(v1, v2, v3, 
                          color = "red", color2 = NA, alpha = 1,
                          fill = TRUE, col.mesh = if (fill) NA else color,
	                  smooth = 0,  material = "default") {
    if (missing(v2) || missing(v3)) {
        if (missing(v2) && missing(v3))
  	    v <- unzipTriangleMatrix(v1)
        else if (missing(v3))
            v <- ve2t(list(vb = v1, ib = v2))
        else stop("unknown form of triangle specification")
	v1 <- v$v1
	v2 <- v$v2
	v3 <- v$v3
    }
    structure(list(v1 = v1, v2 = v2, v3 = v3,
                   color = color, color2 = color2, fill = fill,
                   material = material, col.mesh = col.mesh, alpha = alpha,
                   smooth = smooth),
              class = "Triangles3D")
}

is.Triangles3D <- function(x) identical(class(x), "Triangles3D")

updateTriangles <- function(triangles, color, color2, alpha, fill, col.mesh,
                            material, smooth) {
    if (! missing(color)) triangles$color <- color
    if (! missing(color2)) triangles$color2 <- color2
    if (! missing(fill)) triangles$fill <- fill
    if (! missing(col.mesh)) triangles$col.mesh <- col.mesh
    if (! missing(material)) triangles$material <- material
    if (! missing(alpha)) triangles$alpha <- alpha
    if (! missing(smooth)) triangles$smooth <- smooth
    triangles
}

#**** This assumes comparable scaling of dimensions
#**** 5 is the largest exponent for S that will work; smaller is OK
t2ve <- function (triangles) 
{
    vb <- rbind(triangles$v1, triangles$v2, triangles$v3)
    vbmin <- min(vb)
    vbmax <- max(vb)
    S <- 10^5
    score <- function(v, d) floor(as.vector(v %*% d))
    scale <- function(v) (1 - 1 / S) * (v - vbmin) / (vbmax - vbmin)
    vbs <- scale(vb)
    d <- c(1, S, S^2)
    scores <- score(vbs, d)
    vb <- vb[! duplicated(scores),]
    vbs <- scale(vb)
    scores <- score(vbs, d)
    ib <- rbind(match(score(scale(triangles$v1), d), scores),
                match(score(scale(triangles$v2), d), scores),
                match(score(scale(triangles$v3), d), scores))
    list(vb = t(vb), ib = ib)
}

ve2t <- function(ve) {
    list (v1 = t(ve$vb[,ve$ib[1,]]),
          v2 = t(ve$vb[,ve$ib[2,]]),
          v3 = t(ve$vb[,ve$ib[3,]]))
}

unzipTriangleMatrix <- function(tris) {
    if (ncol(tris) != 3)
        stop("triangle matrix must have three columns.")
    if (nrow(tris) %% 3 != 0)
        stop("number of rows in triangle matrix must be divisible by 3")
    n <- nrow(tris) / 3
    list(v1 = tris[3 * (1 : n) - 2,],
         v2 = tris[3 * (1 : n) - 1,],
         v3 = tris[3 * (1 : n),])
}

zipTriangles <- function(tris) {
    n <- nrow(tris$v1)
    if (nrow(tris$v2) != n || nrow(tris$v3) != n)
        stop("vertex arrays must have the same number of rows")
    v <- matrix(0, nrow = 3 * n, ncol = 3)
    v[3 * (1 : n) - 2,] <- tris$v1
    v[3 * (1 : n) - 1,] <- tris$v2
    v[3 * (1 : n),] <- tris$v3
    v
}

colorTriangles <- function(triangles) {
    if (is.function(triangles$color) || is.function(triangles$color2)) {
        v <- (triangles$v1 + triangles$v2 + triangles$v3) / 3
        if (is.function(triangles$color))
            triangles$color <- triangles$color(v[,1], v[,2], v[,3])
        if (is.function(triangles$color2))
            triangles$color2 <- triangles$color2(v[,1], v[,2], v[,3])
        if (is.function(triangles$col.mesh))
            triangles$col.mesh <- triangles$col.mesh(v[,1], v[,2], v[,3])
    }
    triangles
}

colorScene <- function(scene) {
    if (is.Triangles3D(scene))
        colorTriangles(scene)
    else lapply(scene, colorTriangles)
}

## **** better to make new triangles including only requested components?
canonicalizeAndMergeScene <- function(scene, ...) {
    which <- list(...)
    if (is.Triangles3D(scene)) {
        n.tri <- nrow(scene$v1)
        for (n in which)
            if (length(scene[[n]]) != n.tri)
                scene[[n]] <- rep(scene[[n]], length = n.tri)
        scene
    }
    else {
        scene <- lapply(scene, canonicalizeAndMergeScene, ...)
        x <- scene[[1]]
        x$v1 <- do.call(rbind, lapply(scene, function(x) x$v1))
        x$v2 <- do.call(rbind, lapply(scene, function(x) x$v2))
        x$v3 <- do.call(rbind, lapply(scene, function(x) x$v3))
        for (n in which)
            x[[n]] <- do.call(c, lapply(scene, function(x) x[[n]]))
        x
    }
}

expandTriangleGrid <- function(x, y) {
    nx <- length(x) - 1
    ny <- length(y) - 1
    A <- c(0, 0)
    B <- c(1, 0)
    C <- c(1, 1)
    D <- c(0, 1)
    g <- expand.grid(x = 1 : nx, y = 1 : ny)
    even <- (g$x + g$y) %% 2 == 0
    gx11 <- ifelse(even, g$x + A[1], g$x + A[1])
    gy11 <- ifelse(even, g$y + A[2], g$y + A[2])
    gx12 <- ifelse(even, g$x + A[1], g$x + B[1])
    gy12 <- ifelse(even, g$y + A[2], g$y + B[2])
    i1 <- rbind(cbind(gx11, gy11), cbind(gx12, gy12))
    gx21 <- ifelse(even, g$x + B[1], g$x + B[1])
    gy21 <- ifelse(even, g$y + B[2], g$y + B[2])
    gx22 <- ifelse(even, g$x + C[1], g$x + C[1])
    gy22 <- ifelse(even, g$y + C[2], g$y + C[2])
    i2 <- rbind(cbind(gx21, gy21), cbind(gx22, gy22))
    gx31 <- ifelse(even, g$x + C[1], g$x + D[1])
    gy31 <- ifelse(even, g$y + C[2], g$y + D[2])
    gx32 <- ifelse(even, g$x + D[1], g$x + D[1])
    gy32 <- ifelse(even, g$y + D[2], g$y + D[2])
    i3 <- rbind(cbind(gx31, gy31), cbind(gx32, gy32))
    v1 <- cbind(x[i1[,1]], y[i1[,2]])
    v2 <- cbind(x[i2[,1]], y[i2[,2]])
    v3 <- cbind(x[i3[,1]], y[i3[,2]])
    list(v1 = v1, v2 = v2, v3 = v3)
}

## adapted from lattice ltransform3dto3d
trans3dto3d <- function (x, R.mat) {
    if (length(x) == 0)
        return(x)
    val <- R.mat %*% rbind(t(x), 1)
    val[1, ] <- val[1, ]/val[4, ]
    val[2, ] <- val[2, ]/val[4, ]
    val[3, ] <- val[3, ]/val[4, ]
    t(val[1:3, , drop = FALSE])
}

transformTriangles <- function(triangles, R) {
    tr <- function(v) trans3dto3d(v, R)
    triangles$v1 <- tr(triangles$v1)
    triangles$v2 <- tr(triangles$v2)
    triangles$v3 <- tr(triangles$v3)
    triangles
}

transformScene <- function(scene, rot.mat) {
    if (is.Triangles3D(scene))
        transformTriangles(scene, rot.mat)
    else lapply(scene, transformTriangles, rot.mat)
}

translateTriangles <- function(triangles, x = 0, y = 0, z = 0) {
    M <- diag(4)
    M[1:3,4] <- c(x, y, z)
    transformTriangles(triangles, M)
}

scaleTriangles <- function(triangles, x = 1, y = x, z = x) {
    M <- diag(c(x, y, z, 1))
    transformTriangles(triangles, M)
}

## triangleNormals computes the normal vectors to a collection of
## triangles as the vector crossprocuct of the direction from v1 to v2
## and the direction from v2 to v3.  The result is an n by 3 matrix of
## unit representing the n unit normal vectors.

triangleNormals <- function(triangles) {
   x <- triangles$v2 - triangles$v1
   y <- triangles$v3 - triangles$v2
   z <- cbind(x[,2]*y[,3] - x[,3]*y[,2],
              x[,3]*y[,1] - x[,1]*y[,3],
              x[,1]*y[,2] - x[,2]*y[,1])
   z / sqrt(rowSums(z^2))
}

# adapted from lattice ltransform3dMatrix
trans3dMat <- function (screen, P = diag(4)) {
    givens4 <- function(i, j, gamma) {
        T <- diag(4)
        cgamma <- cos(gamma)
        sgamma <- sin(gamma)
        T[c(i,j),c(i,j)] <- matrix(c(cgamma, sgamma, -sgamma, cgamma), 2, 2)
        T
    }
    screen.names <- names(screen)
    for (i in seq(along = screen.names)) {
        if (screen.names[i] == "x")
            P <- givens4(2, 3, screen[[i]] * pi/180) %*% P
        else if (screen.names[i] == "y")
            P <- givens4(1, 3, -screen[[i]] * pi/180) %*% P #**** whi negative?
        else if (screen.names[i] == "z")
            P <- givens4(1, 2, screen[[i]] * pi/180) %*% P
    }
    P
}

makeViewTransform <- function(ranges, scale, aspect, screen, R.mat) {
    m <- c(mean(ranges$xlim), mean(ranges$ylim), mean(ranges$zlim))
    s <- 0.5 * c(diff(ranges$xlim), diff(ranges$ylim), diff(ranges$zlim))
    if (! scale) s <- rep(max(s), 3)
    else s <- s / c(1, aspect)
    A <- diag(1 / c(s, 1))
    A[1:3, 4] <- -m / s
    trans3dMat(screen, R.mat %*% A)
}

trianglesRanges <- function(triangles, xlim, ylim, zlim) {
    v1 <- triangles$v1
    v2 <- triangles$v2
    v3 <- triangles$v3
    if (is.null(xlim)) xlim <- range(v1[,1], v2[,1], v3[,1], na.rm = TRUE)
    if (is.null(ylim)) ylim <- range(v1[,2], v2[,2], v3[,2], na.rm = TRUE)
    if (is.null(zlim)) zlim <- range(v1[,3], v2[,3], v3[,3], na.rm = TRUE)
    list(xlim = xlim, ylim = ylim, zlim = zlim)
}

sceneRanges <- function(scene, xlim, ylim, zlim) {
    if (is.Triangles3D(scene))
        trianglesRanges(scene, xlim, ylim, zlim)
    else {
        ranges <- lapply(scene, trianglesRanges, xlim, ylim, zlim)
        list(xlim = range(sapply(ranges,function(x) x$xlim)),
             ylim = range(sapply(ranges,function(x) x$ylim)),
             zlim = range(sapply(ranges,function(x) x$zlim)))
    }
}

addTrianglesPerspective <- function(triangles, distance) {
    pt <- function(v) {
        v[, 1] <- v[, 1] / (1 / distance - v[, 3])
        v[, 2] <- v[, 2] / (1 / distance - v[, 3])
        v
    }
    triangles$v1 <- pt(triangles$v1)
    triangles$v2 <- pt(triangles$v2)
    triangles$v3 <- pt(triangles$v3)
    triangles
}

addPerspective <- function(scene, distance) {
    if (is.Triangles3D(scene))
        addTrianglesPerspective(scene, distance)
    else lapply(scene, addTrianglesPerspective, distance)
}

screenRange <- function(v1, v2, v3)
    range(v1[,1:2], v2[,1:2], v3[,1:2], na.rm = TRUE)

vertexTriangles <- function(ve) {
    n.vert <- ncol(ve$vb)
    ib <- ve$ib
    vt <- function(i) which(ib[1,] == i | ib[2,] == i | ib[3,] == i)
    lapply(1 : n.vert, vt)
}

# faster version
vertexTriangles <- function(ve) {
    n.vert <- ncol(ve$vb)
    val <- vector("list", n.vert)
    ib <- ve$ib
    for (i in 1 : ncol(ib)) {
        val[[ib[1,i]]] <- c(val[[ib[1,i]]], i)
        val[[ib[2,i]]] <- c(val[[ib[2,i]]], i)
        val[[ib[3,i]]] <- c(val[[ib[3,i]]], i)
    }
    val
}


vertexNormals <- function(vt, N) {
    vn <- function(tris) {
        z <- apply(N[tris,,drop = FALSE], 2, mean, na.rm = TRUE);
        z <- z / sqrt(sum(z^2))
        if (any(is.na(z))) c(1,0,0) else z
    }
    t(sapply(vt, vn))
}

# faster version
vertexNormals <- function(vt, N) {
    val <- matrix(0, nrow = length(vt), ncol = 3)
    for (i in seq(along = vt)) {
        Ni <- N[vt[[i]],,drop = FALSE]
        Ni1 <- Ni[,1]
        Ni2 <- Ni[,2]
        Ni3 <- Ni[,3]
        z1 <- if (any(is.na(Ni1))) mean(Ni1, na.rm = TRUE)
              else sum(Ni1) / length(Ni1)
        z2 <- if (any(is.na(Ni2))) mean(Ni2, na.rm = TRUE)
              else sum(Ni2) / length(Ni2)
        z3 <- if (any(is.na(Ni3))) mean(Ni3, na.rm = TRUE)
              else sum(Ni3) / length(Ni3)
        z <- c(z1, z2, z3)
        z <- z / sqrt(sum(z^2))
        val[i,] <- if (any(is.na(z))) c(1,0,0) else z
    }
    val
}

interpolateVertexNormals <- function(VN, ib) {
    z <- (VN[ib[1,],] + VN[ib[2,],] + VN[ib[3,],]) / 3
    z / sqrt(rowSums(z^2))
}

vertexColors <- function(vt, col) {
    C <- t(col2rgb(col))
    val <- matrix(0, nrow = length(vt), ncol = 3)
    for (i in seq(along = vt)) {
        vti <- vt[[i]]
        nti <- length(vti)
        Ci <- C[vti,,drop = FALSE]
        Ci1 <- Ci[,1]
        Ci2 <- Ci[,2]
        Ci3 <- Ci[,3]
        val[i,] <- c(sum(Ci1), sum(Ci2), sum(Ci3)) / nti
    }
    val
}

interpolateVertexColors <- function(VC, ib) {
    TC <- (VC[ib[1,],] + VC[ib[2,],] + VC[ib[3,],]) / 3
    rgb(TC[,1], TC[,2], TC[,3], max=255)
}

triangleEdges <- function(vb, ib) {
    edges <- cbind(ib[c(1,2),], ib[c(2,3),], ib[c(3,1),])
    swap <- edges[1,] > edges[2,]
    edges[,swap] <- edges[2:1,swap]
    edges[,! duplicated(edges, MARGIN = 2)]
}

# faster version
triangleEdges <- function(vb, ib) {
    n.vert <- ncol(vb)
    edges <- cbind(ib[c(1,2),], ib[c(2,3),], ib[c(3,1),])
    swap <- edges[1,] > edges[2,]
    edges[,swap] <- edges[2:1,swap]
    score <- as.vector(c(1 + n.vert, 1) %*% edges)
    edges[,! duplicated(score)]
}

triangleMidTriangles <- function(vb, ib, VN) {
    n.vert <- ncol(vb)
    edges <- triangleEdges(vb, ib)
    vb <- (vb[,edges[1,]] + vb[,edges[2,]]) / 2
    d <- c(1 + n.vert, 1)
    scores <- as.vector(d %*% edges)
    mpi <- function(a, b) {
        s <- d[1] * pmin(a, b) + d[2] * pmax(a, b)
        match(s, scores)
    }
    mpi1 <- mpi(ib[1,], ib[2,])
    mpi2 <- mpi(ib[2,], ib[3,])
    mpi3 <- mpi(ib[3,], ib[1,])
    ib <- rbind(mpi1, mpi2, mpi3)
    z <- VN[edges[1,],] + VN[edges[2,],]
    z <- z / sqrt(rowSums(z^2))
    list(vb = vb, ib = ib, VN = z)
}

## surfaceTriangles creates a set of triangles for a grid specified by x,
## y and function falues computed with f if f is a function or taken
## from f if f is a matrix.

surfaceTriangles <- function(x, y, f, 
                             color = "red", color2 = NA,  alpha = 1,
                             fill = TRUE, col.mesh = if (fill) NA else color,
                             smooth = 0, material = "default") {
    if (is.function(f))
        ff <- function(ix, iy) f(x[ix], y[iy])
    else
        ff <- function(ix, iy) f[ix + length(x) * (iy - 1)]
    i <- expandTriangleGrid(1 : length(x), 1 : length(y))
    i1 <- i$v1
    i2 <- i$v2
    i3 <- i$v3
    v1 <- cbind(x[i1[,1]], y[i1[,2]], ff(i1[,1], i1[,2]))
    v2 <- cbind(x[i2[,1]], y[i2[,2]], ff(i2[,1], i2[,2]))
    v3 <- cbind(x[i3[,1]], y[i3[,2]], ff(i3[,1], i3[,2]))
    na1 <- is.na(v1[,1]) | is.na(v1[,2]) | is.na(v1[,3])
    na2 <- is.na(v2[,1]) | is.na(v2[,2]) | is.na(v2[,3])
    na3 <- is.na(v3[,1]) | is.na(v3[,2]) | is.na(v3[,3])
    nna <- ! (na1 | na2 | na3)
    makeTriangles(v1[nna,], v2[nna,], v3[nna,],
                  color = color, color2 = color2, fill = fill, smooth = smooth,
                  material = material, col.mesh = col.mesh, alpha = alpha)
}

