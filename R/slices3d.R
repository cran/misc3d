if (! exists("gray.colors"))
    gray.colors <- function (n, start = 0.3, end = 0.9, gamma = 2.2)
        gray(seq(from = start^gamma, to = end^gamma, length = n)^(1/gamma))

vslice <- function(vol, which, k, tpt = 1) {
    if (length(dim(vol)) == 4)
        switch(which,
               x = vol[k,,,tpt],
               y = vol[,k,,tpt],
               z = vol[,,k,tpt])
     else
        switch(which,
               x = vol[k,,],
               y = vol[,k,],
               z = vol[,,k])
}

slices3d <- function(vol, scale = 0.8, col=gray.colors(512), cross = TRUE) {
    if (! require(tkrplot)) stop("tkrplot is required.");
    r <- range(vol,na.rm = TRUE)
    d <- dim(vol)
    dn <- c("x", "y", "z", "t")
    tt <- tktoplevel()
    bb <- c(round(d[1:3]) / 2, 1)
    bbv <- lapply(bb, tclVar)
    mkimg <- function(which) {
        switch(which,
               x = { i <- 1; j <- 2; k <- 3 },
               y = { i <- 2; j <- 1; k <- 3 },
               z = { i <- 3; j <- 1; k <- 2 })
        f <- function() {
            opar = par(mar=c(0,0,0,0))
            on.exit(par(opar))
            image(vslice(vol, which, bb[i],bb[4]),col=col, zlim = r)
            lines(rep(bb[j]/d[j],100), seq(0,1,len=100))
            lines(seq(0,1,len=100), rep(bb[k]/d[k],100))
        }
        tkrplot(tt, f, hscale = 0.8, vscale = 0.8)
    }
    mkscale <- function(i) {
        f <- function(...) {
            b <- as.numeric(tclvalue(bbv[[i]]))
            if (b != bb[i]) {
                bb[i] <<- b
                if (cross || i == 4)
                    for (j in 1:3) tkrreplot(img[[j]])
                else  tkrreplot(img[[i]])
                tkconfigure(l2, text=bb[i])
            }
        }
        fr <- tkframe(tt)
        s <- tkscale(fr, command=f, from=1, to=d[i], resolution=1,
                variable=bbv[[i]], showvalue=FALSE, orient="horiz")
        l1 <- tklabel(fr, text = dn[i])
        l2 <- tklabel(fr, text = bb[i])
        tkgrid(l1, s, l2)
        fr
    }
    s <- lapply(1:3, mkscale)
    img <- lapply(c("x", "y", "z"), mkimg)
    tkgrid(img[[1]], img[[2]])
    tkgrid(s[[1]],s[[2]])
    tkgrid(img[[3]])
    if (length(d) == 4 && d[4] > 1)
        tkgrid(s[[3]], mkscale(4))
    else tkgrid(s[[3]])
    environment()
}
