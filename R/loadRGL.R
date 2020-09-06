loadRGL <- function() {
    if (! suppressWarnings(requireNamespace("rgl",quietly=TRUE)))
        stop("rgl is mot available")
}


