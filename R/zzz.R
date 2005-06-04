.First.lib <- function(lib, pkg){ 
	    have.rgl <- "package:rgl" %in% search() 
    	    if(!have.rgl) require(rgl) 
        } 

