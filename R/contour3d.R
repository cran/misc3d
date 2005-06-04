contour3d <- function(f, level, 
                      x = 1:dim(f)[1], y = 1:dim(f)[2], z = 1:dim(f)[3], 
                      add = FALSE, ...){ 
    
        PreProcessing <- local({            
		implode <- function(x)    
        		sum(2^(0:7) * x) + 1    
		explode <- function(x)    
        		floor(((x - 1) %% 2^(1:8))/2^(0:7))    
		flip <- function(x)    
        		implode(ifelse(explode(x) == 1, 0, 1))    

                BasicRotation <- matrix(data=c(1,2,3,4,5,6,7,8,5,6,2,1,8,7,3,4,8,7,6,5,4,3,2,1,    
                           4,3,7,8,1,2,6,5,2,6,7,3,1,5,8,4,6,5,8,7,2,1,4,3,    
                           5,1,4,8,6,2,3,7,4,1,2,3,8,5,6,7,3,4,1,2,7,8,5,6,    
                           2,3,4,1,6,7,8,5,6,7,3,2,5,8,4,1,7,8,4,3,6,5,1,2,    
                           8,5,1,4,7,6,2,3,7,3,2,6,8,4,1,5,4,8,5,1,3,7,6,2,    
                           3,2,6,7,4,1,5,8,2,1,5,6,3,4,8,7,1,4,8,5,2,3,7,6,    
                           1,5,6,2,4,8,7,3,5,8,7,6,1,4,3,2,8,4,3,7,5,1,2,6,    
                           3,7,8,4,2,6,5,1,7,6,5,8,3,2,1,4,6,2,1,5,7,3,4,8),     
                           ncol=8, byrow=TRUE)    
                   
               CaseRotation <- matrix(data=c(1,24,2,19,2,17,3,17,2,24,4,24,3,24,6,10,2,15,3,19,    
                             4,17,6,9,3,9,6,8,6,1,9,23,2,20,3,18,4,7,6,16,5,24,    
                             7,5,7,24,12,9,4,20,6,22,8,24,10,24,7,9,15,24,13,20,    
                             6,20,2,21,4,6,3,16,6,4,4,16,8,23,6,14,10,23,5,21,    
                             7,10,7,16,15,9,7,2,13,8,12,23,6,6,3,6,6,17,6,18,    
                             9,18,7,4,13,17,15,18,6,13,7,6,12,16,13,18,6,2,11,24,    
                             7,3,7,12,3,12,2,23,5,23,4,23,7,1,3,14,7,14,6,21,    
                             15,23,4,15,7,19,8,19,13,23,6,11,12,17,10,19,6,23,4,12,    
                             7,18,8,22,13,16,7,13,11,23,13,21,7,15,8,21,13,22,14,24,    
                             8,15,13,11,7,7,8,12,4,22,3,23,7,23,6,24,12,18,6,7,    
                             13,19,9,24,6,19,7,21,11,18,13,24,7,20,15,16,7,22,6,15,    
                             3,22,6,3,15,17,10,22,6,12,12,24,7,11,6,5,3,15,13,10,    
                             7,8,8,20,4,9,7,17,5,22,4,18,2,22,2,22,4,18,5,22,    
                             7,17,4,9,8,20,7,8,13,10,3,15,6,5,7,11,12,24,6,12,    
                             10,22,15,17,6,3,3,22,6,15,7,22,15,16,7,20,13,24,11,18,    
                             7,21,6,19,9,24,13,19,6,7,12,18,6,24,7,23,3,23,4,22,    
                             8,12,7,7,13,11,8,15,14,24,13,22,8,21,7,15,13,21,11,23,    
                             7,13,13,16,8,22,7,18,4,12,6,23,10,19,12,17,6,11,13,23,    
                             8,19,7,19,4,15,15,23,6,21,7,14,3,14,7,1,4,23,5,23,    
                             2,23,3,12,7,12,7,3,11,24,6,2,13,18,12,16,7,6,6,13,    
                             15,18,13,17,7,4,9,18,6,18,6,17,3,6,6,6,12,23,13,8,    
                             7,2,15,9,7,16,7,10,5,21,10,23,6,14,8,23,4,16,6,4,    
                             3,16,4,6,2,21,6,20,13,20,15,24,7,9,10,24,8,24,6,22,    
                             4,20,12,9,7,24,7,5,5,24,6,16,4,7,3,18,2,20,9,23,    
                             6,1,6,8,3,9,6,9,4,17,3,19,2,15,6,10,3,24,4,24,    
                             2,24,3,17,2,17,2,19,1,24),ncol=2,byrow=TRUE)    
            
                CaseRotationFlip <- cbind(CaseRotation, c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1,    
                                  1,1,-1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,-1,-1,1,1,    
                                  1,1,1,1,1,-1,1,1,1,-1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,1,1,1,    
                                 -1,1,-1,-1,-1,1,1,1,1,1,1,1,-1,1,1,-1,-1,1,-1,-1,-1,1,1,1,1,    
                                  1,-1,1,-1,1,1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,    
                                 -1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,-1,1,1,-1,-1,1,1,1,1,1,-1,    
                                 -1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,1,-1,1,1,-1,-1,1,-1,-1,-1,-1,-1,-1,    
                                 -1,1,1,1,-1,1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,1,1,-1,1,-1,-1,-1,    
                                  1,-1,-1,-1,-1,-1,-1,-1,1,1,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,    
                                  1,1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,    
                                 -1,-1,-1,-1,-1,-1))    
                EdgePoints <- matrix(data=c(1,1,2,2,2,3,3,3,4,4,4,1,5,5,6,6,6,7,    
                            7,7,8,8,8,5,9,1,5,10,2,6,11,3,7,12,4,8,13,9,9),ncol=3, byrow=TRUE)    
    
                       
                    
                BasicEdges <- list(c(1,4,9),    
                   c(2,4,9,10),    
                   c(1,4,5,6,9,10),    
                   c(1,4,6,7,9,11),    
                   c(1,4,10,11,12),    
                   c(2,4,6,7,9,10,11),    
                   c(1,2,5,6,7,8,9,10,11,13),    
                   c(9,10,11,12),    
                   c(1,2,7,8,9,11),    
                   c(1,2,3,4,5,6,7,8,13),    
                   c(1,2,6,7,9,12),    
                   c(1,4,5,8,9,10,11,12,13),    
                   c(1,2,3,4,5,6,7,8,9,10,11,12),    
                   c(1,4,7,8,10,11))     
    
                  
                                    
                EdgeSequence1 <- list(    
                           c(1,2,3),    
                           c(1,2,4,2,3,4),    
                           list(c(1,2,5,3,4,6),    
                                c(1,2,6,2,4,6,2,3,4,2,5,3)),    
                           list(c(1,2,5,3,4,6),    
                                c(1,3,5,3,4,5,2,5,4,1,2,6,1,6,3,2,4,6)),    
                           c(1,3,2,2,3,5,3,4,5),    
                           list(c(1,2,6,2,5,6,3,4,7),    
                                c(1,2,7,2,4,7,2,5,4,3,4,5,3,5,6,1,6,3,1,3,7,7,3,4),    
                                c(1,2,7,2,4,7,2,5,4,3,4,5,3,5,6)),    
                           list(c(1,8,2,3,7,6,4,5,9),    
                                c(1,8,2,3,9,4,3,7,9,5,9,7,5,7,6),    
                                c(3,7,6,1,8,4,1,4,5,1,9,2,1,5,9),    
                                c(4,5,9,1,7,6,1,6,2,2,6,3,2,3,8),    
                                c(1,10,2,2,10,9,5,9,10,5,10,6,6,10,7,3,7,10,3,10,4,4,10,8,1,8,10),    
                                c(1,10,2,1,10,7,6,10,7,5,10,6,5,9,10,4,9,10,3,10,4,3,8,10,2,10,8),    
                                c(5,9,10,2,10,9,1,10,2,1,7,10,6,10,7,3,10,6,3,8,10,4,10,8,4,5,10),    
                                c(1,7,6,1,6,9,1,9,2,5,9,6,3,8,4),    
                                c(1,7,8,3,8,7,3,7,6,3,6,5,3,5,4,4,5,9,4,9,8,2,8,9,1,8,2)),    
                           c(1,2,3,1,3,4),    
                           c(1,2,6,1,3,5,1,6,3,3,4,5),    
                           list(c(1,4,5,2,6,3,3,6,7,4,8,5),    
                                c(1,4,3,1,3,2,3,4,8,3,8,7,5,8,6,6,8,7,1,2,6,1,6,5),    
                                c(1,2,9,1,9,5,5,9,8,4,8,9,3,4,9,3,9,7,6,7,9,2,6,9),    
                                c(5,9,6,1,9,5,1,4,9,4,8,9,7,9,8,3,9,7,2,9,3,2,6,9),    
                                c(1,2,5,2,6,5,3,4,8,3,8,7)),    
                           c(1,2,3,1,3,6,1,6,5,3,4,6),    
                           list(c(1,6,2,2,6,8,3,5,4,6,7,8),    
                                c(1,5,3,1,3,6,3,6,7,3,7,4,2,4,8,2,5,4,1,5,2,4,7,8),    
                                c(1,9,2,2,9,5,3,5,9,4,9,8,3,9,4,7,8,9,6,7,9,1,9,6),    
                                c(4,9,5,1,5,9,1,9,2,2,9,8,7,8,9,6,9,7,3,6,9,3,9,4),    
                                c(1,5,2,3,8,4,3,6,8,6,7,8)),    
                            c(1,4,9,2,11,3,5,6,10,7,8,12),    
                            c(1,4,2,1,6,4,1,5,6,3,6,4))     
                      
                EdgeSequence2 <- list(    
                           c(1,3,2),    
                           c(1,4,2,2,4,3),    
                           list(c(1,5,2,3,6,4),    
                                c(1,6,2,2,6,4,2,4,3,2,3,5)),    
                           list(c(1,5,2,3,6,4),    
                                c(1,5,3,3,5,4,2,4,5,1,6,2,1,3,6,2,6,4)),    
                           c(1,2,3,2,5,3,3,5,4),    
                           list(c(1,6,2,2,6,5,3,7,4),    
                                c(1,7,2,2,7,4,2,4,5,3,5,4,3,6,5,1,3,6,1,7,3,7,4,3),    
                                c(1,7,2,2,7,4,2,4,5,3,5,4,3,6,5)),    
                           list(c(1,2,8,3,6,7,4,9,5),    
                                c(1,2,8,3,4,9,3,9,7,5,7,9,5,6,7),    
                                c(3,6,7,1,4,8,1,5,4,1,2,9,1,9,5),    
                                c(4,9,5,1,6,7,1,2,6,2,3,6,2,8,3),    
                                c(1,2,10,2,9,10,5,9,10,5,6,10,6,7,10,3,10,7,3,4,10,4,8,10,1,10,8),    
                                c(1,2,10,1,7,10,6,7,10,5,6,10,5,10,9,4,10,9,3,4,10,3,10,8,2,8,10),    
                                c(5,10,9,2,9,10,1,2,10,1,10,7,6,7,10,3,6,10,3,10,8,4,8,10,4,10,5),    
                                c(1,6,7,1,9,6,1,2,9,5,6,9,3,4,8),    
                                c(1,8,7,3,7,8,3,6,7,3,5,6,3,4,5,4,9,5,4,8,9,2,9,8,1,2,8)),    
                           c(1,3,2,1,4,3),    
                           c(1,6,2,1,5,3,1,3,6,3,5,4),    
                           list(c(1,5,4,2,3,6,3,7,6,4,5,8),    
                                c(1,3,4,1,2,3,3,8,4,3,7,8,5,6,8,6,7,8,1,6,2,1,5,6),    
                                c(1,9,2,1,5,9,5,8,9,4,9,8,3,9,4,3,7,9,6,9,7,2,9,6),    
                                c(5,6,9,1,5,9,1,9,4,4,9,8,7,8,9,3,7,9,2,3,9,2,9,6),    
                                c(1,5,2,2,5,6,3,8,4,3,7,8)),    
                           c(1,3,2,1,6,3,1,5,6,3,6,4),    
                           list(c(1,2,6,2,8,6,3,4,5,6,8,7),    
                                c(1,3,5,1,6,3,3,7,6,3,4,7,2,8,4,2,4,5,1,2,5,4,8,7),    
                                c(1,2,9,2,5,9,3,9,5,4,8,9,3,4,9,7,9,8,6,9,7,1,6,9),    
                                c(4,5,9,1,9,5,1,2,9,2,8,9,7,9,8,6,7,9,3,9,6,3,4,9),    
                                c(1,2,5,3,4,8,3,8,6,6,8,7)),    
                           c(1,9,4,2,3,11,5,10,6,7,12,8),    
                           c(1,2,4,1,4,6,1,6,5,3,4,6))     
    
                    
                GetEdges <- local({    
                        Edges <-            
                        apply(CaseRotationFlip[-c(1,256),], 1, function(x) {    
                                case <- x[1]    
                                rotation <- x[2]    
                                map <- rep(0,8)    
                                for(i in 1:8){    
                                   temp <- as.integer(BasicRotation[rotation,][i])    
                                   map[temp] <- i}    
                                sapply(BasicEdges[[case-1]], function(x){    
                                         if (x!=13){    
                                           EndP1 <- EdgePoints[x,2]    
                                           EndP2 <- EdgePoints[x,3]    
                                           newEdge <- EdgePoints[(EdgePoints[,2]==map[EndP1]    
                                                                      &EdgePoints[,3]==map[EndP2])|    
                                                                     (EdgePoints[,3]==map[EndP1]    
                                                                      &EdgePoints[,2]==map[EndP2]),][1]    
                                          }    
                                          else  newEdge <- 13    
                                          newEdge})    
                        })       
    
                        Case <- cbind(seq(1:256), CaseRotationFlip[,c(1,3)])    
                        Edges <- apply(Case[-c(1,256),], 1, function(x){    
                                case <- x[2]-1          
                                EdgeNum <- x[1]-1    
                                
                                if (x[3]==1)      
                                    sapply(EdgeSequence1[[case]], function(x) Edges[[EdgeNum]][x])     
                                    else sapply(EdgeSequence2[[case]], function(x) Edges[[EdgeNum]][x])     
                         })    
                        Edges    
                })    
    
                BasicFace <- list(c(0),c(0),c(1),c(7),c(0),c(2,7),c(1,2,6,7),c(0),c(0),c(5,6,7),c(0),    
                              c(1,4,7),c(0),c(0))    
                FacePoints <- matrix(data=c(seq(1,6),1,2,4,1,1,5,6,7,7,8,3,7,2,3,3,4,2,6,5,6,8,5,4,8),ncol=5)    
                FacePoints <- cbind(FacePoints, apply(FacePoints[,2:5],1,prod))    
    
                 
                GetFaces <- local({    
                        Faces <-     
                        apply(CaseRotationFlip[-c(1,256),], 1, function(x) {    
                                 case <- x[1]    
                                 rotation <- x[2]    
                                 map <- rep(0,8)    
                                 for(i in 1:8){    
                                     temp <- as.integer(BasicRotation[rotation,][i])    
                                     map[temp] <- i}    
                                 sapply(BasicFace[[case-1]], function(x){    
                                           EndP <- rep(0,4)    
                                           if (x==0) newFace <- 0    
                                           else if (x==7) newFace <- 7    
                                           else {    
                                              for (i in 1:4){    
                                              point <- FacePoints[x,i+1]    
                                              EndP[i] <- map[point]    
                                              }     
                                              newFace<- FacePoints[FacePoints[,6]==prod(EndP[1:4]),][1]    
                                            }    
                                            newFace})    
                          })     
    
                    
                        for (i in 1:254){    
                            for(j in 1:length(Faces[[i]])){    
                            x <- Faces[[i]][j]    
                            if (x!=0){    
                                index <- explode(i+1)    
                                if (sum(index) > 4)     
                                        index <- ifelse(index==0,1,0)    
                                if (x!=7 && index[FacePoints[x,2]]==0)  
                                        Faces[[i]][j] <- -Faces[[i]][j]  
                                if (x==7){ 
                                        tcase <- CaseRotationFlip[i+1,1]-1 
                                        if ((tcase==6 || tcase==10 ||tcase==12) 
                                             && !(index[1]+index[7]==2) && !(index[3]+index[5]==2)) 
                                          Faces[[i]][j] <- -Faces[[i]][j]  
                                        if (tcase==7 && !(index[1]+index[7]==0) && !(index[3]+index[5]==0))  
                                                Faces[[i]][j] <- -Faces[[i]][j]  
                                }   
                        }}}     
                        Faces    
                })                
                                        
                list(Edges=GetEdges, Faces=GetFaces, EdgePoints=EdgePoints, FacePoints=FacePoints,CRF=CaseRotationFlip)    
        })          
    
    

	implode.vec <- function(x1,x2,x3,x4,x5,x6,x7,x8)    
                x1+2*x2+2^2*x3+2^3*x4+2^4*x5+2^5*x6+2^6*x7+2^7*x8 + 1     	

        fgrid <- function(fun, x, y, z) {  
          g <- expand.grid(x = x, y = y, z = z) 
          array(fun(g$x, g$y, g$z), c(length(x), length(y), length(z))) 
           
        }  
   
	levCells <- function(v, level) { 
    		nx <- dim(v)[1] 
    		ny <- dim(v)[2] 
    		nz <- dim(v)[3] 
    		val <- vector("list", nz - 1) 
    		type <- vector("list", nz - 1) 
    		i <- 1:(nx - 1) 
    		j <- 1:(ny - 1) 
    		v1 <- v[,,1,drop=TRUE] 
    		vv1 <- ifelse(v1 > level, 1, 0) 
		 
        	ttt1 <- vv1[i,j] + 2 * vv1[i+1,j] + 4 * vv1[i+1,j+1] + 8 * vv1[i,j+1] 
    		for (k in 1 : (nz - 1)) { 
        		v2 <- v[,,k + 1,drop=TRUE] 
        		vv2 <- ifelse(v2 > level, 1, 0) 
        		ttt2 <- vv2[i,j] + 2 * vv2[i+1,j] + 4 * vv2[i+1,j+1] + 8 * vv2[i,j+1] 
        		ttt <- ttt1 + 16 * ttt2 
        		iii <- ttt > 0 & ttt < 255 
        		val[[k]] <- which(iii) + (nx - 1) * (ny - 1) * (k - 1) 
        		type[[k]] <- as.integer(ttt[iii]) 
        		v1 <- v2 
        		vv1 <- vv2 
        		ttt1 <- ttt2 
    		} 
    		v <- unlist(val) 
    		i <- as.integer((v - 1) %% (nx - 1) + 1) 
    		j <- as.integer(((v - 1) %/% (nx - 1)) %% (ny - 1) + 1) 
   		k <- as.integer((v - 1) %/% ((nx - 1) * (ny - 1)) + 1)
                t <- unlist(type)
                NAS <- which(is.na(t))
                if (length(NAS) > 0) t <- t[-NAS]
    		list(i = i, j = j, k = k, t = t) 
	} 
 
	            
        GetInfo <- function(cube.1){    
               index <- matrix(c(0,1,1,0,0,1,1,0,    
                          0,0,1,1,0,0,1,1,    
                          0,0,0,0,1,1,1,1),    
                        nrow=8) 
               if (interface==1)
                 ax.inc <- c(x[2]-x[1],y[2]-y[1],z[2]-z[1]) 
               else  ax.inc <- c(1,1,1)    
               ver.inc <- t(apply(index,1, function(x) x*ax.inc))    
               cube.co <- kronecker(rep(1,nrow(cube.1)),ver.inc)+kronecker(cube.1,rep(1,8)) 
               if (interface==1) 
                 value <- f(cube.co[,1], cube.co[,2], cube.co[,3])-level   
               else   
                 value <- apply(cube.co, 1, function(x) vol[x[1], x[2], x[3]]) - level             
               if(interface==1) 
                 information <- cbind(cube.co, value)
               else information <- cbind(x[cube.co[,1]],y[cube.co[,2]],z[cube.co[,3]], value)
               information                                  
       }    
    
        GetPoints<-function(edge,p1){ 
                    x1 <- EP[edge,2] 
                    x2 <- EP[edge,3] 
                    c((1-floor(x1/9))*information[p1+x1-1,1]+floor(x1/9)*information[p1,1], 
                      (1-floor(x1/9))*information[p1+x2-1,1]+floor(x1/9)*information[p1+1,1], 
                      (1-floor(x1/9))*information[p1+x1-1,2]+floor(x1/9)*information[p1+1,2], 
                      (1-floor(x1/9))*information[p1+x2-1,2]+floor(x1/9)*information[p1+2,2], 
                      (1-floor(x1/9))*information[p1+x1-1,3]+floor(x1/9)*information[p1+1,3], 
                      (1-floor(x1/9))*information[p1+x2-1,3]+floor(x1/9)*information[p1+5,3], 
                      (1-floor(x1/9))*information[p1+x1-1,4]+floor(x1/9)*(0*information[p1+1,3]+1), 
                      (1-floor(x1/9))*information[p1+x2-1,4]+floor(x1/9)*(0*information[p1+1,3]-1)) 
        }        	     
     
        CalPoint <- function(x1,x2,y1,y2,z1,z2,v1,v2){    
                       
                     s <- v1 / (v1-v2) 
                     x <- x1+s*(x2-x1)    
                     y <- y1+s*(y2-y1)    
                     z <- z1+s*(z2-z1)    
                    
                     c(x,y,z)     
        }    
    
        FaceNo7 <- function(faces, p1){ 
                index <- ifelse(faces > 0, 1, -1) 
        	faces <- abs(faces) 
                e1 <- FP[faces,2] 
                e2 <- FP[faces,3] 
		e3 <- FP[faces,4] 
                e4 <- FP[faces,5] 
		A <- information[p1+e1-1,4] 
                B <- information[p1+e2-1,4]  
                C <- information[p1+e3-1,4] 
		D <- information[p1+e4-1,4]  
                index <- index*ifelse (A*B-C*D > 0, 1, -1) 
                index <- ifelse(index==1, 1, 0) 
                index 
        }   
     
        Face7 <- function(faces, p1){ 
                index <- ifelse(faces > 0, 1, -1) 
        	A0 <- information[p1,4];   B0 <- information[p1+3,4]  
                C0 <- information[p1+2,4]; D0 <- information[p1+1,4]    
                A1 <- information[p1+4,4]; B1 <- information[p1+7,4]  
                C1 <- information[p1+6,4]; D1 <- information[p1+5,4]    
                a <- (A1 - A0)*(C1 - C0) - (B1 - B0)*(D1 - D0)    
                b <- C0*(A1 - A0) + A0*(C1 - C0) - D0*(B1 - B0) - B0*(D1 - D0)    
                c <- A0*C0 - B0*D0 
                tmax <- -b/(2*a) 
                maximum <- a*tmax^2 + b*tmax + c 
                maximum <- ifelse(maximum=="NaN",-1,maximum) 
                cond1 <- ifelse (a < 0, 1 ,0) 
                cond2 <- ifelse (tmax > 0, 1 ,0) 
                cond3 <- ifelse (tmax < 1, 1, 0) 
                cond4 <- ifelse (maximum >0, 1, 0) 
                totalcond <- cond1 * cond2 * cond3 * cond4 
                index <- index*ifelse(totalcond==1, 1, -1) 
                index <- ifelse(index==1, 1, 0) 
        }	 
 
        Render <- function(edges,p1){ 
                if (typeof(edges)=="list"){ 
                  count <- sapply(edges, function(x) length(x)) 
                  edges <- cbind(unlist(edges), rep(p1,count)) 
                } 
                else{ 
                  count <- nrow(edges) 
                  edges <- cbind(as.vector(t(edges)), rep(p1,each=count)) 
                }   
                 
                information <- GetPoints(edges[,1],edges[,2]) 
                information <- matrix(information,ncol=8) 
                information <- CalPoint(information[,1],information[,2],information[,3],information[,4], 
                                information[,5],information[,6],information[,7],information[,8]) 
        	information <- matrix(information,ncol=3) 
        	rgl.triangles(information[,1], 
                      information[,3], 
                      -information[,2],...)   
        }	 
        Render1 <- function(edges){ 
                  if (is.vector(edges)) edges <- matrix(edges, ncol=length(edges)) 
                  p1 <- edges[,1] 
                  count <- ncol(edges)-1 
                  edges <- cbind(as.vector(t(edges[,-1])), rep(p1,each=count)) 
                 
                 
                  information <- GetPoints(edges[,1],edges[,2]) 
                  information <- matrix(information,ncol=8) 
                  information <- CalPoint(information[,1],information[,2],information[,3],information[,4], 
                                information[,5],information[,6],information[,7],information[,8]) 
        	        	information <- matrix(information,ncol=3) 
                  rgl.triangles(information[,1], 
                      information[,3], 
                      -information[,2],...)   
       }	 
     
              
        Faces <- PreProcessing$Faces    
        Edges <- PreProcessing$Edges  
        EP <- PreProcessing$EdgePoints    
        FP <- PreProcessing$FacePoints     
        CR <- PreProcessing$CRF    
 
         
        if(typeof(f)=="closure") interface <- 1 
        else if (is.array(f) && length(dim(f)) == 3) {
          interface <- 2    
        }     
        else stop("vol has to be a function or a 3-dimensional array") 
 
        if (interface==1) {
           NAS <- unique(c(which(is.na(x)),which(is.na(y)),which(is.na(z))))
            if (length(NAS) > 0){
                x <- x[-NAS]
                y <- y[-NAS]
                z <- z[-NAS]
            }  
            vol <- fgrid(f, x, y, z)
        }
        else vol <- f 

        if  (interface==2){ 
            nx <- dim(vol)[1] 
            ny <- dim(vol)[2] 
            nz <- dim(vol)[3] 
            if (length(x) != nx || length(y) != ny || length(z) != nz) 
               stop("dimensions of f do not match x, y, or z")  
        }
        
	if (!add) rgl.clear() 
 
        v <- levCells(vol, level) 
	tcase <- CR[v$t+1,1]-1 

        R <- which(tcase %in% c(1,2,5,8,9,11,13,14)) 
        if (length(R) > 0){ 
                if (interface==1)
                  cube.1 <- cbind(x[v$i[R]],y[v$j[R]],z[v$k[R]]) 
                else    cube.1 <- cbind(v$i[R],v$j[R],v$k[R]) 
                information <- GetInfo(cube.1) 
		information <- rbind(information,rep(0,4)) 
		p1 <- (1:length(R)-1)*8+1 
                cases <- v$t[R] 
                edges <- Edges[cases]  
                Render(edges,p1) 
        } 
 	 
        R <- which(tcase==3)          
        if(length(R) > 0){ 
                if (interface==1)
                  cube.1 <- cbind(x[v$i[R]],y[v$j[R]],z[v$k[R]]) 
                else    cube.1 <- cbind(v$i[R],v$j[R],v$k[R]) 
		information <- GetInfo(cube.1) 
		information <- rbind(information,rep(0,4)) 
		p1 <- (1:length(R)-1)*8+1 
        	cases <- v$t[R] 
		faces <- unlist(Faces[cases])   
                index <- FaceNo7(faces, p1)              
                edges <- matrix(unlist(Edges[cases]), ncol=18, byrow=TRUE) 
		edges <- cbind(edges, p1,index) 
                edges1 <- edges[which(index==0),c(19,1:6)] 
                edges2 <- edges[which(index==1),c(19,7:18)] 
                 
		if (length(edges1)> 0) Render1(edges1) 
                if (length(edges2)> 0) Render1(edges2) 
        }         
        
	R <- which(tcase==4)          
        if(length(R) > 0){ 
                if (interface==1)
                  cube.1 <- cbind(x[v$i[R]],y[v$j[R]],z[v$k[R]]) 
                else    cube.1 <- cbind(v$i[R],v$j[R],v$k[R]) 
		information <- GetInfo(cube.1) 
		information <- rbind(information,rep(0,4)) 
		p1 <- (1:length(R)-1)*8+1 
        	cases <- v$t[R] 
		faces <- unlist(Faces[cases])   
		index <- Face7(faces, p1) 
                edges <- matrix(unlist(Edges[cases]), ncol=24, byrow=TRUE) 
                edges <- cbind(edges, p1,index) 
                edges1 <- edges[which(index==0),c(25,1:6)] 
                edges2 <- edges[which(index==1),c(25,7:18)] 
                if (length(edges1)> 0) Render1(edges1) 
		if (length(edges2)> 0) Render1(edges2) 
	}         
       	 
        R <- which(tcase==6)          
        if(length(R) > 0){ 
		if (interface==1)
                  cube.1 <- cbind(x[v$i[R]],y[v$j[R]],z[v$k[R]]) 
                else    cube.1 <- cbind(v$i[R],v$j[R],v$k[R]) 
       		information <- GetInfo(cube.1) 
		information <- rbind(information,rep(0,4)) 
		p1 <- (1:length(R)-1)*8+1 
        	cases <- v$t[R] 
		faces <- matrix(unlist(Faces[cases]),ncol=2, byrow=TRUE) 
                index1 <- FaceNo7(faces[,1],p1) 
                index2 <- Face7(faces[,2],p1) 
                index <- index1 + 2*index2 
		edges <- matrix(unlist(Edges[cases]), ncol=48, byrow=TRUE) 
                edges <- cbind(edges, p1,index) 
                edges1 <- edges[which(index==0),c(49,1:9)] 
                edges2 <- edges[which(index==2),c(49,10:33)] 
                edges3 <- edges[which(index==1),c(49,34:48)] 
                edges4 <- edges[which(index==3),c(49,34:48)] 
                if (length(edges1)> 0) 	Render1(edges1) 
                if (length(edges2)> 0) 	Render1(edges2) 
                if (length(edges3)> 0)	Render1(edges3) 
                if (length(edges4)> 0)  Render1(edges4) 
	} 
        
        R <- which(tcase==10)          
        if(length(R) > 0){ 
                if (interface==1)
                  cube.1 <- cbind(x[v$i[R]],y[v$j[R]],z[v$k[R]]) 
                else    cube.1 <- cbind(v$i[R],v$j[R],v$k[R]) 
		information <- GetInfo(cube.1) 
		information <- rbind(information,rep(0,4)) 
		p1 <- (1:length(R)-1)*8+1 
        	cases <- v$t[R] 
		faces <- matrix(unlist(Faces[cases]),ncol=3, byrow=TRUE) 
                index1 <- FaceNo7(faces[,1],p1) 
                index2 <- FaceNo7(faces[,2],p1) 
                index3 <- Face7(faces[,3],p1) 
                index <- index1 + 2*index2 + 4*index3 
		edges <- matrix(unlist(Edges[cases]), ncol=96, byrow=TRUE) 
		edges <- cbind(edges, p1,index) 
                edges1 <- edges[which(index==0),c(97,1:12)] 
                edges2 <- edges[which(index==4),c(97,13:36)] 
                edges3 <- edges[which(index==1),c(97,37:60)] 
                edges4 <- edges[which(index==5),c(97,37:60)] 
                edges5 <- edges[which(index==2),c(97,61:84)] 
                edges6 <- edges[which(index==6),c(97,61:84)] 
                edges7 <- edges[which(index==3),c(97,85:96)] 
                edges8 <- edges[which(index==3),c(97,85:96)] 
  
                if (length(edges1)> 0) Render1(edges1) 
                if (length(edges2)> 0) Render1(edges2) 
                if (length(edges3)> 0) Render1(edges3) 
                if (length(edges4)> 0) Render1(edges4) 
		if (length(edges5)> 0) Render1(edges5) 
                if (length(edges6)> 0) Render1(edges6) 
                if (length(edges7)> 0) Render1(edges7) 
                if (length(edges8)> 0) Render1(edges8) 
	} 
 
	R <- which(tcase==7)          
        if(length(R) > 0){ 
                if (interface==1)
                  cube.1 <- cbind(x[v$i[R]],y[v$j[R]],z[v$k[R]]) 
                else    cube.1 <- cbind(v$i[R],v$j[R],v$k[R]) 
		information <- GetInfo(cube.1) 
		information <- rbind(information,rep(0,4)) 
		p1 <- (1:length(R)-1)*8+1 
        	cases <- v$t[R]		 
        	faces <- matrix(unlist(Faces[cases]),ncol=4, byrow=TRUE) 
                index1 <- FaceNo7(faces[,1],p1) 
                index2 <- FaceNo7(faces[,2],p1) 
		index3 <- FaceNo7(faces[,3],p1) 
                index4 <- Face7(faces[,4],p1) 
                index <- index1 + 2*index2 + 4*index3 + 8*index4 
		edges <- matrix(unlist(Edges[cases]), ncol=177, byrow=TRUE) 
                edges <- cbind(edges, p1,index) 
                edges1 <- edges[which(index==0),c(178,1:9)] 
                edges2 <- edges[which(index==8),c(178,1:9)] 
                edges3 <- edges[which(index==4),c(178,10:24)] 
                edges4 <- edges[which(index==12),c(178,10:24)] 
                edges5 <- edges[which(index==2),c(178,25:39)] 
                edges6 <- edges[which(index==10),c(178,25:39)] 
                edges7 <- edges[which(index==1),c(178,40:54)] 
                edges8 <- edges[which(index==9),c(178,40:54)] 
                edges9 <- edges[which(index==6),c(178,55:81)] 
                edges10 <- edges[which(index==14),c(178,55:81)] 
                edges11 <- edges[which(index==5),c(178,82:108)] 
		edges12 <- edges[which(index==13),c(178,82:108)] 
                edges13 <- edges[which(index==3),c(178,109:135)] 
		edges14 <- edges[which(index==11),c(178,109:135)] 
                edges15 <- edges[which(index==15),c(178,136:150)] 
 		edges16 <- edges[which(index==7),c(178,151:177)] 
                 
                if (length(edges1)> 0) Render1(edges1) 
                if (length(edges2)> 0) Render1(edges2) 
                if (length(edges3)> 0) Render1(edges3) 
                if (length(edges4)> 0) Render1(edges4) 
		if (length(edges5)> 0) Render1(edges5) 
                if (length(edges6)> 0) Render1(edges6) 
                if (length(edges7)> 0) Render1(edges7) 
                if (length(edges8)> 0) Render1(edges8) 
		if (length(edges9)> 0) Render1(edges9) 
                if (length(edges10)> 0) Render1(edges10) 
                if (length(edges11)> 0) Render1(edges11) 
                if (length(edges12)> 0) Render1(edges12) 
		if (length(edges13)> 0) Render1(edges13) 
                if (length(edges14)> 0) Render1(edges14) 
                if (length(edges15)> 0) Render1(edges15) 
                if (length(edges16)> 0) Render1(edges16) 
	} 
         
        R <- which(tcase==12)          
        if(length(R) > 0){ 
                if (interface==1)
                  cube.1 <- cbind(x[v$i[R]],y[v$j[R]],z[v$k[R]]) 
                else    cube.1 <- cbind(v$i[R],v$j[R],v$k[R]) 
		information <- GetInfo(cube.1) 
                information <- rbind(information,rep(0,4)) 
		p1 <- (1:length(R)-1)*8+1 
        	cases <- v$t[R] 
        	faces <- matrix(unlist(Faces[cases]),ncol=3, byrow=TRUE) 
                  
                index1 <- FaceNo7(faces[,1],p1) 
                index2 <- FaceNo7(faces[,2],p1) 
		index3 <- Face7(faces[,3],p1) 
                index <- index1 + 2*index2 + 4*index3 
		edges <- matrix(unlist(Edges[cases]), ncol=96, byrow=TRUE) 
                edges <- cbind(edges, p1,index) 
                edges1 <- edges[which(index==0),c(97,1:12)] 
                edges2 <- edges[which(index==4),c(97,13:36)] 
                edges3 <- edges[which(index==2),c(97,37:60)] 
                edges4 <- edges[which(index==6),c(97,37:60)] 
                edges5 <- edges[which(index==1),c(97,61:84)] 
                edges6 <- edges[which(index==5),c(97,61:84)] 
                edges7 <- edges[which(index==3),c(97,85:96)] 
                edges8 <- edges[which(index==7),c(97,85:96)] 
                 
                if (length(edges1)> 0) Render1(edges1) 
                if (length(edges2)> 0) Render1(edges2) 
                if (length(edges3)> 0) Render1(edges3) 
                if (length(edges4)> 0) Render1(edges4) 
                if (length(edges5)> 0) Render1(edges5) 
                if (length(edges6)> 0) Render1(edges6) 
                if (length(edges7)> 0) Render1(edges7) 
                if (length(edges8)> 0) Render1(edges8)  
        } 
                   
    }
