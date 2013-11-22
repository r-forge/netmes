#This file belongs to
#netmes: NETMES, <http://www.bioconductor.org/packages/release/Software.html>
#This R package allows local network based assesments for
##inferring regulatory networks from expression data.
## Copyright (C) August 2009 Frank Emmert-Streib and Gokmen Altay 
##<v@bio-complexity.com>
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU LESSER GENERAL PUBLIC LICENSE
## as published by the Free Software Foundation; either version 3
## of the License, or any later version.
##
## This program is distributed WITHOUT ANY WARRANTY; 
## You can get a copy of the GNU LESSER GENERAL PUBLIC LICENSE
## from
## http://www.r-project.org/Licenses/LGPL-3

motifcount3 <- function(net, rmat, MImat){ 

  # net: directed network (edges 0/1)
  # G: number of nodes (genes)

  G <- dim(net)[1]  # number of genes is assigned inside the local function

  
  moty <- c(0,0,0,0,0)  # five different motifs
  MImoty <- c()
  motyprob <- new.env()  # prob to observe motifs
  for(i in 1:5){
    assign(as.character(i), c(), envir=motyprob)
  }
  motyMI <- new.env()  # prob to observe motifs
  for(i in 1:5){
    assign(as.character(i), c(), envir=motyMI)
  }

  emotnodes <- new.env()
  for(i in 1:5){
    assign(as.character(i), new.env(), envir = emotnodes)
  }
  
  
  for(i in 1:G){
    ind1 <- which(net[i,] != 0)  # outgoing edges
    L1f <- length(ind1)
    ind1 <- append(ind1, which(net[,i] != 0), after = length(ind1)) # append incoming edges

    L1 <- length(ind1)
    if(L1 > 0){
      for(j in 1:L1){
        n2 <- ind1[j]  # second node in the motif
        ind2 <- which(net[n2,] != 0)  # outgoing edges
        ind2 <- setdiff(ind2, i)
        L2f <- length(ind2)
        ind2 <- append(ind2, which(net[,n2] != 0), after = length(ind2)) # append incoming edges
        if(is.element(i, ind2)){
          ind <- which(ind2 == i)
          ind2 <- ind2[-ind]
        }
        
        L2 <- length(ind2)

        if(L2 > 0){
          for(k in 1:L2){
            n3 <- ind2[k] # third node in the motif

            # check what type of the motif (number)
            mtype <- 0
            if(net[i,n3] != 0 | net[n3,i] != 0){ # three edges
              if((j <= L1f & k <= L2f & net[n3,i] != 0) | (j > L1f & k > L2f & net[i,n3] != 0)){ # type 5 (cycle)
                moty[5] <- moty[5] + 1
                mtype <- 5
              }
              else{  # type 4
                moty[4] <- moty[4] + 1
                mtype <- 4
              }
            }
            else{
              if(j <= L1f & k <= L2f){ # type 1 (chain)
                moty[1] <- moty[1] + 1
                mtype <- 1
              }
              if(j <= L1f & k > L2f){ # type 2 (collider)
                moty[2] <- moty[2] + 1
                mtype <- 2
              }
              if(j > L1f & k <= L2f){ # type 3 (fork)
                moty[3] <- moty[3] + 1
                mtype <- 3
              }
              if(j > L1f & k > L2f){ # type 1 (chain)
                moty[1] <- moty[1] + 1
                mtype <- 1
              }
            }

            if(mtype > 0){
              prob <- motifrate(rmat, net, c(i,n2,n3))
              aux <- get(as.character(mtype), envir = motyprob)
              aux <- append(aux, prob, after = length(aux))
              assign(as.character(mtype), aux, envir = motyprob)

              if(mtype == 1){
                assign(as.character(moty[mtype]), c(i,n2,n3), envir = emotnodes$"1")
              }
              if(mtype == 2){
                assign(as.character(moty[mtype]), c(i,n2,n3), envir = emotnodes$"2")
              }
              if(mtype == 3){
                assign(as.character(moty[mtype]), c(i,n2,n3), envir = emotnodes$"3")
              }
              if(mtype == 4){
                assign(as.character(moty[mtype]), c(i,n2,n3), envir = emotnodes$"4")
              }
              if(mtype == 5){
                assign(as.character(moty[mtype]), c(i,n2,n3), envir = emotnodes$"5")
              }


              if(mtype == 2){
                cat("start \n")
                cat( c(i,n2,n3), "\n")
                cat( motifMI(MImat, c(i,n2,n3)), "\n") 
              }
              
              MI <- motifMI(MImat, c(i,n2,n3))
              aux <- get(as.character(mtype), envir = motyMI)
              aux <- append(aux, MI, after = length(aux))
              assign(as.character(mtype), aux, envir = motyMI)
            }
              
            
          }
        }
      }
    }
    
  }

  # remove symmetries
  moty[1:3] <- moty[1:3]/2
  moty[4:5] <- moty[4:5]/6

  res <- new.env() 
  assign("moty", moty, envir = res)         # number of motifs
  assign("motyprob", motyprob, envir = res) # prob to observe a motif
  assign("motyMI", motyMI, envir = res) # average MI of motifs
  assign("emotnodes", emotnodes, envir = res) # node indices for motifs

res
}


