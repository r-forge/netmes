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


evalnetworks <- function(E, G, net, infilepath){


  tnet <- abs(net) + abs(t(net)) # true network (make symmetric)
  tnet <- 1*(tnet > 0)           # in case there were bi-directional connections

  cc <- 0
  Fscores <- c()
  Ithreshold <- c()
  mat <- matrix(0, ncol = G, nrow = G)
  matgALL <- matrix(0, ncol = G, nrow = G)
  for(i in 1:E){
    infile <- paste(infilepath, as.character(i), "/mim.exp", sep = "")
    matf <- readMI(infile, G)
    m <- max(matf)
    if(m == Inf){  # filter problem values
      print("problem")
                                        #break
    }
    else{
      cc <- cc + 1
      mat <- mat + matf

      res <- validate(inet = matf, tnet = tnet)
      #options(warn = 1)
      fs <- fscores(res, beta = 1)
      ind <- which( fs == max(fs))[1]  # use the first
      Ithr <- res[ind,1]
      
      matg <- (matf >= Ithr)*1
      matgALL <- matgALL + matg
      Fscores <- append(Fscores, fs[ind], after = length(Fscores))
      Ithreshold <- append(Ithreshold, Ithr, after = length(Ithreshold))
      
    }
    
    
  #if(i == 11)break
  }

  mat <- mat/cc  # averaged mutual information matrix


  ratematgALL <- matgALL/cc # rate matrix
  indedges <- which(abs(net) > 0, arr.ind = TRUE) 
  TPReALL <- ratematgALL[indedges]
  indedges <- which(abs(net) == 0, arr.ind = TRUE) 
  FNReALL <- ratematgALL[indedges]

  
  evalnetres <- new.env()  # prepare output
  assign("aveMImat", mat, envir=evalnetres)
  assign("rmat", ratematgALL, envir=evalnetres)
  assign("tpre", TPReALL, envir=evalnetres)
  assign("fnre", FNReALL, envir=evalnetres)
  assign("Fscores", Fscores, envir=evalnetres)
  assign("Ithr", Ithreshold, envir=evalnetres)
  assign("net", net, envir=evalnetres)
  
evalnetres
}



