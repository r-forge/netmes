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

TPRall<- function(res){ #object of evalnetworks or evalnetworknames

net <- res$net

#MImat <- mat  # MI matrix used (not normalized)
MImat <- res$aveMImat  # MI matrix used (not normalized)  

ratematgALL <- res$rmat
mat <- res$aveMImat

G <- nrow(net)


Ne <- sum(abs(net))
MI <- vector(mode = "numeric", length = Ne)
TPReALL <- vector(mode = "numeric", length = Ne)
TPReALLac <- c()
TPReALLre <- c()

# MI - mutual information for edges
cc <- 0
indbad <- new.env()
for(i in 1:G){
  aux <- which(net[i,] != 0)
  if(length(aux) > 0){
    for(j in 1:length(aux)){
      cc <- cc + 1
      MI[cc] <- MImat[i,aux[j]]
      TPReALL[cc] <- ratematgALL[i,aux[j]]
      if(TPReALL[cc] < 0.05){
         assign(as.character(eval(cc)), c(i,aux[j]), envir = indbad)
      }
      
      if(net[i,aux[j]] == 1){
        TPReALLac <- append(TPReALLac, ratematgALL[i,aux[j]], after = length(TPReALLac))
      }
      if(net[i,aux[j]] == -1){
        TPReALLre <- append(TPReALLre, ratematgALL[i,aux[j]], after = length(TPReALLre))
      }
    }
  }
}

 TPRs <- new.env()  # prepare output
  assign("TPReALL", TPReALL, envir=TPRs)
  assign("TPReALLac", TPReALLac, envir=TPRs)
  assign("TPReALLre", TPReALLre, envir=TPRs)

TPRs
}



