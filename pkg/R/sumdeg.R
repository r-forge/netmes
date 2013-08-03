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

sumdeg <- function(net){ #non-symmetric true adjanjency matrix

 Ne <- sum(abs(net))

G <- nrow(net)
# sum of degrees
Sm <- vector(mode = "numeric", length = Ne)
cc <- 0
for(i in 1:G){
  aux <- which(net[i,] != 0)
  if(length(aux) > 0){
    for(j in 1:length(aux)){
      cc <- cc + 1
      Sm[cc] <- sum(abs(net[i,])) + sum(abs(net[,aux[j]])) #+ sum(abs(net[aux[j],])) # + EB[cc] 
      #Sm[cc] <- sum(abs(net[,i])) + sum(abs(net[,aux[j]]))
      #Sm[cc] <- - sum(EBmat[i,])  + sum(EBmat[,i]) + sum(EBmat[aux[j],]) - sum(EBmat[,aux[j]])  - EB[cc]
      if(Sm[cc] > 4000)print(c(i, aux[j]))  #????? why 4000
    }
  }
}

Sm
}



