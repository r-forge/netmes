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

motifrate <- function(rmat, net, n){ # prob to observe a motif
  prob <- 0.0
  if((net[n[1], n[2]] != 0) | (net[n[2], n[1]] != 0)){
    prob <- prob + rmat[n[1], n[2]]  # tpr
  }
  else{
    prob <- prob + (1 - rmat[n[1], n[2]]) # tnr
  }
  if((net[n[2], n[3]] != 0) | (net[n[3], n[2]] != 0)){
    prob <- prob + rmat[n[2], n[3]]  # tpr
  }
  else{
    prob <- prob + (1 - rmat[n[2], n[3]]) # tnr
  }
  if((net[n[3], n[1]] != 0) | (net[n[1], n[3]] != 0)){
    prob <- prob + rmat[n[3], n[1]]  # tpr
  }
  else{
    prob <- prob + (1 - rmat[n[3], n[1]]) # tnr
  }

prob <- prob/3

prob
}



