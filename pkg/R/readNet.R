#This file belongs to
#netmes: NETMES, <http://www.bioconductor.org/packages/release/Software.html>
#This R package allows local network based assesments for
##inferring regulatory networks from expression data.
## Copyright (C) August 2009 Frank Emmert-Streib and Gokmen Altay 
##<f.emmert-streib@qub.ac.uk>
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU LESSER GENERAL PUBLIC LICENSE
## as published by the Free Software Foundation; either version 3
## of the License, or any later version.
##
## This program is distributed WITHOUT ANY WARRANTY; 
## You can get a copy of the GNU LESSER GENERAL PUBLIC LICENSE
## from
## http://www.r-project.org/Licenses/LGPL-3

readNet <- function(infile, infilepath){

  file.adr <- paste(infilepath, "1/mim.exp", sep = "");
  cc <-read.table(file.adr)
  G <- dim(cc)[1]  # number of genes is assigned inside the local function

  dat <- matrix(scan(file = infile, what = "character", skip = 0), ncol = 3, byrow = TRUE)
  N <- dim(dat)[1]

  net <- matrix(0, nrow = G, ncol = G)
  for(i in 1:N){
    aux <- dat[i,]
    aux1 <- as.numeric(unlist(strsplit(aux[1], "gene_"))[2]) + 1
    aux2 <- as.numeric(unlist(strsplit(aux[3], "gene_"))[2]) + 1
    if(aux[2] == "ac"){
      net[aux1, aux2] <- 1
    }
    else{
      net[aux1, aux2] <- -1
    }
  }

  diag(net) <- c(1:G)*0  # delete self-loops if existing
  
net
}


