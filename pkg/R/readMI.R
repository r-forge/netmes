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


readMI <- function(infile, infilepath){

  file.adr <- paste(infilepath, "1/mim.exp", sep = "");
  cc <-read.table(file.adr)
  G <- dim(cc)[1]  # number of genes is assigned inside the local function

  names <- scan(file = infile, what = "character", skip = 0, nlines = 1)
  N <- length(names)
  dat <- matrix(scan(file = infile, what = "character", skip = 1), ncol = (N+1), nrow = N, byrow = TRUE)
  
  code <- vector(mode = "numeric", length = G)
  decode <- vector(mode = "numeric", length = G)
  for(i in 0:(G-1)){
    if(i < 10){
      mask <- paste("gene_00", as.character(i), sep = "")
    }
    else{
      mask <- paste("gene_0", as.character(i), sep = "")
    }
    ind <- regexpr(mask, names)
    ind <- which(ind == 1)
    code[(i+1)] <- ind
    decode[ind] <- i+1
  }
  
  mat <- matrix(as.numeric(dat[code,2:(G+1)]), ncol = G, nrow = G)
  matnew <- matrix(0, ncol = G, nrow = G)
  for(i in 1:G){
    ind <- which(mat[i,] > 0)
    val <- mat[i,ind]
    newind <- decode[ind]
    matnew[i,newind] <- val 
  }

  
matnew

}


