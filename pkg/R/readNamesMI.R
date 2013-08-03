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

readNamesMI <- function(infile, G, namesSif){ 

  names <- scan(file = infile, what = "character", skip = 0, nlines = 1)
  N <- length(names)
  dat <- matrix(scan(file = infile, what = "character", skip = 1), ncol = (N+1), nrow = N, byrow = TRUE)

  mat <- matrix(0, ncol = G, nrow = G)
  code <- vector(mode = "numeric", length = G)
  for(i in 1:N){  # sort rows
    aux <- dat[i,]
    ind <- which(namesSif == aux[1])
    code[ind] <- which(names == aux[1])
    mat[ind,] <- as.numeric(aux[2:(N+1)])
  }
  matnew <- mat[,code]  # sort columns
  
matnew

}

