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

motifMI <- function(MImat, n){ # average mutual information 
  
  MI <- (MImat[n[1],n[2]] + MImat[n[2],n[3]] + MImat[n[3],n[1]])/3

MI
}

