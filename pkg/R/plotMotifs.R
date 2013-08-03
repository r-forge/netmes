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

plotMotifs <- function(res){

  moty <- res$moty
  motyprob <- res$motyprob
  motyMI <- res$motyMI

  aux <- paste("& $\\#m$ &", moty[1], "&", moty[2], "&", moty[3], "&", moty[4], "&", moty[5], "\\", sep = " ")
  cat(aux, "\n", sep = "\\")

  mprob <- vector(mode = "numeric", length = 5)
  sdprob <- vector(mode = "numeric", length = 5)
  for(i in 1:5){
    aux <- get(as.character(i), envir = motyprob)
    if(length(aux) > 0){
      mprob[i] <- mean(aux)
      sdprob[i] <- sd(aux)
    }
  }

  mprob <- format(mprob, digits = 3)
  sdprob <- format(sdprob, digits = 2)
  aux <- paste("& $\\bar{p}$ &", mprob[1], "&", mprob[2], "&",mprob[3], "&",mprob[4], "&",mprob[5], "\\", sep = " ")
  cat(aux, "\n", sep = "\\")

  aux <- paste("& $\\sigma(\\bar{p})$", "&", sdprob[1], "&", sdprob[2], "&",sdprob[3], "&",sdprob[4], "&",sdprob[5], "\\", sep = " ")
  cat(aux, "\n", sep = "\\")

  


  mMI <- vector(mode = "numeric", length = 5)
  sdMI <- vector(mode = "numeric", length = 5)
  for(i in 1:5){
    aux <- get(as.character(i), envir = motyMI)
    if(length(aux) > 0){
      mMI[i] <- mean(aux)
      sdMI[i] <- sd(aux)
    }
  }
  mMI <- format(mMI, digits = 2)
  sdMI <- format(sdMI, digits = 3)
  aux <- paste("& $\\bar{I}$ &", mMI[1], "&", mMI[2], "&",mMI[3], "&",mMI[4], "&",mMI[5], "\\", sep = " ")
  cat(aux, "\n", sep = "\\")

  aux <- paste("& $\\sigma(\\bar{I})$ &", sdMI[1], "&", sdMI[2], "&",sdMI[3], "&",sdMI[4], "&",sdMI[5], "\\", sep = " ")
  cat(aux, "\n", sep = "\\")
  
}


