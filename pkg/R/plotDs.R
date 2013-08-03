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

plotDs<- function(TPReALL, Sm, fig=TRUE){ #object of evalnetworks of evalnetworknames

L1 <- min(Sm) # Sm integer values
L2 <- max(Sm)
gem <- vector(mode = "numeric", length = (L2 - L1 + 1))
ges <- vector(mode = "numeric", length = (L2 - L1 + 1))
cc <- 1
for(i in L1:L2){
  inds <- which(Sm == i)
  if(length(inds) > 0){
    gem[cc] <- mean(TPReALL[inds])
    ges[cc] <- sd(TPReALL[inds])
  }
  cc <- cc + 1
}



if(fig==TRUE)
{

postscript(file="Ds.eps", horiz=FALSE, paper = "special", family = "Helvetica", encoding = "TeXtext.enc", width = 10.0, height = 10.0)   #%%%%%%%%%%%%%%%%%%%
par(mar = c(5, 5, 4, 2) )    # margin for plots
plot(c(L1:L2), gem, xlab = "Ds", ylab = "TPR", cex.lab = 2.2, cex.axis = 1.9, ylim = c(0,1))
errbar(c(L1:L2), gem, gem + ges, gem - ges, add = TRUE, xlab = "Ds", ylab = "TPR")
dev.off()

}

 ds <- new.env()  # prepare output
  assign("gem", gem, envir=ds)
  assign("ges", ges, envir=ds)


ds
}



