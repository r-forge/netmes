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

visnet <- function(net, names, mat, Sm=0){

  # net: true network (binary)
  # mat: contains TPRe
  # Sm: edge specific measure

  
  g1 <- graph.adjacency(abs(net), mode = "directed")

  if(length(names) == 1){
    names <- V(g1) + 1
  } # else use names
  V(g1)$names <- names
  V(g1)$label <- names

if(Sm != 0) g1 <- set.edge.attribute(g1, "label", value = Sm )

  g1 <- set.edge.attribute(g1, "label.dist", value = 0.9 )


  
#get.edge.attribute(g1, "label.dist")

  
  
  ce <- c()
  inde <- E(g1)
  for(i in 1:(length(E(g1))  ) ){
    m <- as.character(print(inde[i]))
    ind <- as.numeric(strsplit(m[2], "->")[[1]]) 
    
    if(mat[ind[1], ind[2]] <= 0.25){
      ce <- append(ce, "red", after = length(ce))
    }
    if( (mat[ind[1], ind[2]] <= 0.5) & (mat[ind[1], ind[2]] > 0.25) ){
      ce <- append(ce, "green", after = length(ce))
    }
    if( (mat[ind[1], ind[2]] <= 0.75) & (mat[ind[1], ind[2]] > 0.5) ){
      ce <- append(ce, "blue", after = length(ce))
    }
    if( mat[ind[1], ind[2]] > 0.75 ){
      ce <- append(ce, "black", after = length(ce))
    }
   
  
  }


  g1 <- set.edge.attribute(g1, "color", index = E(g1), ce )
  tkplot(g1, canvas.width = 1000, canvas.height = 1000)

}



