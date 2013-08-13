\name{visnet}
\alias{visnet}

\title{Plotting the inferability of a network with colored edges}
\usage{visnet(net, names, mat, Sm)}
\arguments{
  \item{net}{Adjacency matrix of genes (not symmetric).}
  \item{names}{Names of the genes.}
  \item{mat}{Averaged mutual information matrix for the inferred edges of the overall experiments.}
  \item{Sm}{ Edge specific measure, which is the sum of outgoing edges of the regulating gene and the incoming edges of the regulated gene for an edge between them. Do not enter this parameter if you do not want it to be labelled on each edge in the network.}

}
\value{
      \code{visnet} returns graphical representation of genes with the edge specific measure labeled on each edge and the edges are coloured with four colour (red, green, blue and black from worst to best respectively) regardin the average mutual information of the edge.. 
}
\description{
	  This function takes the interaction file, names of genes, and averaged mutual information matrix for the inferred edges as input and outputs a coloured edge graph for the true gene network.  
}
\details{
      It uses igraph software to plot the directed and coloured true gene network.
}
\references{
Gokmen Altay and Frank Emmert-Streib. Revealing differences in gene network inference algorithms on the network level by ensemble methods. Bioinformatics, 2010, 26(14),1738-1744.

Frank Emmert-Streib, Gokmen Altay. Local network-based measures to assess the inferability of different regulatory networks.  IET Syst. Biol. 2010, Volume 4, Issue 4, p.277-288.
}
\seealso{\code{\link{readSifNet}},   \code{\link{evalnetworks}},    \code{\link{evalnetworksNames}}}
\examples{

# G <- 100 # number of genes in the network
# E <- 5 # number of experiments (data sets)

# infile <- "./netmes/data/networkDAG.sif"

# net <- readNet(infile, G)   # true network (+1 edge/ activator - -1 edges/repr.)

# infilepath <- "./netmes/data/syn/"

#This is an example path. Change this path wrt the path of data in your computer.

# res <- evalnetworks(E, G, net, infilepath)

# Sm <- sumdeg(net)

# visnet(net, c("no"), res$rmat, Sm)

# Sm is optional if wished to be labelled on edges.
# Adjust the plot by selecting "Layout" -> "Fruchterman-Reingold"

# visnet(net, resEC$names, res$rmat, Sm)  # for biological networks-resEC$names containes gene names
# This can be used with evalnetworksNames funtion.


}
\keyword{misc}
