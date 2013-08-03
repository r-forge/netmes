\name{readNet}
\alias{readNet}

\title{Fetching true network from .sif file}
\usage{readNet(infile, G)}
\arguments{
  \item{infile}{ The path for the .sif file that has true interactions of the genes.}
  \item{G}{Number of genes.}
}
\value{
      \code{readNet} returns the genes' relation matrix where activator edges, repressor edges and absence of edges are represented by 1, -1 and 0 respectively. 
}
\description{
	  This function takes the sythetic .sif file, which has interactions of genes, as input and converts it into an adjanceny matrix. 
}
\details{
      It assigns 1 for 'ac' and	and -1 for 're' and then fills by zeros for the rest of the matrix. This function is for synthetic networks as the gene names are expilicit as numbered. Note that the matrix is not symmetric. An example of .sif file can be obtained from SynTRen where the underlying network has gene labels as 'numbers'.
}
\references{
Gokmen Altay and Frank Emmert-Streib. Revealing differences in gene network inference algorithms on the network level by ensemble methods. Bioinformatics, 2010, 26(14),1738-1744.

Frank Emmert-Streib, Gokmen Altay. Local network-based measures to assess the inferability of different regulatory networks.  IET Syst. Biol. 2010, Volume 4, Issue 4, p.277-288.

}
\seealso{}
\examples{

# infile <- "./netmes/data/networkDAG.sif"
# infile <- "/home/zeyneb/Desktop/netmes/data/networkDAG.sif"
infile <- "/extdata/networkDAG.sif"

G <- 100 # number of genes in the network
net <- readNet(infile, G)   
tnet <- abs(net) + abs(t(net)) # true network (make symmetric)
tnet <- 1*(tnet > 0)           # in case there were bi-directional connections

}
\keyword{misc}
