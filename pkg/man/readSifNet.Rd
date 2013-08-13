\name{readSifNet}
\alias{readSifNet}

\title{Fetching true biological network from .sif file}
\usage{readSifNet(infilenet, G)}
\arguments{
  \item{infilenet}{ The path for the .sif file that has true interactions of the genes.}
  \item{G}{Number of genes. }
}
\value{
      \code{readSifNet} returns an object of class '"res"' which is a list with components: 
  \item{net}{the genes' relation matrix where activator edges, repressor edges and absence of edges are represented by 1, -1 and 0 respectively.}
  \item{names}{the gene names corresponding to the row and column names of the adjanceny matrix.}

}
\description{
	  This function is for biological networks where the names of the genes are arbitrary and needs to be derived along with the interaction of the genes. It provides the adjanceny matrix and and the corresponding gene names of the matrix.
}
\details{
      It assigns 1 for 'ac' and	and -1 for 're' and then fills by zeros the rest of the matrix. Note that the matrix is not symmetric.
}
\references{
Gokmen Altay and Frank Emmert-Streib. Revealing differences in gene network inference algorithms on the network level by ensemble methods. Bioinformatics, 2010, 26(14),1738-1744.

Frank Emmert-Streib, Gokmen Altay. Local network-based measures to assess the inferability of different regulatory networks.  IET Syst. Biol. 2010, Volume 4, Issue 4, p.277-288.
}
\seealso{\code{\link{readSifNet}},   \code{\link{evalnetworks}},    \code{\link{evalnetworksNames}}}
\examples{

# infilenet <- "./netmes/data/networkEcoli.sif"

# G <- 100 # number of genes in the network
# resEC <- readSifNet(infilenet, G)   # true network (+1 edge/ activator - -1 edges/repr.
# net <- resEC$net
# tnet <- abs(resEC$net) + abs(t(resEC$net)) # true network (make symmetric)
# tnet <- 1*(tnet > 0)           # in case there were bi-directional connections

}
\keyword{misc}
