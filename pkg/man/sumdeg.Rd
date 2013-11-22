\name{sumdeg}
\alias{sumdeg}

\title{Sum of degree calculation for each node}
\usage{sumdeg(net)}
\arguments{
  \item{net}{adjacency matrix of genes (not symmetric).}
  
}
\value{
   \code{sumdeg} returns sum of degrees of noes as Sm 


}
\description{
	  This function takes the adjacency matrix, and calculate sum of degrees.
}
\details{
      It only considers all edges.
}
\references{
Gokmen Altay and Frank Emmert-Streib. Revealing differences in gene network inference algorithms on the network level by ensemble methods. Bioinformatics, 2010, 26(14),1738-1744.

Frank Emmert-Streib, Gokmen Altay. Local network-based measures to assess the inferability of different regulatory networks.  IET Syst. Biol. 2010, Volume 4, Issue 4, p.277-288.
}
\seealso{\code{\link{plotMotifs}},   \code{\link{evalnetworks}},    \code{\link{evalnetworksNames}}}
\examples{

# E <- 5 # number of experiments (data sets)

# infile <- "./netmes/data/networkDAG.sif"
# infilepath <- "./netmes/data/syn/"

# net <- readNet(infile, infilepath) # true network (+1 edge/ activator - -1 edges/repr.)

# Sm <- sumdeg(net)


}
\keyword{misc}



