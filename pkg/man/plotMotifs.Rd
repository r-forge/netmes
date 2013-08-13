\name{plotMotifs}
\alias{plotMotifs}

\title{Print to screen in latex format of motif statistics}
\usage{plotMotifs(res)}
\arguments{
  \item{res}{an object of class '"res"' which is a list with components: # *** .}

}
\value{
      \code{plotMotifs} returns the table in latex format.
}
\description{
	This function takes the object res (that includes all the explored statistics) and outputs the table in latex format that includes standard deviation, mean and the sum of edges.# ***
}
\details{

}
\references{
Gokmen Altay and Frank Emmert-Streib. Revealing differences in gene network inference algorithms on the network level by ensemble methods. Bioinformatics, 2010, 26(14),1738-1744.

Frank Emmert-Streib, Gokmen Altay. Local network-based measures to assess the inferability of different regulatory networks.  IET Syst. Biol. 2010, Volume 4, Issue 4, p.277-288.
}
\seealso{\code{\link{motifcount3}},    \code{\link{evalnetworksNames}}}
\examples{

# G <- 100 # number of genes in the network
# E <- 5 # number of experiments (data sets)

# infile <- "./netmes/data/networkDAG.sif"

# net <- readNet(infile, G)   # true network (+1 edge/ activator - -1 edges/repr.)

# infilepath <- "./netmes/data/syn/"

# This is an example path. Change this path wrt the path of data in your computer.

# res <- evalnetworks(E, G, net, infilepath)

# res2 <- motifcount3(net, res$rmat, res$aveMImat, G)

# plotMotifs(res2)


}
\keyword{misc}
