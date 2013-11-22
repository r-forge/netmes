\name{TPRall}
\alias{TPRall}

\title{TPR for all of the re or ac types of edges}
\usage{TPRall(res)}
\arguments{
  \item{res}{Object that is the result of evalnetworks or evalnetworknames.}
  
}
\value{
      \code{TPRall} returns an object of class '"TPRs"' which is a list with components: 
  \item{TPReALL}{TPR of all edges.}
  \item{TPReALLac}{TPR of ac types of edges.}
  \item{TPReALLre}{TPR of re types of edges.}

}
\description{
	  This function takes the object that is the result of evalnetworks or evalnetworknames. It results the True positive rate (TPR) of edges for respressor and activator types.
}
\details{
      It only considers ac and re edges.
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
# This is an example path. Change this path wrt the path of data in your computer.

# net <- readNet(infile, infilepath)   # true network (+1 edge/ activator - -1 edges/repr.)


# res <- evalnetworks(E, net, infilepath)

# TPRs <- TPRall(res)  #also save the fig

}
\keyword{misc}



