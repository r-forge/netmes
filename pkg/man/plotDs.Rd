\name{plotDs}
\alias{plotDs}

\title{TPR for all of the re or ac types of edges}
\usage{plotDs(TPReALL, Sm, fig=TRUE)}
\arguments{
  \item{TPReALL}{Object that is the result of TPRall function.}
  \item{Sm}{Sum of degrees of nodes.}
  \item{fig}{ By default it saves the figures of the results. If not wanted then set it FALSE.}

}
\value{
      \code{plotDs} returns an object of class '"ds"' which is a list with components: 
  \item{gem}{Mean values of the degree of nodes.}
  \item{ges}{Standard deviation values of the degree of nodes.}

  \item{Ds.eps}{(If fig is set as TRUE or by default) Saves the figure Ds.eps to the current directory where x-axis sum of degree values, Ds, and y-axis is for corresponding TPR values.}


}
\description{
	  This function takes the object that is the result of TPRall function and Sm and computes Ds values and plot it.
}
\details{
     The influence of the degree of nodes on the inferability of edges is computed and plotted as Ds.eps.
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

# net <- readNet(infile, infilepath) # true network (+1 edge/ activator - -1 edges/repr.)

# Sm <- sumdeg(net)

# res <- evalnetworks(E, net, infilepath)

# TPRs <- TPRall(res)  #also save the fig

# Ds <-plotDs(TPRs$TPReALL, Sm)

}
\keyword{misc}



