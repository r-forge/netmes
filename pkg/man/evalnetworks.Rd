\name{evalnetworks}
\alias{evalnetworks}

\title{Computing all the statistics from data of sythetic network}
\usage{evalnetworks(E, G, net, infilepath)}
\arguments{
  \item{E}{number of experiments.}
  \item{G}{Number of genes. }
  \item{net}{Adjacency matrix of genes (not symmetric).}
  \item{infilepath}{The path where all the data files kept.}

}
\value{
      \code{evalnetworks} returns an object of class '"res"' which is a list with components: 
  \item{aveMImat}{averaged mutual information matrix of the overall experiments.}
  \item{rmat}{averaged mutual information matrix for the inferred edges of the overall experiments.}
  \item{tpre}{All the averaged mutual information values for the inferred edges.}
  \item{fnre}{All the averaged mutual information values for the non-inferred edges.}
  \item{Fscores}{Maximum F-score values for all the experiments.}
  \item{Ithr}{The threshold mutual information values, which gives maximum F-score, for all the experiments.}
  \item{net}{Adjacency matrix of genes (not symmetric).}


}
\description{
	  This function takes the adjacency matrix, the data sets folder path and number of mutual information steps for the therhold mutual information. It returns averaged mutual information matrix, averaged mutual information matrix for the inferred edges, the averaged mutual information values for the inferred edges, the averaged mutual information values for the non-inferred edges, maximum F-scores, threshold mutual information values that gives maximum F-score, and also the same input adjacency matrix of genes.      
}
\details{
      It reads all the mutual information matrices from the folders and take the average and then derive the other parameters from it.
}
\references{
Gokmen Altay and Frank Emmert-Streib. Revealing differences in gene network inference algorithms on the network level by ensemble methods. Bioinformatics, 2010, 26(14),1738-1744.

Frank Emmert-Streib, Gokmen Altay. Local network-based measures to assess the inferability of different regulatory networks.  IET Syst. Biol. 2010, Volume 4, Issue 4, p.277-288.

}
\seealso{ \code{\link{readNet}} }

\examples{

G <- 100 # number of genes in the network
E <- 5 # number of experiments (data sets)

# infile <- "./netmes/data/networkDAG.sif"
# infile <- "/home/zeyneb/Desktop/netmes/data/networkDAG.sif"
infile <- "/extdata/networkDAG.sif"

net <- readNet(infile, G)   # true network (+1 edge/ activator - -1 edges/repr.)

# infilepath <- "./netmes/data/syn/"
# infilepath <- "/home/zeyneb/Desktop/netmes/data/syn/" 
infilepath <- "/extdata/syn/"

#This is an example path. Change this path wrt the path of data in your computer.

res <- evalnetworks(E, G, net, infilepath)

#boxplot(res$Fscores) #optional

#save(res, file = "resDAG")  # is used for saving boxplots 

}
\keyword{misc}



