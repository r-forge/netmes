\name{evalnetworksNames}
\alias{evalnetworksNames}

\title{Computing all the statistics from data of biological network}
\usage{evalnetworksNames(E, net, infilepath, namesSif)}
\arguments{
  \item{E}{number of experiments.}
  \item{net}{Adjacency matrix of genes (not symmetric).}
  \item{infilepath}{The path of the directory, where all the data files kept.}
  \item{namesSif}{names of the genes corresponding to the adjacency matrix.}

}
\value{
      \code{evalnetworksNames} returns an object of class '"res"' which is a list with components: 
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
      It reads all the mutual information matrices from the folders and take the average and then derive the other parameters from it. This function is designed for biological networks where the gene names need to be derived along with the gene interactions from the files.
}
\references{
Gokmen Altay and Frank Emmert-Streib. Revealing differences in gene network inference algorithms on the network level by ensemble methods. Bioinformatics, 2010, 26(14),1738-1744.

Frank Emmert-Streib, Gokmen Altay. Local network-based measures to assess the inferability of different regulatory networks.  IET Syst. Biol. 2010, Volume 4, Issue 4, p.277-288.

}
\seealso{\code{\link{readSifNet}},   \code{\link{readNet}},    \code{\link{evalnetworksNames}}}
\examples{

# E <- 5 # number of experiments

# infilenet <- "./netmes/data/networkEcoli.sif"
# infilepath <- "./netmes/data/bio/"
# This is an example path. Change this path wrt the path of data in your computer.

# resEC <- readSifNet(infilenet, infilepath)   # true network (+1 edge/ activator - -1 edges/repr.
# net <- resEC$net

# res <- evalnetworksNames(E, resEC$net, infilepath, resEC$names)

# boxplot(res$Fscores)

# save(res, file = "res_EColi")
# if desired for later analysis.

}
\keyword{misc}



