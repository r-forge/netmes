\name{motifcount3}
\alias{motifcount3}

\title{Motif statistics for regulatory networks}
\usage{motifcount3(net, rmat, MImat, G)}
\arguments{
  \item{net}{adjacency matrix of genes (not symmetric).}
  \item{rmat}{averaged mutual information matrix for the inferred edges of the overall experiments. }
  \item{MImat}{averaged mutual information matrix of the overall experiments.}
  \item{G}{ number of genes in the network.}
}
\value{
      \code{motifcount3} returns an object of class '"res"' which is a list with components: 
  \item{moty}{number of motifs.}
  \item{motyprob}{probability to observe a motif.}
  \item{motyMI}{average mutual information of motifs.}
  \item{emotnodes}{ node indices for motifs.}

}
\description{
	  This function takes the adjacency matrix, averaged mutual information matrix for the inferred edges, averaged mutual information matrix and number of genes as input. It returns   

 the data sets folder path and number of mutual information steps for the therhold mutual information. Considering 3 gene case, it computes number of motifs, probability to observe a motif, average mutual information of motifs and  node indices for motifs.
}
\details{
      It only consider 3 genes case and search for all the possible motifs with 3 genes.
}
\references{
Gokmen Altay and Frank Emmert-Streib. Revealing differences in gene network inference algorithms on the network level by ensemble methods. Bioinformatics, 2010, 26(14),1738-1744.

Frank Emmert-Streib, Gokmen Altay. Local network-based measures to assess the inferability of different regulatory networks.  IET Syst. Biol. 2010, Volume 4, Issue 4, p.277-288.
}
\seealso{\code{\link{plotMotifs}},   \code{\link{evalnetworks}},    \code{\link{evalnetworksNames}}}
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

res2 <- motifcount3(net, res$rmat, res$aveMImat, G)

plotMotifs(res2)

# if True Reconstruction Rate wished to be plotted
envm <- res2$motyprob
y <- c()
xg <- c()
i <- 2
aux <- get(as.character(i), envir = envm)
y <- append(y, aux, after = length(y))
xg <- append(xg, rep(1, length(aux)), after = length(xg))
boxplot(y ~ xg, names = c("DAG"), ylab = expression(p), cex.lab = 1.7, cex.axis = 1.6, col = c("darkred") )



}
\keyword{misc}



