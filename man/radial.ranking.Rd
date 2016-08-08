\name{radial.ranking}

\alias{radial.ranking}

\title{Radial Ranking of MST}

\description{Rank vertices in an object of class \code{igraph} (see package 
\code{\link{igraph}} for the definition of class \code{igraph}) that consists 
of a minimum spanning tree (MST) or the union of multiple MSTs radially 
such that vertices with higher depth and distance from the centroid are 
given higher ranks.}

\usage{radial.ranking(object)}

\arguments{
  \item{object}{object of class \code{igraph} that consists of a minimum 
spanning tree or the union of multiple spanning trees.}
}

\details{Rank nodes in an object of class \code{igraph} (see package 
\code{\link{igraph}}) that consists of a minimum spanning tree (MST) or 
the union of multiple MSTs radially. The MST is rooted at the node of 
smallest geodesic distance (centroid) and nodes with largest depths 
from the root are assigned higher ranks. Hence, ranks are increasing 
radially from the root of the MST (Friedman and Rafsky 1979).}

\value{Numeric vector giving the radial node ranks in the MST or union of MSTs.}

\references{
Rahmatallah Y., Emmert-Streib F. and Glazko G. (2012) Gene set analysis for 
self-contained tests: complex null and specific alternative hypotheses. 
Bioinformatics \bold{28}, 3073--3080.

Friedman J. and Rafsky L. (1979) Multivariate generalization of the 
Wald-Wolfowitz and Smirnov two-sample tests. Ann. Stat. \bold{7}, 697--717.
}

\author{Yasir Rahmatallah and Galina Glazko}

\seealso{\code{\link{HDP.ranking}}, \code{\link{RKStest}}, \code{\link{RMDtest}}.}

\examples{
## generate random data using normal distribution
## generate 20 features in 20 samples
object <- matrix(rnorm(400),20,20)
objt <- aperm(object, c(2,1))
## calculate the weight matrix
Wmat <- as.matrix(dist(objt, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
## create a weighted undirectional graph from the weight matrix
gr <- graph_from_adjacency_matrix(Wmat, weighted = TRUE, mode = "undirected")
## find the minimum spanning tree
MST <- mst(gr)
radial.ranks <- radial.ranking(MST)
}

\keyword{arith}
