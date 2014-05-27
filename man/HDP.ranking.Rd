\name{HDP.ranking}

\alias{HDP.ranking}

\title{High Directed Preorder Ranking of MST}

\description{Rank nodes in an object of class \code{igraph} (see package 
\code{\link{igraph}} for the definition of class \code{igraph}) containing a 
minimum spanning tree (MST) according to the High Directed Preorder traversal 
of the tree.}

\usage{HDP.ranking(mst)}

\arguments{
  \item{mst}{object of class \code{igraph} containing the minimum 
spanning tree.}
}

\details{Rank nodes in an object of class \code{igraph} (see package 
\code{\link{igraph}}) containing a minimum spanning tree (MST). The MST is 
rooted at a node with the largest geodesic distance and the rest of the nodes 
are ranked according to the high directed preorder (HDP) traversal of the tree 
(Friedman and Rafsky 1979).}

\value{Numeric vector giving the node ranks according to HDP traversal 
of the MST.}

\references{
Rahmatallah Y., Emmert-Streib F. and Glazko G. (2012) Gene set analysis for 
self-contained tests: complex null and specific alternative hypotheses. 
Bioinformatics \bold{28}, 3073--3080.

Friedman J. and Rafsky L. (1979) Multivariate generalization of the 
Wald-Wolfowitz and Smirnov two-sample tests. Ann. Stat. \bold{7}, 697--717.
}

\author{Yasir Rahmatallah and Galina Glazko}

\note{This function does not work properly if there is any node in the MST 
with more than 26 links. However, this situation is almost impossible for a 
dataset composed of a few hundreds or less of samples.}

\seealso{\code{\link{KStest}}.}

\examples{
## generate random data using normal distribution
## generate 20 features in 20 samples
object <- matrix(rnorm(400),20,20)
objt <- aperm(object, c(2,1))
## calculate the weight matrix
Wmat <- as.matrix(dist(objt, method = "euclidean", diag = TRUE, upper = TRUE, p = 2))
## create a weighted undirectional graph from the weight matrix
gr <- graph.adjacency(Wmat, weighted = TRUE, mode = "undirected")
## find the minimum spanning tree
mst <- minimum.spanning.tree(gr)
ranks <- HDP.ranking(mst)
plot(mst)
}

\keyword{arith}
