\name{findMST2.PPI}

\alias{findMST2.PPI}

\title{Union of the First and Second Minimum Spanning Trees for PPI Networks}

\description{Find the union of the first and second minimum spanning trees for 
protein-protein interaction (PPI) networks.}

\usage{findMST2.PPI(object, return.MST2only=TRUE)}

\arguments{
  \item{object}{an object of class \code{igraph} representing the PPI network.}
  \item{return.MST2only}{logical. If \code{TRUE} (default), an object of 
class \code{igraph} containing the MST2 is returned. If \code{FALSE}, a 
list of length three containing objects of class \code{igraph} is returned. 
The first and second objects are the first and second MSTs, respectively. 
The third is the union of the first and second, MST2.}
}

\details{This function produces the union of the first and second minimum 
spanning trees (MSTs) as an \code{igraph} object (check package 
\code{\link{igraph}} for details). It can as well return the first and 
second minimum spanning trees when \code{return.MST2only} is \code{FALSE}.

For the graph \eqn{G(V,E)} where V is the set of vertices and E is the set of 
edges, the first MST is defined as the acyclic subset \eqn{T_{1} \subseteq E} 
that connects all vertices in V and whose total length
\eqn{\sum_{i,j \in T_{1}} d(v_{i},v_{j})} is minimal 
(Rahmatallah et. al. 2014). The second MST is defined as the MST of the 
reduced graph \eqn{G(V,E-T_{1})}. The union of the first and second MSTs is 
denoted as MST2.

It was shown in Zybailov et. al. 2016 that MST2 can be informative as a 
graphical visualization tool in deciphering the properties of protein-protein 
interaction (PPI) networks by highlighting the minimum set of essential 
interactions among proteins. Most influential proteins with many interactions 
tend to occupy central position and have relatively high connectivity degree 
in the MST2 because the shortest paths connecting the vertices of the first 
and second MSTs tend to pass through the verteces corresponding to these proteins. 
In contrast, proteins with few interactions most likely occupy non-central 
positions in the MST2 and have a degree of 2.
}

\value{
If \code{return.MST2only=TRUE} (default), function \code{findMST2.PPI} 
returns an object of class \code{igraph} representing the MST2. If 
\code{return.MST2only=FALSE}, function \code{findMST2.PPI} returns a list 
of length 3 with the following components:
\item{MST2}{an object of class \code{igraph} containing the union of the 
first and second MSTs.}
\item{first.mst}{an object of class \code{igraph} containing the first MST.}
\item{second.mst}{an object of class \code{igraph} containing the second MST.}
}

\references{
Zybailov B., Byrd A., Glazko G., Rahmatallah Y. and Raney K. (2016) 
Protein-protein interaction analysis for functional characterization of 
helicases. Methods, \bold{108}, 56--64.
}

\author{Yasir Rahmatallah and Galina Glazko}

\seealso{\code{\link{GSNCAtest}}, \code{\link{plotMST2.pathway}}.}

\examples{
## generate a random undirected graph with power-law 
## distribution degree where minimum degree is 4 and 
## maximum degree is 100
set.seed(123)
degs <- sample(c(4:100), 100, replace=TRUE, prob=c(4:100)^-2)
if(floor(sum(degs)/2) != (sum(degs)/2)) degs[1] <- degs[1] + 1
randomGraph <- sample_degseq(degs, method="vl")
## find MST2 of the random graph and highlight vertices 
## with degree greater than 10 with red color
mst2.ppi <- findMST2.PPI(object=randomGraph, return.MST2only=TRUE)
degs <- degree(mst2.ppi)
ind <- which(degs > 10)
V(mst2.ppi)$color <- "yellow"
V(mst2.ppi)$color[ind] <- "red"
}

\keyword{dplot}
\keyword{graphs}
