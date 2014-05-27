\name{findMST2}

\alias{findMST2}

\title{Union of the First and Second Minimum Spanning Trees}

\description{Find the union of the first and second minimum spanning trees.}

\usage{findMST2(object, cor.method="pearson", min.sd=1e-3, return.MST2only=FALSE)}

\arguments{
  \item{object}{a numeric matrix with columns and rows respectively 
corresponding to samples and features.}
  \item{cor.method}{a character string indicating which correlation 
coefficient is to be computed.
Possible values are \dQuote{\code{pearson}} (default), \dQuote{\code{spearman}} 
and \dQuote{\code{kendall}}.}
  \item{min.sd}{the minimum allowed standard deviation for any feature. 
If any feature has a standard deviation smaller than \code{min.sd} the 
execution stops and an error message is returned.}
  \item{return.MST2only}{logical. If \code{FALSE} (default) a list of length 
three containing objects of class 
\code{igraph} is returned. The first and second are the first and second MSTs 
respectively. The third is the union of the first and second, MST2. If 
\code{TRUE}, an object of class \code{igraph} containing the MST2 is returned.}
}

\details{This function produces the union of the first and second minimum 
spanning trees (MSTs) as an object of class \code{igraph} (check package 
\code{\link{igraph}} for details). It can as well return the first and 
second minimum spanning trees when \code{return.MST2only} is \code{FALSE} 
(default). It starts by calculating the correlation (coexpression) matrix and 
using it to obtain a weighting matrix for a complete graph using the equation 
\eqn{w_{ij} = 1 - |r_{ij}|} where \eqn{r_{ij}} is the correlation between 
features \eqn{i} and \eqn{j} and \eqn{w_{ij}} is the weight of the link between 
vertices (nodes) \eqn{i} and \eqn{j} in the graph \eqn{G(V,E)}.

For the graph \eqn{G(V,E)} where V is the set of vertices and E is the set of 
edges, the first MST is defined as the acyclic subset \eqn{T_{1} \subseteq E} 
that connects all vertices in V and whose total length
\eqn{\sum_{i,j \in T_{1}} d(v_{i},v_{j})} is minimal 
(Rahmatallah et. al. 2014). The second MST is defined as the MST of the 
reduced graph \eqn{G(V,E-T_{1})}. The union of the first and second MSTs is 
denoted as MST2.

It was shown in Rahmatallah et. al. 2014 that MST2 can be used as a graphical 
visualization tool to highlight the most highly correlated genes in the 
correlation network. A gene that is highly correlated with all the other genes 
tends to occupy a central position and has a relatively high degree in the MST2 
because the shortest paths connecting the vertices of the first and second MSTs 
tend to pass through the vertex corresponding to this gene. In contrast, a gene 
with low intergene correlations most likely occupies a non-central position in 
the MST2 and has a degree of 2.

In rare cases, a feature may have a constant or nearly constant level across 
the samples. This results in a zero or a tiny standard deviation. Such case 
produces an error in command \code{cor} used to compute the correlations 
between features. To avoid this situation, standard deviations are checked in 
advance and if any is found below the minimum limit \code{min.sd} 
(default is \code{1e-3}), the execution stops and an error message is returned 
indicating the the number of feature causing the problem (if only one the 
index of that feature is given too). 
}

\value{
If \code{return.MST2only} is \code{TRUE}, \code{findMST2} returns an object 
of class \code{igraph} containing the MST2. If \code{return.MST2only} is 
\code{FALSE} (default), \code{findMST2} returns a list with the following 
components:
\item{MST2}{an object of class \code{igraph} containing the union of the 
first MST and second MST.}
\item{first.mst}{an object of class \code{igraph} containing the first MST.}
\item{second.mst}{an object of class \code{igraph} containing the second MST.}
}

\references{
Rahmatallah Y., Emmert-Streib F. and Glazko G. (2014) Gene sets net 
correlations analysis (GSNCA): a multivariate differential coexpression test 
for gene sets. Bioinformatics \bold{30}, 360--368.
}

\author{Yasir Rahmatallah and Galina Glazko}

\seealso{\code{\link{GSNCAtest}}, \code{\link{plotMST2.pathway}}.}

\examples{
## generate a dataset of 20 features and 20 samples
## use multivariate normal distribution with different covariance matrices
library(MASS)
ngenes <- 20
nsamples <- 20
zero_vector <- array(0,c(1,ngenes))
## create a covariance matrix with high off-diagonal elements
## for the first 5 features and low for the remaining 15 features
cov_mtrx <- diag(ngenes)
cov_mtrx[!diag(ngenes)] <- 0.1
mask <- diag(ngenes/4)
mask[!diag(ngenes/4)] <- 0.6
cov_mtrx[1:(ngenes/4),1:(ngenes/4)] <- mask
gp <- mvrnorm(nsamples, zero_vector, cov_mtrx)
dataset <- aperm(gp, c(2,1))
## findMST2 returns a list of length 3
## trees[[1]] is an object of class igraph containing the MST2
trees <- findMST2(dataset)
}

\keyword{dplot}
