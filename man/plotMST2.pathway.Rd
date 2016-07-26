\name{plotMST2.pathway}

\alias{plotMST2.pathway}

\title{Plot MST2 for a pathway in two conditions}

\description{This is a wrapper function which uses function 
\code{\link{findMST2}} to find the union of the first and second minimum 
spanning trees (or MST2) of the correlation network for a feature set (pathway) 
under two conditions. It plots the MST2 of the correlation network of the 
feature set under both conditions side-by-side and highlights hub nodes to 
facilitate a visual comparison.}

\usage{plotMST2.pathway(object, group, name=NULL, cor.method="pearson", 
min.sd=1e-3, legend.size=1, leg.x=-0.8, leg.y=1.5, return.weights=FALSE, 
group1.name="Group 1", group2.name="Group 2", label.size=1, 
label.color="black", label.dist=0.5, vertex.size=8, vertex.label.font=1, 
edge.width=1)}

\arguments{
  \item{object}{a numeric matrix with columns and rows respectively 
corresponding to samples and features. Gene names are provided to this 
function as the rownames of this matrix.}
  \item{group}{a numeric vector indicating group associations for samples. 
Possible values are 1 and 2.}
  \item{name}{an optional character string giving the name of the feature set 
(gene set). If given, the name will be displayed at the top of the plot.}
  \item{cor.method}{a character string indicating which correlation 
coefficient is to be computed. Possible values are \dQuote{\code{pearson}} 
(default), \dQuote{\code{spearman}} and \dQuote{\code{kendall}}. Default 
value is \dQuote{\code{pearson}}.}
  \item{min.sd}{a numeric value indicating the minimum allowed standard 
deviation for any feature. If any feature has a standard deviation 
smaller than \code{min.sd} then the execution stops and an error message 
is returned. Default value is 1e-3.}
  \item{legend.size}{an optional numeric value controlling the relative font 
size of the legend to the default font size. Default is 1.}
  \item{leg.x}{a numeric value indicating the amount of horizontal 
shift of the legend box to allow better positioning in the plot.}
  \item{leg.y}{a numeric value indicating the amount of vertical 
shift of the legend box to allow better positioning in the plot.}
  \item{return.weights}{logical. Default value is FALSE. If the weight 
factors aasigned to the genes by the GSNCA method are desired, setting 
this parameter to TRUE returns the weight factors in a matrix with 2 columns 
(for class 1 and class 2) and number of rows equal to the number of genes 
in the gene set. If the \code{rownames} of \code{object} are provided, then 
they will be used as \code{rownames} for the returned matrix. If the 
\code{rownames} of \code{object} are abscent, node labels will be set to 
\code{as.character(c(1:nrow(object)))}.}
  \item{group1.name}{an optional character string to be presented as the 
given name for class 1 in the plot. Default value is \dQuote{\code{Group 1}}}
  \item{group2.name}{an optional character string to be presented as the 
given name for class 2 in the plot. Default value is \dQuote{\code{Group 2}}}
  \item{label.size}{a numeric value passed to argument vertex.label.cex in 
command plot.igraph to specify the vertex label size. Default value is 1.}
  \item{label.color}{a character string specifying the color of vertex 
labels. Default value is \dQuote{\code{black}}.}
  \item{label.dist}{a numeric value passed to argument vertex.label.dist in 
command plot.igraph to specify the distance between vertex labels and the 
centers of vertices. Default value is 0.5.}
  \item{vertex.size}{a numeric value passed to argument vertex.size in 
command plot.igraph to specify the vertex size. Default value is 8.}
  \item{vertex.label.font}{a numeric value passed to argument 
vertex.label.font in command plot.igraph to specify the used font type. 
Default value is 1.}
  \item{edge.width}{a numeric value passed to argument edge.width in 
command plot.igraph to specify the edge width in the plot.}
}

\details{This is a wrapper plotting function for the convenience of users. It 
uses function \code{\link{findMST2}} to find the union of the first and second 
minimum spanning trees (or MST2) of the correlation network for a feature set 
(pathway) under two conditions and plots them side-by-side. It also lists the 
hub nodes and their weight factors (\eqn{w}) under each condition (see 
Rahmatallah et. al. 2014 for details). The range in which weight factors fall 
is indicated by the node colors defined in the legend. Weight factor have 
values mostly ranging between 0.5 (low coexpression) and 1.5 (high 
coexpression). To allow the users more control over plotting parameters and 
to present different feature sets appropriately, two optional arguments were 
introduced: \code{legend.size} and \code{label.size}. Node lables will be the 
names of the features in the set, i.e. \code{rownames(object)}. If the 
\code{rownames} attribute is not set for \code{object}, node labels will be 
set to \code{as.character(c(1:nrow(object)))}.

The weight factors, inferred from the Gene Sets Net Correlations Analysis 
(GSNCA) method (see \code{\link{GSNCAtest}}), correlate to some extent with 
genes centralities in the MST2: genes with large weights are placed near the 
center of the MST2, and genes with small weights are placed on the periphery 
(Rahmatallah et. al. 2014). Adopting network terminology, a gene with the 
largest weight is a hub gene, coexpressed with most of the other genes in a 
pathway (see \code{\link{findMST2}}). Therefore, MST2 is a convenient 
graphical visualization tool to examine the pathways tested by the GSNCA 
method (see \code{\link{GSNCAtest}}).

The correlation (coexpression) network is obtained using the weight matrix 
\eqn{W} with elements \eqn{w_{ij} = 1 - |r_{ij}|} where \eqn{r_{ij}} is the 
correlation between features \eqn{i} and \eqn{j} and \eqn{w_{ij}} is the weight 
of the link between vertices (nodes) \eqn{i} and \eqn{j} in the network. The 
correlation coefficient used is indicated by the argument \code{cor.method} 
with three possible values: \dQuote{\code{pearson}} (default), 
\dQuote{\code{spearman}} and \dQuote{\code{kendall}}.

In some cases (especially for RNA-Seq count data), a feature (or more) may 
have a constant or nearly constant level across the samples in one or both 
conditions. This results in a zero or a tiny standard deviation. Such case 
produces an error in command \code{cor} used to compute the correlation 
coefficients between features. To avoid this situation, standard deviations 
are checked in advance and if any is found below the minimum limit 
\code{min.sd} (default is \code{1e-3}), the execution stops and an error 
message is returned indicating the number of feature causing the problem 
(if only one the index of that feature is given too).}

\references{
Rahmatallah Y., Emmert-Streib F. and Glazko G. (2014) Gene sets net 
correlations analysis (GSNCA): a multivariate differential coexpression test 
for gene sets. Bioinformatics \bold{30}, 360--368.
}

\author{Yasir Rahmatallah and Galina Glazko}

\note{This function is suitable for a feature set of roughly 80 features or 
less. It works for feature sets with larger number of features but the 
placements of nodes and their labels in the plot will be too crowded 
for a useful visual presentation.}

\seealso{\code{\link{findMST2}}, \code{\link{GSNCAtest}}.}

\examples{
## generate a feature set of length 20 in two conditions
## each condition has 20 samples
## use multivariate normal distribution with different covariance matrices
library(MASS)
ngenes <- 20
nsamples <- 40
zero_vector <- array(0,c(1,ngenes))
## create a covariance matrix with low off-diagonal elements
cov_mtrx1 <- diag(ngenes)
cov_mtrx1[!diag(ngenes)] <- 0.1
## create a covariance matrix with high off-diagonal elements
## for the first 5 features and low for the rest 15 features
cov_mtrx2 <- diag(ngenes)
cov_mtrx2[!diag(ngenes)] <- 0.1
mask <- diag(ngenes/4)
mask[!diag(ngenes/4)] <- 0.6
cov_mtrx2[1:(ngenes/4),1:(ngenes/4)] <- mask
gp1 <- mvrnorm((nsamples/2), zero_vector, cov_mtrx1)
gp2 <- mvrnorm((nsamples/2), zero_vector, cov_mtrx2)
gp <- rbind(gp1,gp2)
dataset <- aperm(gp, c(2,1))
## first 20 samples belong to group 1
## second 20 samples belong to group 2
## since rowname(object)=NULL, node labels will be automatically 
## set to as.character(c(1:nrow(object))) 
plotMST2.pathway(object=dataset, group=c(rep(1,20),rep(2,20)), 
name="Example Pathway") 
}

\keyword{dplot}
\keyword{graphs}
