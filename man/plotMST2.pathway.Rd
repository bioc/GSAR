\name{plotMST2.pathway}

\alias{plotMST2.pathway}

\title{Plot MST2 for a pathway in two conditions}

\description{This is a wrapper function which uses function 
\code{\link{findMST2}} to find the union of the first and second minimum 
spanning trees (or MST2) of the correlation network for a feature set (pathway) 
under two conditions. It plots the MST2 of the correlation network of the 
feature set under both conditions side-by-side and highlights hub nodes to 
facilitate a visual comparison.}

\usage{plotMST2.pathway(object, group, name=NULL, legend.size=1, label.size=1, 
cor.method="pearson", min.sd=1e-3)}

\arguments{
  \item{object}{a numeric matrix with columns and rows respectively 
corresponding to samples and features.}
  \item{group}{a numeric vector indicating group associations for samples. 
Possible values are 1 and 2.}
  \item{name}{an optional character string giving the name of the feature set 
(pathway). If given, the name will be displayed at the top of the plot.}
  \item{legend.size}{an optional numeric value controlling the relative font 
size of the legend to the default font size. Default is 1.}
  \item{label.size}{an optional numeric value controlling the relative font 
size of vertex labels to the default font size. Default is 1.}
  \item{cor.method}{a character string indicating which correlation 
coefficient is to be computed. Possible values are \dQuote{\code{pearson}} 
(default), \dQuote{\code{spearman}} and \dQuote{\code{kendall}}.}
  \item{min.sd}{the minimum allowed standard deviation for any feature. If 
any feature has a standard deviation smaller than \code{min.sd} the execution 
stops and an error message is returned.}
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

\keyword{misc}
