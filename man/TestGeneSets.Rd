\name{TestGeneSets}

\alias{TestGeneSets}

\title{Test a List of Gene Sets Using a Specific Statistical Method}

\description{A wrapper function that invokes a specific statistical method 
from the ones available in package GSAR (see Rahmatallah et. al. 2014, 
Rahmatallah et. al. 2012 and Friedman and Rafsky 1979 for details) to 
test a list of gene sets in a sequential order and returns results in 
a list object.}

\usage{TestGeneSets(object, group, geneSets=NULL, min.size=10, max.size=500, 
test=NULL, nperm=1000, mst.order=1, pvalue.only=TRUE)}

\arguments{
  \item{object}{a numeric matrix with columns and rows respectively 
corresponding to samples and features.}
  \item{group}{a numeric vector indicating group associations for samples. 
Possible values are 1 and 2.}
  \item{geneSets}{a list of character vectors providing the identifiers of 
features to be considered in each gene set.}
  \item{min.size}{a numeric value indicating the minimum allowed gene set 
size. Default value is 10.}
  \item{max.size}{a numeric value indicating the maximum allowed gene set 
size. Default value is 500.}
  \item{test}{a character parameter indicating which statistical method 
to use for testing the gene sets. Must be one of \dQuote{\code{GSNCAtest}}, 
\dQuote{\code{WWtest}}, \dQuote{\code{KStest}}, \dQuote{\code{MDtest}}, 
\dQuote{\code{RKStest}}, \dQuote{\code{RMDtest}}.}
  \item{nperm}{number of permutations used to estimate the null distribution 
of the test statistic. If not given, a default value 1000 is used.}
  \item{mst.order}{numeric value to indicate the consideration of the union 
of the first \code{mst.order} MSTs when \dQuote{\code{RKStest}} or 
\dQuote{\code{RMDtest}} are used. Default value is 1. Maximum allowed 
value is 5.}
  \item{pvalue.only}{logical. If \code{TRUE} (default), the p-value is 
returned. If \code{FALSE} a list of length three containing the observed 
statistic, the vector of permuted statistics, and the p-value is returned.}
}

\details{This is a wrapper function that facilitates the use of any 
statistical method in package GSAR for multiple gene sets that are 
provided in a list object. The function filters out any gene that is 
abscent in the considered data (input parameter \code{object}) from 
the gene sets and discard any set that is too small in size (has less 
than \code{min.size} genes) or too large (has more than 
\code{max.size} genes). The function performs the specified method 
for all the remaining gene sets in a sequential order and return 
results in a list object.}

\value{
A list object of length equals the length of the provided gene set list. 
When \code{pvalue.only=TRUE} (default), each item in the returned list 
by function \code{TestGeneSets} consists of a numeric p-value indicating 
the attained significance level obtained by the specified method. When 
\code{pvalue.only=FALSE}, each item in the returned list is a list of 
length 3 with the following components:
\item{statistic}{the value of the observed test statistic.}
\item{perm.stat}{numeric vector of the resulting test statistic for 
\code{nperm} random permutations of sample labels.}
\item{p.value}{p-value indicating the attained significance level.}
}

\references{
Rahmatallah Y., Emmert-Streib F. and Glazko G. (2014) Gene sets net 
correlations analysis (GSNCA): a multivariate differential coexpression test 
for gene sets. Bioinformatics \bold{30}, 360--368.

Rahmatallah Y., Emmert-Streib F. and Glazko G. (2012) Gene set analysis for 
self-contained tests: complex null and specific alternative hypotheses. 
Bioinformatics \bold{28}, 3073--3080.

Friedman J. and Rafsky L. (1979) Multivariate generalization of the 
Wald-Wolfowitz and Smirnov two-sample tests. Ann. Stat. \bold{7}, 697--717.
}

\author{Yasir Rahmatallah and Galina Glazko}

\seealso{\code{\link{GSNCAtest}}, \code{\link{WWtest}}, \code{\link{MDtest}}, 
\code{\link{KStest}}, \code{\link{RKStest}}, \code{\link{RMDtest}}.}

\examples{
## generate a feature set of size 50 in two conditions
## where each condition has 20 samples
## use multivariate normal distribution
library(MASS)
ngenes <- 50
nsamples <- 40
## let the mean vector have zeros of length 50 for both conditions
zero_vector <- array(0,c(1,ngenes))

## set the covariance matrix to be an identity matrix for both conditions
cov_mtrx <- diag(ngenes)
gp <- mvrnorm(nsamples, zero_vector, cov_mtrx)
## apply a mean shift of 5 to the first 10 features under condition 1
gp[1:20,1:10] <- gp[1:20,1:10] + 5
dataset <- aperm(gp, c(2,1))
## assign a unique identifier to each gene
rownames(dataset) <- as.character(c(1:ngenes))
## first 20 samples belong to condition 1
## second 20 samples belong to condition 2
sample.labels <- c(rep(1,20),rep(2,20))
## construct 3 named gene sets such that they respectively consist of 
## genes 1 to 20, 11 to 40, and 31 to 50. Notice that gene sets 
## can have intersections and can be of different sizes 
## Sine only the first 10 genes have a significant difference between 
## the two conditions the only the first gene set (set1) returns a 
## small p-value when KStest is selected
geneSets <- list("set1"=as.character(c(1:20)), "set2"=as.character(c(11:40)), 
"set3"=as.character(c(31:40)))
results <- TestGeneSets(object=dataset, group=sample.labels, 
geneSets=geneSets, test="KStest")
}

\keyword{multivariate}
\keyword{nonparametric}
