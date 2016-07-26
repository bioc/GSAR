\name{KStest}

\alias{KStest}

\title{Multivariate Kolmogorov-Smirnov Test of Means}

\description{Performs two-sample nonparametric multivariate test of means 
based on the minimum spanning tree (MST) and Kolmogorov-Smirnov statistic. 
It tests the null hypothesis that a set of features has the
same mean in two conditions versus different means.}

\usage{KStest(object, group, nperm=1000, pvalue.only=TRUE)}

\arguments{
  \item{object}{a numeric matrix with columns and rows respectively 
corresponding to samples and features.}
  \item{group}{a numeric vector indicating group associations for samples. 
Possible values are 1 and 2.}
  \item{nperm}{number of permutations used to estimate the null distribution 
of the test statistic. If not given, a default value 1000 is used.}
  \item{pvalue.only}{logical. If \code{TRUE} (default), the p-value is 
returned. If \code{FALSE} a list of length three containing the observed 
statistic, the vector of permuted statistics, and the p-value is returned.}
}

\details{This function tests the null hypothesis that a set of features has 
no shift between two conditions. It performs a two-sample nonparametric 
multivariate test based on the minimum spanning tree (MST) and 
Kolmogorov-Smirnov statistic as proposed by Friedman and Rafsky (1979). 
The MST of the weighted undirectional graph created from the samples is found. 
The nodes of the MST are ranked based on their position in the MST. The MST 
is rooted at the node with largest geodisic distance (rank 1) and then nodes 
are ranked in the High Directed Preorder (HDP) traversal of the tree 
(Rahmatallah et. al. 2012). The quantity \eqn{d_i = (r_i / n_1) - (s_i / n_2)} 
is calculated where \eqn{r_i}(\eqn{s_i}) is the number of nodes (samples) 
from condition 1(2) which ranked lower than \eqn{i}, \eqn{1 \le i \le N} and 
\eqn{N} is the total number of samples. The Kolmogorov-Smirnov statistic is 
given by the maximum absolute difference 
\eqn{D = \sqrt{\frac{n_{1}n_{2}}{n_{1}+n_{2}}} max|d_i|}. The performance of this 
test under different alternative hypotheses was thoroughly examind 
in Rahmatallah et. al. (2012). The null distribution of the test statistic 
is estimated by permuting sample labels \code{nperm} times and calculating 
the test statistic for each. P-value is calculated as
 
\deqn{p.value = \frac{\sum_{k=1}^{nperm} I \left[ D_{k} \geq D_{obs} \right] + 1}{nperm + 1}}

where \eqn{D_{k}} is the test statistic for permutation \code{k}, \eqn{D_{obs}} is the 
observed test statistic, and \code{I} is the indicator function.
}

\value{
When \code{pvalue.only=TRUE} (default), function \code{KStest} returns 
the p-value indicating the attained significance level. When 
\code{pvalue.only=FALSE}, function \code{KStest} produces a list of 
length 3 with the following components:
\item{statistic}{the value of the observed test statistic.}
\item{perm.stat}{numeric vector of the resulting test statistic for 
\code{nperm} random permutations of sample labels.}
\item{p.value}{p-value indicating the attained significance level.}
}

\references{
Rahmatallah Y., Emmert-Streib F. and Glazko G. (2012) Gene set analysis for 
self-contained tests: complex null and specific alternative hypotheses. 
Bioinformatics \bold{28}, 3073--3080.

Friedman J. and Rafsky L. (1979) Multivariate generalization of the 
Wald-Wolfowitz and Smirnov two-sample tests. Ann. Stat. \bold{7}, 697--717.
}

\author{Yasir Rahmatallah and Galina Glazko}

\note{This function invokes function \code{\link{HDP.ranking}} which does not 
work properly if there is any node in the MST with more than 26 links. 
However, this situation is almost impossible for a dataset composed of a few 
hundreds or less of samples.}

\seealso{\code{\link{MDtest}}, \code{\link{WWtest}}, \code{\link{RKStest}}, 
\code{\link{RMDtest}}, \code{\link{HDP.ranking}}.}

\examples{
## generate a feature set of length 20 in two conditions
## each condition has 20 samples
## use multivariate normal distribution
library(MASS)
ngenes <- 20
nsamples <- 40
## let the mean vector have zeros of length 20 in both conditions
zero_vector <- array(0,c(1,ngenes))
## set the covariance matrix to be an identity matrix for both conditions
cov_mtrx <- diag(ngenes)
gp <- mvrnorm(nsamples, zero_vector, cov_mtrx)
## apply a mean shift of 3 to half of the features under condition 1
gp[1:20,1:10] <- gp[1:20,1:10] + 3
dataset <- aperm(gp, c(2,1))
## first 20 samples belong to condition 1
## second 20 samples belong to condition 2
pvalue <- KStest(object=dataset, group=c(rep(1,20),rep(2,20)))
}

\keyword{multivariate}
\keyword{nonparametric}
