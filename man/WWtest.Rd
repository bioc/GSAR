\name{WWtest}

\alias{WWtest}

\title{Multivariate Generalization of the Wald-Wolfowitz Runs Test}

\description{Performs two-sample nonparametric multivariate generalization of 
the Wald-Wolfowitz runs test based on the minimum spanning tree (MST). It 
tests the alternative hypothesis that a set of features has different 
distributions in two conditions against the null hypothesis of having the 
same distribution.}

\usage{WWtest(object, group, nperm=1000, pvalue.only=TRUE)}

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

\details{This function tests the alternative hypothesis that a set of features 
has different distributions in two conditions against the null hypothesis of 
having the same distribution. It performs the two-sample nonparametric 
multivariate generalization of the Wald-Wolfowitz runs test based on the 
minimum spanning tree (MST) as proposed by Friedman and Rafsky (1979). The 
performance of this test under different alternative hypotheses was thoroughly 
examind in Rahmatallah et. al. (2012). The null distribution of the test 
statistic is estimated by permuting sample labels \code{nperm} times and 
calculating the test statistic for each. P-value is calculated as

\deqn{p.value = \frac{\sum_{k=1}^{nperm} I \left[ W_{k} \leq W_{obs} \right] + 1}{nperm + 1}}

where \eqn{W_{k}} is the test statistic for permutation \code{k}, \eqn{W_{obs}} is the 
observed test statistic, and \code{I} is the indicator function.}

\value{
When \code{pvalue.only=TRUE} (default), function \code{WWtest} returns the 
p-value indicating the attained significance level. When \code{pvalue.only=FALSE}, 
function \code{WWtest} produces a list of length 3 with the following components:
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
\author{
Yasir Rahmatallah and Galina Glazko
}

\seealso{\code{\link{KStest}}, \code{\link{RKStest}}, \code{\link{MDtest}}, 
\code{\link{RMDtest}}.}

\examples{
## generate a feature set of length 20 in two conditions
## each condition has 20 samples
## use multivariate normal distribution
library(MASS)
ngenes <- 20
nsamples <- 40
## let the mean vector have zeros of length 20 for condition 1
zero_vector <- array(0,c(1,ngenes))
## let the mean vector have 2s of length 20 for condition 2
mu_vector <- array(2,c(1,ngenes))
## set the covariance matrix to be an identity matrix
cov_mtrx <- diag(ngenes)
gp1 <- mvrnorm((nsamples/2), zero_vector, cov_mtrx)
gp2 <- mvrnorm((nsamples/2), mu_vector, cov_mtrx)
## combine the data of two conditions into one dataset
gp <- rbind(gp1,gp2)
dataset <- aperm(gp, c(2,1))
## first 20 samples belong to group 1
## second 20 samples belong to group 2
pvalue <- WWtest(object=dataset, group=c(rep(1,20),rep(2,20))) 
}

\keyword{multivariate}
\keyword{nonparametric}
