\name{ADtest}

\alias{ADtest}

\title{Multivariate Anderson-Darling Test of Means}

\description{Performs two-sample nonparametric multivariate test of means 
based on the minimum spanning tree (MST) and Anderson-Darling statistic. 
It tests the null hypothesis that a set of features has the
same mean in two conditions versus different means.}

\usage{ADtest(object, group, nperm=1000, pvalue.only=TRUE)}

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
Anderson-Darling statistic as proposed in Rahmatallah and Glazko (2024). 
The MST of the weighted undirectional graph created from the samples is found. 
The nodes of the MST are ranked based on their position in the MST. The MST 
is rooted at the node with largest geodisic distance and then nodes 
are ranked in the High Directed Preorder (HDP) traversal of the tree 
(Rahmatallah et. al. 2012). The Anderson-Darling statistic can be defined as

\deqn{A = \frac{1}{n_{1} n_{2}} \sum_{i=1}^{N-1} \frac{(r_{i} N - i n_{1})^{2}}{i (N-i)}} 

where \eqn{r_i} is the number of nodes (samples) from condition 1 which ranked lower 
than \eqn{i}, \eqn{1 \le i \le N}, \eqn{N} is the total number of samples, and 
\eqn{n_{1}} and \eqn{n_{2}} are respectively the number of samples in groups 1 and 2. 
The Anderson-Darling statistic put more emphasis on the tails of the deviation 
between the empirical cumulative distribution functions (CDFs) of samples of two 
phenotypes in the MST. The performance of this test under different alternative 
hypotheses was examind in Rahmatallah and Glazko (2024). The null distribution of 
the test statistic is estimated by permuting sample labels \code{nperm} times and 
calculating the test statistic for each. P-value is calculated as
 
\deqn{p.value = \frac{\sum_{k=1}^{nperm} I \left[ A_{k} \geq A_{obs} \right] + 1}{nperm + 1}}

where \eqn{A_{k}} is the test statistic for permutation \code{k}, \eqn{A_{obs}} is the 
observed test statistic, and \code{I} is the indicator function.
}

\value{
When \code{pvalue.only=TRUE} (default), function \code{ADtest} returns 
the p-value indicating the attained significance level. When 
\code{pvalue.only=FALSE}, function \code{ADtest} produces a list of 
length 3 with the following components:
\item{statistic}{the value of the observed test statistic.}
\item{perm.stat}{numeric vector of the resulting test statistic for 
\code{nperm} random permutations of sample labels.}
\item{p.value}{p-value indicating the attained significance level.}
}

\references{
Rahmatallah Y. and Glazko G. (2024) Gene Set Analysis: improving data 
interpretability with new differential variance tests. 09 September 2024, 
PREPRINT (Version 1) available at Research Square 
[https://doi.org/10.21203/rs.3.rs-4888767/v1].

Rahmatallah Y., Emmert-Streib F. and Glazko G. (2012) Gene set analysis for 
self-contained tests: complex null and specific alternative hypotheses. 
Bioinformatics \bold{28}, 3073--3080.

}

\author{Yasir Rahmatallah and Galina Glazko}

\note{This function invokes function \code{\link{HDP.ranking}} which does not 
work properly if there is any node in the MST with more than 26 links. 
However, this situation is almost impossible for a dataset composed of a few 
hundreds or less of samples.}

\seealso{\code{\link{KStest}}, \code{\link{MDtest}}, \code{\link{CVMtest}}, 
\code{\link{RADtest}}, \code{\link{RKStest}}, \code{\link{RMDtest}}, 
\code{\link{RCVMtest}}, \code{\link{HDP.ranking}}.}

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
pvalue <- ADtest(object=dataset, group=c(rep(1,20),rep(2,20)))
}

\keyword{multivariate}
\keyword{nonparametric}
