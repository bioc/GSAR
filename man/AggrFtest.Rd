\name{AggrFtest}

\alias{AggrFtest}

\title{Aggregated F-Test of Variance Using Fisher's Probability Combining Method}

\description{Performs two-sample nonparametric test of variance. The univariate 
F-test is used for every gene in the gene set and the resulted p-values are 
aggregated together using Fisher's probability combining method and used as 
the test statistic. The null distribution of the test statistic is estimated by 
permuting sample labels and calculating the test statistic for a large number 
of times. This statistic tests the null hypothesis that none of the genes shows 
significant difference in variance between two conditions against the alternative 
hypothesis that at least one gene shows significant difference in variance 
between two conditions according to the F-test.}

\usage{AggrFtest(object, group, nperm=1000, pvalue.only=TRUE)}

\arguments{
  \item{object}{a numeric matrix with columns and rows respectively 
corresponding to samples and features (genes).}
  \item{group}{a numeric vector indicating group associations for samples. 
Possible values are 1 and 2.}
  \item{nperm}{a numeric value indicating the number of permutations 
used to estimate the null distribution of the test statistic. If not 
given, a default value 1000 is used.}
  \item{pvalue.only}{logical. If \code{TRUE} (default), the p-value is 
returned. If \code{FALSE} a list of length three containing the observed 
statistic, the vector of permuted statistics, and the p-value is returned.}
}

\details{This function tests the null hypothesis that none of the genes 
in a gene set shows a significant difference in variance between two 
conditions according to the F-test against the alternative hypothesis 
that at least one gene shows significant difference in variance according 
to the F-test. It performs a two-sample nonparametric test of variance by 
using the univariate F-test for every gene in a set, adjust for multiple 
testing using the Benjamini and Hochberg method (also known as FDR) as 
shown in Benjamini and Hochberg (1995), and then aggregates the obtained 
adjusted p-values using Fisher's probability combining method to get a 
test statistic (T) for the gene set

\deqn{T = -2 \sum_{i=1}^{p} \log_{e} (p_{i})}

where \eqn{p_{i}} is the adjusted p-value of the univariate F-test for 
gene \code{i}. The null distribution of the test statistic is estimated 
by permuting sample labels \code{nperm} times and calculating the test 
statistic T for each. P-value is calculated as
 
\deqn{p.value = \frac{\sum_{k=1}^{nperm} I \left[ T_{k} \geq T_{obs} \right] + 1}{nperm + 1}}

where \eqn{T_{k}} is the test statistic for permutation \code{k}, \eqn{T_{obs}} is the 
observed test statistic, and \code{I} is the indicator function.
}

\value{
When \code{pvalue.only=TRUE} (default), function \code{AggrFtest} returns 
the p-value indicating the attained significance level. When 
\code{pvalue.only=FALSE}, function \code{AggrFtest} produces a list of 
length 3 with the following components:
\item{statistic}{the value of the observed test statistic.}
\item{perm.stat}{numeric vector of the resulting test statistic for 
\code{nperm} random permutations of sample labels.}
\item{p.value}{p-value indicating the attained significance level.}
}

\references{
Benjamini Y. and Hochberg Y. (1995) Controlling the false discovery rate: a 
practical and powerful approach to multiple testing. Journal of the Royal 
Statistical Society Series B \bold{57}, 289--300.
}

\author{Yasir Rahmatallah and Galina Glazko}

\seealso{\code{\link{RKStest}}, \code{\link{RMDtest}}, 
\code{\link{RADtest}}, \code{\link{RCVMtest}}.}

\examples{
## generate a feature set of length 20 in two conditions
## each condition has 20 samples
## use multivariate normal distribution
library(MASS)
ngenes <- 20
nsamples <- 40
## let the mean vector have zeros of length 20 for both conditions
zero_vector <- array(0,c(1,ngenes))
## set the covariance matrix to be an identity matrix for condition 1
cov_mtrx <- diag(ngenes)
gp1 <- mvrnorm((nsamples/2), zero_vector, cov_mtrx)
## set some scale difference in the covariance matrix for condition 2
cov_mtrx <- cov_mtrx*3
gp2 <- mvrnorm((nsamples/2), zero_vector, cov_mtrx)
## combine the data of two conditions into one dataset
gp <- rbind(gp1,gp2)
dataset <- aperm(gp, c(2,1))
## first 20 samples belong to group 1
## second 20 samples belong to group 2
pvalue <- AggrFtest(object=dataset, group=c(rep(1,20),rep(2,20)))
}

\keyword{nonparametric}
