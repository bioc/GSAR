\name{GSNCAtest}

\alias{GSNCAtest}

\title{Gene Sets Net Correlations Analysis}

\description{Performs Gene Sets Net Correlation Analysis (GSNCA) test to 
detect differentially coexpressed gene sets.}

\usage{GSNCAtest(object, group, nperm=1000, cor.method="pearson", check.sd=TRUE, 
min.sd=1e-3, max.skip=10)}

\arguments{
  \item{object}{a numeric matrix with columns and rows respectively 
corresponding to samples and features.}
  \item{group}{a numeric vector indicating group associations for samples. 
Possible values are 1 and 2.}
  \item{nperm}{number of permutations used to estimate the null distribution 
of the test statistic.
If not given, a default value 1000 is used.}
  \item{cor.method}{a character string indicating which correlation 
coefficient is to be computed. Possible values are \dQuote{\code{pearson}} 
(default), \dQuote{\code{spearman}} and \dQuote{\code{kendall}}.}
  \item{check.sd}{logical. Should the standard deviations of features checked 
for small values before the intergene correlations are computed? Default is 
TRUE (recommended).}
  \item{min.sd}{the minimum allowed standard deviation for any feature. If 
any feature has a standard deviation smaller than \code{min.sd} the execution 
stops and an error message is returned.}
  \item{max.skip}{maximum number of skipped random permutations which yield 
any feature with a standard deviation less than \code{min.sd}.}
}

\details{This function performs the Gene Sets Net Correlations Analysis 
(GSNCA), a two-sample nonparametric multivariate differential coexpression test 
that accounts for the correlation structure between features (genes).
The test assigns weight factors to features under one condition and adjust 
these weights simultaneously such that equality is achieved between each 
feature's weight and the sum of its weighted absolute correlations with other 
features in the feature set. The problem is solved as an eigenvector problem 
with a unique solution (see Rahmatallah et. al. 2014 for details). The test 
statistic \eqn{w_{GSNCA}} is given by the first norm between the scaled 
weight vectors \eqn{w^{(1)}} and \eqn{w^{(2)}} (each vector is multiplied 
by its norm) between two conditions

\deqn{w_{GSNCA} = \sum_{i=1}^{p} | w_{i}^{(1)} - w_{i}^{(2)} |}

This test statistic tests the null hypothesis that \eqn{w_{GSNCA}=0} against 
the alternative that \eqn{w_{GSNCA}} does not equal to zero. 
The performance of this test was thoroughly examind in 
Rahmatallah et. al. (2014). The null distribution of the test statistic is 
estimated by permuting sample labels \code{nperm} times and calculating the 
test statistic for each. P-value is calculated as 
\code{p.value}=(\code{b} + 1)/(\code{nperm} + 1) where \code{b} is the number 
of permutations giving a more extreme statistic than the observed test 
statistic.

In the case of RNA-seq count data, some non-expressed genes may have zero 
counts across the samples under one or two conditions. Such situation 
results in zero or tiny standard deviation for one or more features. 
Such case produces an error in command \code{cor} used to compute the 
correlation coefficients between features. 
To avoid this situation, standard deviations are checked in advance when 
\code{check.sd} is \code{TRUE} (default) and if any is found below the minimum 
limit \code{min.sd} (default is \code{1e-3}), the execution stops and an error 
message is returned indicating the number of feature causing the problem 
(if only one the index of that feature is given too). 
If a feature has nearly a constant level for some samples under both 
conditions, permuting sample labels may group such samples under one condition 
and produce a standard deviation smaller than \code{min.sd}. To allow the test 
to skip such permutations without causing excessive delay, we set an upper 
limit for the number of allowed skips by the argument \code{max.skip} 
(default is 10). If the upper limit is exceeded, an error message is 
returned. Allowing this skipping may or may not solve the issue depending on 
the proportion of samples causing the problem in the feature set.

If the user is certain that the tested feature sets contain no feature with 
nearly equal levels over many samples (such as the case with microarrays), 
the checking stage for tiny standard deviations can be skipped by setting 
\code{check.sd} to \code{FALSE} in order to reduce the execution time.}

\value{
GSNCAtest produces a list with the following components:
\item{statistic}{the value of the observed test statistic.}
\item{perm.stat}{numeric vector of the resulting test statistic for 
\code{nperm} random permutations of sample labels.}
\item{p.value}{p-value indicating the attained significance level.}
}

\references{
Rahmatallah Y., Emmert-Streib F. and Glazko G. (2014) Gene sets net 
correlations analysis (GSNCA): a multivariate differential coexpression test 
for gene sets. Bioinformatics \bold{30}, 360--368.
}

\author{Yasir Rahmatallah and Galina Glazko}

\seealso{\code{\link{findMST2}}, \code{\link{plotMST2.pathway}}.}

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
result <- GSNCAtest(object=dataset, group=c(rep(1,20),rep(2,20))) 
pvalue <- result$p.value
}

\keyword{multivariate}
\keyword{nonparametric}
