\name{GSAR-package}
\alias{GSAR-package}
\alias{GSAR}
\docType{package}
\title{Gene Set Analysis in R}

\description{Package GSAR provides a set of statistical methods for 
self-contained gene set analysis. It consists of two-sample multivariate 
nonparametric statistical methods to test a null hypothesis 
against specific alternative hypotheses, such as differences in shift 
(functions \code{\link{KStest}} and \code{\link{MDtest}}), scale (functions 
\code{\link{RKStest}}, \code{\link{RMDtest}}, and \code{\link{AggrFtest}}) 
or correlation structure (function \code{\link{GSNCAtest}}) between two 
conditions. It also offers a graphical visualization tool for correlation 
networks to examine the change in the net correlation structure of a gene 
set between two conditions (function \code{\link{plotMST2.pathway}}). 
The visualization scheme is based on the minimum spanning trees (MSTs). 
Function \code{\link{findMST2}} is used to find the unioin of the first 
and second MSTs. The same tool works as well for protein-protein 
interaction (PPI) networks to highlight the most essential interactions 
among proteins and reveal fine network structure as was already shown in 
Zybailov et. al. 2016. Function \code{\link{findMST2.PPI}} is used to 
find the unioin of the first and second MSTs of PPI networks.
Some of the methods available in this package were 
proposed in Rahmatallah et. al. 2014 and Friedman and Rafsky 1979. The 
performance of different methods available in this package was thoroughly 
tested using simulated data and microarray datasets in 
Rahmatallah et. al. 2012 and Rahmatallah et. al. 2014. These methods can 
also be applied to RNA-Seq count data given that proper normalization is 
used. Proper normalization must take into account both the within-sample 
differences (mainly gene length) and between-samples differences 
(library size or sequencing depth). However, because the count data often 
follows the negative binomial distribution, special attention should be 
paid to applying the variance tests (\code{\link{RKStest}}, 
\code{\link{RMDtest}}, and \code{\link{AggrFtest}}). The variance of the 
negative binomial distribution is proportional to it's mean and 
multivariate tests of variance designed specifically for RNA-seq count 
data are virtually unavailable. The performance of variance tests in this 
package with count data highly depends on the used normalization and 
remains currenly under-explored.}

\author{
Yasir Rahmatallah <yrahmatallah@uams.edu>, 
Galina Glazko <gvglazko@uams.edu>

Maintainer: Yasir Rahmatallah <yrahmatallah@uams.edu>, 
Galina Glazko <gvglazko@uams.edu>
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

Zybailov B., Byrd A., Glazko G., Rahmatallah Y. and Raney K. (2016) 
Protein-protein interaction analysis for functional characterization of 
helicases. Methods, \bold{108}, 56--64.
}

\keyword{ package }

\seealso{\code{\link[igraph:igraph-package]{igraph}}.}
