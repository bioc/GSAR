\name{GSAR-package}
\alias{GSAR-package}
\alias{GSAR}
\docType{package}
\title{Gene Set Analysis in R}

\description{Package GSAR provides a set of statistical methods for 
self-contained gene set analysis. It consists of two-sample multivariate 
nonparametric statistical methods to test a null hypothesis 
against specific alternative hypotheses, such as differences in shift 
(function \code{\link{KStest}}), scale (function \code{\link{RKStest}}) or 
correlation structure (function \code{\link{GSNCAtest}}) between two conditions. 
It also offers a graphical visualization tool for correlation networks to 
examine the change in the net correlation structure of a gene set between two 
conditions based on the minimum spanning trees 
(function \code{\link{plotMST2.pathway}}). 
The methods available in this package were proposed in Rahmatallah et. al. 2014 
and Friedman and Rafsky 1979. The performance of different methods available 
in this package was thoroughly tested using simulated data and microarray 
datasets in Rahmatallah et. al. 2012 and Rahmatallah et. al. 2014. These 
methods (except \code{\link{RKStest}}) can also be applied to RNA-Seq count 
data given that proper normalization is used. Proper normalization must take 
into account both the within-sample differences (mainly gene length) and 
between-samples differences (library size or sequencing depth). 
A recent work examining the performance of the tests in package GSAR using 
simulated and real RNA-Seq dataset was submitted for publication and will be 
cited in future versions of this package.}

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
}

\keyword{ package }

\seealso{\code{\link[igraph:igraph-package]{igraph}}.}
