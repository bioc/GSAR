\name{GSAR-package}
\alias{GSAR-package}
\alias{GSAR}
\docType{package}
\title{Gene Set Analysis in R}

\description{Package GSAR provides a set of statistical methods for 
self-contained gene set analysis. It consists of two-sample multivariate 
nonparametric statistical methods to test a null hypothesis 
against specific alternative hypotheses, such as differences in shift 
(functions \code{\link{KStest}}, \code{\link{MDtest}}, \code{\link{ADtest}}, 
and \code{\link{CVMtest}}), scale (functions \code{\link{RKStest}}, 
\code{\link{RMDtest}}, \code{\link{RADtest}}, \code{\link{RCVMtest}}, and 
\code{\link{AggrFtest}}) or correlation structure 
(function \code{\link{GSNCAtest}}) between two conditions. It also 
offers a graphical visualization tool for correlation 
networks to examine the change in the net correlation structure of a gene 
set between two conditions (function \code{\link{plotMST2.pathway}}). 
The visualization scheme is based on the minimum spanning trees (MSTs). 
Function \code{\link{findMST2}} is used to find the unioin of the first 
and second MSTs. The same tool works as well for protein-protein 
interaction (PPI) networks to highlight the most essential interactions 
among proteins and reveal fine network structure as was already shown in 
Zybailov et. al. 2016. Function \code{\link{findMST2.PPI}} is used to 
find the unioin of the first and second MSTs of PPI networks.
Gene set analysis methods available in this package were 
proposed in Rahmatallah et. al. 2012, Rahmatallah et. al. 2014, and 
Rahmatallah and Glazko 2024. The performance of different methods was 
tested using simulation and real gene expression data. These methods can 
be applied to RNA-Seq count data given that proper normalization is 
used. Proper normalization must take into account both the within-sample 
differences (mainly gene length) and between-samples differences 
(library size). However, because the count data often 
follows the negative binomial distribution, special attention should be 
paid to applying the variance tests (\code{\link{RKStest}}, 
\code{\link{RMDtest}}, \code{\link{RADtest}}, \code{\link{RCVMtest}}, 
and \code{\link{AggrFtest}}). The variance of the 
negative binomial distribution is proportional to it's mean and 
multivariate tests of variance designed specifically for RNA-seq count 
data are under-explored.}

\author{
Yasir Rahmatallah <yrahmatallah@uams.edu>, 
Galina Glazko <gvglazko@uams.edu>

Maintainer: Yasir Rahmatallah <yrahmatallah@uams.edu>, 
Galina Glazko <gvglazko@uams.edu>
}

\references{
Rahmatallah Y. and Glazko G. (2024) Gene Set Analysis: improving data 
interpretability with new differential variance tests. 09 September 2024, 
PREPRINT (Version 1) available at Research Square 
[https://doi.org/10.21203/rs.3.rs-4888767/v1].

Rahmatallah Y., Emmert-Streib F. and Glazko G. (2014) Gene sets net 
correlations analysis (GSNCA): a multivariate differential coexpression test 
for gene sets. Bioinformatics \bold{30}, 360--368.

Rahmatallah Y., Emmert-Streib F. and Glazko G. (2012) Gene set analysis for 
self-contained tests: complex null and specific alternative hypotheses. 
Bioinformatics \bold{28}, 3073--3080.

Zybailov B., Byrd A., Glazko G., Rahmatallah Y. and Raney K. (2016) 
Protein-protein interaction analysis for functional characterization of 
helicases. Methods, \bold{108}, 56--64.
}

\keyword{ package }

\seealso{\code{\link{igraph}}.}
