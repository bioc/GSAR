\name{GSARnews}
\title{GSAR News}
\encoding{UTF-8}

\section{Version 1.40.0}{\itemize{

\item
New function ADtest is introduced. It implements a nonparametric multivariate test of means based on the HDP ranking of samples in the MST and the Anderson-Darling statistic.

\item
New function RADtest is introduced. It implements a nonparametric multivariate test of variance based on the radial ranking of samples in the MST and the Anderson-Darling statistic.

\item
New function CVMtest is introduced. It implements a nonparametric multivariate test of means based on the HDP ranking of samples in the MST and the Cramer-Von Mises statistic.

\item
New function RCVMtest is introduced. It implements a nonparametric multivariate test of variance based on the radial ranking of samples in the MST and the Cramer-Von Mises statistic.
}}

\section{Version 1.8.0}{\itemize{

\item
New function MDtest is introduced. It implements a nonparametric multivariate test of means based on sample ranking in the MST similar to function KStest, but the test statistic is the mean deviation between the CDFs of two conditions.

\item
New function RMDtest is introduced. It implements a nonparametric multivariate test of variance based on sample ranking in the MST similar to function RKStest, but the test statistic is the mean deviation between the CDFs of two conditions.

\item
New function AggrFtest is introduced. It implements a nonparametric test of variance by aggregating the univariate p-values obtained by the F-test using Fisher's probability combining method. It test the hypothesis that all genes in a gene set show no significant difference in variance between two conditions against the alternative hypothesis that at least one gene in the gene set shows significant difference in variance between two conditions.

\item
New function findMST2.PPI is introduced. It finds the union of the first and second MSTs similar to function findMST2, but it accepts an object of class igraph as input rather that a matrix of gene expression data. The input igraph object represents a protein-protein interaction (PPI) network that can be binary or weighted, directed or undirected.

\item
New wrapper function TestGeneSets is introduced. It performs a specific statistical method from the ones available in package GSAR for multiple gene sets. The gene sets are provided as a list of character vectors where each entry has the feature (gene) identifiers in a single gene set.

\item
New argument pvalue.only added to all available statistical methods in the package. When pvalue.only=TRUE (default), each statistical method returns the p-value only. When pvalue.only=FALSE, each statistical method returns a list of length 3 consisting of the observed statistic, vector of permuted statistics, and p-value.

\item
New arguments leg.x, leg.y, group1.name, group2.name, label.color, label.dist, vertex.size, vertex.label.font, and edge.width added to function plotMST2.pathway to allow more flexibility in generating plots. The values of most of these arguments are passed to function plot.igraph.
}}

\section{Version 1.4.0}{\itemize{

\item
New argument return.weights added to function plotMST2.pathway. If return.weights=TRUE, the weight vectors found by GSNCA for the genes under two classes are returned as a 2-column matrix.
}}

\section{Version 1.2.2}{\itemize{

\item
Bug fix in function HDP.ranking to accomodate the changes in package igraph (version 1.0.1).

\item
Additional steps added to the code of the second case study in the vignette to filter C2 gene sets properly.

\item
Minor changes in the vignette.
}}

\section{Version 1.0.0}{\itemize{

\item
The package provides two-sample nonparametric multivariate statistical methods to test specific alternative hypotheses 
against a null hypothesis.

\item
GSAR depends on package igraph to handle graphs in objects of class \code{igraph} and uses some functions too.

\item
New capabilities and future changes will be reported in subsequent versions.
}}
