\name{GSARnews}
\title{GSAR News}
\encoding{UTF-8}

\section{Version 1.4.0}{\itemize{

\item
New argument return.weights added to function plotMST2.pathway. If return.weights=TRUE, the weight vectors found by GSNCA for the genes under two classes are returned as a 2-column matrix.

\item
A few changes in the notations and the codes of the case studies used in the vignette.
}}

\section{Version 1.2.1}{\itemize{

\item
Bug fix in function HDP.ranking to accomodate the changes in package igraph (version 1.0.1).

\item
Additional steps added to the code of the second case study to filter C2 gene sets properly.
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
