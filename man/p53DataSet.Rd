\name{p53DataSet}

\alias{p53DataSet}

\docType{data}

\title{p53 Dataset of the NCI-60 Cell Lines}

\description{A matrix of gene expression profiles for a processed version of 
the p53 dataset obtained from the NCI-60 cell lines using the hgu95av2 
microarray platform.}

\usage{data(p53DataSet)}

\format{
A matrix of 8655 rows and 50 columns where rows correspond to genes and 
columns correspond to samples. Gene symbol identifiers are used for rows. 
Column names indicate the class of the samples (wild type p53 or mutated p53) 
with the first 17 column names starting with \code{WT1} and ending with 
\code{WT17} and next 33 column names starting with \code{MUT1} and ending with 
\code{MUT33}.
}

\details{
p53 is a major tumor suppressor protein. The p53 dataset comprises 50 samples 
of NCI-60 cell lines differentiated based on the status of the TP53 gene: 17 
cell lines carrying wild type (WT) TP53 and 33 cell lines carrying mutated 
(MUT) TP53 (Olivier et. al. 2002, Subramanian et. al. 2005). Transcriptional 
profiles obtained from microarrays of platform hgu95av2 were obtained from 
the available datasets at the GSEA Broad Institute's website.

Probe level intensities were quantile normalized and transformed to the log 
scale using log2(1 + intensity). Probes originally had Affymetrix identifiers 
which are mapped to unique gene symbol identifiers. Probes without mapping 
to entrez and gene symbol identifiers were discarded. Probes with duplicate 
intensities were assessed and the probe with the largest absolute value of 
t-statistic between WT and MUT conditions was selected as the gene match. 
Genes were assigned gene symbol identifiers and columns were assigned names 
indicating weither they belong to WT or MUT condition. The columns were 
sorted such that the first 17 columns are WT samples and the next 33 
columns are the MUT samples. \code{p53DataSet} was used in the analysis 
presented in Rahmatallah et. al. 2014.
}

\source{
Broad Institute (http://www.broadinstitute.org/gsea/datasets.jsp)
}

\references{
Rahmatallah Y., Emmert-Streib F. and Glazko G. (2014) Gene sets net 
correlations analysis (GSNCA): a multivariate differential coexpression 
test for gene sets. Bioinformatics \bold{30}, 360--368.

Subramanian A., Tamayo P., Mootha V., Mukherjee S., Ebert B., Gillette M., 
Paulovich A., Pomeroy S., Golub T., Lander E. and Mesirov J. (2005) Gene set 
enrichment analysis: A knowledge-based approach for interpreting genome-wide 
expression profiles. Proc. Natl. Acad. Sci. \bold{102}, 15545--15550.

Olivier M., Eeles R., Hollstein M., Khan M., Harris C. and Hainaut P. (2002) 
The IARC TP53 database: new online mutation analysis and recommendations to 
users. Hum. Mutat. \bold{19}, 607--614.
}

\examples{
data(p53DataSet)
dim(p53DataSet)
}

\keyword{datasets}
