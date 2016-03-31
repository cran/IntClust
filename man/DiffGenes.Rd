\name{DiffGenes}
\alias{DiffGenes}
\title{
Differential gene expressions for multiple results
}
\description{
The function \code{DiffGenes} will, given the output of a certain method, 
look for genes that are differentially expressed for each cluster
by applying the limma function to that cluster and compare it to
all other clusters simultaneously. If a list of outputs of several
methods is provided, DiffGenes will perform the limma function
for each method. 
}
\usage{
DiffGenes(List,Selection=NULL, GeneExpr = NULL, nrclusters = NULL, method = "limma", 
sign = 0.05, topG = NULL, fusionsLog = TRUE, WeightClust = TRUE, 
names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
A list of the clustering outputs to be compared. The first element
of the list will be used as the reference in \code{ReorderToReference}.
}
 \item{Selection}{
If differential gene expression should be investigated for a specific selection of compounds, this selection can be provided here.
Selection can be of the type "character" (names of the compounds) or "numeric" (the number of specific cluster).
 }
  \item{GeneExpr}{
The gene expression matrix or ExpressionSet of the objects. The rows should correspond with the genes.
}
  \item{nrclusters}{
Optional. The number of clusters to cut the dendrogram in. The number of clusters should not be
specified if the interest lies only in a specific selection of compounds which is known by name. 
Otherwise, it is required.
}
  \item{method}{
The method to applied to look for DE genes. For now, only the limma
method
is available
}
  \item{sign}{
The significance level to be handled.
}
  \item{topG}{
Overrules sign. The number of top genes to be shown.
}
  \item{fusionsLog}{
To be handed to \code{ReorderToReference}.
}
  \item{WeightClust}{
To be handed to \code{ReorderToReference}.
}
  \item{names}{
Optional. Names of the methods.
}
}
\details{
The function rearranges the clusters of the methods to a 
reference method such that a comparison is made easier. 
Given a list of methods, it calls upon \code{ReorderToReference} to
rearrange the number of clusters according to the first
element of the list which will be used as the reference.
}
\value{
The returned value is a list with an element per method. Each element
contains a list per cluster with the following elements:
\item{Compounds}{A list with the elements LeadCpds (the compounds of interest) and 
OrderedCpds (all compounds in the order of the clustering result)}
\item{Genes}{A list with the elements TopDE (a table with information on the top
genes) and AllDE (a table with information on all genes)}
}
\references{
SMYTH, G. K. (2004). Linear models and empirical Bayes methods for assessing
differential expression in microarray experiments. Statistical Applications
in Genetics and Molecular Biology. 3(1).
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
data(geneMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)

L=list(MCF7_T ,MCF7_F)

MCF7_FT_DE = DiffGenes(L,GeneExpr=geneMat,nrclusters=7,method="limma",
sign=0.05,topG=10,fusionsLog=TRUE,WeightClust=TRUE)

}
\keyword{Differential Gene Expression }
\keyword{ limma }
