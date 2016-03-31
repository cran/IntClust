\name{PreparePathway}
\alias{PreparePathway}
\title{
Preparing a data set for pathway analysis
}
\description{
The functions for pathway analysis in this package can also work on 
results of the integrated data functions. However, a differential
gene expression needs to be conducted to perform pathway analysis.
The function \code{PreparePathway} checks if the necessary elements
are present in the data structures and if not, the elements such as
p-values are created. It is an internal function to all pathway analysis
functions but can be used separately as well.
}
\usage{
PreparePathway(Object, GeneExpr, topG, sign)
}
\arguments{
  \item{Object}{
A list with at least an element with the name "Compounds" such that the function knows
 which compounds to test for differential gene expression. If the elements "Genes"
 and "pvalsgenes" are present as well, these will be collected and the gene expression is
 not analyzed.
}
  \item{GeneExpr}{
The gene expression matrix or ExpressionSet of the objects. The rows should correspond with the genes.
}
  \item{topG}{
Overrules sign. The number of top genes to be returned in the result.
If not specified, only the significant genes are shown.
}
  \item{sign}{
The significance level to be handled.
}
}
\value{
The returned value is a list with three elements:
 \item{pvalsgenses }{This is a list with that contains a vector of raw
 p-values for every group of tested compounds.}
 \item{Compounds}{This is a list with that contains another list per group
 of tested compounds. Every list contains the lead compounds and the ordered
 compounds.}
  \item{Genes }{This is a list with that contains contains another list per group
 of tested compounds. Every list contains two data frames, one with information
 on the top genes and one with information on all genes.}
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(geneMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)


L1=list(MCF7_F)

Comps1=FindCluster(L1,nrclusters=7,select=c(1,1))
Comps2=FindCluster(L1,nrclusters=7,select=c(1,2))
Comps3=FindCluster(L1,nrclusters=7,select=c(1,3))


L2=list()

L2$'Cluster 1'$Compounds$LeadCpds=Comps1
L2$'Cluster 2'$Compounds$LeadCpds=Comps2
L2$'Cluster 3'$Compounds$LeadCpds=Comps2

MCF7_PreparePaths=PreparePathway(Object=L2,GeneExpr=geneMat,topG=NULL,sign=0.05)
str(MCF7_PreparePaths)
}
