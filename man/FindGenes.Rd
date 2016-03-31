\name{FindGenes}
\alias{FindGenes}
\title{
Investigates whether genes are differential expressed in multiple clusters
}
\description{
Due to the shifting of compounds over the clusters for the different methods,
it is possible that the same gene is found significant for a different cluster
in another method. These can be tracked with the \code{FindGenes} function. Per method
and per cluster, it will take note of the genes found significant and investigate
if these were also find for another cluster in another method.
}
\usage{
FindGenes(DataLimma, names = NULL)
}
\arguments{
  \item{DataLimma}{
Preferably an output of the \code{DiffGenes} function. If not, an ID element of the
top genes must be present for each cluster of each method specified in the data structure.
}
  \item{names}{
Optional. Names of the methods.
}
}
\value{
The returned value is a list with an element per cluster and per cluster one
for every gene.  Per gene, a vector is given which contain the methods for
which the gene was found. If the cluster is changed compared to the reference
method of DataLimma,  this is indicated with an underscore.
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

MCF7_DiffGenes_FandT10=DiffGenes(list(MCF7_F,MCF7_T),Selection=NULL,GeneExpr=geneMat,
nrclusters=7,method="limma",sign=0.05,top=10,fusionsLog = TRUE, WeightClust = TRUE, 
names = NULL)

MCF7_SharedGenes=FindGenes(DataLimma=MCF7_DiffGenes_FandT10,names=c("FP","TP"))

}
