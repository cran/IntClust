\name{FindCluster}
\alias{FindCluster}
\title{
Find a selection of compounds in the output of \code{ReorderToReference}
}
\description{
\code{FindCluster} selects the compounds belonging to a cluster 
after the results of the methods have been rearranged by the \code{ReorderToReference}.
}
\usage{
FindCluster(List, nrclusters=NULL, select = c(1, 1), fusionsLog = TRUE, 
WeightClust = TRUE, names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
A list of the clustering outputs to be compared. The first element of the
list will be used as the reference in \code{ReorderToReference}.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in.
}
  \item{select}{
The row (the method) and the number of the cluster to select.
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
\value{
A character vector containing the names of the compounds in the
selected cluster.
}

\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)

L=list(MCF7_F,MCF7_T)
names=c("FP","TP")

Comps=FindCluster(L,nrclusters=7,select=c(1,4))
Comps

}