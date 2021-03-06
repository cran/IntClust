% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{FindCluster}
\alias{FindCluster}
\title{Find a selection of objects in the output of \code{ReorderToReference}}
\usage{
FindCluster(List, nrclusters = NULL, select = c(1, 1), fusionsLog = TRUE,
  weightclust = TRUE, names = NULL)
}
\arguments{
\item{List}{A list of the clustering outputs to be compared. The first
element of the list will be used as the reference in
\code{ReorderToReference}.}

\item{nrclusters}{The number of clusters to cut the dendrogram in. Default is NULL.}

\item{select}{The row (the method) and the number of the cluster to select. Default is c(1,1).}

\item{fusionsLog}{Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE}

\item{weightclust}{Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.}

\item{names}{Optional. Names of the methods. Default is NULL.}
}
\value{
A character vector containing the names of the objects in the
selected cluster.
}
\description{
\code{FindCluster} selects the objects belonging to a cluster after the
results of the methods have been rearranged by the
\code{ReorderToReference}.
}
\examples{

data(fingerprintMat)
data(targetMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)

L=list(MCF7_F,MCF7_T)
names=c("FP","TP")

Comps=FindCluster(List=L,nrclusters=7,select=c(1,4))
Comps
}
