\name{SharedComps}
\alias{SharedComps}
\title{
Intersection of clusters over multiple methods
}
\description{
The \code{SharedComps} function is an easy way to select the compounds that are shared over clusters of different methods.
}
\usage{
SharedComps(List,nrclusters=NULL,fusionsLog=FALSE,WeightClust=FALSE,names=NULL)
}
\arguments{
  \item{List}{
A list of clustering outputs or the output of the \code{DiffGenes} function. The first element of the list
will be used as a reference in \code{ReorderToReference}.
}
\item{nrclusters}{
If List is the output several clustering methods, it has to be provided in how many clusters to cut the dendrograms in.
}
  \item{fusionsLog}{
To be handed to \code{ReorderToReference}.
}
  \item{WeightClust}{
To be handed to \code{ReorderToReference}.
}
\item{names}{
Names of the methods or clusters.
}
}
\value{
A vector containing the shared compounds of all listed elements.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
data(geneMat)
data(GeneInfo)
data(ListGO)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)

L=list(MCF7_F,MCF7_T)
names=c('FP','TP')

Comps=SharedComps(List=L,nrclusters=7,fusionsLog=FALSE,WeightClust=FALSE,names=names)

}
