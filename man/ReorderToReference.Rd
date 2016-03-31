\name{ReorderToReference}
\alias{ReorderToReference}
\title{
Order the outputs of the clustering methods against a reference
}
\description{
When multiple methods are performed on a data set, it is interesting to compare
their results. However, a comparison is not easily done since different
methods leads to a different ordering of the objects. The \code{ReorderToReference} 
rearranges the cluster to a reference method.
}
\usage{
ReorderToReference(List, nrclusters = NULL, fusionsLog = FALSE, WeightClust = 
FALSE, names = NULL)
}
\arguments{
  \item{List}{
A list of clustering outputs to be compared. The first element of the
list will be used as the reference.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in.
}
  \item{fusionsLog}{
Logical indicator for the fusion of clusters.
}
  \item{WeightClust}{
Optional. To be used for the outputs of CEC, WeightedClust or WeightedSimCLust.
Then only the result of the Clust element is considered.
}
  \item{names}{
Optional. Names of the methods.
}
}
\details{
It is interesting to compare the results of the methods described in the methodology. All methods
result in a dendrogram which is cut into a specific number of clusters with the \code{cutree} function.
This results in an numbering of cluster based on the ordering of the names in the data and not on
the order in which they are grouped into clusters. However, different methods lead to different
clusters and it is possible that cluster $1$ of one method will not be the cluster that has the most in
common with cluster 1 of another method. This makes comparisons rather difficult. Therefore the
ReorderToReference function was written which takes one method as a reference and rearranges the
cluster numbers of the other methods to this reference such that clusters are appointed to that cluster
they have the most in common with. The result of this function is a matrix of which the columns are
in the order of the clustering of the compounds of the referenced method and the rows represent the
methods. Each cell contains the number of the cluster the compound is in for that method compared
to the method used as a reference. This function is applied in the functions \code{SimilarityMeasure},
\code{DiffGenes}, \code{Pathways} and \code{ComparePlot}. It is a possibility that 2 or more clusters
are fused together compared to the reference method. If this is true, the function will alert the user
and will ask to put the parameter fusionsLog to true. Since \code{ReorderToReference} is often used as an internal
function, also for visualization, it will print out how many more colors should be specified for those
clusters that did not find a suitable match. This can be due to fusion or complete segregation of its 
compounds into other clusters.
}
\value{
A matrix of which the cells indicate to what cluster the compounds belong to according to the
rearranged methods. 
}
\author{
Marijke Van Moerbeke
}
\note{
The \code{ReorderToReference} function was optimized for the situations presented by the data
sets at hand. It is noted that the function might fail in a particular situation which results
in a infinite loop.
}
\examples{
data(fingerprintMat)
data(targetMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_ADC=ADC(list(fingerprintMat,targetMat),distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward")

L=list(MCF7_F,MCF7_ADC,MCF7_T)
names=c("FP","ADC","TP")

MCF7_Matrix=ReorderToReference(L,nrclusters = 7, fusionsLog = FALSE, WeightClust = 
FALSE, names = names)

}

