\name{SimilarityMeasure}
\alias{SimilarityMeasure}
\title{
A measure of similarity for the outputs of the different methods
}
\description{
The function \code{SimilarityMeasure} computes the similarity of the methods. 
Given a list of outputs as input, the first element will be seen as the reference. 
Function \code{MatrixFunction} is called upon and the cluster numbers are 
rearranged according to the reference. Per method, \code{SimilarityMeasure} 
investigates which objects have the same cluster number in reference and said method.
This number is divided by the total number of objects and used as a similarity measure.
}
\usage{
SimilarityMeasure(List, nrclusters = NULL, fusionsLog = TRUE, 
WeightClust = TRUE, names = NULL)
}
\arguments{
  \item{List}{
A list of clustering outputs to be compared. The first element of the
list will be used as the reference in \code{ReorderToReference}.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in.~~
}
  \item{fusionsLog}{
To be handed to MatrixFunction.
}
  \item{WeightClust}{
To be handed to MatrixFunction.
}
  \item{names}{
Optional. Names of the methods.
}
}
\value{
A vector of similarity measures, one for each method given as input.
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

MCF7_SimFandT=SimilarityMeasure(L,nrclusters=7,fusionsLog=TRUE,WeightClust=TRUE,
names=names)

}
\keyword{ Similarity Meausure }

