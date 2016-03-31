\name{CharacteristicFeatures}
\alias{CharacteristicFeatures}
\title{
Determining the characteristic features of a cluster
}
\description{
The function \code{CharacteristicFeatures} requires as input a list of one or multiple clustering results.
It is capable of selecting the binary features which determine a cluster with the help of the fisher's exact test.
}
\usage{
CharacteristicFeatures(List,Selection=NULL,BinData,ContData = NULL,
Datanames=NULL,nrclusters=NULL,sign=0.05,topC=NULL,fusionsLog=TRUE,
WeightClust=TRUE,names=NULL)
}
\arguments{
  \item{List}{
A list of the clustering outputs to be compared. The first element
of the list will be used as the reference in \code{ReorderToReference}.
}
 \item{Selection}{
If differential gene expression should be investigated for a specific selection of compounds, this selection can be provided here.
Selection can be of the type "character" (names of the compounds) or "numeric" (the number of specific cluster).
 }
  \item{BinData}{
A list of the binary feature data matrices. These will be evaluated with the fisher's extact test.
}
  \item{ContData}{
A list of continuous data sets of the compounds. These will be evaluated with the t-test.
}
  \item{Datanames}{
A vector with the names of the binary data matrices.
}
  \item{nrclusters}{
Optional. The number of clusters to cut the dendrogram in. The number of clusters should not be
specified if the interest lies only in a specific selection of compounds which is known by name. 
Otherwise, it is required.
}
  \item{sign}{
The significance level to be handled.
}
  \item{topC}{
Overrules sign. The number of features to display for each cluster. 
If not specified, only the significant genes are shown.
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
\item{Characteristics}{A list with an element per defined binary data matrix in BinData and continuous data in ContData. Each element is again a list with the elements TopFeat (a table with information on the top
features) and AllFeat (a table with information on all features)}
}
\author{
Marijke Van Moerbeke
}
\examples{
\dontrun{
data(fingerprintMat)
data(targetMat)
data(geneMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)

L=list(MCF7_T ,MCF7_F)

MCF7_Char=CharacteristicFeatures(List=L,Selection=NULL,BinData=list(fingerprintMat,
targetMat),Datanames=c("F","T"),nrclusters=7,topC=NULL,sign=0.05,fusionsLog=TRUE,WeightClust=TRUE,
names=c("F","T"))
}
}
