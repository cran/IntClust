\name{CompareInteractive}
\alias{CompareInteractive}
\title{
Interactive comparison of clustering results for a specific cluster or method.
}
\description{
A visual comparison of all methods is handy to see which compounds will always
cluster together independent of the applied methods. The function \code{CompareInteractive} 
plots the comparison over the specified methods. A cluster or method can than be identified by
clicking and is plotted separately against the single source or other specified methods.
}
\usage{

CompareInteractive(ListM,ListS,nrclusters=NULL,cols=NULL,fusionsLogM
=FALSE,fusionsLogS=FALSE,WeightClustM=FALSE,WeightClustS=FALSE,
namesM=NULL,namesS=NULL,marginsM=c(2,2.5,2,2.5),marginsS=c(8,2.5,2,2.5)
,Interactive=TRUE,N=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ListM}{
A list of the multiple source clustering or other methods to be compared and from which a cluster or 
method will be identified. The first element of the list will be used as the reference in
\code{ReorderToReference}.
}
  \item{ListS}{
A list of the single source clustering or other methods the identified result will be
compared to. The first element of the list will be used as the reference in
\code{ReorderToReference}.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in.
}
  \item{cols}{
The hex codes of the colors to be used.
}
  \item{fusionsLogM}{
The fusionsLog parameter for the elements in ListM. To be handed to \code{ReorderToReference}.
}
  \item{fusionsLogS}{
The fusionslog parameter for the elements in ListS. To be handed to \code{ReorderToReference}.
}
  \item{WeightClustM}{
The WeightClust parameter for the elements in ListM. To be handed to \code{ReorderToReference}.
}
  \item{WeightClustS}{
The WeightClust parameter for the elements in ListS. To be handed to \code{ReorderToReference}.
}
  \item{namesS}{
Optional. Names of the single source clusterings to be used as labels for the columns.
}
  \item{namesM}{
Optional. Names of the multiple source clusterings to be used as labels for the columns.
}
  \item{marginsM}{
Optional. Margins to be used for the plot for the elements is ListM after the identification.
}
  \item{marginsS}{
Optional. Margins to be used for the plot for the elements is ListS after the identification.
}
  \item{Interactive}{
Optional. Do you want an interactive plot? Defaults to TRUE, if not the function provides
the same as \code{ComparePlot} for the elements in ListM.
}
  \item{N}{
The number of methods/clusters you want to identify.
}
  \item{\dots}{
Other options which can be given to the \code{color2D.matplot} function.
}
}
\details{
This function relies on \code{ComparePlot} to plot the results. 
}
\value{
The returned value is a plot of the comparison of the elements of ListM. On this plot multiple clusters and/or methods can 
be identified. Click on a cluster of a specific method to see how that cluster of that method compares to the elements in ListS.
Click left next to a row to identify a all cluster of a specific method. A new plotting window will appear for every identification.
}
\author{
Marijke Van Moerbeke
}
\examples{
\dontrun{
data(fingerprintMat)
data(targetMat)
data(Colors2)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)

L=list(fingerprintMat,targetMat)

MCF7_W=WeightedClust(L,type="data",distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,
method=NULL,weight=seq(1,0,-0.1),WeightClust=0.5,clust="agnes",linkage="ward",StopRange=FALSE)

ListM=list(MCF7_W)
namesM=c(seq(1.0,0.0,-0.1))

ListS=list(MCF7_F,MCF7_T)
namesS=c("FP","TP")

CompareInteractive(ListM,ListS,nrclusters=7,cols=Colors2,fusionsLogM=FALSE,
fusionsLogS=FALSE,WeightClustM=FALSE,WeightClustS=TRUE,namesM,namesS,
marginsM=c(2,2.5,2,2.5),marginsS=c(8,2.5,2,2.5),Interactive=TRUE,N=1)
}
}
\seealso{
\code{\link{ComparePlot}}
}
