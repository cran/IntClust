\name{LabelPlot}
\alias{LabelPlot}
\title{
Coloring specific leaves of a dendrogram
}
\description{
Just as the function \code{ClusterCols}, \code{LabelCols} as its own plotting
function \code{LabelPlot} which plots the dendrogram.
}
\usage{
LabelPlot(Data, Sel1, Sel2 = NULL, col1, col2 = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
The result of a method which contains the dendrogram to be colored.
}
  \item{Sel1}{
The selection of objects to be colored.
}
  \item{Sel2}{
An optional second selection to be colored.
}
  \item{col1}{
The color for the first selection.
}
  \item{col2}{
The color for the optional second selection.
}
  \item{\dots}{
Other options which can be given to the plot function.
}
}

\value{
A plot of the dendrogram of which the leaves of the selection(s)
are colored.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)

ClustF_6=cutree(MCF7_F$Clust,6)

SelF=rownames(fingerprintMat)[ClustF_6==6]
SelF

LabelPlot(MCF7_F,Sel1=SelF,Sel2=NULL,col1='darkorchid')

}
\keyword{ Dendrogram }

