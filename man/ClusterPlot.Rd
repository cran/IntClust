\name{ClusterPlot}
\alias{ClusterPlot}
\title{
Plot a dendrogram with leaves colored by a result of choice
}
\description{
The above described function \code{ClusterCols} is used in the function
Clusterplot which actually plots the dendrogram made by ClusterCols.
Further, given the outputs of any other functions, it is capable of 
selecting the elements needed for ClusterCols.
}
\usage{
ClusterPlot(Data1, Data2=NULL, nrclusters = NULL, cols = NULL,plottype="new",
location=NULL,ColorComps = NULL, ...)
}
\arguments{
  \item{Data1}{
The resulting clustering of a method which contains the dendrogram to
be colored.
}
  \item{Data2}{
Optional. The resulting clustering of another method , i.e. the resulting clustering
on which the colors should be based.
}
  \item{nrclusters}{
Optional. The number of clusters to cut the dendrogram in. If not specified
the dendrogram will be drawn without colors to discern the different clusters.
}
  \item{cols}{
The colors for the clusters if nrclusters is specified.
}
  \item{plottype}{
Should be one of "pdf","new" or "sweave". If "pdf", a location should be provided in
"location" and the figure is saved there. If "new" a new graphic device is opened and if "sweave",
the figure is made compatible to appear in a sweave or knitr document.
}
\item{location}{
If plottype is "pdf", a location should be provided in "location" and the figure is saved there. 
}
\item{ColorComps}{
If only a specific set of compounds need to be highlighted, this can be specified here. 
The compounds should be given in a character vector. If specified, all other compound labels will be colored black.
}
  \item{\dots}{
Other options which can be given to the plot function.
}
}
\details{
This function relies on the internal \code{ClusterCols} function.
}
\value{
A plot of the dendrogram of the first clustering result with colored leaves.
If a second clustering result is given in Data2, the colors are based on this clustering result.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
data(Colors1)

MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)

ClusterPlot(MCF7_T ,nrclusters=7,cols=Colors1,plottype="new",location=NULL,
main="Clustering on Target Predictions: Dendrogram",ylim=c(-0.1,1.8))
}
\keyword{ Dendrogram }
