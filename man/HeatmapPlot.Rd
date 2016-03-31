\name{HeatmapPlot}
\alias{HeatmapPlot}
\title{
Comparing two clustering results with a heatmap
}
\description{
The \code{HeatmapCols} function calculates the distance between two outputs of
clustering methods and plots the resulting heatmap. The function heatmap.2 is
called upon to make the actual plot of the heatmap. It is noted that for thi
s function the number of colors should be one more than the number of clusters
to color the so called zero cells in the distance matrix.
}
\usage{
HeatmapPlot(Data1, Data2, names = NULL, nrclusters = NULL, 
cols = NULL,plottype="new",location=NULL)
}
\arguments{
  \item{Data1}{
The resulting clustering of method 1.
}
  \item{Data2}{
The resulting clustering of method 2.
}
  \item{names}{
The names of the objects in the data sets.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in.
}
  \item{cols}{
The colors to be used for the clusters.
}
  \item{plottype}{
Should be one of "pdf","new" or "sweave". If "pdf", a location should be provided in
"location" and the figure is saved there. If "new" a new graphic device is opened and if "sweave",
the figure is made compatible to appear in a sweave or knitr document, i.e. no new device is
opened and the plot appears in the current device or document.
}
  \item{location}{
If plottype is "pdf", a location should be provided in "location" and the figure is saved there. 
}
}
\details{
Another way to compare to methods is via an adaptation of heatmaps. The input of this function is
the resulting clustering (the Clust element of the list) of two methods and can be seen as: method
1 versus method 2. The dendrograms are cut into a specific number of clusters. Each cluster
of method 2 and its members are given a distinct color represented by a number. These are the
clusters to which a comparison is made. A matrix is set up of which the columns are determined
by the ordering of clustering of method 2 and the rows by the ordering of method 1. Every column
represent one object just as every row and every column represent the color of its cluster. A function
visits every cell of the matrix. If the objects represented by the cell are still together in a cluster, the
color of the column is passed to the cell. This creates the distance matrix which can be given to the
HeatmapCols function to create the heatmap.
}
\value{
A heatmap based on the distance matrix created by the function with the dendrogram
of method 2 on top of the plot and the one from method 1 on the left. The names
of the compounds are depicted on the bottom in the order of clustering of method
2 and on the right by the ordering of method 1. Vertically the cluster of method
2 can be seen while horizontally those of method 1 are portrayed.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
data(Colors2)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
clust="agnes",linkage="ward",gap=FALSE,maxK=55)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
clust="agnes",linkage="ward",gap=FALSE,maxK=55)


L=list(MCF7_F,MCF7_T)
names=c("FP","TP")

HeatmapPlot(MCF7_T,MCF7_F,names=rownames(fingerprintMat),nrclusters=7,cols=Colors2,plottype="new",
location=NULL)

}
\keyword{ Heatmap }

