\name{HeatmapSelection}
\alias{HeatmapSelection}
\title{
A function to select a group of compounds via the similarity heatmap.
}
\description{
The function \code{HeatmapSelection} plots the similarity values
between compounds. The plot is similar to the one produced by 
\code{SimilarityHeatmap} but without the dendrograms on the sides.
The function is rather explorative and experimental and is to be
used with some caution. By clicking in the plot, the user can select 
a group of compounds of interest. See more in \code{Details}.
}
\usage{
HeatmapSelection(Data,type=c("data","dist","clust","sim"),
distmeasure="tanimoto",normalize=FALSE,method="Q",cutoff=NULL,
percentile=FALSE,dendrogram=NULL,width=7,height=7)
}
\arguments{
  \item{Data}{
The data of which a heatmap should be drawn. 
}
  \item{type}{
The type of data. Data can either be the data itself ("data"), the outcome 
of a clustering method ("clust'), a distance matrix ("dist") or a 
similarity matrix ("sim").
}
\item{distmeasure}{
If type is "data", a distance measure for the clustering should be specified.
}
\item{normalize}{
Logical. If type is "data", it can be specified whether the data should be normalized.
}
\item{method}{
If type is "data" and normalize is TRUE, a method for normalization should be specified.
See \code{Normalization}.
}
  \item{cutoff}{
Optional. If a cutoff value is specified, all values lower are put to zero while all
other values are kept. This helps to highlight the most similar compounds.
}
  \item{percentile}{
Logical. The cutoff value can be a percentile. If one want the cutoff value to be the 90th
percentile of the data, one should specify cutoff = 0.90 and percentile = TRUE.
}
  \item{dendrogram}{
Optional. If the clustering results of the data is already available and should not be recalculated,
this results can be provided here. Otherwise, it will be calculated given the data. This is necessary
to have the compounds in their order of clustering on the plot.
}
\item{width}{
The width of the plot to be made. This can be adjusted since the default size might not show a clear picture.
}
\item{height}{
The height of the plot to be made. This can be adjusted since the default size might not show a clear picture.
}
}
\details{
A similarity heatmap is created in the same way as in \code{SimilarityHeatmap}. The user is now free
 to select two points on the heatmap. It is advised that these two points are in opposite corners of a square 
 that indicates a high similarity among the compounds. The points do not have to be the exact corners of the
group of interest, a little deviation is allowed as rows and columns of the selected subset of the matrix with 
sum equal to 1 are filtered out. A sum equal to one, implies that the compound is only similar to itself.

The function is meant to be explorative but is experimental. The goal was to make the selection of interesting compounds
easier as sometimes the labels of the dendrograms are too distorted to be read. If the figure is exported to
a pdf file with an appropriate width and height, the labels can be become readable again.
}
\value{
A heatmap with the names of the compounds on the right and bottom.
Once points are selected, it will return the names of the compounds that are in the 
selected square provided that these show similarity among each other.
}
\author{
Marijke Van Moerbeke
}
\examples{
\dontrun{
data(fingerprintMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55)

HeatmapSelection(Data=MCF7_F$DistM,type="dist",cutoff=0.90,percentile=TRUE,
dendrogram=MCF7_F,width=7,height=7)
}
}
\keyword{Heatmap}
