\name{SelectnrClusters}
\alias{SelectnrClusters}
\title{
Determines an optimal number of clusters based on silhouette widths
}
\description{
The function \code{SelectnrClusters} determines an optimal optimal number of clusters 
based by calculating silhouettes widths for a sequence of clusters. See "Details" for a more elaborate
description.
}
\usage{
SelectnrClusters(List,type=c("data","dist","pam"),distmeasure=c("tanimoto","tanimoto")
,normalize=FALSE,method=NULL,nrclusters = seq(5, 25, 1),names=NULL,StopRange=FALSE,
plottype="new",location=NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
A list of matrices of the same type. It is assumed the rows are corresponding
with the objects.
}
  \item{type}{
Type indicates whether the provided matrices in "List" are either data matrices, distance
matrices or clustering results obtained with \code{pam} of the \pkg{cluster} package.
Type should be one of "data", "dist" or"pam".		
}
  \item{distmeasure}{
A character vector with the distance measure for each data matrix.
Should be one of "tanimoto", "euclidean", "jaccard","hamming".
}
  \item{normalize}{
Logical. Indicates whether to normalize the distance matrices or not.
This is recommended if different distance types are used. More details
on standardization in \code{Normalization}.
}
  \item{method}{
A method of normalization. Should be one of "Quantile","Fisher-Yates",
"standardize","Range" or any of the first letters of these names.
}
  \item{nrclusters}{
A sequence of numbers of clusters to cut the dendrogram in.
}
  \item{names}{
The labels to give to the elements in List.
}
  \item{StopRange}{
Logical. Indicates whether the distance matrices with values not between zero and one should be standardized to have so.
If FALSE the range normalization is performed. See \code{Normalization}. If TRUE, the distance matrices are not changed.
This is recommended if different types of data are used such that these are comparable.
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
If the object provided in List are data or distance matrices clustering around medoids is performed with the
\code{pam} function of the \pkg{cluster} package. Of the obtained pam objects, average silhouette widths are retrieved.
A silhouette width represents how well an object lies in its current cluster. Values around one are an indication 
of an appropriate clustering while values around zero show that the object might as well lie in the neighbouring cluster.
The average silhouette width is a measure of how tightly grouped the data is.  This is performed for every number of cluster
for every object provided in List. Then the average is taken for every number of clusters over the provided objects. This results
in one average value per number of clusters. The number width the maximal average silhouette width is 
chosen as the optimal number of clusters.
}
\value{
A plots are made showing the average silhouette widths of the provided objects for each number 
of clusters. Further, a list with two elements is returned:
\item{Silhouette_Widths}{A data frame with the silhouette widths for each object and the average silhouette widths per number of clusters}
\item{Optimal_Nr_of_CLusters}{The determined optimal number of cluster }
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)

List=list(fingerprintMat,targetMat)

NrClusters=SelectnrClusters(List=List,type="data",distmeasure=c("tanimoto",
"tanimoto"),nrclusters=seq(5,10),normalize=FALSE,method=NULL,names=c("FP","TP"),
StopRange=FALSE,plottype="new",location=NULL)

NrClusters
}

