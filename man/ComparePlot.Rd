\name{ComparePlot}
\alias{ComparePlot}
\title{
Comparison of clustering results over multiple results
}
\description{
A visual comparison of all methods is handy to see which compounds will always
cluster together independent of the applied methods. To this aid the function 
\code{ComparePlot} has been written.
}
\usage{
ComparePlot(List, nrclusters = NULL, cols = NULL, fusionsLog = FALSE, 
WeightClust = FALSE, names = NULL,margins = c(8.1, 3.1, 3.1, 4.1),
plottype="new",location=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
A list of the outputs from the methods to be compared. The first element of the
list will be used as the reference in \code{ReorderToReference}.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in.
}
  \item{cols}{
The hex codes of the colors to be used.
}
  \item{fusionsLog}{
To be handed to \code{ReorderToReference}.
}
  \item{WeightClust}{
To be handed to \code{ReorderToReference}.
}
  \item{names}{
Optional. Names of the methods to be used as labels for the columns.
}
  \item{margins}{
Optional. Margins to be used for the plot.
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
  \item{\dots}{
Other options which can be given to the \code{color2D.matplot} function.
}
}
\details{
This function makes use of the functions \code{ReorderToReference} and \code{Colorsnames}. 
Given a list with the outputs of several methods, the first step is to call upon 
\code{ ReorderToReference} and to produce a matrix of which the columns are ordered according
to the ordering of the objects of the first method in the list. Each cell represent the
number of the cluster the object belongs to for a specific method indicated by the rows.
The clusters are arranged in such a way that these correspond to that one cluster of the
referenced method that they have the most in common with. The function \code{color2D.matplot}
produces a plot of this matrix but needs a vector indicating the names of the colors to be used.
This is where \code{ColorsNames} comes in. A vector of the color names of the output of the
\code{ReorderToReference} is created and handed to \code{color2D.matplot}. It is optional to adjust
the margins of the plot and to give a vector with the names of the methods which will be used
as labels for the rows in the plot. The labels for the columns are the names of the object in
the order of clustering of the referenced method. Further, the similarity measures of the methods
compared to the reference will be computed and shown on the right side of the plot.
}
\value{
A plot which translates the matrix output of the function \code{ReorderToReference} in which the columns
represent the objects in the ordering the referenced method and the rows the outputs of the 
given methods. Each cluster is given a distinct color. This way it can be easily observed which
objects will cluster together. The labels on the right side of the plot are the similarity
measures computed by \code{SimilarityMeasure}.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
data(Colors2)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)

L=list(MCF7_F,MCF7_T)
names=c("FP","TP")

ComparePlot(L,nrclusters=7,cols=Colors2,fusionsLog=TRUE,WeightClust=TRUE,names=names,
margins=c(9.1,4.1,4.1,4.1),plottype="new",location=NULL)
}
\seealso{
\code{\link{ReorderToReference}},\code{\link{ColorsNames}}
}
