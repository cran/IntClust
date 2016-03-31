\name{BoxPlotDistance}
\alias{BoxPlotDistance}
\title{
Box plots of one distance matrix categorized against another distance matrix.
}
\description{
Given two distance matrices, the function categorizes one distance matrix and
produces a box plot from the other distance matrix against the created 
categories. The option is available to choose one of the plots or to have both
plots. The function also works on outputs from ADEC and CEC functions which do 
not have distance matrices but incidence matrices.
}
\usage{
BoxPlotDistance(Data1, Data2, type=c('data','dist','clusters'), distmeasure="tanimoto",
normalize=FALSE,method=NULL,lab1, lab2, limits1=NULL, limits2=NULL, plot = 1,
StopRange=FALSE,plottype="new",location=NULL)
}
\arguments{
  \item{Data1}{
The first data matrix, cluster outcome or distance matrix to be plotted.
}
  \item{Data2}{
The second data matrix, cluster outcome or distance matrix to be plotted.
}
  \item{type}{
Type indicates the kind of data provided as input. Should be one of "data", "cluster" or "distance".
The type "cluster" is used if the data is the output of one of the integrated data cluster
functions of the package.
}
  \item{distmeasure}{
Choice of metric for the dissimilarity matrix (character) and should only be specified if type is "data".
Should be one of "tanimoto", "euclidean", "jaccard","hamming".
}
  \item{normalize}{
Logical. Indicates whether to normalize the distance matrices or not.
This is recommended if different distance types are used. More details
on normalization in \code{Normalization}.
}
  \item{method}{
A method of normalization. Should be one of "Quantile","Fisher-Yates",
"standardize","Range" or any of the first letters of these names.
}
  \item{lab1}{
The label to plot for Data1.
}
  \item{lab2}{
The label to plot for Data2.
}
  \item{limits1}{
The limits for the categories of Data1.
}
  \item{limits2}{
The limits for the categories of Data2.
}
  \item{plot}{
The type of plots:
1 - Plot the values of Data1 versus the categories of Data2.
2 - Plot the values of Data2 versus the categories of Data1.
3 - Plot both types 1 and 2.
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
\value{
One/multiple box plots.
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


BoxPlotDistance(MCF7_F,MCF7_T,type="cluster",lab1="F",lab2="T",limits1=c(0.3,0.7),
limits2=c(0.3,0.7),plot=1,StopRange=FALSE,plottype="new", location=NULL)
}
\keyword{ Box Plot}