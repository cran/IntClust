\name{BinFeaturesPlot}
\alias{BinFeaturesPlot}
\title{
Plot of a selection of features
}
\description{
The function \code{BinFeaturesPlot} plots the binary data matrix for a 
selection of features. It is possible to separate between compounds
of interest and the other compounds. This is a visualization to see which
characteristics are (not) expressed in a specific cluster.
}
\usage{
BinFeaturesPlot(LeadCpds, OrderLab, Features, Data, ColorLab, nrclusters = NULL, 
cols = NULL, name = c("FP"),colors1 = c("gray90", "blue"), colors2 = c("gray90",
"green"),margins=c(5.5,3.5,0.5,5.5),plottype="new",location=NULL)
}
\arguments{
  \item{LeadCpds}{
A character vector containing the compounds one wants to separate 
from the others.
}
  \item{Features}{
A character vector containing the selection of features to be plotted.
}
  \item{Data}{
The data matrix the features are derived from.
}
  \item{OrderLab}{
Optional. If the compounds are to set in a specific order of a specific clustering.
}
  \item{ColorLab}{
Optional. The clustering result that determines the clusters of the labels of the
objects in the plot.
}
  \item{nrclusters}{
Optional. The number of clusters to consider if ColorLab is specified.
}
  \item{cols}{
The colors for the clusters of the labels of the objects as determined by ColorLab.
}
  \item{name}{
The overall name to give the features.
}
  \item{colors1}{
 Colors to indicate the present and abscence of features of the compounds not in LeadComps.
  }
  \item{colors2}{
 Colors to indicate the present and abscence of features of the LeadComps.
 
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
}
\value{
A plot indicating the values of the features of the LeadCpds in green 
and those of the others in blue.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
data(geneMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55)
\dontrun{
MCF7_Interactive=ChooseFeatures(Interactive=TRUE,LeadCpds=NULL,ClusterResult=MCF7_T,
ClusterColors=MCF7_F,BinData=list(fingerprintMat),Datanames=c("FP"),geneMat,
topChar = 20, topG = 20,nrclusters=7,N=1)

Lead=MCF7_Interactive$"Choice 1"$Compounds$LeadCpds
Feat=MCF7_Interactive$"Choice 1"$Characteristic$FP

BinFeaturesPlot(LeadCpds=Lead,Features=Feat,Data=fingerprintMat,OrderLab=MCF7_F,ColorLab=MCF7_F,
nrclusters=7,cols=Colors1,name=c("FP"),margins=c(5.5,3.5,0.5,5.5),plottype="new",location=NULL)
}
}
\keyword{Features Plot}
