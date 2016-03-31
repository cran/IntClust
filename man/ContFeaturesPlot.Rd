\name{ContFeaturesPlot}
\alias{ContFeaturesPlot}
\title{
Plot of continuous features
}
\description{
The function \code{BioassayPlot} plots the bioassay values for the compounds.
It is possible to separate between compounds of interest and the other
compounds. This is a visualization to see how bio-assays differ from cluster 
to cluster. 
}
\usage{
ContFeaturesPlot(LeadCpds,Data,nrclusters,OrderLab,ColorLab=NULL,cols=NULL,
ylab="bio-assays",AddLegend=TRUE,margins=c(5.5,3.5,0.5,8.7),plottype="new",location=NULL)
}
\arguments{
  \item{LeadCpds}{
A character vector containing the compounds one wants to separate 
from the others.
}
  \item{Data}{
The bio-assay data matrix.
}
  \item{nrclusters}{
Optional. The number of clusters to consider if ColorLab is specified.
}
  \item{OrderLab}{
Optional. If the compounds are to set in a specific order of a specific method.
}
  \item{ColorLab}{
The clustering result that determines the color of the labels of the
objects in the plot. If NULL, the labels are black.
}
  \item{cols}{
The colors for the labels of the objects.
}
  \item{ylab}{
The lable of the y-axis.
}
  \item{AddLegend}{
Logical. Indicates whether a legend should be added to the plot.  
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
A plot in which the values of the bio-assays for the LeadCpds are separeted
from the others.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(Colors2)
Comps=c("Cpd1", "Cpd2", "Cpd3", "Cpd4", "Cpd5")

Data=matrix(sample(15, size = 50*5, replace = TRUE), nrow = 50, ncol = 5)
colnames(Data)=colnames(Data, do.NULL = FALSE, prefix = "col")
rownames(Data)=rownames(Data, do.NULL = FALSE, prefix = "row")
for(i in 1:50){
	rownames(Data)[i]=paste("Cpd",i,sep="")
}

ContFeaturesPlot(LeadCpds=Comps,OrderLab=rownames(Data),ColorLab=NULL,Data=Data,
nrclusters=7,cols=Colors2,ylab="bio-assays",AddLegend=TRUE,margins=c(5.5,3.5,0.5,8.7),
plottype="new",location=NULL)



}
\keyword{Features Plot}
