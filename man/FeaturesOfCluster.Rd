\name{FeaturesOfCluster}
\alias{FeaturesOfCluster}
\title{
Lists all features present in a selected cluster of compounds
}
\description{
The function \code{FeaturesOfCluster} lists the number of features compounds of
the cluster have in common. A threshold can be set selecting among how many
compounds of the cluster the features should be shared. An optional plot of
the features is available.
}
\usage{
FeaturesOfCluster(LeadCpds,Data,Threshold=1,Plot=TRUE,plottype="new",location=NULL)
}
\arguments{
  \item{LeadCpds}{
A character vector containing the compounds one wants to investigate in
terms of features.
}
  \item{Data}{
The data matrix the features are derived from.
}
  \item{Threshold}{
The number of compounds the features at least should be shared amongst. Default is set to 
1 implying that the features should be present in at least one of the compounds
specified in leadCpds.
}
  \item{Plot}{
 Logical. Indicates whether or not a BinFeaturesPlot should be set up for the
 selectcion of compounds and discovered features.

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
and those of the others in blue. It list all features which are present in at least
the threshold number of compounds. By including all other compounds as well, one can see
whether features are common in the compounds or rather specific for the cluster. 

Further, it returns a list with 2 items. The first indicates the number of shared features among the 
compounds. This provides an overview of which compounds are more similar than others. The second 
item is a character vector of the plotted features such that these can be retrieved for further
investigation.
}
\author{
Marijke Van Moerbeke
}
\examples{
\dontrun{
data(fingerprintMat)

Lead=rownames(fingerprintMat)[1:5]

FeaturesOfCluster(LeadCpds=Lead,Data=fingerprintMat,
Threshold=1,Plot=TRUE,plottype="new",location=NULL)
}
}
\keyword{Features Plot}
