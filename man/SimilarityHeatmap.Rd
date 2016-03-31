\name{SimilarityHeatmap}
\alias{SimilarityHeatmap}
\title{
A heatmap of similarity values between compounds
}
\description{
The function \code{SimilarityHeatmap} plots the similarity values
between compounds. The darker the shade, the more similar compounds
are. The option is available to set a cutoff value to highlight the most
similar compounds.
}
\usage{
SimilarityHeatmap(Data,type=c("data","clust","sim","dist"),
distmeasure="tanimoto",normalize=FALSE,method="Q",cutoff=NULL,
percentile=FALSE,plottype="new",location=NULL)

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
If data is of type "clust", the distance matrix is extracted from the result and transformed 
to a similarity matrix. Possibly a range normalization is performed. If data is of type "dist",
it is also transformed to a similarity matrix and cluster is performed on the distances. If data
is of type "sim", the data is tranformed to a distance matrix on which clustering is performed.
Once the similarity mattrix is obtained, the cutoff value is applied and a heatmap is drawn. If
no cutoff value is desired, one can leave the default NULL specification.

}
\value{
A heatmap with the names of the compounds on the right and bottom and
a dendrogram of the clustering at the left and top.
}
\author{
Marijke Van Moerbeke
}
\examples{
\dontrun{
data(fingerprintMat)


MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55)

SimilarityHeatmap(Data=MCF7_F,type="clust",cutoff=0.90,percentile=TRUE)
SimilarityHeatmap(Data=MCF7_F,type="clust",cutoff=0.75,percentile=FALSE)

}
}

\keyword{Heatmap}
