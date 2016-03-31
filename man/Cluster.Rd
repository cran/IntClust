\name{Cluster}
\alias{Cluster}
\title{
Perform clustering on a single data source
}
\description{
The function \code{Cluster} was written to perform clustering on
a single source of information, i.e one data matrix. The option is
available to compute the gap statistic to determine the optimal
number of clusters. 
}
\usage{
Cluster(Data,type=c("data","dist"), distmeasure = "tanimoto",
normalize=FALSE,method=NULL, clust = "agnes", linkage ="ward",alpha=0.625
,gap = TRUE,maxK = 50,StopRange=FALSE)
}
\arguments{
  \item{Data}{
 A matrix containing the data. It is assumed the rows are corresponding
 with the objects.
}
	\item{type}{
Type indicates whether the provided matrix in "Data" is either a data or a distance
matrix obtained from the data. If type="dist" the calculation of the distance
matrix is skipped. Type should be one of "data" or "dist".	
}
  \item{distmeasure}{
Choice of metric for the dissimilarity matrix (character).
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
  \item{clust}{
Choice of clustering function (character). Defaults to "agnes".
}
  \item{linkage}{
Choice of inter group dissimilarity (character). Defaults to "ward".
}
\item{alpha}{
The parameter alpha to be used in the "flexible" linkage of the agnes function. 
Defaults to 0.625 and is only used if the linkage is set to "flexible"
}
  \item{gap}{
Logical. Indicator if gap statistics should be computed. 
Setting to $FALSE$ will greatly reduce the computation time.
}
  \item{maxK}{
The maximum number of clusters to be considered during the gap. 
}
  \item{StopRange}{
Logical. Indicates whether the distance matrices with values not between zero and one should be standardized to have so.
If FALSE the range normalization is performed. See \code{Normalization}. If TRUE, the distance matrices are not changed.
This is recommended if different types of data are used such that these are comparable.
}
}
\details{
The gap statistic is determined by the criteria described by the cluster package:
firstSEmax, globalSEmax, firstmax,globalmax, Tibs2001SEmax. The number of 
iterations is set to a default of 500. The implemented distances to be used for
the dissimilarity matrix are jaccard, tanimoto and euclidean. The jaccard distances
were computed with the \code{dist.binary(\ldots,method=1)} function in the ade4
package and the euclidean ones with the \code{daisy} function in again the cluster
package. The Tanimoto distances were implemented manually. 
}
\value{
The returned value is a list with two elements:
\item{DistM}{The distance matrix of the data matrix}
\item{Clust}{The resulting clustering}
If the gap option was indicated to be true, another 3 elements are joined to the list.
Clust\_gap contains the output from the function to compute the gap statistics and
gapdata is a subset of this output. Both can be used to make plots to visualize the 
gap statistic. The final component is k which is a matrix containing the optimal number
of clusters determined by each criterion mentioned earlier.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",alpha=0.625,gap=FALSE,maxK=55
,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",alpha=0.625,gap=FALSE,maxK=55
,StopRange=FALSE)

}
\note{
For now, the only option is to carry out agglomerative
hierarchical clustering as it was implemented
in the \code{agnes} function in the cluster package.
}
\keyword{ Clustering }
