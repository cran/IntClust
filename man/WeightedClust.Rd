\name{WeightedClust}
\alias{WeightedClust}
\title{
Weighted clustering
}
\description{
Weighted clustering is performed with the function \code{WeightedClust}. Given a
list of the data matrices, a dissimilarity matrix is computed of each with
the provided distance measures. These matrices are then combined resulting
in a weighted dissimilarity matrix. Hierarchical clustering is performed
on this weighted combination with the agnes function and the ward link
}
\usage{
WeightedClust(List,type=c("data","dist","clusters"),
distmeasure = c("tanimoto", "tanimoto"),normalize=FALSE,method=NULL,
weight = seq(1, 0, -0.1), WeightClust = 0.5, clust="agnes",
linkage = "ward",alpha=0.625,StopRange=FALSE)
}
\arguments{
  \item{List}{
A list of matrices of the same type. It is assumed the rows are corresponding
with the objects.
}
  \item{type}{
Type indicates whether the provided matrices in "List" are either data matrices, distance
matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
matrices is skipped and if type="clusters" the single source clustering is skipped.
Type should be one of "data", "dist" or"clusters".			
}
  \item{distmeasure}{
A vector of the distance measures to be used on each data matrix.
Should be of "tanimoto", "euclidean", "jaccard","hamming".
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
 \item{weight}{
 Optional. A list of different weight combinations for the data sets in List.
 If NULL, the weights are determined to be rqual for each data set.
 It is further possible to fix weights for some data matrices and to
 let it vary randomly for the remaining data sets. An example is provided in the details.
}
  \item{WeightClust}{
A weight for which the result will be put aside of the other results. This was
done for comparative reason and easy access.
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
  \item{StopRange}{
Logical. Indicates whether the distance matrices with values not between zero and one should be standardized to have so.
If FALSE the range normalization is performed. See \code{Normalization}. If TRUE, the distance matrices are not changed.
This is recommended if different types of data are used such that these are comparable.
}
}
\details{
The weight combinations should be provided as elements in a list. For three data
matrices an example could be: weights=list(c(0.5,0.2,0.3),c(0.1,0.5,0.4)). To provide
a fixed weight for some data sets and let it vary randomly for others, the element "x" 
indicates a free parameter. An example is weights=list(c(0.7,"x","x")). The weight 0.7 
is now fixed for the first data matrix while the remaining 0.3 weight will be divided over
the other two data sets. This implies that every combination of the sequence from 0 to 0.3
with steps of 0.1 will be reported and clustering will be performed for each.
}
\value{
The returned value is a list of four elements:
\item{DistM}{A list with the distance matrix for each
data structure}
\item{WeightedDist }{A list with the weighted distance matrices}
\item{Results}{The hierarchical clustering result for each element in 
IncidenceComb}
\item{Clust}{The result for the weight specified in Clustweight}
The value has class 'Weighted'
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
L=list(fingerprintMat,targetMat)

MCF7_Weighted=WeightedClust(L,type="data", distmeasure=c("tanimoto","tanimoto"),
normalize=FALSE,method=NULL,weight=seq(1,0,-0.1),WeightClust=0.5,clust="agnes",linkage="ward"
,alpha=0.625,StopRange=FALSE)

}
\note{
For now, only hierarchical clustering with the \code{agnes} 
function is implemented.
}

\keyword{ Weighted Clustering }
