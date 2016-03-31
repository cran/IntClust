\name{WonM}
\alias{WonM}
\title{
Weighting on Membership
}
\description{
Weighting on membership is performed with the \code{WonM} function. The first step is to compute the
appropriate distance matrices for each data source and to use these for hierarchical clustering. This
is executed with the agnes function and the ward link. The user may specify a range of values for
the number of clusters to cut the resulting dendrograms in. For each value of number of clusters,
an incidence matrix is computed and these are added for each data source separately. Eventually,
the sums of the incidence matrices are joined together as well, resulting in one consensus matrix.
Hierarchical clustering is performed on the consensus matrix to obtain the final clustering result.
}
\usage{
WonM(List,type=c("data","dist","clusters"), distmeasure = c("tanimoto",
"tanimoto"),normalize=FALSE,method=NULL,nrclusters = seq(5, 25, 1), clust =
"agnes", linkage=c("flexible","flexible"),alpha=0.625,StopRange=FALSE)
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
Should be of "tanimoto", "euclidean", "jaccard","hamming"..
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
A sequence of the number of clusters to cut the dendrogram in.
}
  \item{clust}{
Choice of clustering function (character). Defaults to "agnes".
}
  \item{linkage}{
A vector with the choice of inter group dissimilarity (character) for each data set.
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
\value{
The returned value is list with four elements:
\item{DistM}{A list with the distance matrix for each
data structure}
\item{ClustSep }{The hierarchical clustering result on each
data set}
\item{Consensus }{The computed consensus matrix over
all data sources}
\item{Clust}{The resulting clustering}
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
L=list(fingerprintMat,targetMat)

MCF7_WonM=WonM(L,type="data",distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,
method=NULL,nrclusters=seq(5,25),clust="agnes",linkage=c("flexible","flexible"),
alpha=0.625,StopRange=FALSE)
}
\note{
For now, only hierarchical clustering with the \code{agnes} 
function is implemented.
}
\keyword{ Weighting on Membership }
\keyword{ Integrative Clustering }
