\name{SNFb}
\alias{SNFb}

\title{
Similarity Network Fusion - version b
}
\description{
Function \code{SNFb} performs SNF but first determines the subsets
of neighbours and then normalization is performed on the neighbours only.
The function is based on the functions \code{affinityMatrix} and \code{snf}
from the SNFtool package.
}
\usage{
SNFb(List,type=c("data","dist","clusters"),distmeasure = c("tanimoto",
"tanimoto"),normalize=FALSE,method=NULL,NN = 20, mu = 0.5,T = 20, clust =
"agnes", linkage = "ward",alpha=0.625,StopRange=FALSE)
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
  \item{method}{
A method of normalization. Should be one of "Quantile","Fisher-Yates",
"standardize","Range" or any of the first letters of these names.
}
  \item{normalize}{
Logical. Indicates whether to normalize the distance matrices or not.
This is recommended if different distance types are used. More details
on standardization in \code{Normalization}.
}
  \item{NN}{
The number of neighbours to be used in the procedure. 
}
  \item{mu}{
The parameter epsilon. The value is recommended to be between 0.3 and 0.8.
}
  \item{T}{
The number of iterations.
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
\value{
The returned value is a list with two elements:
 \item{FusedM }{The fused similarity matrix}
  \item{DistM }{The distance matrix computed by
 subtracting FusedM from one}
 \item{Clust}{The resulting clustering}

}
\references{
WANG, B., MEZLINI, M. A., DEMIR, F., FIUME, M., TU, Z., BRUDNO, M.,
HAIBE-KAINS, B., GOLDENBERG, A. (2014). Similarity Network Fusion
for aggregating data types on a genomic scale. Nature. 11(3) pp. 
333-337.
WANG, B., MEZLINI, M. A., DEMIR, F., FIUME, M., TU, Z., BRUDNO, M.,
HAIBE-KAINS, B., GOLDENBERG, A. (2014). SNFtool: Similarity Network Fusion.
R package version 2.2
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
L=list(fingerprintMat,targetMat)

MCF7_SNFb=SNFb(list(fingerprintMat,targetMat),type="data",distmeasure=c("tanimoto",
"tanimoto"),normalize=FALSE,method=NULL,NN=10,mu=0.5,T=20,clust="agnes",linkage="ward"
,alpha=0.625,StopRange=FALSE)
}
\note{
For now, only hierarchical clustering with the \code{agnes} 
function link is implemented.
}

\seealso{
\code{\link{SNF}},\code{\link{SNFa}},\code{\link{SNFc}}
}
\keyword{ Similarity Network Fusion}
\keyword{ Integrative Clustering }
