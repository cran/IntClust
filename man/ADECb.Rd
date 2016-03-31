\name{ADECb}
\alias{ADECb}
\title{
Aggregated Data Ensemble Clustering - version b
}
\description{
Function \code{ADECb} performs aggregated data ensemble clustering in which
in every iteration the total number of features are used in the clustering 
procedure. However, the function is capable of cutting the resulting 
dendrogram several times, each time into a different number of cluster.
}
\usage{
ADECb(List, distmeasure = "tanimoto",normalize=FALSE,method=NULL,
nrclusters = seq(5, 25, 1),  clust = "agnes", linkage = "ward",
alpha=0.625)
}

\arguments{
  \item{List}{
A list of data matrices of the same type. It is assumed the 
rows are corresponding with the objects.
}
  \item{distmeasure}{
The distance measure to be used on the fused data matrix (character).
Should be one of "tanimoto", "euclidean", "jaccard","hamming".
}
  \item{normalize}{
Logical. Indicates whether to normalize the distance matrices or not.
This is recommended if different distance types are used. More details
on normalization in  \code{Normalization}.
}
  \item{method}{
A method of normalization. Should be one of "Quantile","Fisher-Yates",
"standardize","Range" or any of the first letters of these names.
}
  \item{nrclusters}{
A sequence of numbers of clusters to cut the dendrogram in.
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
}
\details{
ADECb starts with the merging of the data matrices into one larger data matrix. 
Then, ensemble clustering is performed on the fused data. This comes down to 
repeatedly applying hierarchical clustering. All features will be used in every iteration. 
Variation is inserted by not splitting the dendrogram a single time into one specific 
number of clusters but multiple times and for a range of numbers of clusters.
More information can be found in Fodeh et al. (2013).
}
\value{
The returned value is a list with the following three elements.
\item{AllData}{Fused data matrix of the data matrices}
\item{S}{The resulting co-association matrix}
\item{Clust}{The resulting clustering}
The value has class 'ADEC'. The Clust element will be of interest for further
applications.
}
\references{
FODEH, J. S., BRANDT, C., LUONG, B. T., HADDAD, A., SCHULTZ, M., MURPHY, T.,
KRAUTHAMMER, M. (2013). Complementary Ensemble Clustering of Biomedical Data.
J Biomed Inform. 46(3) pp.436-443.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
L=list(fingerprintMat,targetMat)
MCF7_ADECb=ADECb(L,distmeasure="tanimoto",normalize=FALSE,method=NULL,
nrclusters=seq(5,25),clust="agnes",linkage="ward",alpha=0.625)
}
\note{
For now, only hierarchical clustering with the \code{agnes} 
function is implemented.
}
\seealso{
\code{\link{ADEC}},\code{\link{ADECa}},\code{\link{ADECc}}
}
\keyword{ Aggregated Data Ensemble Clustering  }
\keyword{ Integrative Clustering }