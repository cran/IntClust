\name{ADC}
\alias{ADC}
\title{
Aggregated Data Clustering
}
\description{
In order to perform aggregated data clustering, the \code{ADClust}
function was written. The data matrices are aggregated into one
and hierarchical clustering is performed.
}
\usage{
ADC(List, distmeasure = "tanimoto",normalize=FALSE,method=NULL,clust = "agnes",
linkage = "ward",alpha=0.625)
}
\arguments{
  \item{List}{
A list of data matrices of the same type. It is assumed the 
rows are corresponding with the objects.
}
  \item{distmeasure}{
Choice of metric for the dissimilarity matrix (character).
Should be one of "tanimoto", "euclidean", "jaccard", "hamming".
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
In order to perform aggregated data clustering, the \code{ADC}
function was written. A list of data matrices of the same type 
(continuous or binary) is required as input which are combined 
into a single (larger) matrix. Hierarchical clustering is performed
with the \code{agnes} function and the ward link on the resulting
data matrix and an applicable distance measure is indicated by the
user.
}
\value{
The returned value is a list with the following three elements.
\item{AllData}{Fused data matrix of the data matrices}
\item{DistM}{The distance matrix computed from the AllData element}
\item{Clust}{The resulting clustering}
The value has class 'ADC'. The Clust element will be of interest for
 further applications.
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
MCF7_ADC=ADC(L,distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",
linkage="ward",alpha=0.625)
}
\note{
For now, only hierarchical clustering with the \code{agnes} 
function is implemented.
}
\keyword{ Aggregated Data Clustering }
\keyword{ Integrative Clustering }
