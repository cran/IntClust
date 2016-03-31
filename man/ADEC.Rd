\name{ADEC}
\alias{ADEC}
\title{
Aggregated Data Ensemble Clustering 
}
\description{
Function \code{ADEC} performs which the functions ADECa, ADECb and ADECc
is specified by the user.
}
\usage{
ADEC(List, distmeasure = "tanimoto",normalize=FALSE,method=NULL, t = 10,
r = NULL, nrclusters = NULL, clust = "agnes", linkage = "ward",alpha=0.625
,ResampleFeatures=TRUE)
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
on normalization in \code{Normalization}.
}
  \item{method}{
A method of normalization. Should be one of "Quantile","Fisher-Yates",
"standardize","Range" or any of the first letters of these names.
}
  \item{t}{
The number of iterations.
}
  \item{r}{
Optional. The number of features to take for the random sample.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in. 
If a sequence is specified either ADECb or ADECc is performed. 
A fixed number of clusters defaults to ADECa
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
\item{ResampleFeatures}{
Logical. Whether the features should be resamples. If TRUE, either ADECa or ADECc is performed.
}
}
\details{
See the details of ADECa, ADECb and ADEDc for more information.
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
MCF7_ADECa=ADEC(L,distmeasure="tanimoto",normalize=FALSE,method=NULL,t=25,r=NULL,
nrclusters=7,clust="agnes",linkage="ward",alpha=0.625,ResampleFeatures=TRUE)
}
\note{
For now, only hierarchical clustering with the \code{agnes} 
function implemented.
}
\seealso{
\code{\link{ADECa}}, \code{\link{ADECb}},\code{\link{ADECc}}
}
\keyword{ Aggregated Data Ensemble Clustering  }
\keyword{ Integrative Clustering }
