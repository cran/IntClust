\name{CECa}
\alias{CECa}
\title{
Complementary Ensemble Clustering - version a
}
\description{
Function \code{CECa} performs complementary ensemble clustering in
which in every iteration the number of random samples taken is 
randomly set between m/2 and m-1 with m the total number of features.
The number of features to sample can also be specified by the user.
}
\usage{
CECa(List, distmeasure = c("tanimoto", "tanimoto"),normalize=FALSE,method=NULL,
t = 10, r = NULL, nrclusters = NULL, weight = NULL, clust = "agnes",
linkage=c("flexible","flexible"),alpha=0.625, WeightClust = 0.5,StopRange=FALSE)
}
\arguments{
  \item{List}{
A list of data matrices. It is assumed the  rows are corresponding
 with the objects.
}
  \item{distmeasure}{
A character vector with the distance measure for each data matrix.
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
A vector of the number of clusters to cut the dendrogram of each clustering result in.
}
 \item{weight}{
 Optional. A list of different weight combinations for the data sets in List.
 If NULL, the weights are determined to be equal for each data set.
 It is further possible to fix weights for some data matrices and to
 let it vary randomly for the remaining data sets. An example is provided in the details.
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
  \item{WeightClust}{
A weight for which the result will be put aside of the other results.
This was done for comparative reason and easy access.
}
  \item{StopRange}{
Logical. Indicates whether the distance matrices with values not between zero and one should be standardized to have so.
If FALSE the range normalization is performed. See \code{Normalization}. If TRUE, the distance matrices are not changed.
This is recommended if different types of data are used such that these are comparable.
}
}
\details{
Ensemble clustering is performed on each data matrix. This comes down to 
repeatedly applying hierarchical clustering. A random sample of features is taken 
in each application. Afterwards the incidence matrices are combined in 
a weighted sum and hierarchical clustering is performed once more.
More information can be found in Fodeh et al. (2013).

The weight combinations should be provided as elements in a list. For three data
matrices an example could be: weights=list(c(0.5,0.2,0.3),c(0.1,0.5,0.4)). To provide
a fixed weight for some data sets and let it vary randomly for others, the element "x" 
indicates a free parameter. An example is weights=list(c(0.7,"x","x")). The weight 0.7 
is now fixed for the first data matrix while the remaining 0.3 weight will be divided over
the other two data sets. This implies that every combination of the sequence from 0 to 0.3
with steps of 0.1 will be reported and clustering will be performed for each.
}
\value{
The returned value is a list with the following four elements.
\item{Incidence}{The summed incidence matrices for each data matrix}
\item{IncidenceComb}{The co-association matrix after a weighted sum of the
elements of Incidence for each weight}
\item{Results}{The hierarchical clustering result for each element in 
IncidenceComb}
\item{Clust}{The result for the weight specified in Clustweight}
The value has class 'CEC'
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

MCF7_CECa=CECa(List=L,distmeasure=c("tanimoto","tanimoto"),
normalize=FALSE,method=NULL,t=25,r=NULL,nrclusters=c(7,7),
clust="agnes",linkage=c("flexible","flexible"),alpha=0.625,StopRange=FALSE)
}
\note{
For now, only hierarchical clustering with the \code{agnes} 
function is implemented.
}
\seealso{
\code{\link{CEC}},\code{\link{CECb}},\code{\link{CECc}}
}
\keyword{ Complementary Ensemble Clustering }
\keyword{ Integrative Clustering }
