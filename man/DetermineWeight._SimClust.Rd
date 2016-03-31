\name{DetermineWeight_SimClust}
\alias{DetermineWeight_SimClust}
\title{
Determines an optimal weight for weighted clustering by similarity weighted clustering.
}
\description{
The function \code{DetermineWeight_SimClust} determines an optimal weight for performing weighted
similarity clustering on by applying similarity clustering. For each given weight, is each
separate clustering compared to the clustering on a weighted dissimilarity matrix and a
Jaccard coefficient is calculated. The ratio of the Jaccard coefficients closets to one
indicates an optimal weight.
}
\usage{
DetermineWeight_SimClust(List, type = c("data","dist","clusters"),
weight=seq(0, 1, by = 0.01),nrclusters = NULL, distmeasure = c("tanimoto",
"tanimoto"),normalize=FALSE,method=NULL,clust = "agnes",linkage=c("ward","ward")
,alpha=0.625,gap = FALSE, maxK = 50, names = c("B", "FP"),StopRange=FALSE,
plottype="new",location=NULL)

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
  \item{weight}{
 Optional. A list of different weight combinations for the data sets in List.
 If NULL, the weights are determined to be rqual for each data set.
 It is further possible to fix weights for some data matrices and to
 let it vary randomly for the remaining data sets. An example is provided in the details.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in. This is necessary for the
computation of the Jaccard coefficient.
}
  \item{distmeasure}{
A character vector with the distance measure for each data matrix.
Should be one of "tanimoto", "euclidean", "jaccard","hamming".
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
  \item{gap}{
Logical. Whether or not to calculate the gap statistic in the clustering
on each data matrix separately. Only if type="data".
}
  \item{maxK}{
The maximal number of clusters to consider in calculating the gap
statistic. Only if type="data".
}
  \item{names}{
The labels to give to the elements in List.
}
  \item{StopRange}{
Logical. Indicates whether the distance matrices with values not between zero and one should be standardized to have so.
If FALSE the range normalization is performed. See \code{Normalization}. If TRUE, the distance matrices are not changed.
This is recommended if different types of data are used such that these are comparable.
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
If the type of List is data, an hierarchical clustering is performed on each
data matrix separately. After obtaining clustering results for the two data
matrices, the distance matrices are extracted. If these are not calculated
with the same distance measure, they are normalized to be in the same range.
For each weight, a weighted linear combination of the distance matrices is
taken and hierarchical clustering is performed once again. The resulting
clustering is compared to each of the separate clustering results and a 
Jaccard coefficient is computed. The ratio of the Jaccard coefficients
closets to one, indicates an optimal weight. A plot of all the ratios is
produced with an extra indication for the optimal weight. 

The weight combinations should be provided as elements in a list. For three data
matrices an example could be: weights=list(c(0.5,0.2,0.3),c(0.1,0.5,0.4)). To provide
a fixed weight for some data sets and let it vary randomly for others, the element "x" 
indicates a free parameter. An example is weights=list(c(0.7,"x","x")). The weight 0.7 
is now fixed for the first data matrix while the remaining 0.3 weight will be divided over
the other two data sets. This implies that every combination of the sequence from 0 to 0.3
with steps of 0.1 will be reported and clustering will be performed for each.
}
\value{
The returned value is a list with three elements:
\item{ClustSep}{The result of \code{Cluster} for each single element of List}
\item{Result}{A data frame with the Jaccard coefficients and their ratios 
for each weight}
\item{Weight}{The optimal weight}
}
\references{
PERUALILA-TAN, N., SHKEDY, Z., TALLOEN, W., GOEHLMANN, H. W. H.,
QSTAR Consortium, VAN MOERBEKE, M., KASIM, A., (in press). Weighted-Similarity Based
Clustering of Chemical Structure and Bioactivity Data in Early Drug Discovery. Journal of
Bioinformatics and Computational Biology. 
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",alpha=0.625,gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",alpha=0.625,gap=FALSE,maxK=55,StopRange=FALSE)

L=list(MCF7_F,MCF7_T)

MCF7_Weight=DetermineWeight_SimClust(List=L,type="clusters",weight=seq(0,1,by=0.01),nrclusters=7,
distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,method=NULL,clust="agnes",
linkage=c("flexible","flexible"),alpha=0.625,gap=FALSE,maxK=50,names=c("FP","TP"),StopRange=FALSE,
plottype="new",location=NULL)

}
\seealso{
\code{\link{WeightedSimClust}}
}
\keyword{ Weighted Similarity Clustering }
\keyword{ Integrated Data Clustering }
