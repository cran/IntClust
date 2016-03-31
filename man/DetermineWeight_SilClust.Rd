\name{DetermineWeight_SilClust}
\alias{DetermineWeight_SilClust}
\title{
Determines an optimal weight for weighted clustering by silhouettes widths.
}
\description{
The function \code{DetermineWeight_SilClust} determines an optimal weight for weighted
similarity clustering by calculating silhouettes widths. See "Details" for a more elaborate
description.
}
\usage{
DetermineWeight_SilClust(List,type=c("data","dist","clusters"),weight=seq(0,1,by=0.01),
distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,method=NULL,nrclusters=NULL,
names=NULL,nboot=1000,StopRange=FALSE,plottype="new",location=NULL)
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
 If NULL, the weight is sequence from 0 to 1  with steps of 0.1 and a result is produced 
 for each weight. It is further possible to fix weights for some data matrices and to
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
  \item{names}{
The labels to give to the elements in List.
}
  \item{nboot}{
Number of bootstraps to be run.
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
For each given weight, a linear combination of the distance matrices of the single data sources is 
obtained. For these distance matrices, medoid clustering with nrclusters is set up by the \code{pam} function of the \pkg{cluster}
and the silhouette widths are retrieved. These widths indicates how well an object fits in its current cluster. Values around one
indicate an appropriate cluster. The silhouette widths are regressed in function of the cluster membership determined by the objects. First, in function
of the cluster membership determined by the weighted combination. Then, also in function of the cluster membership determined by the single 
source clustering. The regression function is fit by the \code{lm} function and the \code{r.squared} value is retrieved. The\code{r.squared} value
indicates how much of the variance of the silhouette widths is explained by the membership. Optimally this value is high.

Next, a statistic is determined. Suppose that RWW is the \code{r.squared} retrieved from regressing the weighted silhouette widths versus the weighted cluster
membership and RWX the \code{r.squared} retrieved from regressing the weighted silhouette widths versus the cluster membership determined by data X. 
If M is total number of data sources, than statistic is obtained as:
\deqn{Stat=abs(M*RWW-\sum{RWX})}

The lower the statistical value, the better the weighted clustering is explained by the single data sources. The goal is to obtain the weights for which this
value is minimized. Via bootstrapping a p-value is obtained for every statistic.

The weight combinations should be provided as elements in a list. For three data
matrices an example could be: weights=list(c(0.5,0.2,0.3),c(0.1,0.5,0.4)). To provide
a fixed weight for some data sets and let it vary randomly for others, the element "x" 
indicates a free parameter. An example is weights=list(c(0.7,"x","x")). The weight 0.7 
is now fixed for the first data matrix while the remaining 0.3 weight will be divided over
the other two data sets. This implies that every combination of the sequence from 0 to 0.3
with steps of 0.1 will be reported and clustering will be performed for each.
}
\value{
Two plots are made: one of the statistical values versus the weights and one of the p-values versus the weights.
Further, a list with two elements is returned:
\item{Result}{A data frame with the statistic for each weight combination}
\item{Weight}{The optimal weight}
}
\author{
Marijke Van Moerbeke
}
\examples{
\dontrun{
data(fingerprintMat)
data(targetMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)

L=list(MCF7_F,MCF7_T)

MC7_Weight=DetermineWeight_SilClust(List=L,type="clusters",weight=seq(0,1,by=0.01),nrclusters=7,
distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,method=NULL,names=c("FP","TP"),nboot=100,
StopRange=FALSE,plottype="new",location=NULL)
}
}
\keyword{ Weighted Clustering }
\keyword{ Integrated Data Clustering }
