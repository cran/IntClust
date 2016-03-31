\name{CompareSilCluster}
\alias{CompareSilCluster}
\title{
Compares medoid clustering results based on silhouette widths
}
\description{
The function \code{CompareSilCluster} compares the results of two medoid clusterings.
The null hypothesis is that the clustering is identical. A test statistic is calcluated and a p-value 
obtained with bootstrapping. See "Details" for a more elaborate description.
}
\usage{
CompareSilCluster(List,type=c("data","dist"),distmeasure=c("tanimoto",
"tanimoto"),normalize=FALSE,method=NULL,nrclusters=NULL,names=NULL,
nboot=1000,StopRange=FALSE,plottype="new",location=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
A list of matrices of the same type. It is assumed the rows are corresponding
with the objects.
}
  \item{type}{
Type indicates whether the provided matrices in "List" are either data matrices, distance
matrices.	
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
  \item{nrclusters}{
The number of clusters to cut the dendrogram in. This is necessary for the
computation of the Jaccard coefficient.
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
For the data or distance matrices in List, medoid clustering with nrclusters is set up by the \code{pam} function of the \pkg{cluster}
and the silhouette widths are retrieved. These widths indicates how well an object fits in its current cluster. Values around one
indicate an appropriate cluster while values around zero indicate that the object might as well lie in its neighbouring cluster.
The silhouette widths are than regressed in function of the cluster membership of the objects. First the widths are modeled according
to the cluster membership of object these were derived from. Next, these are modeled in function of the membership determined by the 
other object. The regression function is fit by the \code{lm} function and the \code{r.squared} value is retrieved. The\code{r.squared} value
indicates how much of the variance of the silhouette widths is explained by the membership. Optimally this value is high.

Next, a statistic is determined. Suppose that RXX is the \code{r.squared} retrieved from regressing the silhouette widths of object X versus the corresponding cluster
membership of object X and RXY the \code{r.squared} retrieved from regressing the silhouette widths of object X versus the cluster membership determined by object Y and vice versa. 
The statistic is obtained as:
\deqn{Stat=abs(\sum{RXX}-\sum{RXY})}


The lower the statistical value, the better the clustering is explained by the sources. Via bootstrapping a p-value is obtained.
}
\value{
A plots are made of the density of the statistic under the null hypotheses. The p-value is also indicated on this plot.
Further, a list with two elements is returned:
\item{Observed Statistic}{The observed statistical value}
\item{P-Value}{The P-value of the obtained statistic retrieved after bootstrapping}
}
\author{
Marijke Van Moerbeke
}
\examples{
\dontrun{
data(fingerprintMat)
data(targetMat)

List=list(fingerprintMat,targetMat)

Comparison=CompareSilCluster(List=List,type="data",
distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,method=NULL,
nrclusters=7,names=NULL,nboot=1000,StopRange=FALSE,plottype="new",location=NULL)

Comparison
}
}

