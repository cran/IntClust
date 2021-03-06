% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{CVAA}
\alias{CVAA}
\title{Cumulative voting-based aggregation algorithm}
\usage{
CVAA(Reference = NULL, nrclustersR = 7, List, typeL = c("data", "dist",
  "clust"), distmeasure = c("tanimoto", "tanimoto"), normalize = c(FALSE,
  FALSE), method = c(NULL, NULL), clust = "agnes", linkage = c("flexible",
  "flexible"), alpha = 0.625, nrclusters = c(7, 7), gap = FALSE,
  maxK = 15, votingMethod = c("CVAA", "W-CVAA"), optimalk = nrclustersR)
}
\arguments{
\item{Reference}{The reference structure to be updated.}

\item{nrclustersR}{The number of clusters present in the reference structure. Default is 7.}

\item{List}{A list of data matrices. It is assumed the rows are corresponding with the objects.}

\item{typeL}{indicates whether the provided matrices in "List" are either data matrices, distance
matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
matrices is skipped and if type="clusters" the single source clustering is skipped.
Type should be one of "data", "dist" or "clusters".}

\item{distmeasure}{A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").}

\item{normalize}{Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}.}

\item{method}{A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.}

\item{clust}{Choice of clustering function (character). Defaults to "agnes".}

\item{linkage}{Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.}

\item{alpha}{The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible"}

\item{nrclusters}{The number of clusters to divide each individual dendrogram in. Default is c(7,7) for two data sets.}

\item{gap}{Logical. Whether the optimal number of clusters should be determined with the gap statistic. Defaults to FALSE.}

\item{maxK}{The maximal number of clusters to investigate in the gap statistic. Default is 15.}

\item{votingMethod}{The method to be performed: "CVAA","W-CVAA".}

\item{optimalk}{An estimate of the final optimal number of clusters. Default is nrclustersR.}
}
\value{
The returned value is a list of two elements:
\item{DistM}{A NULL object}
\item{Clust}{The resulting clustering}
The value has class 'Ensemble'.
}
\description{
The \code{CVAA} includes the ensemble clustering methods CVAA and W-CVAA which are voting-based consensus methods.
}
\details{
\insertCite{Saeed2012}{IntClust} describe the Cumulative Voting-based Aggregation Algorithm (CVAA) and introduce a 
variation Weighted Cumulative Voting-based Aggregation Algorithm (W-CVAA, \insertCite{Saeed2014}{IntClust}). In the CVAA algorithm, one data 
partitioning is chosen as the reference partition. In a first step each partition is relabelled with respect to the reference partition in 
search of an optimal relabelling. In the next step a consensus partition is obtained. The W-CVAA algorithm is similar but appoints weights 
to each individual partition. The weights are based on the mutual information of the partition measured by the Shannon entropy.
}
\examples{
data(fingerprintMat)
data(targetMat)
L=list(fingerprintMat,targetMat)

MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)

MCF7_CVAA=CVAA(Reference=MCF7_T,nrclustersR=7,List=L,typeL="data",
distmeasure=c("tanimoto", "tanimoto"),normalize=c(FALSE,FALSE),method=
c(NULL,NULL),clust="agnes",linkage = c("flexible","flexible"),alpha=0.625,
nrclusters=c(7,7),gap = FALSE, maxK = 15,votingMethod="CVAA",optimalk=7)
}
\references{
\insertRef{Saeed2012}{IntClust} \insertRef{Saeed2014}{IntClust}
}
