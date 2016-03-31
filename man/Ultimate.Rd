\name{Ultimate}
\alias{Ultimate}
\title{
Function that performs any aggregated data function
}
\description{
The function \code{Ultimate} has the ability to perform multiple of the methods
listed above simultaneously. The only necessary input are the data matrices
and specification of the options. First, clustering is based on each data matrix
separately after which the specified integrative analysis methods are conducted.
A plot comparing the results is made automatically with ComparePlot. If weights
are involved in the method, a comparison plot of the results for these weights
is made as well.
}
\usage{
Ultimate(List,type=c("data","dist","clusters"),distmeasure,normalize=FALSE,method=NULL,
StopRange=FALSE,NN = 20,mu = 0.5, T = 20, t = 10, r = NULL, nrclusters = NULL,
nrclusterssep = c(7, 7),nrclustersseq = NULL, weight = NULL, Clustweight = 0.5, 
clust = "agnes", linkage=c("ward","ward"),alpha=0.625, gap = FALSE, maxK = 50,
 IntClust = c("ADC", "ADECa", "ADECb","ADECc", "WonM", "CECa", "CECb", "CECc",
  "WeightedClust", "WeightedSim", "SNFa", "SNFb", "SNFc"), fusionsLog = TRUE,
 WeightClust= TRUE, PlotCompare = FALSE, cols = NULL, ...)  
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
A list of matrices of the same type. It is assumed that the rows are corresponding
to the objects.
}
  \item{type}{
Type indicates whether the provided matrices in "List" are either data or distance
matrices obtained from the data. If type="dist" the calculation of the distance
matrices is skipped and the methods ADC ,ADECa, ADECb, ADECc, CECa, CECb and
CECc are not performed. Type should be one of "data" or "dist".		
}
  \item{distmeasure}{
A vector of the distance measures to be used on each data matrix.
}
  \item{normalize}{
Logical. Indicates whether to normalize the distance matrices or not.
This is recommended if different distance types are used. More details
on standardization in \code{Distance}.
}
  \item{method}{
A method of normalization. Should be one of "Quantile","Fisher-Yates",
"standardize","Range" or any of the first letters of these names.
}
  \item{StopRange}{
Logical. Indicates whether the distance matrices with values not between zero and one should be standardized to have so.
If FALSE the range normalization is performed. See \code{Normalization}. If TRUE, the distance matrices are not changed.
This is recommended if different types of data are used such that these are comparable.
}
  \item{NN}{
The number of neighbours to be used in SNF.
}
  \item{mu}{
A parameter in SNF.
}
  \item{T}{
The number of iterations in SNF.
}
  \item{t}{
The number of iterations in ADEC and CEC.
}
  \item{r}{
Optional. The number of features to take for the random sample in ADEC and
CEC.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in for ADEC and the plot.
}
  \item{nrclusterssep}{
Optional. Vector of the number of clusters to cut the dendrogram in of each data
source. If NULL, the value of nclusters is used for each.
}
  \item{nrclustersseq}{
The sequence of number of clusters to cut the dendrogram in for ADECb, CECb
and WonM.
}
  \item{weight}{
The weights to be used in CEC and WeightedClust.
}
  \item{Clustweight}{
Optional. To be used for the outputs of CEC or WeightedClust. Then only the
result of the Clust element is considered.
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
Logical. Indicator if gap statistics should be computed. 
Setting to $FALSE$ will greatly reduce the computation time.
}
  \item{maxK}{
The maximum number of clusters to be considered during the gap.
}
  \item{IntClust}{
Specification of the methods to be applied.
}
  \item{fusionsLog}{
To be handed to MatrixFunction.
}
  \item{WeightClust}{
To be handed to MatrixFunction.
}
  \item{PlotCompare}{
Logical. Should the plot over the methods and weight be produced?
}
  \item{cols}{
Color scheme to be used in the plots.
}
  \item{\dots}{
Options to be given to ComparePlot.
}
}
\value{
The output of \code{Ultimate} is a list . The first element contains
the results of the clustering of the first data source and the
last element on the second data source . In between are the results of the
integrative methods.
}
\author{
Marijke Van Moerbeke
}
\examples{
\dontrun{
data(fingerprintMat)
data(targetMat)
data(Colors2)
L=list(fingerprintMat,targetMat)

MCF7_All=Ultimate(L,type="data",distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,
method=NULL,StopRange=FALSE,NN=20,alpha=0.5,T=20,t=25,r=NULL,nrclusters=7,
nrclustersseq=c(5,25,1),weight=seq(1,0,-0.1),Clustweight=0.5,clust="agnes",
linkage=c("ward","ward"),alpha=0.625,gap=FALSE,IntClust=c("ADC","ADECa","ADECb",
"ADECc","WonM","CECa","CECb","CECc","WeightedClust","WeightedSim",
"SNFa","SNFb","SNFc"),fusionsLog=TRUE,WeightClust=TRUE,PlotCompare=TRUE,
cols=Colors2)
}
}