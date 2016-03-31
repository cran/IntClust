\name{ChooseCluster}
\alias{ChooseCluster}
\title{
Interactive plot to determine DE Genes and DE features for a specific cluster
}
\description{
If desired, the function produced a dendrogram of a clustering results. One or
multiple cluster can be indicated by a mouse click. From these clusters DE genes
and characteristic features are determined. It is also possible to provide the
compounds of interest without producing the plot.
}
\usage{
ChooseCluster(Interactive = TRUE, LeadCpds=NULL, ClusterResult, ColorLab 
= NULL,BinData=NULL,ContData=NULL,Datanames = c("FP"), GeneExpr, topChar = 20, topG = 20, 
sign = 0.05, nrclusters = NULL, cols = NULL, N = 1)

}
\arguments{
  \item{Interactive}{
Logical. Produce plot or not. Defaults to TRUE.
}
  \item{LeadCpds}{
A list of the compounds of the clusters of interest. If Interactive=TRUE, 
these are determined by the mouse-click and it defaults to NULL.
}
  \item{ClusterResult}{
The output of one of the aggregated cluster functions,
The clustering result of interest.
}
  \item{ColorLab}{
The clustering result the dendrogram should be colored after as in 
\code{ClusterPlot}. It is the output of one of the clustering functions.
}
 \item{BinData}{
A list of the binary feature data matrices. These will be evaluated with the fisher's extact test.
}
  \item{ContData}{
A list of continuous data sets of the compounds. These will be evaluated with the t-test.
}
  \item{Datanames}{
A character vector of the labels to give to the matrices in BinData.
}
  \item{GeneExpr}{
A gene expression matrix, may also be an ExpressionSet. The rows should correspond with the genes.
}
  \item{topChar}{
The number of top characteristics to return. If NULL, only the significant 
characteristics are saved.
}
  \item{topG}{
The number of top genes to return. If NULL, only the significant 
genes are saved.
}
  \item{sign}{
The significance level.
}
  \item{nrclusters}{
Optional. The number of clusters to cut the dendrogram in. If NULL, the dendrogram
will be plotted without colors to discern the different clusters.
}
  \item{cols}{
The colors to use in the dendrogram.
}
  \item{N}{
The number of clusters one wants to identify by a mouse click.
}
}
\details{
The DE genes are determined by testing for significance of the specified cluster
versus all other compounds combined. This is performed by the limma function. The
binary features are evaluated with the fisher exact test while the continuous
features are tested with the t-test. Multiplicity correction is included.
}
\value{
The returned value is a list with one element per cluster of interest indicated by the prefix "Choice". 
This element is again a list with the following three elements:
\item{Compounds}{A list with the elements LeadCpds (the compounds of interest) and 
OrderedCpds (all compounds in the order of the clustering result)}
\item{Characteristics}{The found (top) characteristics of the feature data}
\item{Genes}{A list with the elements TopDE (a table with information on the top
genes) and AllDE (a table with information on all genes)}
}
\author{
Marijke Van Moerbeke
}
\examples{
\dontrun{
data(fingerprintMat)
data(targetMat)
data(geneMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)

MCF7_Interactive=ChooseCluster(Interactive=TRUE,LeadCpds=NULL,ClusterResult=MCF7_T,
ColorLab=MCF7_F,BinData=list(fingerprintMat),Datanames=c("FP"),geneMat,
topChar = 20, topG = 20,nrclusters=7,N=1)
}

}
\keyword{ Differential expression }
\keyword{ Interactive plot }
