\name{PathwayAnalysis}
\alias{PathwayAnalysis}
\title{
Pathway Analysis
}
\description{
The \code{PathwayAnalysis} function combines the functions \code{PathwaysIter} and \code{Geneset.intersect} such that only one function should be called.
}
\usage{
PathwayAnalysis(List, Selection=NULL, GeneExpr = NULL, nrclusters = NULL,
method = c("limma", "MLP"), GeneInfo = NULL, geneSetSource = "GOBP",
topP = NULL, topG = NULL, GENESET = NULL, sign = 0.05, niter = 10, 
fusionsLog = TRUE, WeightClust = TRUE, names = NULL,seperatetables=FALSE,
separatepvals=FALSE)

}
\arguments{
  \item{List}{
A list of clustering outputs or output of the\code{DiffGenes} function. The first element of the
list will be used as the reference in \code{ReorderToReference}. The output of \code{ChooseFeatures} is also accepted.
}
 \item{Selection}{
If pathway analysis should be conducted for a specific selection of compounds, this selection can be provided here.
Selection can be of the type "character" (names of the compounds) or "numeric" (the number of specific cluster).
 }
  \item{GeneExpr}{
The gene expression matrix of the objects. The rows should correspond with the genes.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in.
}
  \item{method}{
The method to applied to look for DE genes. For now, only the limma method
is available.
}
  \item{GeneInfo}{
A data frame with at least the columns ENTREZID and SYMBOL. This is necessary to connect the
symbolic names of the genes with their EntrezID in the correct order. The order of the gene is here not
in the order of the rownames of the gene expression matrix but in the order of their significance.
}
  \item{geneSetSource}{
The source for the getGeneSets function, defaults to "GOBP".
}
  \item{topP}{
Overrules sign. The number of pathways to display for each cluster. 
If not specified, only the significant genes are shown.
}
  \item{topG}{
 Overrules sign. The number of top genes to be returned in the result.
 If not specified, only the significant genes are shown.
}
  \item{GENESET}{
Optional. Can provide own candidate gene sets.
}
  \item{sign}{
The significance level to be handled.
}
  \item{niter}{
The number of times to perform pathway analysis.
}
  \item{fusionsLog}{
To be handed to ReorderToReference.
}
  \item{WeightClust}{
To be handed to ReorderToReference.
}
  \item{names}{
Optional. Names of the methods.
}
  \item{seperatetables}{
Logical. If TRUE, a separate element is created per cluster. containing
the pathways for each iteration.
}
  \item{separatepvals}{
Logical. If TRUE, the p-values of the each iteration of each pathway
in the intersection is given. If FALSE, only the mean p-value is provided.
}
}
\value{
The output is a list with an element per method. For each method, it
is portrayed per cluster which pathways belong to the intersection over
all iterations and their corresponding mean p-values.
}
\author{
Marijke Van Moerbeke
}
\examples{
\dontrun{
data(fingerprintMat)
data(targetMat)
data(geneMat)
data(GeneInfo)
data(GS)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)

L=list(MCF7_F,MCF7_T)
names=c('FP','TP')

MCF7_PathsFandT=PathwaysAnalysis(L, GeneExpr = geneMat, nrclusters = 7, method = c("limma", 
"MLP"), GeneInfo = GeneInfo, geneSetSource = "GOBP", topP = NULL, 
topG = NULL, GENESET = GS, sign = 0.05,niter=2,fusionsLog = TRUE, WeightClust = TRUE, 
 names =names,seperatetables=FALSE,separatepvals=FALSE)
}
}
\keyword{ Pathway Analysis}

