\name{Pathways}
\alias{Pathways}
\title{
Pathway analysis for multiple clustering results
}
\description{
A pathway analysis per the cluster per method is conducted.
}
\usage{
Pathways(List, Selection=NULL, GeneExpr = NULL, nrclusters = NULL, method = 
c("limma", "MLP"),GeneInfo = NULL, geneSetSource = "GOBP", topP = NULL,
topG = NULL,GENESET = NULL, sign = 0.05, fusionsLog = TRUE, WeightClust = TRUE, 
names = NULL)
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
The gene expression matrix or ExpressionSet of the objects. The rows should correspond with the genes.
}
  \item{nrclusters}{
Optional. The number of clusters to cut the dendrogram in. The number of clusters should not be
specified if the interest lies only in a specific selection of compounds which is known by name. 
Otherwise, it is required.
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
  \item{fusionsLog}{
To be handed to ReorderToReference.
}
  \item{WeightClust}{
To be handed to ReorderToReference.
}
  \item{names}{
Optional. Names of the methods.
}
}
\details{
After finding differently expressed genes, it can be investigated whether pathways are related to
those genes. This can be done with the help of the function \code{Pathways} which makes use of the 
\code{MLP} function of the MLP package. Given the output of a method, the cutree function is performed
which results into a specific number of clusters. For each cluster, the limma method is performed
comparing this cluster to the other clusters. This to obtain the necessary p-values of the genes.
These are used as the input for the \code{MLP} function to find interesting pathways. By default 
the candidate gene sets are determined by the \code{AnnotateEntrezIDtoGO} function. The default
source will be GOBP, but this can be altered.  Further, it is also possible to provide own candidate
gene sets in the form of a list of pathway categories in which each component contains a vector of
Entrez Gene identifiers related to that particular pathway. The default values for the minimum and
maximum number of genes in a gene set for it to be considered were used. For MLP this is 
respectively 5 and 100. If a list of outputs of several methods is provided as data input, the
cluster numbers are rearranged according to a reference method. The first method is taken as the 
reference and ReorderToReference is applied to get the correct ordering. When the clusters haven been
re-appointed, the pathway analysis as described above is performed for each cluster of each method.
}
\value{
The returned value is a list with an element per cluster per method. This element is again a
list with the following four elements:
\item{Compounds}{A list with the elements LeadCpds (the compounds of interest) and 
OrderedCpds (all compounds in the order of the clustering result)}
\item{Characteristics}{The found (top) characteristics of the feauture data}
\item{Genes}{A list with the elements TopDE (a table with information on the top
genes) and AllDE (a table with information on all genes)}
\item{Pathways}{A list with the element ranked.genesets.table which is a data frame containing
the genesets, their p-values and their descriptions. The second element is nr.genesets and
contains the used and total number of genesets.}
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

MCF7_PathsFandT=Pathways(L, GeneExpr = geneMat, nrclusters = 7, method = c("limma", 
"MLP"), GeneInfo = GeneInfo, geneSetSource = "GOBP", topP = NULL, 
topG = NULL, GENESET = GS, sign = 0.05,fusionsLog = TRUE, WeightClust = TRUE, 
 names =names)
 }
}
\seealso{
\code{\link{PathwaysIter}}
}
\keyword{ Pathway Analysis }

