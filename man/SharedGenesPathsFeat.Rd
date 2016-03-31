\name{SharedGenesPathsFeat}
\alias{SharedGenesPathsFeat}
\title{
Intersection of genes and pathways over multiple methods
}
\description{
It is interesting to investigate exactly which and how many differently expressed genes, pathways and characteristics
are shared by the clusters over the different methods. The function \code{SharedGenesPathsFeat} will provide this
information. Given the outputs of the \code{DiffGenes}, the \code{Geneset.intersect} function and/or \code{CharacteristicFeatures},
it investigates how many genes, pathways and/or characteristics are expressed by each cluster per method, how many
of these are shared over the methods and which ones are shared including their respective p-values
of each method and a mean p-value. This is very handy to look into the shared genes and pathways
of clusters that share many objects but also of those that only share only a few. Further, the result
also includes the number of compounds per cluster per method and how many of these are shared
over the methods. The input can also be focused for a specific selection of compounds or a specific cluster.
}
\usage{
SharedGenesPathsFeat(DataLimma = NULL, DataMLP = NULL, DataFeat=NULL,
names = NULL, Selection=FALSE)
}
\arguments{
  \item{DataLimma}{
Optional. The output of a \code{DiffGenes} function.
}
  \item{DataMLP}{
Optional. The output of \code{Geneset.intersect} function.
}
  \item{DataFeat}{
Optional. The output of \code{CharacteristicFeatures} function.
}
  \item{names}{
Optional. Names of the methods or "Selection" if it only considers a selection of compounds.
}
  \item{Selection}{
Logical. Do the results concern only a selection of compounds or a specific cluster? 
If yes, then Selection should be put to TRUE. Otherwise all compounds and clusters 
are considered.
}
}
\value{
The result of the \code{SharedGenesPathsFeat} function is a list with two elements. The first element Table is a table
indicating how many genes, pathways and/or characteristics were found to be differentially expressed and how
many of these are shared. The table also contains the number of compounds shared between the
clusters of the different methods. The second element Which is another list with a component per
cluster. Each component consists of four vectors: SharedComps indicating which objects were shared
across the methods, SharedGenes represents the shared genes, SharedPaths shows the
shared pathways and SharedFeat the shared features.
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
data(ListGO)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)

L=list(MCF7_F,MCF7_T)
names=c('FP','TP')

MCF7_Paths_FandT=PathwaysIter(L,GeneExpr=geneMat,nrclusters=7,method=c("limma", "MLP"),
ENTREZID=GeneInfo[,1],geneSetSource = "GOBP",top=NULL,topG=NULL,GENESET=ListGO,sign=0.05,
niter=2,fusionsLog=TRUE,WeightClust=TRUE,names=c("FP","TP"))

MCF7_Paths_intersection=Geneset.intersect(MCF7_Paths_FandT,0.05,names=names,
seperatetables=FALSE,separatepvals=FALSE)

MCF7_DiffGenes_FandT10=DiffGenes(list(MCF7_F,MCF7_T),geneMat,nrclusters=7,"limma",0.05,top=10)

MCF7_Char=CharacteristicFeatures(list(MCF7_F,MCF7_T),Selection=NULL,BinData=list(fingerprintMat,
targetMat),Datanames=c("F","T"),nrclusters=7,top=NULL,sign=0.05,fusionsLog=TRUE,WeightClust=TRUE,
names=c("F","T"))

MCF7_Shared=SharedGenesPathsFeat(DataLimma=MCF7_DiffGenes_FandT10,DataMLP=
MCF7_Paths_intersection,DataFeat=MCF7_Char)
str(MCF7_Shared)
}
}


