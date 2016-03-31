\name{ProfilePlot}
\alias{ProfilePlot}
\title{
Plotting gene profiles
}
\description{
In \code{ProfilePlot}, the gene profiles of the significant genes for a specific
cluster are shown on 1 plot. Therefore, each gene is normalized by subtracting its
the mean.
}
\usage{
ProfilePlot(Genes, Comps, GeneExpr = NULL, 
Raw = FALSE, OrderLab = NULL, ColorLab = NULL, nrclusters = NULL,
cols = NULL, AddLegend = TRUE,  margins = c(8.1, 4.1, 1.1, 6.5),
extra = 5,plottype="new",location=NULL, ...)
}
\arguments{
  \item{Genes}{
The genes to be plotted.
}
  \item{Comps}{
The objects to be plotted or to be separated from the other objects.
}
  \item{GeneExpr}{
The gene expression matrix or ExpressionSet of the objects.
}
  \item{Raw}{
Logical. Should raw p-values be plotted?
}
  \item{OrderLab}{
Optional. If the compounds are to set in a specific order of a specific method.
}
  \item{ColorLab}{
The clustering result that determines the color of the labels of the
objects in the plot.
}
  \item{nrclusters}{
Optional. The number of clusters to cut the dendrogram in.
}
  \item{cols}{
Optional. The color to use for the objects in Clusters for each method.
}
  \item{AddLegend}{
Optional. Whether a legend of the colors should be added to the plot.
}
  \item{margins}{
Optional. Margins to be used for the plot.
}
  \item{extra}{
The space between the plot and the legend.
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
  \item{\dots}{
Optional parameter to be handed to the plot function.
}
}
\value{
A plot which contains multiple gene profiles. A distinction is made between the values for the
objects in Comps and the others.
}
\author{
Marijke Van Moerbeke
}
\seealso{
 \code{\link{ProfilePlot}}
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
GeneInfo=GeneInfo,geneSetSource = "GOBP",top=NULL,topG=NULL,GENESET=ListGO,sign=0.05,
niter=2,fusionsLog=TRUE,WeightClust=TRUE,names=c("FP","TP"))

MCF7_Paths_intersection=Geneset.intersect(MCF7_Paths_FandT,0.05,names=names,
seperatetables=FALSE,separatepvals=FALSE)

MCF7_DiffGenes_FandT10=DiffGenes(list(MCF7_F,MCF7_T),geneMat,nrclusters=7,"limma",0.05,top=10)

MCF7_Shared10=Shared(DataLimma=MCF7_DiffGenes_FandT10,DataMLP=MCF7_Paths_intersection)

Comps=SharedComps(list(MCF7_DiffGenes_FandT10$`Method 1`$"Cluster 1",
MCF7_DiffGenes_FandT10$`Method 2`$"Cluster 1"))

MCF7_SharedGenes=FindGenes(DataLimma=MCF7_DiffGenes_FandT10,names=c("FP","TP"))

Genes=names(MCF7_SharedGenes[[1]])[-c(2,4,5)]

ListC=list(MCF7_DiffGenes_FandT10[[1]][[1]]$Compounds$LeadCpds,
MCF7_DiffGenes_FandT10[[2]][[1]]$Compounds$LeadCpds)

colscl=ColorPalette(colors=c("red","green","purple","brown","blue","orange"),ncols=9)

ProfilePlot(Genes=Genes,Comps=Comps,GeneExpr=geneMat,Raw=FALSE,OrderLab=MCF7_F,
ColorLab=NULL,nrcluster=7,Clusters=ListC,cols=colscl,AddLegend=TRUE,
usedgenes=Genes,margins=c(8.1,4.1,1.1,6.5),plottype="new",location=NULL,cex=0.75)
}
}
\keyword{ Gene Profile }

