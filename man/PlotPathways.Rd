\name{PlotPathways}
\alias{PlotPathways}
\title{
A GO plot of a pathway analysis output.
}
\description{
The \code{PlotPathways} function takes an output of the \code{PathwayAnalysis} fucntion and plots a GO graph with the
help of the \code{plotGOgraph} function of the MLP package.
}
\usage{
PlotPathways(Pathways,nRow=5,main=NULL,plottype="new",location=NULL)
}
\arguments{
  \item{Pathways}{
One element of the output list returned by \code{PathwayAnalysis} or \code{Geneset.intersect}.
}
  \item{nRow}{
Number of GO IDs for which to produce the plot
}
  \item{main}{
Title of the plot.
}
  \item{plottype}{
Should be one of "pdf","new" or "sweave". If "pdf", a location should be provided in
"location" and the figure is saved there. If "new" a new graphic device is opened and if "sweave",
the figure is made compatible to appear in a sweave or knitr document.
}
\item{location}{
If plottype is "pdf", a location should be provided in "location" and the figure is saved there. 
}
}
\value{
The output is a GO graph.
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

MCF7_PathsFandT=PathwayAnalysis(L, GeneExpr = geneMat, nrclusters = 7, method = c("limma", 
"MLP"), ENTREZID = GeneInfo[, 1], geneSetSource = "GOBP", topP = NULL, 
topG = NULL, GENESET = GS, sign = 0.05,niter=2,fusionsLog = TRUE, WeightClust = TRUE, 
 names =names,seperatetables=FALSE,separatepvals=FALSE)
 
PlotPathways(MCF7_PathsFandT$FP$Pathways,nRow=5,main=NULL)
}
}
\keyword{ Pathway Analysis}

