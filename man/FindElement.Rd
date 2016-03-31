\name{FindElement}
\alias{FindElement}
\title{
Find an element in a data structure
}
\description{
The function \code{FindElement} is used internally in the
\code{PreparePathway} function but might come in handy for other uses
as well. Given the name of an object, the function searches for that
object in the data structure and extracts it. When multiple objects 
have the same name, all are extracted.
}
\usage{
FindElement(What, Object, Element = list())
}
\arguments{
  \item{What}{
A character string indicating which object to look for.
}
  \item{Object}{
The data structure to look into. Only the classes data frame
and list are supported.
}
  \item{Element}{
Not to be specified by the user.
}
}
\value{
The returned value is a list with an element for each object
found. The element contains everything the object contained
in the original data structure.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
data(geneMat)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)

MCF7_DiffGenes_FandT10=DiffGenes(list(MCF7_F,MCF7_T),Selection=NULL,GeneExpr=geneMat,
nrclusters=7,method="limma",sign=0.05,top=10,fusionsLog = TRUE, WeightClust = TRUE, 
names = NULL)

Find=FindElement('TopDE',MCF7_DiffGenes_FandT10)
}
