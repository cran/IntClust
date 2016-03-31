\name{TrackCluster}
\alias{TrackCluster}
\title{
Follow a cluster over multiple methods
}
\description{
It is often desired to track a specific selection of object over the
different methods and/or weights. This can be done with the 
\code{ClusterDistribution}. For every method, it is tracked where the objects
of the selections are situated.
}
\usage{
TrackCluster(List, Selection, nrclusters=NULL, followMaxComps = FALSE, 
followClust = TRUE, fusionsLog = TRUE, WeightClust = TRUE, names = NULL
, SelectionPlot = TRUE, Table = TRUE, CompleteSelectionPlot = FALSE, 
ClusterPlot=FALSE,cols=NULL,legendposx=0.5,legendposy=2.4,plottype="new",location=NULL)
}
\arguments{
  \item{List}{
A list of the clustering outputs.The first element of the
list will be used as the reference in \code{ReorderToReference}.
}
  \item{Selection}{
The selection of objects to follow or a specific cluster number.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in.
}
  \item{followMaxComps}{
Logical for plot. Whether to follow the maximum of objects.
}
  \item{followClust}{
Logical for plot. Whether to follow the specific cluster.
}
  \item{fusionsLog}{
To be handed to \code{ReorderToReference}.
}
  \item{WeightClust}{
To be handed to \code{ReorderToReference}.
}
  \item{names}{
Optional. Names of the methods.
}
  \item{SelectionPlot}{
Logical. Should a plot be produced. Depending on followMaxComps and followClust
it focuses on the maximum of compounds or a cluster. It will not be indicated to which cluster compounds moved.
}
  \item{Table}{
Logical. Should a table with the compounds per method and the shared compounds
be produced?
}
  \item{CompleteSelectionPlot}{
Logical. Should the complete distribution of the selection be plotted? This implies that it will be indicated to which cluster compounds will move.
}
  \item{ClusterPlot}{
Logical. Plot of specific cluster.
}
  \item{cols}{
The colors used for the different clusters.
}
  \item{legendposx}{
The x-coordinate of the legend on all plots. 
 }
 \item{legendposy}{
The y-coordinate of the legend on all plots. 
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
\details{
The result is provided with extra information as which compounds of the
original selection can be found in this cluster and which are extra. 
Further, plots of the distribution of the compounds can be produced. 
One plot follows the complete distribution of the cluster while 
another one focuses on either the maximum number of compounds or a 
specific cluster, whatever is specified. It are the number of compounds 
that are plotted and the first element indicated the number of compounds
in the selection. A table can be produced as well, that separates the 
objects that are shared over all methods from those extra in the original
selection and extra for the other methods. The \code{ReorderToReference} is applied
to make sure that the clusters are comparable over the methods.

The function is experimental and might not work in specific cases. Please let us know such that we can improve its functionality.
}
\value{
The returned value is a list with an element for every method. This element
is another list with the following elements:
\item{Selection}{The selection of compounds to follow}
\item{nr.clusters}{the number of clusters the selection is divided over}
\item{nr.min.max.together }{the minimum and maximum number of compounds
found together}
\item{perc.min.max.together}{minimum and maximum percentage of compounds
found together}
\item{AllClusters}{A list with an element per cluster that contains at least one 
of the compounds in Selection. The list contains the cluster number, the complete
cluster, the objects that originally could be found in this cluster and which
object were joined extra to it.}
Depending on whether followMaxComps or followClust is specified, the
cluster of interest is mentioned separately as well for easy access. If the 
option was specified to create a table, this can be found under the Table
element. Each plot that was specified to be created is plotted in a new window
in the graphics console.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
data(Colors1)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)

L=list(MCF7_F,MCF7_T)
names=c("FP","TP")

Comps=FindCluster(L,nrclusters=7,select=c(1,4))
Comps

CompsFPAll=TrackCluster(List=L,Selection=Comps,nrclusters=7,followMaxComps=TRUE,
followClust=FALSE,fusionsLog=TRUE,WeightClust=TRUE,names=names,SelectionPlot=TRUE,
Table=TRUE,CompleteSelectionPlot=TRUE,cols=Colors1,plottype="new")

}
