### R code from vignette source 'IntClust_Vignette.Rnw'

###################################################
### code chunk number 1: Install (eval = FALSE)
###################################################
## install.packages("IntClust")
## library(IntClust)
## data(fingerprintMat)
## data(targetMat)


###################################################
### code chunk number 2: SingleClust (eval = FALSE)
###################################################
## MCF7_F <- Cluster(Data=fingerprintMat,type="data",distmeasure="tanimoto",
##           normalize=FALSE,method=NULL,clust="agnes",linkage="flexible",gap=FALSE) 
## 		
## MCF7_T <- Cluster(Data=targetMat,type="data",distmeasure="tanimoto",
##           normalize=FALSE,method=NULL,clust="agnes",linkage="flexible",gap=FALSE)


###################################################
### code chunk number 3: NrClusters (eval = FALSE)
###################################################
## List=list(fingerprintMat,targetMat)
## NrClusters=SelectnrClusters(List=List,type="data",distmeasure=c("tanimoto",
##            "tanimoto"),nrclusters=seq(5,20),normalize=c(FALSE,FALSE),
##             names=c("FP","TP"))


###################################################
### code chunk number 4: Colours (eval = FALSE)
###################################################
## Colours <- ColorPalette(colors=c("chocolate","firebrick2","darkgoldenrod2",
##            "darkgreen","blue2","darkorchid3","deeppink"),ncols=7) 


###################################################
### code chunk number 5: Fig1 (eval = FALSE)
###################################################
## ClusterPlot(Data1=MCF7_F,nrclusters=7,cols=Colours,main="Clustering on 
##            Fingerprints: Dendrogram",ylim=c(-0.1,1.8))
## 
## ClusterPlot(Data1=MCF7_T,nrclusters=7,cols=Colours,colorComps=NULL,main="Clustering on 
##             Targets: Dendrogram",ylim=c(-0.1,2.5))


###################################################
### code chunk number 6: ADC (eval = FALSE)
###################################################
## L=list(fingerprintMat,targetMat)
## MCF7_ADC=ADC(List=L,distmeasure="tanimoto",normalize=FALSE,clust="agnes",
##          linkage="flexible")


###################################################
### code chunk number 7: Weighted (eval = FALSE)
###################################################
## L=list(fingerprintMat,targetMat)
## MCF7_Weighted=WeightedClust(L,type="data",distmeasure=c("tanimoto","tanimoto"),
##               normalize=c(FALSE,FALSE),weight=seq(0,1,0.1),weightclust=0.5,
##               StopRange=FALSE)


###################################################
### code chunk number 8: WonM (eval = FALSE)
###################################################
## L=list(fingerprintMat,targetMat)
## MCF7_WonM=WonM(List=L,type="data",distmeasure=c("tanimoto","tanimoto"),
##           normalize=c(FALSE,FALSE),nrclusters=seq(5,25),linkage=
##           c("flexible","flexible"))


###################################################
### code chunk number 9: ComparePlot (eval = FALSE)
###################################################
## L=list(MCF7_F,MCF7_ADC,MCF7_WonM,MCF7_Weighted,MCF7_T)
## N=c("FP","ADC","WonM",paste("Weight",seq(1,0,-0.1),sep=" "),"TP")
## ComparePlot(L,nrclusters=7,cols=Colours,fusionsLog=TRUE,weightclust=FALSE,names=N,
## margins=c(9.1,4.1,4.1,4.1),plottype="new",location=NULL)


###################################################
### code chunk number 10: FindCluster (eval = FALSE)
###################################################
## Comps=FindCluster(List=L,nrclusters=7,select=c(1,6))
## Comps 


###################################################
### code chunk number 11: Track (eval = FALSE)
###################################################
## Tracking=TrackCluster(List=L,Selection=Comps,nrclusters=7,followMaxComps=FALSE,
## 		followClust=TRUE,fusionsLog=TRUE,weightclust=FALSE,
## 		names=N,selectionPlot=FALSE,table=FALSE,legendposy=2.4,
## 		completeSelectionPlot=TRUE,cols=Colours,plottype="sweave",
## 		location=NULL)


###################################################
### code chunk number 12: Features (eval = FALSE)
###################################################
## MCF7_Feat=ChooseCluster(Interactive=FALSE,leadCpds=list(Comps),clusterResult=MCF7_F,
##           colorLab=MCF7_F,binData=list(fingerprintMat,targetMat),datanames=
##           c("FP","TP"),topChar = 20,topG = 20)


###################################################
### code chunk number 13: FeatPlot (eval = FALSE)
###################################################
## BinFeaturesPlot_SingleData(leadCpds=Comps,orderLab=MCF7_F,features=MCF7_Feat
##                            $Characteristics$FP$TopFeat$Names,data=fingerprintMat,
##                            colorLab=MCF7_F,nrclusters=7,cols=Colours,name=c("FP"))
##                  
## BinFeaturesPlot_SingleData(leadCpds=Comps,orderLab=MCF7_F,features=MCF7_Feat
##                            $Characteristics$TP$TopFeat$Names,data=targetMat,
##                            colorLab=MCF7_F,nrclusters=7,cols=Colours,name=c("TP"))


###################################################
### code chunk number 14: DEGenes (eval = FALSE)
###################################################
## data(geneMat)
## MCF7_Genes=DiffGenes(List=NULL,Selection=Comps,geneExpr=geneMat,method="limma",
##            sign=0.05,topG=10)
## Genes=MCF7_Genes$Selection$Genes$TopDE$ID				  


###################################################
### code chunk number 15: GenePlot (eval = FALSE)
###################################################
## ProfilePlot(Genes=Genes[1:5],Comps=Comps,geneExpr=geneMat,raw=FALSE,
##             order=MCF7_F,color=MCF7_F,nrclusters=7,cols=Colours,
##             addLegend=TRUE,margins=c(8.1,4.1,1.1,6.5),plottype="sweave",
##             location=NULL)


###################################################
### code chunk number 16: Pathways (eval = FALSE)
###################################################
## data(GeneInfo)
## data(GS)
## L=list(MCF7_Genes)
## 
## MCF7_Paths=PathwayAnalysis(List=L,Selection=Comps,geneExpr=geneMat,
##            method = c("limma","MLP"),geneInfo=GeneInfo,
##            geneSetSource="GOBP",topP = NULL,topG=NULL,GENESET=GS,
##            sign = 0.05,niter=2)
## 
## PlotPathways(MCF7_Paths$Selection$Pathways)


###################################################
### code chunk number 17: sessionInfo
###################################################
toLatex(sessionInfo())


