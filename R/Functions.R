## THE REIDS PACKAGE ##


## IMPORTS ##
#' @import Biobase 
#' @import cluster
#' @import stats
#' @import graphics
#' @import grDevices
#' @import ggplot2
#' @import circlize
#' @import data.table
#' @importFrom Rdpack reprompt
#' @importFrom igraph graph.adjacency as_data_frame minimum.spanning.tree plot.igraph V E
#' @importFrom data.table rbindlist
#' @importFrom limma topTable lmFit eBayes
#' @importFrom ade4 acm.disjonctif dist.binary
#' @importFrom gridExtra grid.arrange
#' @importFrom plyr join_all try_default
#' @importFrom gplots heatmap.2
#' @importFrom gtools permute permutations
#' @importFrom plotrix color2D.matplot
#' @importFrom e1071 hamming.distance
#' @importFrom pls stdize
#' @importFrom utils combn
#' @importFrom FactoMineR MCA
#' @importFrom analogue distance
#' @importFrom lsa cosine
#' @importFrom SNFtool affinityMatrix SNF


## FUNCTIONS ##

## Multi-Source Clustering Procedures ##

# Direct Clustering Methods

#' @title Aggregated data clustering
#' 
#' @description Aggregated Data Clustering (ADC) is a direct clustering multi-source technique. ADC merges the columns of all data sets into a single large data set on which a final clustering is performed.
#' 
#' @export
#' @param List A list of data matrices of the same type. It is assumed the rows are corresponding with the objects.
#' @param distmeasure Choice of metric for the dissimilarity matrix (character). Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to "tanimoto".
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to FALSE. This is recommended if different distance types are used. More details on normalization in \code{Normalization}.
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is NULL.
#' @param clust Choice of clustering function (character). Defaults to "agnes". 
#' @param linkage Choice of inter group dissimilarity (character). Defaults to "flexible".
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible".
#' @return The returned value is a list with the following three elements.
#' 	\item{AllData}{Fused data matrix of the data matrices}
#' \item{DistM}{The distance matrix computed from the AllData element}
#' 	\item{Clust}{The resulting clustering}
#' 	The value has class "ADC". The Clust element will be of interest for further applications.
#' @details In order to perform aggregated data clustering, the \code{ADC} function was written. A list of data matrices of the same type (continuous or binary) is required as input which are combined into a single (larger) matrix. Hierarchical clustering is performed
#' 	with the \code{agnes} function and the ward link on the resulting data matrix and an applicable distance measure is indicated by the user.
#' @references 
#' \insertRef{Fodeh2013}{IntClust}
#' @examples
#' data(fingerprintMat)
#' data(targetMat)
#' L=list(fingerprintMat,targetMat)
#' MCF7_ADC=ADC(List=L,distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",
#' linkage="flexible",alpha=0.625)
ADC<-function(List,distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="flexible",alpha=0.625){
	
	#Checking required data types and methods:
	if(class(List) != "list"){
		stop("Data must be of type list")
	}
	
	
	#Fuse variables into 1 data Matrix
	
	OrderNames=rownames(List[[1]])
	for(i in 1:length(List)){
		List[[i]]=List[[i]][OrderNames,]
	}
	
	AllData<-NULL
	for (i in 1:length(List)){
		if(i==1){
			AllData=List[[1]]
		}
		else{
			AllData=cbind(AllData,List[[i]])
		}
	}
	
	#Compute Distance Matrix on AllData
	
	AllDataDist=Distance(AllData,distmeasure,normalize,method)
	
	#Perform hierarchical clustering with ward link on distance matrix
	
	HClust = cluster::agnes(AllDataDist,diss=TRUE,method=linkage,par.method=alpha)		
	
	
	out=list(AllData=AllData,DistM=AllDataDist,Clust=HClust)
	attr(out,'method')<-'ADC'	
	return(out)
	
}


#' @title Aggregated data ensemble clustering
#' 
#' @description Aggregated Data Ensemble Clustering (ADEC) is a direct clustering multi-source technique. ADEC is an iterative procedure which starts with the merging of the data sets. In each iteration, a random sample of the features is selected and/or a resulting dendrogram is divided into k clusters for a range of values of k.
#' 
#' @export
#' @param List A list of data matrices of the same type. It is assumed the rows are corresponding with the objects.
#' @param distmeasure Choice of metric for the dissimilarity matrix (character). Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to "tanimoto".
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to FALSE. This is recommended if different distance types are used. More details on normalization in \code{Normalization}.
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is NULL.
#' @param t The number of iterations. Defaults to 10.
#' @param r The number of features to take for the random sample. If NULL (default), all features are considered.
#' @param nrclusters A sequence of numbers of clusters to cut the dendrogram in. If NULL (default), the function stops.
#' @param clust Choice of clustering function (character). Defaults to "agnes". 
#' @param linkage Choice of inter group dissimilarity (character). Defaults to "flexible".
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible".
#' @return The returned value is a list with the following three elements.
#' \item{AllData}{Fused data matrix of the data matrices}
#' \item{DistM}{The resulting co-association matrix}
#' \item{Clust}{The resulting clustering}
#' The value has class 'ADEC'. The Clust element will be of interest for further applications.
#' @details If r is specified and nrclusters is a fixed number, only a random sampling of the features will be performed for the t iterations (ADECa). If r is NULL and the nrclusters is a sequence, the clustering is performedon all features and the dendrogam is divided into clusters for the values of nrclusters (ADECb). If both r is specified and nrclusters is a sequence, the combination is performed (ADECc).
#' After every iteration, either be random sampling, multiple divisions of the dendrogram or both, an incidence matrix is set up. All incidence matrices are summed and represent the distance matrix on which a final clustering is performed. 
#' @references 
#' \insertRef{Fodeh2013}{IntClust}
#' @examples
#' data(fingerprintMat)
#' data(targetMat)
#' L=list(fingerprintMat,targetMat)
#' MCF7_ADEC=ADEC(List=L,distmeasure="tanimoto",normalize=FALSE,method=NULL,t=100, 
#' r=100,nrcluster=seq(1,10,1),clust="agnes",linkage="flexible",alpha=0.625)
ADEC<-function(List,distmeasure="tanimoto",normalize=FALSE,method=NULL,t=10,r=NULL,nrclusters=NULL,clust="agnes",linkage="flexible",alpha=0.625){
	
	if(class(List) != "list"){
		stop("Data must be of type lists")
	}
	
	if(is.null(nrclusters)){
		stop("Give a number of cluters to cut the dendrogram into.")
	}
	
	if(!is.null(r)&length(nrclusters)==1){
		message("Performing a random sampling of the features with a fixed number of clusters.")
	}else if(is.null(r)&length(nrclusters)>1){
		message("Dividing the dendrogram in k clusters for a range of values of k.")
		t=1
	}else if(!is.null(r)&length(nrclusters)>1){
		message("Performing a random sampling of the features and dividing the dendrogram in k clusters for a range of values of k.")
	}
	else{
		stop("Specify r and/or nrclusters in order to perform an ADEC method.")
	}
	
	#Fuse A1 and A2 into 1 Data Matrix
	
	OrderNames=rownames(List[[1]])
	for(i in 1:length(List)){
		List[[i]]=List[[i]][OrderNames,]
	}
	
	AllData<-NULL
	for (i in 1:length(List)){
		if(i==1){
			AllData=List[[1]]
		}
		else{
			AllData=cbind(AllData,List[[i]])
		}
	}
	
	
	#take random sample of features
	
	nc=ncol(AllData)
	evenn=function(x){if(x%%2!=0) x=x-1 else x}
	nc=evenn(nc)
	
	
	#Put up Incidence matrix
	Incidence=matrix(0,dim(List[[i]])[1],dim(List[[i]])[1])
	rownames(Incidence)=rownames(AllData)
	colnames(Incidence)=rownames(AllData)
	
	
	#Repeat for t iterations
	
	
	for(g in 1:t){
		#message(g)
		
		#if r is not fixed: changes per iteration. Need 1 value for r.
		if(is.null(r)){
			r=ncol(AllData)
		}
		
		#take random sample:
		ZeroPresent=TRUE
		while(ZeroPresent){
			temp1=sample(ncol(AllData),r,replace=FALSE)
			
			A_prime=AllData[,temp1]
			
			if(all(rowSums(A_prime)!=0)){
				ZeroPresent=FALSE
			}
			
		}
		#Step 2: apply hierarchical clustering on A_prime  + cut tree into nrclusters
		
		DistM=Distance(A_prime,distmeasure,normalize,method)
		
		HClust_A_prime=cluster::agnes(DistM,diss=TRUE,method=linkage,par.method=alpha)
		
		
		for(k in 1:length(nrclusters)){
			#message(k)
			Temp=stats::cutree(HClust_A_prime,nrclusters[k])	
			MembersofClust=matrix(1,dim(List[[1]])[1],dim(List[[1]])[1])
			
			for(l in 1:length(Temp)){
				label=Temp[l]
				sameclust=which(Temp==label)
				MembersofClust[l,sameclust]=0
			}
			Incidence=Incidence+MembersofClust
		}	
		
	}
	
	Clust=cluster::agnes(Incidence,diss=TRUE,method="ward",par.method=alpha)
	
	out=list(AllData=AllData,DistM=Incidence,Clust=Clust)
	attr(out,'method')<-'ADEC'	
	return(out)
	
}

# Similarity-based approaches

#' @title Weighted clustering
#' 
#' @description Weighted Clustering (Weighted) is a similarity-based multi-source clustering technique. Weighted clustering is performed with the function \code{WeightedClust}. Given a
#' list of the data matrices, a dissimilarity matrix is computed of each with
#' the provided distance measures. These matrices are then combined resulting
#' in a weighted dissimilarity matrix. Hierarchical clustering is performed
#' on this weighted combination with the agnes function and the ward link
#'  
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clusters".	
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}. 
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param StopRange Logical. Indicates whether the distance matrices with values not between zero and one should be standardized to have values between zero and one.
#' If FALSE the range normalization is performed. See \code{Normalization}. If TRUE, the distance matrices are not changed.
#' This is recommended if different types of data are used such that these are comparable. Default is FALSE.
#' @param weight Optional. A list of different weight combinations for the data sets in List.
#' If NULL, the weights are determined to be equal for each data set.
#' It is further possible to fix weights for some data matrices and to
#'  let it vary randomly for the remaining data sets. Defaults to seq(1,0,-0.1).  An example is provided in the details.
#' @param weightclust A weight for which the result will be put aside of the other results. This was done for comparative reason and easy access. Defaults to 0.5 (two data sets)
#' @param clust Choice of clustering function (character). Defaults to "agnes". 
#' @param linkage Choice of inter group dissimilarity (character) for the final clustering. Defaults to "ward".
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible".
#' @return The returned value is a list of four elements:
#' \item{DistM}{A list with the distance matrix for each data structure}
#' \item{WeightedDist }{A list with the weighted distance matrices}
#' \item{Results}{The hierarchical clustering result for each element in WeightedDist}
#' \item{Clust}{The result for the weight specified in Clustweight}
#' The value has class 'Weighted'.
#' @details The weight combinations should be provided as elements in a list. For three data
#' matrices an example could be: weights=list(c(0.5,0.2,0.3),c(0.1,0.5,0.4)). To provide
#' a fixed weight for some data sets and let it vary randomly for others, the element "x" 
#' indicates a free parameter. An example is weights=list(c(0.7,"x","x")). The weight 0.7 
#' is now fixed for the first data matrix while the remaining 0.3 weight will be divided over
#' the other two data sets. This implies that every combination of the sequence from 0 to 0.3
#' with steps of 0.1 will be reported and clustering will be performed for each.
#' @references 
#' \insertRef{PerualilaTan2016}{IntClust}
#' @examples
#' data(fingerprintMat)
#' data(targetMat)
#' L=list(fingerprintMat,targetMat)
#' 
#' MCF7_Weighted=WeightedClust(List=L,type="data",distmeasure=c("tanimoto","tanimoto"),
#' normalize=c(FALSE,FALSE),method=c(NULL,NULL),StopRange=FALSE,weight=seq(1,0,-0.1),
#' weightclust=0.5,clust="agnes",linkage="ward",alpha=0.625)
WeightedClust <- function(List,type=c("data","dist","clusters"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),StopRange=FALSE,weight=seq(1,0,-0.1),weightclust=0.5,clust="agnes",linkage="ward",alpha=0.625){ # weight = weight to data1
	
	
	#Step 1: compute distance matrices:
	type<-match.arg(type)
	
	CheckDist<-function(Dist,StopRange){
		if(StopRange==FALSE & !(0<=min(Dist) & max(Dist)<=1)){
			#message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
			Dist=Normalization(Dist,method="Range")
		}
		else{
			Dist=Dist
		}
	}
	
	
	if(type=="data"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,,drop=FALSE]
		}
		Dist=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i],normalize[i],method[i]))
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
	}
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		Dist=List
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
	}
	else{
		Dist=lapply(seq(length(List)),function(i) return(List[[i]]$Dist))
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		OrderNames=rownames(DistM[[1]])
		for(i in 1:length(DistM)){
			DistM[[i]]=DistM[[i]][OrderNames,OrderNames]
		}
	}
	
	#Step 2: Weighted linear combination of the distance matrices:
	if(is.null(weight)){
		equalweights=1/length(List)
		weight=list(rep(equalweights,length(List)))
		
	}
	else if(class(weight)=='list' & length(weight[[1]])!=length(List)){
		stop("Give a weight for each data matrix or specify a sequence of weights")
	}
	
	if(class(weight)!="list"){
		condition<-function(l){		
			l=as.numeric(l)
			if( sum(l)==1 ){  #working with characters since with the numeric values of comb or permutations something goes not the way is should: 0.999999999<0.7+0.3<1??
				#return(row.match(l,t1))
				return(l)
			}
			else(return(0))
		}
		
		if(all(seq(1,0,-0.1)!=weight)){
			for(i in 1:length(weight)){
				rest=1-weight[i]
				if(!(rest%in%weight)){
					weight=c(weight,rest)
				}
			}
		}
		
		
		
		
		t1=gtools::permutations(n=length(weight),r=length(List),v=as.character(weight),repeats.allowed = TRUE)
		t2=lapply(seq_len(nrow(t1)), function(i) if(sum(as.numeric(t1[i,]))==1) return(as.numeric(t1[i,])) else return(0)) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
		t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
		t4=t2[which(t3!=0)]
		weight=lapply(seq(length(t4)),function(i) rev(t4[[i]]))
		
	}
	if(class(weight)=="list" & "x" %in% weight[[1]]){ #x indicates a free weight
		newweight=list()
		for(k in 1:length(weight)){
			w=weight[[k]]
			weightsfordata=which(w!="x") #position of the provided weight = position of the data to which the weight is given
			givenweights=as.numeric(w[weightsfordata])
			
			stilltodistribute=1-sum(givenweights)
			
			newweights=seq(stilltodistribute,0,-0.1)
			
			t1=gtools::permutations(n=length(newweights),r=length(List)-length(weightsfordata),v=as.character(newweights),repeats.allowed = TRUE)
			Input1=as.list(seq_len(nrow(t1)))
			Input2=lapply(seq(length(Input1)),function(i) {Input1[[i]][length(Input1[[i]])+1]=stilltodistribute
						return(Input1[[i]])})
			t2=lapply(seq(length(Input2)), FUN=function(i){if(sum(as.numeric(t1[Input2[[i]][1],])+0.00000000000000002775)==Input2[[i]][2]) return(as.numeric(t1[i,])) else return(0)}) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
			t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
			weightsforotherdata=t2[which(t3!=0)]
			
			new=list()
			for(i in 1:length(weightsforotherdata)){
				w1=weightsforotherdata[[i]]
				new[[i]]=rep(0,length(List))
				new[[i]][weightsfordata]=givenweights
				new[[i]][which(new[[i]]==0)]=w1
			}
			
			newweight[k]=new
		}
		
		weight=newweight
	}
	weightedcomb<-function(w,Dist){
		temp=lapply(seq_len(length(Dist)),function(i) w[i]*Dist[[i]])
		temp=Reduce("+",temp)	
		return(temp)
	}
	
	DistClust=NULL
	Clust=NULL
	
	DistM=lapply(seq(length(weight)),function(i) weightedcomb(weight[[i]],Dist=Dist))
	namesweights=c()
	WeightedClust=lapply(seq(length(weight)),function(i) cluster::agnes(DistM[[i]],diss=TRUE,method=linkage,par.method=alpha))
	for(i in 1:length(WeightedClust)){
		namesweights=c(namesweights,paste("Weight",weight[i],sep=" "))
		if(all(weight[[i]]==weightclust)){
			Clust=WeightedClust[[i]]	
			DistClust=DistM[[i]]
		}
	}	
	
	if(is.null(DistClust)){
		DistClust=weightedcomb(weightclust,Dist=Dist)
		Temp=cluster::agnes(DistClust,diss=TRUE,method=linkage,par.method=alpha)
		Clust=Temp
	}
	
	Results=lapply(seq(1,length(WeightedClust)),function(i) return(c("DistM"=DistM[i],"Clust"=WeightedClust[i])))
	names(Results)=namesweights
	
	# return list with objects
	out=list(Dist=Dist,Results=Results,Clust=list("DistM"=DistClust,"Clust"=Clust))
	attr(out,'method')<-'Weighted'
	return(out)
	
}


#' @title Similarity network fusion
#' 
#' @description Similarity Network Fusion (SNF) is a similarity-based multi-source clustering technique. SNF consists of two steps. In the initial step a similarity network is set up for each data matrix. The network is the visualization of the similarity matrix as a weighted graph with the objects as vertices and the pairwise similarities as weights on the edges. In the network-fusion step, each network is iteratively updated with information of the other network which results in more alike networks every time. This eventually converges to a single network.
#' @export
#' @param List A list of data matrices of the same type. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clusters".
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}. 
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param StopRange Logical. Indicates whether the distance matrices with values not between zero and one should be standardized to have so.
#' If FALSE the range normalization is performed. See \code{Normalization}. If TRUE, the distance matrices are not changed.
#' This is recommended if different types of data are used such that these are comparable. Default is FALSE.
#' @param NN The number of neighbours to be used in the procedure. Defaults to 20.
#' @param mu The parameter epsilon. The value is recommended to be between 0.3 and 0.8. Defaults to 0.5.
#' @param T The number of iterations.
#' @param clust Choice of clustering function (character). Defaults to "agnes". 
#' @param linkage Choice of inter group dissimilarity (character) for the final clustering. Defaults to "ward".
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible"
#' @return The returned value is a list with the following three elements.
#' \item{FusedM }{The fused similarity matrix}
#' \item{DistM }{The distance matrix computed by subtracting FusedM from one}
#' \item{Clust}{The resulting clustering}
#' The value has class 'SNF'.
#' @details If r is specified and nrclusters is a fixed number, only a random sampling of the features will be performed for the t iterations (ADECa). If r is NULL and the nrclusters is a sequence, the clustering is performedon all features and the dendrogam is divided into clusters for the values of nrclusters (ADECb). If both r is specified and nrclusters is a sequence, the combination is performed (ADECc).
#' After every iteration, either be random sampling, multiple divisions of the dendrogram or both, an incidence matrix is set up. All incidence matrices are summed and represent the distance matrix on which a final clustering is performed. 
#' @references 
#' \insertRef{Wang2014a}{IntClust}
#' @examples
#' data(fingerprintMat)
#' data(targetMat)
#' L=list(fingerprintMat,targetMat)
#' MCF7_SNF=SNF(List=L,type="data",distmeasure=c("tanimoto","tanimoto"),normalize=
#' c(FALSE,FALSE),method=c(NULL,NULL),StopRange=FALSE,NN=10,mu=0.5,T=20,clust="agnes",
#' linkage="ward",alpha=0.625)
SNF=function(List,type=c("data","dist","clusters"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),StopRange=FALSE,NN=20,mu=0.5,T=20,clust="agnes",linkage="ward",alpha=0.625){
	
	#Checking required data types and methods:
	if(class(List) != "list"){
		stop("Data must be of type list")
	}
	
	if(mu<0.3 | mu >0.8){
		message("Warning: mu is recommended to be between 0.3 and 0.8 for the SNF method. Default is 0.5.")
	}
	
	
	CheckDist<-function(Dist,StopRange){
		if(StopRange==FALSE & !(0<=min(Dist) & max(Dist)<=1)){
			#message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
			Dist=Normalization(Dist,method="Range")
		}
		else{
			Dist=Dist
		}
	}
	
	
	#STEP 1: Distance Matrices
	if(type=="data"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,,drop=FALSE]
		}
		DistM=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i],normalize[i],method[i]))
		DistM=lapply(seq(length(DistM)),function(i) CheckDist(DistM[[i]],StopRange))
	}
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		DistM=List
		DistM=lapply(seq(length(DistM)),function(i) CheckDist(DistM[[i]],StopRange))
	}
	else{
		DistM=lapply(seq(length(List)),function(i) return(List[[i]]$DistM))
		DistM=lapply(seq(length(DistM)),function(i) CheckDist(DistM[[i]],StopRange))
		
		OrderNames=rownames(DistM[[1]])
		for(i in 1:length(DistM)){
			DistM[[i]]=DistM[[i]][OrderNames,OrderNames]
		}
	}
	
	
	#STEP 2: Affinity Matrices
	
	AffM=lapply(seq(length(List)), function(i) SNFtool::affinityMatrix(DistM[[i]], NN, mu))
	
	#STEP 3: Fuse Networks Into 1 Single Network
	
	SNF_FusedM=SNFtool::SNF(AffM, NN, T)
	rownames(SNF_FusedM)=rownames(List[[1]])
	colnames(SNF_FusedM)=rownames(List[[1]])
	Dist=1-SNF_FusedM
	
	#STEP 4: Perform Hierarchical Clustering with WARD Link
	
	HClust = cluster::agnes(Dist,diss=TRUE,method=linkage,par.method=alpha)		
	
	
	#Output= list with the fused matrix and the performed clustering
	out=list(SNF_FusedM=SNF_FusedM,DistM=Dist,Clust=HClust)
	attr(out,'method')<-'SNF'
	return(out)
}

# Graph-based consensus approaches

#' @title Ensemble clustering
#' 
#' @description The \code{EnsembleClustering} includes the ensemble clustering methods CSPA, HGPA and MCLA which are graph-based consensus methods.
#'  
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clusters".	
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}. 
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param clust Choice of clustering function (character). Defaults to "agnes". 
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible,", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible"
#' @param nrclusters The number of clusters to divide each individual dendrogram in. Default is c(7,7) for two data sets.
#' @param gap Logical. Whether the optimal number of clusters should be determined with the gap statistic. Default is FALSE.
#' @param maxK The maximal number of clusters to investigate in the gap statistic. Default is 15.
#' @param ensembleMethod The method to be performed: "CSPA", "HGPA", "MCLA" or "Best".
#' @param waitingtime The time in seconds to wait until the MATLAB results are generated. Defaults to 300.
#' @param file_number The specific file number to be placed as a tag in the file generated by MATLAB. Defaults to 00.
#' @param executable Logical. Whether the MATLAB functions are performed via an executable on the command line (TRUE, only possible for Linux systems) or by calling on MATLAB directly (FALSE). Defaults to FALSE. The files EnsembleClusteringC.m (CSPA), EnsembleClusteringH.m (HGPA), EnsembleClusteringM.m (MCLA) and MetisAlgorithm.m are present in the inst folder to be transformed in executables.
#' @return The returned value is a list of two elements:
#' \item{DistM}{A list with the distance matrix for each data structure}
#' \item{Clust}{The resulting clustering}
#' The value has class 'Ensemble'.
#' @details \insertCite{Strehl2002}{IntClust} introduce three heuristic algorithms to solve the cluster ensemble problem. 
#' Each method starts by transforming the clustering solutions into a single hypergraph in which a hyperedge represents a single cluster.
#' The Cluster-based Similarity Partitioning Algorithm (CSPA) transforms the hypergraph into an overall similarity matrix which entries 
#' represent the fraction of clusterings in which two objects are in the same cluster. The similarity matrix is considered as a 
#' graph and the objects are reclustered with the graph partitioning algorithm METIS \insertCite{Karypis1998}{IntClust}. Hyper-Graph Partitioning 
#' Algorithm (HGPA) partitions the hypergraph directly by cutting a minimal number of hyperedges. It aims to obtain connected components of 
#' approximately the same dimension. The partitioning algorithm is HMetis \insertCite{Karypis1997}{IntClust}. The Meta-CLustering Algorithm (MCLA) 
#' computes a similarity between the hyperedges (clusters) based on the Jaccard index. The resulting similarity matrix is used to build 
#' a meta-graph which is partitioned by the METIS algorithm \insertRef{Karypis1998}{IntClust} into resulting meta-clusters. 
#' The final partition of the objects is obtaining by appointing each object to the meta-cluster to which it is assigned the most. The \code{R} code calls on the MATLAB code provided by \insertCite{Strehl2002}{IntClust}. The MATLAB functions are included in the inst folder and should be located in the working directory. Shell script for the executable can be found in the inst folder as well.
#' @references \insertRef{Strehl2002}{IntClust}  \insertRef{Karypis1997}{IntClust}  \insertRef{Karypis1998}{IntClust}  
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' L=list(fingerprintMat,targetMat)
#' 
#' MCF7_CSPA=EnsembleClustering(List=L,type="data",distmeasure=c("tanimoto",
#' "tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),StopRange=FALSE,
#' clust="agnes",linkage=c("flexible","flexible"),nrclusters=c(7,7),gap=FALSE,
#' maxK=15,ensembleMethod="CSPA",executable=FALSE)
#' 
#' }
EnsembleClustering<-function(List,type=c("data","dist","clust"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust = "agnes", linkage = c("flexible","flexible"),alpha=0.625,nrclusters=7,gap = FALSE, maxK = 15,ensembleMethod=c("CSPA","HGPA","MCLA","Best"),waitingtime=300,file_number=00,executable=FALSE){
	# Step 1: Generate several clustering results on the objects
	# Is the clustering of the objects on each data set enough (sometimes only 2 data sets) or should we produce multiple for each soource?
	
	if(type=="data"){
		
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,,drop=FALSE]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type="data",distmeasure[i],normalize[i],method[i],clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))
		
		Dist=lapply(seq(length(List)),function(i) Clusterings[[i]]$DistM)
	
		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=clusters
			}
		}
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
	}
	
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize=FALSE,method=NULL,clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))
		
		Dist=List

		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=ceiling(mean(clusters))
			}
		}
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
		
		
	}
	else if(type=="clust"){
		
		Clusterings=List
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		if(is.null(nrclusters)){
			stop("Please specify a number of clusters")
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
		
	}
	#### Step 1 Complete
	
	# Step 2: Get data in format to be given to matlab =>  matrix with a row per clustering result
	if(type=="data" | type=="dist"){
		nc=nrow(List[[1]])
	}
	else if(type=="clust"){
		nc=nrow(List[[1]]$DistM)
	}
	
	matlabdata=matrix(unlist(Clusters),ncol=nc,byrow=TRUE)
	utils::write.table(matlabdata,file=paste("matlabdata_",file_number,".csv",sep=""),sep=",",col.names=FALSE,row.names=FALSE)
	
	if(ensembleMethod=="CSPA"){
		if(executable){
			system(paste("./EnsembleClusteringC ",file_number,sep=""),intern=TRUE)
		}
		else{
		matlabcode=c(
				paste("cls = csvread('matlabdata_",file_number,".csv')",sep=""),
				paste("ClusterEnsemble_",file_number,"=cspa(cls,",file_number,")",sep=""),
				paste("csvwrite('ClusterEnsemble_",file_number,".csv', ClusterEnsemble_",file_number,")",sep=""))
		#writeLines(matlabcode, con=paste("EnsembleClusteringC_",file_number,".m",sep=""))
		
		#file.remove(file=paste("EnsembleClusteringC_",file_number,".m",sep=""))
		}
	}	
	if(ensembleMethod=="HGPA"){
		if(executable){
			system(paste("./EnsembleClusteringH ",file_number,sep=""),intern=TRUE)
		}
		else{
		matlabcode=c(
				paste("cls = csvread('matlabdata_",file_number,".csv')",sep=""),
				paste("ClusterEnsemble_",file_number,"=hgpa(cls,",file_number,")",sep=""),
				paste("csvwrite('ClusterEnsemble_",file_number,".csv', ClusterEnsemble_",file_number,")",sep=""))
		}
#		writeLines(matlabcode, con=paste("EnsembleClusteringH_",file_number,".m",sep=""))
		
	}	
	if(ensembleMethod=="MCLA"){
		if(executable){
			system(paste("./EnsembleClusteringM ",file_number,sep=""),intern=TRUE)
		}
		else{
		matlabcode=c(
				paste("cls = csvread('matlabdata_",file_number,".csv')",sep=""),
				paste("ClusterEnsemble_",file_number,"=mcla(cls,",file_number,")",sep=""),
				paste("csvwrite('ClusterEnsemble_",file_number,".csv', ClusterEnsemble_",file_number,")",sep=""))
		}
#		writeLines(matlabcode, con=paste("EnsembleClusteringM_",file_number,".m",sep=""))
		
	}
	if(ensembleMethod=="Best"){
		if(executable){
			system(paste("matlab -nodisplay -r \"run('EnsembleClusteringB ",file_number,"'); exit\"",sep=""),intern=TRUE)
		}
		else{
		matlabcode=c(
				paste("cls = csvread('matlabdata_",file_number,".csv')",sep=""),
				paste("ClusterEnsemble_",file_number,"=clusterensemble(cls)",sep=""),
				paste("csvwrite('ClusterEnsemble_",file_number,".csv', ClusterEnsemble_",file_number,")",sep=""))
		}
#		writeLines(matlabcode, con=paste("EnsembleClusteringB_",file_number,".m",sep=""))	
	}
	
	if(!executable){
		writeLines(matlabcode, con=paste("EnsembleClustering_",file_number,".m",sep=""))
		system(paste("matlab -nodisplay -r \"run('EnsembleClustering_",file_number,".m'); exit\"",sep=""),intern=TRUE)
	}
	Continue=FALSE
	time=0
	while(Continue==FALSE){
		Sys.sleep(15)
		time=time+15
		Continue=file.exists(paste("ClusterEnsemble_",file_number,".csv",sep=""))
		if(time>waitingtime & Continue==FALSE){
			stop(paste("Waited",waitingtime, "seconds for completion of the ensemble clustering procedure. Increase waiting time to continue.",sep=" "))
		}
	}

	ClusterEnsemble <- utils::read.table(paste("ClusterEnsemble_",file_number,".csv",sep=""),sep=",")
	ClusterEnsemble=as.vector(as.matrix(ClusterEnsemble))
	names(ClusterEnsemble)=rownames(List[[1]])
	
	
	DistM=NULL

	
	Clust=list()
	Clust$order=sort(ClusterEnsemble,index=TRUE)$ix
	Clust$order.lab=names(ClusterEnsemble[sort(ClusterEnsemble,index=TRUE)$ix])
	Clust$Clusters=ClusterEnsemble
	
	Out=list("DistM"=DistM,"Clust"=Clust)
	attr(Out,"method")="Ensemble"
	
	file.remove(paste("ClusterEnsemble_",file_number,".csv",sep=""))
	file.remove(paste("matlabdata_",file_number,".csv",sep=""))
	
	if(!executable){
		file.remove(paste("EnsembleClustering_",file_number,".m",sep=""))
	}	

	return(Out)	
}


#' @title Hybrid bipartite graph formulation
#' 
#' @description Hybrid Bipartite Graph Formulation (HBGF) is a graph-based consensus multi-source clustering technique. The method builds a 
#' bipartite graph in which the two types of vertices are represented by the objects on one hand and the clusters of the partitions on the 
#' other hand. An edge is only present between an object vertex and a cluster vertex indicating that the object belongs to that cluster. 
#' The graph can be partitioned with the Spectral clustering \insertCite{Ng2000}{IntClust}.
#'  
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clusters".	
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}. 
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param clust Choice of clustering function (character). Defaults to "agnes". 
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible"
#' @param nrclusters The number of clusters to divide each individual dendrogram in. Default is c(7,7) for two data sets.
#' @param gap Logical. Whether the optimal number of clusters should be determined with the gap statistic. Default is FALSE.
#' @param maxK The maximal number of clusters to investigate in the gap statistic. Default is 15.
#' @param graphPartitioning A character string indicating the preferred graph partitioning algorithm. For now only spectral clustering ("Spec") is implemented. Defaults to "Spec".
#' @param optimalk An estimate of the final optimal number of clusters. Default is 7.
#' @return The returned value is a list of two elements:
#' \item{DistM}{A NULL object}
#' \item{Clust}{The resulting clustering}
#' The value has class 'Ensemble'.
#' @references \insertRef{Fern2004}{IntClust} 
#' @examples
#' data(fingerprintMat)
#' data(targetMat)
#' L=list(fingerprintMat,targetMat)
#' 
#' MCF7_HBGF=HBGF(List=L,type="data",distmeasure=c("tanimoto","tanimoto"),normalize=
#' c(FALSE,FALSE),method=c(NULL,NULL),clust="agnes",linkage = c("flexible",
#' "flexible"),nrclusters=c(7,7),gap = FALSE, maxK = 15,graphPartitioning="Spec",
#' optimalk=7)
HBGF<-function(List,type=c("data","dist","clust"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust = "agnes", linkage = c("flexible","flexible"),alpha=0.625,nrclusters=c(7,7),gap = FALSE, maxK = 15,graphPartitioning="Spec",optimalk=7){
	# Step 1: Generate several clustering results on the objects
	# Is the clustering of the objects on each data set enough (sometimes only 2 data sets) or should we produce multiple for each soource?
	
	if(type=="data"){
		
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type="data",distmeasure[i],normalize[i],method[i],clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))
		
		Dist=lapply(seq(length(List)),function(i) Clusterings[[i]]$DistM)
		
		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=clusters
			}
		}
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
	}
	
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize=FALSE,method=NULL,clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))
		
		Dist=List
		
		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=ceiling(mean(clusters))
			}
		}
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
		
		
	}
	else if(type=="clust"){
		
		Clusterings=List
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		if(is.null(nrclusters)){
			stop("Please specify a number of clusters")
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
		
	}
	#### Step 1 Complete
	
	# Step 2: construct a graph G from the cluster ensemble: based on matlab code from Brodley and Fern (http://web.engr.oregonstate.edu/~xfern/)
	
	ClusterEnsembles=as.data.frame(Clusters)
	rownames(ClusterEnsembles)=rownames(List[[1]])
	for(i in 1:ncol(ClusterEnsembles)){
		colnames(ClusterEnsembles)[i]=paste("Solution_",i,sep="")
	}
	
	
	A=lapply(seq(length(Clusters)),function(i) model.matrix(~as.factor(Clusters[[i]])-1))
	A=Reduce(cbind,A)
	rownames(A)=rownames(List[[1]])
	for(c in 1:ncol(A)){
		colnames(A)[c]=paste("Cluster ",c,sep="")
	}
	
	
	
	#Need W or just the connectivity matrix A? Ferns Algorthym uses A
	if(graphPartitioning=="Spec"){
		D<-diag(sqrt(colSums(A)))
		L<-A%*%solve(D)
		SingularValues=svd(L,nu=optimalk,nv=optimalk)
		U=SingularValues$u
		V=SingularValues$v
		
		U=U/matrix(rep(sqrt(rowSums(U^2)),optimalk),nrow=nrow(A),ncol=optimalk,byrow=FALSE)
		V=V/matrix(rep(sqrt(rowSums(V^2)),optimalk),nrow=ncol(A),ncol=optimalk,byrow=FALSE)
		
		permutated=gtools::permute(1:nrow(A))
		centers=U[permutated[1],,drop=FALSE]
		
		c=rep(0,nrow(A))
		c[permutated[1]]=2*optimalk
		
		#finsih this last for loop
		
		for(j in 2:optimalk){
			
			c=c+abs(U%*%t(centers[j-1,,drop=FALSE]))
			m = which.min(c)
			centers = rbind(centers,U[m,])
			c[m] = 2*optimalk
			
		}
		
		#k-means clustering
		
		clusterid=try(kmeans(x=rbind(U,V),centers=centers,iter.max=200),silent=TRUE)
		if(class(clusterid)=="try-error"){
			clusterid=try(kmeans(x=rbind(U,V),centers=length(centers),iter.max=200),silent=TRUE)
		}
		Clusters=clusterid$cluster[1:nrow(A)]
		names(Clusters)=rownames(A)
		
		clusters=unique(Clusters)
		order=c()
		for(j in clusters){
			order=c(order,which(Clusters==j))
		}
		
		order.lab=as.character(order)
		
	}
	
	Out=list(DistM=NULL,Clust=list(order=order,order.lab=order.lab,Clusters=Clusters))
	attr(Out,"method")="Ensemble"
	return(Out)
	
}


#' @title Clustering aggregation
#' 
#' @description The \code{ClusteringAggregation} includes the ensemble clustering methods Balls, Agglomerative (Aggl.) and Furthest which are graph-based consensus methods.
#'  
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clusters".	
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}. 
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param clust Choice of clustering function (character). Defaults to "agnes". 
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible"
#' @param nrclusters The number of clusters to divide each individual dendrogram in. Default is c(7,7) for two data sets.
#' @param gap Logical. Whether the optimal number of clusters should be determined with the gap statistic. Default is FALSE.
#' @param maxK The maximal number of clusters to investigate in the gap statistic. Default is 15.
#' @param agglMethod The method to be performed: "Balls","Aggl","Furthest" or "LocalSearch".
#' @param improve Logical. If TRUE, a local search is performed to improve the obtained results. Default is TRUE.
#' @param distThresh_B  A distance threshold for the Balls algoritme. Default is 0.5.
#' @param distThresh_A A distance threshold for the Aggl. algoritme. Default is 0.8.
#' @details \insertCite{Gionis2007}{IntClust} propose heuristic algorithms in order to find a solution for the consensus problem. In a first step, a 
#' weighted graph is built from the objects with weights between two vertices determined by the fraction of clusterings that place the two vertices 
#' in different clusters. In a second step, an algorithm searches for the partition that minimizes the total number of disagreements with the given 
#' partitions. The Balls algorithm is an iterative process which finds a cluster for the consensus partition in each iteration. For each object $i$, 
#' all objects at a distance of at most 0.5 are collected and the average distance of this set to the $i$th object of interest is calculated. If the
#'  average distance is less or equal to a parameter $alpha$ the objects form a cluster; otherwise the object forms a singleton. The Agglomerative 
#' (Aggl.) algorithm starts by considering every object as a singleton cluster. Next, the two closest clusters are merged if the average distance 
#' between the clusters is less than 0.5. If there are no two clusters with an average distance smaller than 0.5, the algorithm stops and returns 
#' the created clusters as a solution. The Furthest algorithm starts by placing all objects into a single cluster. In each iteration, the pair of 
#' objects that are the furthest apart are considered as new cluster centers. The remaining objects are appointed to the center that increases the 
#' cost of the partition the least and the new cost is computed. The cost is the sum of the all distances between the obtained partition and the 
#' partitions in the ensemble. The iteration continues until the cost of the new partition is higher than the previous partition. 
#' @return The returned value is a list of two elements:
#' \item{DistM}{A NULL object}
#' \item{Clust}{The resulting clustering}
#' The value has class 'Ensemble'.
#' @references \insertRef{Gionis2007}{IntClust} 
#' @examples
#' data(fingerprintMat)
#' data(targetMat)
#' L=list(fingerprintMat,targetMat)
#' 
#' MCF7_Aggl=ClusteringAggregation(List=L,type="data",distmeasure=c("tanimoto","tanimoto"),
#' normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust="agnes",linkage = c("flexible",
#' "flexible"),alpha=0.625,nrclusters=c(7,7),gap = FALSE, maxK = 15,agglMethod="Aggl",
#' improve=TRUE,distThresh_B=0.5,distThresh_A=0.8)
ClusteringAggregation<-function(List,type=c("data","dist","clust"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust = "agnes", linkage = c("flexible","flexible"),alpha=0.625,nrclusters=c(7,7),gap = FALSE, maxK = 15,agglMethod=c("Balls","Aggl","Furthest","LocalSearch"),improve=TRUE,distThresh_B=0.5,distThresh_A=0.8){
	# Step 1: Generate several clustering results on the objects
	# Is the clustering of the objects on each data set enough (sometimes only 2 data sets) or should we produce multiple for each soource?
	
	if(type=="data"){
		
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type="data",distmeasure[i],normalize[i],method[i],clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))
		
		Dist=lapply(seq(length(List)),function(i) Clusterings[[i]]$DistM)
		
		
		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=clusters
			}
		}
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
	}
	
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize=FALSE,method=NULL,clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))
		
		Dist=List
		
		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=ceiling(mean(clusters))
			}
		}
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
		
		
	}
	else if(type=="clust"){
		
		Clusterings=List
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		if(is.null(nrclusters)){
			stop("Please specify a number of clusters")
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
		
	}
	
	# Step 2: Set up the distance matrix X
	H=list()
	for(i in 1:length(List)){
		Bin=matrix(0,nrow=max(Clusters[[i]]),ncol=length(Clusters[[i]]))
		for(j in 1:nrow(Bin)){
			Bin[j,which(Clusters[[i]]==j)]=1
		}
		H[[i]]=as.data.frame(Bin)	
	}
	T=data.table::rbindlist(H)
	S=t(as.matrix(T))%*%(as.matrix(T))
	S=S/length(List)
	X=1-S
	if(type!="clust"){
		rownames(X)=rownames(List[[1]])
		colnames(X)=rownames(List[[1]])
	}
	else{
		rownames(X)=rownames(List[[1]]$DistM)
		colnames(X)=rownames(List[[1]]$DistM)
	}
	
	# Step 3: The specified algorithm on the distance matrix X
	
	if(agglMethod=="Balls"){
		AggClusters=list()
		SortedVertices=sort(rowSums(X))
		while(length(SortedVertices)>0){
			u=SortedVertices[1]
			B=which(SortedVertices[-c(1)]<=u+0.5)
			if(length(B)>0){
				AvDistBtoU=mean(X[which(rownames(X)==names(u)),which(colnames(X)%in%names(B))])
				
				if(AvDistBtoU<=distThresh_B){
					AggClusters[[length(AggClusters)+1]]=c(names(u),names(B))
					SortedVertices=SortedVertices[-which(names(SortedVertices)%in%c(names(u),names(B)))]
				}
				else{
					AggClusters[[length(AggClusters)+1]]=names(u)
					SortedVertices=SortedVertices[-which(names(SortedVertices)%in%c(names(u)))]
				}
			}
			else{
				AggClusters[[length(AggClusters)+1]]=names(u)
				SortedVertices=SortedVertices[-which(names(SortedVertices)%in%c(names(u)))]
			}
		}
		
		nrclusters=length(AggClusters)
		order.lab=unlist(AggClusters)
		
		if(type!="clust"){
			order=match(order.lab,rownames(List[[1]]))
		}
		else{
			order=match(order.lab,rownames(List[[1]]$DistM))
		}
		
		clusters=sapply(1:nrclusters, function(i) rep(i,length(AggClusters[[i]])))
		Clusters=unlist(clusters)  
		names(Clusters)=order.lab
		
		if(type!="clust"){
			Clusters=Clusters[rownames(List[[1]])]
		}
		else{
			Clusters=Clusters[rownames(List[[1]]$DistM)]
		}
		
		if(improve==TRUE){agglMethod="LocalSearch"}
	}
	
	if(agglMethod=="Aggl"){
		#Continue <- readline("This algorithm might take some time. Do you want to proceed?")  
		#if(Continue=="yes"){
		AvDist<-function(DistMat, pairs,temp_clusters){
			Cluster1=names(temp_clusters)[which(temp_clusters==pairs[1])]
			Cluster2=names(temp_clusters)[which(temp_clusters==pairs[2])]
			AvDist1to2=mean(DistMat[which(rownames(DistMat)%in%Cluster2),which(colnames(DistMat)%in%Cluster1)])
			
			return(AvDist1to2)
			
		}
		
		temp_clusters=seq(1:length(rownames(X)))
		names(temp_clusters)=rownames(X)
		
		Continue=TRUE
#			time=0
#			times=0
		while(Continue==TRUE){
			
#				temp<-proc.time()
			Pairs=combn(unique(temp_clusters),2)
			AvDistances=apply(Pairs,2,AvDist,DistMat=X,temp_clusters=temp_clusters)
			
			if(min(AvDistances)<distThresh_A){
				ChosenPair=Pairs[,which.min(AvDistances)]
				temp_clusters[which(temp_clusters==ChosenPair[2])]=ChosenPair[1]
			}
			else{
				Continue=FALSE				
			}
			
			if(length(unique(temp_clusters))==1){
				Continue=FALSE # all have been put into the same cluster
			}	
			
		}
		
		Clusters=as.numeric(as.factor(temp_clusters))
		names(Clusters)=names(temp_clusters)
		
		clusters=unique(Clusters)
		order=c()
		for(j in clusters){
			order=c(order,which(Clusters==j))
		}
		
		order.lab=as.character(order)
		if(improve==TRUE){agglMethod="LocalSearch"}
		#}
		#else{
		#	stop("Algorithm stopped. No output produced.")
		#}
	}
	
	if(agglMethod=="Furthest"){
		
		clusters=rep(0,length(rownames(X)))
		names(clusters)=rownames(X)
		
		centers=c(rownames(X)[which(X==max(X),arr.ind = TRUE)[1,1]],colnames(X)[which(X==max(X),arr.ind = TRUE)[1,2]])
		for(i in 1:length(centers)){
			clusters[which(names(clusters)%in%centers[i])]=i
		}
		
		for(j in names(clusters)[-c(which(names(clusters)%in%centers))]){
			DistancesToCenters=X[which(rownames(X)==j),which(colnames(X)%in%centers)]
			AssignedCenter=names(DistancesToCenters)[which.min(DistancesToCenters)]
			clusters[which(names(clusters)==j)]=clusters[which(names(clusters)==AssignedCenter)]
			
		}
		
		Part_1=0
		Part_2=0
		
		Pairs=combn(length(rownames(X)),2)
		for(k in 1:ncol(Pairs)){
			Pair=rownames(X)[Pairs[,k]]
			
			if(clusters[which(names(clusters)%in%Pair[1])]==clusters[which(names(clusters)%in%Pair[2])]){
				Part_1=Part_1+X[which(rownames(X)%in%Pair[1]),which(colnames(X)%in%Pair[2])]
			}
			else{
				Part_2=Part_2+(1-X[which(rownames(X)%in%Pair[1]),which(colnames(X)%in%Pair[2])])
			}
			
		}
		
		Cost_clusters=Part_1+Part_2
		Solution=clusters
		
		Continue=TRUE

		while(Continue==TRUE){
		
			AvDistToCenters=c()
			for(i in names(clusters)[-c(which(names(clusters)%in%centers))]){
				DistancesToCenters=X[which(rownames(X)==i),which(colnames(X)%in%centers)]
				AvDistToCenters=c(AvDistToCenters,mean(DistancesToCenters))	
			}
			names(AvDistToCenters)=names(clusters)[-c(which(names(clusters)%in%centers))]
			
			new_center=names(AvDistToCenters)[which.max(AvDistToCenters)]
			
			centers<-c(centers,new_center)
			
			for(i in 1:length(centers)){
				clusters[which(names(clusters)%in%centers[i])]=i
			}
			
			for(j in names(clusters)[-c(which(names(clusters)%in%centers))]){
				DistancesToCenters=X[which(rownames(X)==j),which(colnames(X)%in%centers)]
				AssignedCenter=names(DistancesToCenters)[which.min(DistancesToCenters)]
				clusters[which(names(clusters)==j)]=clusters[which(names(clusters)==AssignedCenter)]		
			}
			
			Part_1=0
			Part_2=0
			
			Pairs=combn(length(rownames(X)),2)
			for(k in 1:ncol(Pairs)){
				Pair=rownames(X)[Pairs[,k]]
				
				if(clusters[which(names(clusters)%in%Pair[1])]==clusters[which(names(clusters)%in%Pair[2])]){
					Part_1=Part_1+X[which(rownames(X)%in%Pair[1]),which(colnames(X)%in%Pair[2])]
				}
				else{
					Part_2=Part_2+(1-X[which(rownames(X)%in%Pair[1]),which(colnames(X)%in%Pair[2])])
				}
				
			}
			
			Cost_clusters_new=Part_1+Part_2
			
			
			if(Cost_clusters_new<Cost_clusters){
				Cost_clusters=Cost_clusters_new
				Solution=clusters
			}
			else{
				Continue=FALSE
				clusters=Solution
			}
			
			
		}
		
		Clusters=clusters
		clusters=unique(Clusters)
		order=c()
		for(j in clusters){
			order=c(order,which(Clusters==j))
		}
		
		order.lab=as.character(order)
		
		if(improve==TRUE){
			agglMethod="LocalSearch"
		}
	}
	
	#print(AgglMethod)
	
	if(agglMethod=="LocalSearch"){
		UpdateClustering=Clusters 
		
		Stay=rep(0,length(UpdateClustering))
		names(Stay)=names(UpdateClustering)
		Moved=c()

		while(length(which(Stay=="Stay"))!=length(UpdateClustering)){  #continue untill all nodes want to stay where they are
		
			for(i in names(UpdateClustering)){   #Pick a node 
				
				PresentCluster=UpdateClustering[which(names(UpdateClustering)==i)]  # Present cluster of the node
				ClusterSizes=table(UpdateClustering)  # Current cluster sizes
				
				Cost_Cj=c()   #Compute the cost of moving node i to any cluster
				for(j in sort(unique(UpdateClustering))){
					Cj=names(UpdateClustering)[which(UpdateClustering==j)]
					M_v_Cj=sum(X[which(rownames(X)==i),which(colnames(X)%in%Cj)])
					Cost_Cj=c(Cost_Cj,M_v_Cj)				
				}
				names(Cost_Cj)=sort(unique(UpdateClustering))
				
				Dist_k=c()  #Compute the cost of moving node i to any cluster or a singleton
				for(k in sort(unique(UpdateClustering))){
					ClusterToMoveTo=k
					Cost_Singl=sum(ClusterSizes[which(names(ClusterSizes)!=ClusterToMoveTo)]-Cost_Cj[which(names(Cost_Cj)!=ClusterToMoveTo)])
					Dist_k=c(Dist_k,Cost_Cj[which(names(ClusterSizes)==ClusterToMoveTo)]+Cost_Singl)
				}	
				PotentialClusters=which(Dist_k==min(Dist_k)) # To which clusters does the move hae minimal cost
				
				
				if(PresentCluster%in%PotentialClusters){
					Stay[which(names(UpdateClustering)==i)]="Stay"	#If present cluster is among the potential ones: stay	
				}
				else{
					UpdateClustering[which(names(UpdateClustering)==i)]=PotentialClusters[1]  #else move and if mover before stay is changed to 0 again
					if(Stay[which(names(UpdateClustering)==i)]=="Stay"){
						Stay[which(names(UpdateClustering)==i)]=0
					}
					Moved=c(Moved,i)
					Moved=unique(Moved)
				}
			}

		}
		
		Clusters=UpdateClustering
		
		clusters=unique(Clusters)
		order=c()
		for(j in clusters){
			order=c(order,which(Clusters==j))
		}
		
		order.lab=as.character(order)
		
	}
	
	Out=list(DistM=NULL,Clust=list(order=order,order.lab=order.lab,Clusters=Clusters))
	attr(Out,"method")="Ensemble"
	return(Out)
}


# Voting-based consensus approaches

#' @title Cumulative voting-based aggregation algorithm
#' 
#' @description The \code{CVAA} includes the ensemble clustering methods CVAA and W-CVAA which are voting-based consensus methods.
#'  
#' @export
#' @param Reference The reference structure to be updated.
#' @param nrclustersR The number of clusters present in the reference structure. Default is 7.
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param typeL indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clusters".	
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}. 
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param clust Choice of clustering function (character). Defaults to "agnes". 
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible"
#' @param nrclusters The number of clusters to divide each individual dendrogram in. Default is c(7,7) for two data sets.
#' @param gap Logical. Whether the optimal number of clusters should be determined with the gap statistic. Defaults to FALSE.
#' @param maxK The maximal number of clusters to investigate in the gap statistic. Default is 15.
#' @param votingMethod The method to be performed: "CVAA","W-CVAA".
#' @param optimalk An estimate of the final optimal number of clusters. Default is nrclustersR.
#' @details \insertCite{Saeed2012}{IntClust} describe the Cumulative Voting-based Aggregation Algorithm (CVAA) and introduce a 
#' variation Weighted Cumulative Voting-based Aggregation Algorithm (W-CVAA, \insertCite{Saeed2014}{IntClust}). In the CVAA algorithm, one data 
#' partitioning is chosen as the reference partition. In a first step each partition is relabelled with respect to the reference partition in 
#' search of an optimal relabelling. In the next step a consensus partition is obtained. The W-CVAA algorithm is similar but appoints weights 
#' to each individual partition. The weights are based on the mutual information of the partition measured by the Shannon entropy. 
#' @return The returned value is a list of two elements:
#' \item{DistM}{A NULL object}
#' \item{Clust}{The resulting clustering}
#' The value has class 'Ensemble'.
#' @references \insertRef{Saeed2012}{IntClust} \insertRef{Saeed2014}{IntClust}  
#' @examples
#' data(fingerprintMat)
#' data(targetMat)
#' L=list(fingerprintMat,targetMat)
#' 
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' MCF7_CVAA=CVAA(Reference=MCF7_T,nrclustersR=7,List=L,typeL="data",
#' distmeasure=c("tanimoto", "tanimoto"),normalize=c(FALSE,FALSE),method=
#' c(NULL,NULL),clust="agnes",linkage = c("flexible","flexible"),alpha=0.625,
#' nrclusters=c(7,7),gap = FALSE, maxK = 15,votingMethod="CVAA",optimalk=7)
CVAA<-function(Reference=NULL,nrclustersR=7,List,typeL=c("data","dist","clust"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust = "agnes", linkage = c("flexible","flexible"),alpha=0.625,nrclusters=c(7,7),gap = FALSE, maxK = 15,votingMethod=c("CVAA","W-CVAA"),optimalk=nrclustersR){
	
	#needs a reference partition: which one to choose? Weighted? Or one of the single source clusterings? Let user decide but suggest Weighted as default
	#Reference can be a "method": List needs to be data or dist
	#if Reference is a clust, list can still be anything
	
	if(typeL=="data"){
		
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type="data",distmeasure[i],normalize[i],method[i],clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))
		
		Dist=lapply(seq(length(List)),function(i) Clusterings[[i]]$DistM)
		
		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=clusters
			}
		}
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
		
		U_all=lapply(seq(length(Clusters)),function(i) model.matrix(~as.factor(Clusters[[i]])-1))
		for(i in 1:length(U_all)){
			colnames(U_all[[i]])=seq(1,ncol(U_all[[i]]))
			rownames(U_all[[i]])=rownames(List[[i]])
		}
		
	}
	
	else if(typeL=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],typeL,distmeasure[i],normalize=FALSE,method=NULL,clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))
		
		Dist=List
	
		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=ceiling(mean(clusters))
			}
		}
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
		U_all=lapply(seq(length(Clusters)),function(i) model.matrix(~as.factor(Clusters[[i]])-1))
		for(i in 1:length(U_all)){
			colnames(U_all[[i]])=seq(1,ncol(U_all[[i]]))
			rownames(U_all[[i]])=rownames(List[[i]])
		}
		
		
	}
	else if(typeL=="clust"){
		
		Clusterings=List
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		if(is.null(nrclusters)){
			stop("Please specify a number of clusters")
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
		U_all=lapply(seq(length(Clusters)),function(i) model.matrix(~as.factor(Clusters[[i]])-1))
		for(i in 1:length(U_all)){
			colnames(U_all[[i]])=seq(1,ncol(U_all[[i]]))
			rownames(U_all[[i]])=rownames(List[[i]]$DistM)
			
		}
		
	}
	
	if(!(is.null(Reference))){
		
		ListNew=list()
		element=0
		
		if(attributes(Reference)$method != "CEC" & attributes(Reference)$method != "Weighted" & attributes(Reference)$method!= "WeightedSim"){
			ResultsClust=list()
			ResultsClust[[1]]=list()
			ResultsClust[[1]][[1]]=Reference
			names(ResultsClust[[1]])[1]="Clust"
			element=element+1					
			ListNew[[element]]=ResultsClust[[1]]
			#attr(ListNew[element],"method")="Weights"
		}
		else if(attributes(Reference)$method=="CEC" | attributes(Reference)$method=="Weighted" | attributes(Reference)$method == "WeightedSim"){
			ResultsClust=list()
			ResultsClust[[1]]=list()
			if(attributes(Reference)$method != "WeightedSim"){
				ResultsClust[[1]][[1]]=Reference$Clust
				names(ResultsClust[[1]])[1]="Clust"				
				element=element+1		
				ListNew[[element]]=ResultsClust[[1]]
				attr(ListNew[[element]]$Clust,"method")="Weights"
				
			}			
			else{
				for (j in 1:length(Reference$Results)){
					ResultsClust[[j]]=list()
					ResultsClust[[j]][[1]]=Reference$Results[[j]]
					names(ResultsClust[[j]])[1]="Clust"
					attr(ResultsClust[[1]],"method")="Weights"
					element=element+1					
					ListNew[[element]]=ResultsClust[[j]]
					ListNew=unlist(ListNew,recursive=FALSE)
				}		
			}		
		}	
		
		Reference=ListNew
		
		if(is.null(nrclustersR)){
			stop("Please specify a number of clusters for the reference partition")
		}
		
		CutTree<-function(Data,nrclusters){
			if(attributes(Data$Clust)$method == "Ensemble"){
				Clusters=Data$Clust$Clust$Clusters
				names(Clusters)=NULL
			}
			else{
				Clusters=cutree(Data$Clust$Clust,k=nrclusters)
			}
			return(Clusters)
		}
		
		Ref_Clusters=CutTree(Reference[[1]],nrclusters=nrclustersR)	
		
		U_0=model.matrix(~as.factor(Ref_Clusters)-1)
		colnames(U_0)=seq(1,ncol(U_0))  #assign reference to U_0
		rownames(U_0)=rownames(Reference[[1]]$Clust$DistM)
		U_Ref=U_0
		
	}
	else{
		stop("No Reference specified")
	}
	
	if(votingMethod=="CVAA"){
		
		if(is.null(Reference)){
			random_part<-sample(1:length(U_all),size=1)
			U_0=U_all[[random_part]]
			U_all=U_all[-random_part]
		}
		
		for(i in 1:length(U_all)){ #for i in 2 to b to:
			
			W_i=solve((t(U_all[[i]])%*%U_all[[i]]))%*%t(U_all[[i]])%*%U_0
			
			V_i=U_all[[i]]%*%W_i
			
			U_0=((i-1)/i)*U_0+(1/i)*V_i			
		}
	}
	
	else if(votingMethod=="W-CVAA"){
		
		#Preparation for the Weights for the partitions
		H_c=c()
		for(i in 1:length(U_all)){
			P_c=(colSums(U_all[[i]]))/nrow(U_all[[i]])
			H_c=c(H_c,-sum(P_c*log(P_c)))
		}
		
		#message("Reordering the clusterings in decreasing order of average amount of information")
		Order=sort(H_c,decreasing=TRUE,index.return=TRUE)$ix
		U_all=U_all[Order]
		
		if(is.null(Reference)){
			#message("Performing the Ada-cVote algorithm by Ayad and Kamel")	
			U_0=U_all[Order[1]]
			U_all=U_all[Order[-c(1)]]
		}
		
		for(i in 1:length(U_all)){ #for i in 2 to b to:
			
			# The weight for partition i
			T_i=(H_c[i])/sum(H_c)
			
			
			W_i=solve((t(U_all[[i]])%*%U_all[[i]]))%*%t(U_all[[i]])%*%U_0
			
			V_i=U_all[[i]]%*%W_i
			
			U_0=((i-1)/i)*U_0+(T_i/i)*V_i			
		}
		
		
	}
	
	U_consensus=U_0	
	
	if(!is.null(optimalk) & optimalk==nrclustersR){
		
		Clusters=apply(U_consensus,1,function(i) which.max(i))
		names(Clusters)=rownames(U_consensus)
		clusters=unique(Clusters)
		order=c()
		for(j in clusters){
			order=c(order,which(Clusters==j))
		}
		
		order.lab=as.character(order)
	}
	
	else if(!is.null(optimalk)){
		
		p_hat_joint_c_x=(1/nrow(U_consensus))*U_consensus
		p_hat_x=1/nrow(U_consensus)
		
		#The JS-Alink Algorithm (applied to get the best partition out of the consensus partion U_consensus).
		#This does not do anything if the optimal number of clusters is equal to the number of clusters of the reference. It will when cluster of the reference are to be merged (less clusters than the reference)
		
		Pairs=combn(ncol(U_consensus),2)
		
		JS=c()
		JS_ALink<-function(Pair,U_Ref,p_hat_joint_c_x){
			p_hat_cl=colSums(U_Ref)[Pair[1]]/nrow(U_Ref)    #the number of objects assigned to cluster l in the reference divided by the total
			p_hat_cm=colSums(U_Ref)[Pair[2]]/nrow(U_Ref)  	#the number of objects assigned to cluster m in the reference divided by the total
			
			p_hat_x_cl=p_hat_joint_c_x[,Pair[1]]/p_hat_cl
			p_hat_x_cm=p_hat_joint_c_x[,Pair[2]]/p_hat_cm
			
			alpha_1=p_hat_cl/(p_hat_cl+p_hat_cm)
			alpha_2=p_hat_cm/(p_hat_cl+p_hat_cm)
			
			
			temp1=alpha_1*p_hat_x_cl+alpha_2*p_hat_x_cm
			temp2=temp1*log(temp1)
			temp2[which(is.na(temp2))]=0
			part1=-sum(temp2)
			
			temp3=p_hat_x_cl*log(p_hat_x_cl)
			temp3[which(is.na(temp3))]=0
			part2=alpha_1*(-sum(temp3))
			
			temp4=p_hat_x_cm*log(p_hat_x_cm)
			temp4[which(is.na(temp4))]=0
			part3=alpha_2*(-sum(temp4))
			
			
			JS_cl_cm=part1-part2-part3
			JS=c(JS,JS_cl_cm)
			
		}
		JS=apply(Pairs,2,function(i) JS_ALink(Pair=i,U_Ref,p_hat_joint_c_x))
		
		Clust=agnes(JS,diss=TRUE,method="average")
		
		if(is.null(optimalk)){
			k_part=cutree(as.hclust(Clust),k=optimalk)
			
			priors=c()
			jointdistr=c()
			for(a in 1:optimalk){
				S_g=which(k_part==a)
				prior=sum(colSums(U_Ref)[S_g]/nrow(U_Ref)) 
				priors=c(priors,prior)
				
				joint_temp=rowSums(p_hat_joint_c_x[,S_g,drop=FALSE]/(colSums(U_Ref)[S_g,drop=FALSE]/nrow(U_Ref)))*colSums(U_Ref)[S_g,drop=FALSE]/nrow(U_Ref)
				jointdistr=cbind(jointdistr,joint_temp)
				
			}
			colnames(jointdistr)=c(1:ncol(jointdistr))
			
			U_hat=jointdistr/p_hat_x
			Clusters=apply(U_consensus,1,function(i) which.max(i))
			names(Clusters)=rownames(U_consensus)
			
			clusters=unique(Clusters)
			order=c()
			for(j in clusters){
				order=c(order,which(Clusters==j))
			}
			
			order.lab=as.character(order)
			
		}
	}
	else if(is.null(optimalk)){
		#message("No optimal number of clusters was specified. The consenus matrix will be returned but an optimal partition was not extracted by the JS-ALink algorithm")
		Out=U_consensus
		return(Out)
		
#		# implemented as interpreted form the text of Ayad (not sure of interpretation)
#		# Relying on the merge component of Clust
		# Not sure what the interpretation is of "lifetime of cluster" : not open for use
#	
#		LifeTime<-function(Merging){
#			lifetimes=c()
#			for(a in 1:nrow(Merging)){
#				Merge=Merging[a,]
#				
#				GetC1C2<-function(Merge,C1=c(),C2=c(),Data){
#					if(sign(Merge)[1]=="-1"){
#						C1=abs(Merge[1])
#					}
#					else{
#						C1=GetC1C2(Merge=Data[Merge[1],],C1=C1,C2=C2,Data=Data[1:Merge[1],])
#						C1=unlist(C1)
#						
#					}
#					if(sign(Merge)[2]=="-1"){
#						C2=abs(Merge[2])
#					}
#					else{
#						C2=GetC1C2(Merge=Data[Merge[2],],C1=C1,C2=C2,Data=Data[1:Merge[2],])
#						C2=unlist(C2)
#					}
#					
#					C1C2=list(C1,C2)
#					return(C1C2)
#					
#				}
#                C1C2=GetC1C2(Merge,Data=Merging[1:a,])    
#				C1=C1C2[[1]]
#				C2=C1C2[[2]]
#				
#				MergedPairs=expand.grid(C1,C2)
#				distances=c()
#				for(b in 1:nrow(MergedPairs)){
#					WhichPair=which(Pairs[1,]==MergedPairs[b,1]&Pairs[2,]==MergedPairs[b,2])
#					distances=c(distances,JS[WhichPair])
#				}
#				lifetimes=c(lifetimes,mean(distances))
#			}
#			
#			return(lifetimes)
#	
#		}
#		LifeTimes<-LifeTime(Merging=Clust$merge)
#		
#		Optimalk=which.max(lifetimes)
#		
		
	}
	
	Out=list(DistM=NULL,Clust=list(order=order,order.lab=order.lab,Clusters=Clusters))
	attr(Out,"method")="Ensemble"
	return(Out)
	
}


#' @title Consensus clustering
#' 
#' @description The \code{ConsensusClustering} includes the ensemble clustering methods IVC, IPVC and IVC which are voting-based consensus methods.
#'  
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clusters".	
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}. 
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param clust Choice of clustering function (character). Defaults to "agnes". 
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible"
#' @param nrclusters The number of clusters to divide each individual dendrogram in. Default is c(7,7) for two data sets.
#' @param gap Logical. Whether the optimal number of clusters should be determined with the gap statistic. Defaults to FALSE.
#' @param maxK The maximal number of clusters to investigate in the gap statistic. Default is 15.
#' @param votingMethod The method to be performed: "IVC", "IPVC,"IVC".
#' @param optimalk An estimate of the final optimal number of clusters. Default is 7.
#' @details \insertCite{Nguyen2007}{IntClust} propose three EM-like consensus clustering algorithms: Iterative Voting Consensus (IVC), 
#' Iterative Probabilistic Voting Consensus (IPVC) and Iterative Pairwise Consensus (IPC). Given a number of clusters $k$, the methods 
#' iteratively compute the cluster centers and reassign each object to the closest center. IVC and IPVC represent the cluster centers 
#' by a vector of the majority votes of the cluster labels of all points belonging to the cluster in each partition. For the reassignment, 
#' IVC uses the Hamming distance to compute the distance between the data points and the cluster centers. IPVC is a 
#' refinement of IVC as the distance function takes into account the proportion that each feature of a point differs from the points in 
#' the cluster. The IPC algorithms is slightly different since the original clusters are built from a similarity matrix which represents 
#' the ratio of the number of partitions in which two objects reside in the same cluster. The distance between a data point and a cluster 
#' center is the average of the similarity values between the data point and the points residing in the cluster. The iteration ends when the
#' consensus partition does not change.
#' @return The returned value is a list of two elements:
#' \item{DistM}{A NULL object}
#' \item{Clust}{The resulting clustering}
#' The value has class 'Ensemble'.
#' @references \insertRef{Nguyen2007}{IntClust}
#' @examples
#' data(fingerprintMat)
#' data(targetMat)
#' L=list(fingerprintMat,targetMat)
#' 
#' MCF7_IVC=ConsensusClustering(List=L,type="data",distmeasure=c("tanimoto", "tanimoto"),
#' normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust="agnes",linkage = c("flexible",
#' "flexible"),alpha=0.625,nrclusters=c(7,7),gap = FALSE, maxK = 15,
#' votingMethod="IVC",optimalk=7)
ConsensusClustering<-function(List,type=c("data","dist","clust"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust = "agnes", linkage = c("flexible","flexible"),alpha=0.625,nrclusters=c(7,7),gap = FALSE, maxK = 15,votingMethod=c("IVC","IPVC","IPC"),optimalk=7){
	
	if(type=="data"){
		
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,,drop=FALSE]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type="data",distmeasure[i],normalize[i],method[i],clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))
		
		Dist=lapply(seq(length(List)),function(i) Clusterings[[i]]$DistM)
		
		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=clusters
			}
		}
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
		
	}
	
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize=FALSE,method=NULL,clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))
		
		Dist=List
	
		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=ceiling(mean(clusters))
			}
		}
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
		
		
	}
	else if(type=="clust"){
		
		Clusterings=List
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		if(is.null(nrclusters)){
			stop("Please specify a number of clusters")
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
		
	}
	
	
	Y=Reduce(cbind,Clusters)
	rownames(Y)=rownames(List[[1]])	
	colnames(Y)=seq(1,ncol(Y))
	
	
	if(votingMethod=="IVC"){
		target=sample(optimalk,nrow(Y),replace=TRUE)
		update=rep(0,nrow(Y))
		Continue=TRUE
		while(Continue){
			
			#representation of each cluster
			P=list()
			y_P=list()
			for(i in unique(target)){
				P_i=which(target==i)	
				y_Pi=apply(Y,2,function(k) as.numeric(names(which.max(table(k[P_i])))))
				P[[i]]=P_i
				y_P[[i]]=y_Pi
				
			}
			
			#update y
			for(y in 1:length(target)){
				distances=c()
				for(a in 1:length(y_P)){
					distances=c(distances,sum(Y[y,] != y_P[[a]])	)
				}
				
				update[y]=which.min(distances)
				
			}
			if(all(update==target)){			
				Continue=FALSE
				#print("here")
			}
			else{
				target=update
			}	
			
		}
		
	}
	else if(votingMethod=="IPVC"){
		target=sample(optimalk,nrow(Y),replace=TRUE)
		update=rep(0,nrow(Y))
		Continue=TRUE
		
		while(Continue){
			
			#representation of each cluster
			P=list()
			n_P=list()
			for(i in unique(target)){
				P_i=which(target==i)	
				P[[i]]=P_i
				n_P[[i]]=length(P_i)
				
			}
			
			#update y
			for(y in 1:length(target)){
				distances=c()
				
				for(a in 1:length(P)){	
					dist=0
					for(b in 1:ncol(Y)){
						dist=dist+sum(Y[y,b] != Y[P[[a]],b])/n_P[[a]]	
					}
					distances=c(distances,dist)
				}
				
				update[y]=which.min(distances)
				
			}
			if(all(update==target)){			
				Continue=FALSE
				#print("here")
			}
			else{
				target=update
			}	
			
		}
		
		
		
	}
	else if(votingMethod=="IPC"){
		
		
		#similarity matrix
		H=list()
		for(i in 1:length(List)){
			Bin=matrix(0,nrow=max(Clusters[[i]]),ncol=length(Clusters[[i]]))
			for(j in 1:nrow(Bin)){
				Bin[j,which(Clusters[[i]]==j)]=1
			}
			H[[i]]=as.data.frame(Bin)	
		}
		T=data.table::rbindlist(H)
		S=t(as.matrix(T))%*%(as.matrix(T))
		S=S/length(List)
		
		if(type!="clust"){
			rownames(S)=rownames(List[[1]])
			colnames(S)=rownames(List[[1]])
		}
		else{
			rownames(S)=rownames(List[[1]]$DistM)
			colnames(S)=rownames(List[[1]]$DistM)
		}
		
		target=sample(optimalk,nrow(Y),replace=TRUE)
		update=rep(0,nrow(Y))
		Continue=TRUE
		
		while(Continue){
			
			#representation of each cluster
			P=list()
			n_P=list()
			for(i in unique(target)){
				P_i=which(target==i)	
				P[[i]]=P_i
				n_P[[i]]=length(P_i)
				
			}
			
			#update x
			for(x in 1:length(target)){
				similarities=c()
				
				for(a in 1:length(P)){	
					sim=sum(S[x,P[[a]]])/n_P[[a]]	
					similarities=c(similarities,sim)
				}
				
				update[x]=which.max(similarities)
			}
			
			
			if(all(update==target)){			
				Continue=FALSE
				#print("here")
			}
			else{
				target=update
			}	
			
		}
		
	}
	
	
	
	Clusters=target
	
	if(type!="clust"){
		names(Clusters)=rownames(List[[1]])
	}
	else{
		names(Clusters)=rownames(List[[1]]$DistM)
	}
	
	order=sort(Clusters)
	order.lab=names(order)
	Out=list(DistM=NULL,Clust=list(order=order,order.lab=order.lab,Clusters=Clusters))
	attr(Out,"method")="Ensemble"
	return(Out)
	
}	

#' @title Evidence accumulation
#' 
#' @description The Evidence Accumulation (EA, \insertCite{Fred2005}{IntClust}) sets up the same similarity matrix and 
#' applies a minimum spanning tree (MST,\insertCite{Prim1957}{IntClust}) algorithm to the corresponding graph. The algorithm will 
#' connect all objects with the shortest past. In order to recover the clusters, weak edges are cut at a threshold value $t$.
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or"clusters".	
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}. 
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param clust Choice of clustering function (character). Defaults to "agnes". 
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible"
#' @param nrclusters The number of clusters to divide each individual dendrogram in. Default is c(7,7) for two data sets.
#' @param gap Logical. Whether the optimal number of clusters should be determined with the gap statistic. Defaults to FALSE.
#' @param maxK The maximal number of clusters to investigate in the gap statistic. Default is 15.
#' @param graphPartitioning The method to be performed: "MTS", "SL", "SL_agnes".  
#' @param t A threshold to cut weak edges. Default is NULL.
#' @return The returned value is a list of two elements:
#' \item{DistM}{A NULL object}
#' \item{Clust}{The resulting clustering}
#' The value has class 'Ensemble'.
#' @references \insertRef{Fred2005}{IntClust} \insertRef{Prim1957}{IntClust}
#' @examples
#' data(fingerprintMat)
#' data(targetMat)
#' L=list(fingerprintMat,targetMat)
#' 
#' MCF7_EA=EvidenceAccumulation(List=L,type="data",distmeasure=c("tanimoto", "tanimoto"),
#' normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust="agnes",linkage = c("flexible",
#' "flexible"),alpha=0.625,nrclusters=c(7,7),gap = FALSE, maxK = 15,graphPartitioning="MTS")
EvidenceAccumulation<-function(List,type=c("data","dist","clust"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust = "agnes", linkage = c("flexible","flexible"),alpha=0.625,nrclusters=c(7,7),gap = FALSE, maxK = 15,graphPartitioning=c("MTS","SL","SL_agnes"),t=NULL){
	
	#needs a reference partition: which one to choose? Weighted? Or one of the single source clusterings? Let user decide but suggest Weighted as default
	#Reference can be a "method": List needs to be data or dist
	#if Reference is a clust, list can still be anything
	
	if(type=="data"){
		
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,,drop=FALSE]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type="data",distmeasure[i],normalize[i],method[i],clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))
		
		Dist=lapply(seq(length(List)),function(i) Clusterings[[i]]$DistM)
		
		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=clusters
			}
		}
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
		
	}
	
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize=FALSE,method=NULL,clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))
		
		Dist=List
		
		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=ceiling(mean(clusters))
			}
		}
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
		
		
	}
	else if(type=="clust"){
		
		Clusterings=List
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		if(is.null(nrclusters)){
			stop("Please specify a number of clusters")
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
		
	}
	
	#similarity matrix
	H=list()
	for(i in 1:length(List)){
		Bin=matrix(0,nrow=max(Clusters[[i]]),ncol=length(Clusters[[i]]))
		for(j in 1:nrow(Bin)){
			Bin[j,which(Clusters[[i]]==j)]=1
		}
		H[[i]]=as.data.frame(Bin)	
	}
	T=data.table::rbindlist(H)
	S=t(as.matrix(T))%*%(as.matrix(T))
	S=S/length(List)
	
	if(type!="clust"){
		rownames(S)=rownames(List[[1]])
		colnames(S)=rownames(List[[1]])
	}
	else{
		rownames(S)=rownames(List[[1]]$DistM)
		colnames(S)=rownames(List[[1]]$DistM)
	}
	
	
	if(graphPartitioning=="MST"){
		
		#are assigning to the first encounter to break ties and get clusters: if changed  to join all if one is in common, all end up in 1 cluster
		
		Graph=igraph::graph.adjacency(adjmatrix=S, mode=c( "undirected"), weighted=TRUE, diag=TRUE,add.colnames=TRUE, add.rownames=NA)
		MST_Graph=igraph::as_data_frame(igraph::minimum.spanning.tree(Graph))
		
		if(!is.null(t)){
			MST_Graph=MST_Graph[-c(which(as.numeric(MST_Graph[,3])<t)),]
		}
		
		Partition=list()
		Partition[[1]]=c(MST_Graph[1,1],MST_Graph[1,2])
		Placed=rep(FALSE,nrow(S))
		Placed[Partition[[1]]]=TRUE
		for(j in 2:nrow(MST_Graph)){
			Edge=as.numeric(MST_Graph[j,c(1:2)])		
			k=1
			
			while(any(Placed[Edge]==FALSE)){
				
				if(k>length(Partition)){
					Partition[[length(Partition)+1]]=Edge
					Placed[Edge]=TRUE
				}	
				
				else if(Edge[1]%in%Partition[[k]] | Edge[2]%in%Partition[[k]]){
					Partition[[k]]=c(unlist(Partition[[k]]),Edge)
					Partition[[k]]=unique(Partition[[k]])
					Placed[Edge]=TRUE
				}
				k=k+1
			}	
			
			
		}	
		
		order=unlist(Partition)
		order.lab=rownames(S)[order]
		t1<-sapply(1:length(Partition),function(i) rep(i,length(Partition[[i]])))
		t2<-sapply(1:length(Partition),function(i) rownames(S)[Partition[[i]]])
		t3<-cbind(unlist(t1),unlist(t2))
		clus=as.numeric(t3[,1])
		names(clus)=t3[,2]
		Clusters=clus[rownames(S)]
		if(any(is.na(Clusters))){
			Clusters[which(is.na(Clusters))]=c(1:length(which(is.na(Clusters))))
			names(Clusters)=rownames(S)
		}
		
		
	}	
	
	else if(graphPartitioning=="SL"){
		if(is.null(t)){
			stop("Specify a treshold t for the Sl algorithm")
		}
		#better to let run for a number of thresholds t (default is 0.5)
		#Not discerning enough if only 2 datasets with 1 clustering each: need more to get more clusters
		
		Clusters=list()
		for(a in 1:nrow(S)){
			if(a==1){
				Clusters[[length(Clusters)+1]]=a
			}
			else{
				for(b in 1:(a-1)){
					
					if(S[a,b]>t){
						
						if(a%in%unlist(Clusters) & b%in%unlist(Clusters)){
							
							Clus1=which(sapply(1:length(Clusters),function(i) a%in%Clusters[[i]]))
							Clus2=which(sapply(1:length(Clusters),function(i) b%in%Clusters[[i]]))
							if(Clus1 != Clus2){
								Clusters[[Clus1]]=unique(c(unlist(Clusters[[Clus1]]),unlist(Clusters[[Clus2]])))
								Clusters[[Clus2]]=NULL
							}	
							
						}	
						else if(a%in%unlist(Clusters) & !(b%in%unlist(Clusters))){
							Clus1=which(sapply(1:length(Clusters),function(i) a%in%Clusters[[i]]))
							Clusters[[Clus1]]=unique(c(unlist(Clusters[[Clus1]]),b))
							
						}
						else if(!(a%in%unlist(Clusters)) & b%in%unlist(Clusters)){
							Clus1=which(sapply(1:length(Clusters),function(i) b%in%Clusters[[i]]))
							if(a==56){
								#print("here")
							}
							Clusters[[Clus1]]=unique(c(unlist(Clusters[[Clus1]]),a))
						}
						else{
							Clusters[[length(Clusters)+1]]=unique(c(a,b))
						}
						
					}
					
				}
			}
		}
		
		Singletons=which(!seq(1,nrow(S))%in%unlist(Clusters))
		if(length(Singletons)>0){
			
			for(l in Singletons){
				Clusters[[length(Clusters)+1]]=l
			}
			
		}
		
		
		order=unlist(Clusters)
		order.lab=rownames(S)[order]
		t1<-sapply(1:length(Clusters),function(i) rep(i,length(Clusters[[i]])))
		t2<-sapply(1:length(Clusters),function(i) rownames(S)[Clusters[[i]]])
		t3<-cbind(unlist(t1),unlist(t2))
		clus=as.numeric(t3[,1])
		names(clus)=t3[,2]
		Clusters=clus[rownames(S)]
		
		#Or perform single linkage on the similarity matrix S
		#Clusters=agnes(1-S,diss=TRUE,method="single")
		
	}
	
	else if(graphPartitioning=="SL_agnes"){
		Clust=agnes(1-S,diss=TRUE,method="single")
		Out=list(DistM=NULL,Clust=list(order=order,order.lab=order.lab,Clusters=Clusters))
		attr(Out,"method")="Single Clustering"
		return(Out)
		
	}		
	
	Out=list(DistM=NULL,Clust=list(order=order,Clusters=Clusters))
	attr(Out,"method")="Ensemble"
	return(Out)
	
}	

#' @title Link based clustering
#' 
#' @description The \code{LinkBasedClustering} includes the ensemble clustering methods cts, srs and asrs which are voting-based consensus methods.
#'  
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or"clusters".	
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}. 
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param clust Choice of clustering function (character). Defaults to "agnes". 
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible"
#' @param nrclusters The number of clusters to divide each individual dendrogram in. Default is c(7,7) for two data sets.
#' @param gap Logical. Whether the optimal number of clusters should be determined with the gap statistic. Defaults to FALSE.
#' @param maxK The maximal number of clusters to investigate in the gap statistic. Default is 15.
#' @param linkBasedMethod The method to be performed: "cts", "srs", "asrs". 
#' @param decayfactor The decay factor to be specified for the methods. Defaults to 0.8.
#' @param niter The number of iterations. Default is 5.
#' @param linkBasedLinkage The linkage to be used in the final clustering. Default is "ward".
#' @param waitingtime The time in seconds to wait until the MATLAB results are generated. Defaults to 300.
#' @param file_number The specific file number to be placed as a tag in the file generated by MATLAB. Defaults to 00.
#' @param executable Logical. Whether the MATLAB functions are performed via an executable on the command line (TRUE, only possible for Linux systems) or by calling on MATLAB directly (FALSE). Defaults to FALSE. The files LinkBasedClusteringcts.m (cts), LinkBasedClusteringsrs.m (srs), LinkBasedClusteringasrs.m (asrs) and and MetisAlgorithm.m are present in the inst folder to be transformed in executables.
#' @details  \insertCite{Iam-on2010}{IntClust} describe three methods for link-based clustering based on a co-association matrix: 
#' Connected-Triple Based Similarity (CTS), SimRank Based Similarity (SRS) and approximate SimRank-based similarity (ASRS). 
#' The methods compute a similarity matrix based on additional information. CTS incorporates information regarding the shared third 
#' link between two objects. SRS works based on the assumption that neighbours are similar if their neighbours are similar as well.
#' The ASRS is introduced as a more efficient implementation of SRS. The \code{R} code calls on the MATLAB code provided by \insertCite{Iam-on2010}{IntClust}. The MATLAB functions are included in the inst folder and should be located in the working directory. Shell script for the executable can be found in the inst folder as well.
#' @return The returned value is a list of two elements:
#' \item{DistM}{The resulting distance matrix}
#' \item{Clust}{The resulting clustering}
#' The value has class 'LinkBased'.
#' @references \insertRef{Iam-on2010}{IntClust}
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' L=list(fingerprintMat,targetMat)
#' 
#' MCF7_cts=LinkBasedClustering(List=L,type="data",distmeasure=c("tanimoto", "tanimoto")
#' ,normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust="agnes",linkage = c("flexible",
#' "flexible"),alpha=0.625,nrclusters=c(7,7),gap = FALSE, maxK = 15,linkBasedMethod="cts",
#' decayfactor=0.8,niter=5,linkBasedLinkage="ward",waitingtime=300,file_number=00)
#' }
LinkBasedClustering<-function(List,type=c("data","dist","clust"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust = "agnes", linkage = c("flexible","flexible"),alpha=0.625,nrclusters=c(7,7),gap = FALSE, maxK = 15,linkBasedMethod=c("cts","srs","asrs"),decayfactor=0.8,niter=5,linkBasedLinkage="ward",waitingtime=300,file_number=00,executable=FALSE){
	# Step 1: Generate several clustering results on the objects
	# Is the clustering of the objects on each data set enough (sometimes only 2 data sets) or should we produce multiple for each soource?
	
	if(type=="data"){
		
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type="data",distmeasure[i],normalize[i],method[i],clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))
		
		Dist=lapply(seq(length(List)),function(i) Clusterings[[i]]$DistM)
		
		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=clusters
			}
		}
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
	}
	
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize=FALSE,method=NULL,clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))
		
		Dist=List

		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=clusters
			}
		}
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
		
		
	}
	else if(type=="clust"){
		
		Clusterings=List
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		if(is.null(nrclusters)){
			stop("Please specify a number of clusters")
		}
		
		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
		
	}
	#### Step 1 Complete
	
	# Step 2: Get data in format to be given to matlab =>  matrix with a row per clustering result
	if(type=="data" | type=="dist"){
		nc=nrow(List[[1]])
	}
	else if(type=="clust"){
		nc=nrow(List[[1]]$DistM)
	}
	matlabdata=t(matrix(unlist(Clusters),ncol=nc,byrow=TRUE))
	utils::write.table(matlabdata,file=paste("matlabdata_",file_number,".csv",sep=""),sep=",",col.names=FALSE,row.names=FALSE)
	
	if(linkBasedMethod=="cts"){
		if(executable){
			system(paste("./LinkBasedClusteringcts ",file_number,sep=""),intern=TRUE)
		}
		else{
			matlabcode=c(
				paste("M = csvread('matlabdata_",file_number,".csv')",sep=""),
				paste("S_",file_number,"=cts(M,",decayfactor,")",sep=""),					
				paste("csvwrite('S_",file_number,".csv', S_",file_number,")",sep=""))
		}
		
	}	
	if(linkBasedMethod=="srs"){
		if(executable){
			system(paste("./LinkBasedClusteringsrs ",file_number,sep=""),intern=TRUE)
		}
		else{
		matlabcode=c(
				paste("M = csvread('matlabdata_",file_number,".csv')",sep=""),
				paste("S_",file_number,"=srs(M,",decayfactor,",",niter,")",sep=""),					
				paste("csvwrite('S_",file_number,".csv', S_",file_number,")",sep=""))
		}
		
	}	
	if(linkBasedMethod=="asrs"){
		if(executable){
			system(paste("./LinkBasedClusteringasrs ",file_number,sep=""),intern=TRUE)
		}
		else{
		matlabcode=c(
				paste("M = csvread('matlabdata_",file_number,".csv')",sep=""),
				paste("S_",file_number,"=asrs(M,",decayfactor,")",sep=""),					
				paste("csvwrite('S_",file_number,".csv', S_",file_number,")",sep=""))
		}
		
	}
	
	if(!executable){
		writeLines(matlabcode, con=paste("LinkBasedClustering_",file_number,".m",sep=""))
		system(paste("matlab -nodisplay -r \"run('LinkBasedClustering_",file_number,".m'); exit\"",sep=""),intern=TRUE)
	}
	Continue=FALSE
	time=0
	while(Continue==FALSE){
		Sys.sleep(15)
		time=time+15
		Continue=file.exists(paste("S_",file_number,".csv",sep=""))
		if(time>waitingtime & Continue==FALSE){
			stop(paste("Waited",waitingtime, "seconds for completion of the ensemble clustering procedure. Increase waiting time to continue.",sep=" "))
		}
	}
	
	
	SimM=utils::read.table(paste("S_",file_number,".csv",sep=""),sep=",")
	DistM=as.matrix(1-SimM)
	rownames(DistM)=rownames(Clusterings[[1]]$DistM)
	colnames(DistM)=rownames(Clusterings[[1]]$DistM)

	LinkBasedCluster=cluster::agnes(DistM,diss=TRUE,method=linkBasedLinkage,par.method=alpha)
	
	Out=list("DistM"=DistM,"Clust"=LinkBasedCluster)
	attr(Out,"method")="LinkBased"
	
	
	file.remove(paste("matlabdata_",file_number,".csv",sep=""))
	file.remove(paste("S_",file_number,".csv",sep=""))
	if(!executable){
		file.remove(paste("LinkBasedClustering_",file_number,".m",sep=""))
	}
	return(Out)	
}


#' @title Complementary ensemble clustering
#' 
#' @description Complementary Ensemble Clustering (CEC) Complementary Ensemble Clustering (CEC, \cite{Fodeh2013}) shows similarities with ADEC. 
#' However, instead of merging the data matrices, ensemble clustering is performedon each data matrix separately. The resulting incidence 
#' matrices for each data sets are combined in a weighted linear equation. The weighted incidence matrix is the input for the final clustering 
#' algorithm. Similarly as ADEC, there are versions depending of the specification of the number of features to sample and the number of clusters.
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}. 
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param t The number of iterations. Defaults to 10.
#' @param r A vector with the number of features to take for the random sample for each element in List. If NULL (default), all features are considered.
#' @param nrclusters A list with a sequence of numbers of clusters to cut the dendrogram in for each element in List. If NULL (default), the function stops.
#' @param weight The weights for the weighted linear combination.
#' @param clust Choice of clustering function (character). Defaults to "agnes". 
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible"
#' @param weightclust A weight for which the result will be put aside of the other results. This was done for comparative reason and easy access.
#' @return The returned value is a list of four elements:
#' \item{DistM}{The resulting incidence matrix}
#' \item{Results}{The hierarchical clustering result for each element in WeightedDist}
#' \item{Clust}{The result for the weight specified in Clustweight}
#' The value has class 'CEC'.
#' @details If r is specified and nrclusters is a fixed number, only a random sampling of the features will be performed for the t iterations (CECa). If r is NULL and the nrclusters is a sequence, the clustering is performedon all features and the dendrogam is divided into clusters for the values of nrclusters (CECb). If both r is specified and nrclusters is a sequence, the combination is performed (CECc).
#' After every iteration, either be random sampling, multiple divisions of the dendrogram or both, an incidence matrix is set up. All incidence matrices are summed and represent the distance matrix on which a final clustering is performed. 
#' @references 
#' \insertRef{Fodeh2013}{IntClust}
#' @examples
#' data(fingerprintMat)
#' data(targetMat)
#' L=list(fingerprintMat,targetMat)
#' 
#' MCF7_CEC=CEC(List=L,distmeasure=c("tanimoto","tanimoto"),normalize=FALSE,method=NULL
#' ,t=100, r=c(100,100), nrclusters=list(seq(2,10,1),seq(2,10,1)),clust="agnes",linkage=
#' c("flexible","flexible"),alpha=0.625,weightclust=0.5)
CEC<-function(List,distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),t=10,r=NULL,nrclusters=NULL,weight=NULL,clust="agnes",linkage=c("flexible","flexible"),alpha=0.625,weightclust=0.5){
	
	if(class(List) != "list"){
		stop("Data must be of type list")
	}
	
	if(is.null(nrclusters)){
		stop("Give a number of clusters to cut the dendrogram into for each data modality.")
	}
	
	if(!is.null(r)&length(nrclusters)==1){
		message("Performing a random sampling of the features with a fixed number of clusters.")
	}else if(is.null(r)&length(nrclusters)>1){
		message("Dividing the dendrogram in k clusters for a range of values of k.")
		t=1
	}else if(!is.null(r)&length(nrclusters)>1){
		message("Performing a random sampling of the features and dividing the dendrogram in k clusters for a range of values of k.")
	}
	else{
		stop("Specify r and/or nrclusters in order to perform an ADEC method.")
	}
	
	#Put all data in the same order
	OrderNames=rownames(List[[1]])
	for(i in 1:length(List)){
		List[[i]]=List[[i]][OrderNames,]
	}
	
	#Put up Incidence matrix for each data modality
	nc=c()
	Incidence=list()
	for (i in 1:length(List)){
		Incidence[[i]]=matrix(0,dim(List[[i]])[1],dim(List[[i]])[1])
		rownames(Incidence[[i]])=rownames(List[[i]])
		colnames(Incidence[[i]])=rownames(List[[i]])
		nc=c(nc,ncol(List[[i]]))
	}
	evenn=function(x){if(x%%2!=0)x=x-1 else x}
	nc=lapply(nc,FUN=evenn)
	nc=unlist(nc)
	
	#Repeat for t iterations
	for(g in 1:t){
		#message(g)		
		if(is.null(r)){
			r=unlist(sapply(List,ncol))
		}
		
		#take random sample:
		A_prime=list()
		for(i in 1:length(r)){
			A=List[[i]]
			temp=sample(ncol(A),r[i],replace=FALSE)
			A_prime[[i]]=A[,temp]
			
			Ok=FALSE
			while(Ok==FALSE){
				if(any(rowSums(A_prime[[i]])==0)){
					temp=sample(ncol(A),r[i],replace=FALSE)
					A_prime[[i]]=A[,temp]
				}
				else{
					Ok=TRUE
				}				
			}			
		}
		
		#protect against zero rows:
		
		
		#Step 2: apply hierarchical clustering on each + cut tree into nrclusters
		
		DistM=lapply(seq(length(A_prime)),function(i) Distance(A_prime[[i]],distmeasure=distmeasure[i],normalize[i],method[i]))
		
		
		HClust_A_prime=lapply(seq(length(DistM)),function(i) cluster::agnes(DistM[[i]],diss=TRUE,method=linkage[i],par.method=alpha))
		
		
		for(o in seq(length(HClust_A_prime))){	
			for(k in 1:length(nrclusters)){
				Temp=stats::cutree(HClust_A_prime[[o]],nrclusters[[o]][k])	
				MembersofClust=matrix(1,dim(List[[o]])[1],dim(List[[o]])[1])
					
				for(l in 1:length(Temp)){
					label=Temp[l]
					sameclust=which(Temp==label)
					MembersofClust[l,sameclust]=0
				}
				Incidence[[o]]=Incidence[[o]]+MembersofClust	
			}
		}
		
	
		
		
	}
	
	if(is.null(weight)){
		equalweights=1/length(List)
		weight=list(rep(equalweights,length(List)))	
	}
	else if(class(weight)=='list' & length(weight[[1]])!=length(List)){
		stop("Give a weight for each data matrix or specify a sequence of weights")
	}
	else{
		#message('The weights are considered to be a sequence, each situation is investigated')
	}
	
	if(class(weight)!="list"){
		condition<-function(l){		
			l=as.numeric(l)
			if( sum(l)==1 ){  #working with characters since with the numeric values of comb or permutations something goes not the way is should: 0.999999999<0.7+0.3<1??
				#return(row.match(l,t1))
				return(l)
			}
			else(return(0))
		}
		
		if(all(seq(1,0,-0.1)!=weight)){
			for(i in 1:length(weight)){
				rest=1-weight[i]
				if(!(rest%in%weight)){
					weight=c(weight,rest)
				}
			}
		}
		
		
		t1=gtools::permutations(n=length(weight),r=length(List),v=as.character(weight),repeats.allowed = TRUE)
		t2=lapply(seq_len(nrow(t1)), function(i) if(sum(as.numeric(t1[i,]))==1) return(as.numeric(t1[i,])) else return(0)) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
		t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
		t4=t2[which(t3!=0)]
		weight=lapply(seq(length(t4)),function(i) rev(t4[[i]]))
		
	}
	
	if(class(weight)=="list" & "x" %in% weight[[1]]){ #x indicates a free weight
		for(i in 1:length(weight)){
			w=weight[[i]]
			weightsfordata=which(w!="x") #position of the provided weight = position of the data to which the weight is given
			givenweights=as.numeric(w[weightsfordata])
			
			stilltodistribute=1-sum(givenweights)
			
			newweights=seq(stilltodistribute,0,-0.1)
			
			t1=gtools::permutations(n=length(newweights),r=length(List)-length(weightsfordata),v=as.character(newweights),repeats.allowed = TRUE)
			Input1=as.list(seq_len(nrow(t1)))
			Input2=lapply(seq(length(Input1)),function(i) {Input1[[i]][length(Input1[[i]])+1]=stilltodistribute
						return(Input1[[i]])})
			t2=lapply(seq(length(Input2)), FUN=function(i){if(sum(as.numeric(t1[Input2[[i]][1],])+0.00000000000000002775)==Input2[[i]][2]) return(as.numeric(t1[i,])) else return(0)}) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
			t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
			weightsforotherdata=t2[which(t3!=0)]
			
			new=list()
			for(i in 1:length(weightsforotherdata)){
				w1=weightsforotherdata[[i]]
				new[[i]]=rep(0,length(List))
				new[[i]][weightsfordata]=givenweights
				new[[i]][which(new[[i]]==0)]=w1
			}
			
			weight=new
		}
	}
	
	weightedcomb<-function(w,Dist){
		temp=lapply(seq_len(length(Dist)),function(i) w[i]*Dist[[i]])
		temp=Reduce("+",temp)	
		return(temp)
	}
	
	IncidenceComb=lapply(weight,weightedcomb,Incidence)
	namesweights=c()	
	CEC=list()
	for (i in 1:length(IncidenceComb)){
		CEC[[i]]=cluster::agnes(IncidenceComb[[i]],diss=TRUE,method="ward")
		namesweights=c(namesweights,paste("Weight",weight[i],sep=" "))
		if(all(weight[[i]]==weightclust)){
			Clust=CEC[i]
			DistClust=IncidenceComb[i]
		}
	}
	
	Results=lapply(seq(1,length(CEC)),function(i) return(c("DistM"=IncidenceComb[i],"Clust"=CEC[i])))
	names(Results)=namesweights
	
	out=list(Incidence=Incidence,Results=Results,Clust=c("DistM"=DistClust,"Clust"=Clust))
	attr(out,'method')<-'CEC'
	return(out)	
}


#' @title Weighting on membership
#' 
#' @description Weighting on Membership (WonM) is similar to CEC as the dendrograms are divided into clusters for a range of values for the
#' number of clustes. However, instead of weighting the sum of the incidences matrices, the final matrix for clustering is the normal sum of 
#' all incidence matrices.
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or"clusters".	
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}. 
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param nrclusters A sequence of numbers of clusters to cut the dendrogram in. Defaults is a sequence of 5 to 25.
#' @param clust Choice of clustering function (character). Defaults to "agnes". 
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible"
#' @return The returned value is a list of two elements:
#' \item{DistM}{The resulting incidence matrix}
#' \item{Clust}{The resulting clusters}
#' The value has class 'WonM'.
#' @examples
#' data(fingerprintMat)
#' data(targetMat)
#' L=list(fingerprintMat,targetMat)
#' 
#' MCF7_WonM=WonM(List=L,type="data",distmeasure=c("tanimoto", "tanimoto"),
#' normalize=c(FALSE,FALSE),method=c(NULL,NULL),nrclusters=seq(5,25,1),
#' clust="agnes",linkage=c("flexible","flexible"),alpha=0.625)
WonM=function(List,type=c("data","dist","clusters"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),nrclusters=seq(5,25,1),clust="agnes",linkage=c("flexible","flexible"),alpha=0.625){
	
	type<-match.arg(type)
	
	
	#Step 1: Distance Matrices
	if(type=="data"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,,drop=FALSE]
		}
		Dist=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i],normalize[i],method[i]))
	}
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		Dist=List
	}
	else{
		Dist=lapply(seq(length(List)),function(i) return(List[[i]]$DistM))
		
		OrderNames=rownames(Dist[[1]])
		for(i in 1:length(Dist)){
			Dist[[i]]=Dist[[i]][OrderNames,OrderNames]
		}
	}
	
	
	#Step 2: perform hierarchical clustering on both distance matrices
	
	HClustering=lapply(seq(length(List)),function(i) cluster::agnes(Dist[[i]],diss=TRUE,method=linkage[i],par.method=alpha))
	
	
	#Step 3: cut the dendrograms into a range of K values
	
	#Give 0 to pair belonging together, give 1 to a pair not belonging together : ==> Distances created otherwise similarities.
	ClusterMembers<-function(HClust,nrclusters){
		Temp=lapply(seq(length(nrclusters)),function(i) stats::cutree(HClust,nrclusters[i]))		
		CM=lapply(seq(length(nrclusters)),function(i) matrix(0,dim(List[[1]])[1],dim(List[[1]])[1]))
		
		clusters<-function(temp,cm){
			for(l in 1:length(temp)){
				label=temp[l]
				sameclust=which(temp==label)
				cm[l,sameclust]=1		
			}
			return(cm)
		}
		
		CM2=lapply(seq(length(nrclusters)),function(i) clusters(temp=Temp[[i]],cm=CM[[i]]))
		Consensus2=Reduce("+",CM2)
		Consensus2=Consensus2/length(nrclusters)
		return(Consensus2)
		
		
	}
	
	Consensus=lapply(seq(length(List)), function(i) ClusterMembers(HClustering[[i]],nrclusters))
	
	OverallConsensus=Reduce("+",Consensus)	
	OverallConsensus=OverallConsensus/length(List)
	OverallConsensus=as.matrix(OverallConsensus)
	rownames(OverallConsensus)=rownames(Dist[[1]])
	colnames(OverallConsensus)=rownames(Dist[[1]])
	OverallConsensusD=1-OverallConsensus
	OverallClusteringR=cluster::agnes(OverallConsensusD,diss=TRUE,method="ward")
	
	out=list(DistM=OverallConsensusD,Clust=OverallClusteringR)
	attr(out,'method')<-'WonM'
	return(out)
}


#' @title Multi-source ABC clustering
#' 
#' @description The Aggregating Bundles of Clusters (ABC, \insertCite{Amaratunga2008}{IntClust}) was originally developed for a single 
#' gene expression data. We extend this method to incorporate multiple data sets of any source. Multi-Source ABC (M-ABC) is an iterative 
#' algorithm in which for each iteration a random sample of objects and features is taken of each data set. A clustering algorithm is run 
#' on each subset and an incidence matrix $C$ is set up by dividing the resulting dendrogram in $k$ clusters. After $r$ iterations, 
#' all incidence matrices are summed and divided by number of times two objects were selected simultaneously. This similarity value is 
#' transformed into a dissimilarity measure expressing the number of times the objects are not clustered together when both are selected. 
#' The obtained matrix is used a input into a clustering algorithm.
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param transpose Logical, whether the data should be transposed to have the ABC orginal format of rows being the variables and columns the samples. Defaults to TRUE.
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param weighting Logical value indicating whether the rows should be weighted in the resampling. Default is c(FALSE,FALSE) for two data sets.
#' @param stat The statistic to be used in weighing the rows. Currently the F-statistic, Coefficient of Variation, Double Bump statistic, and Variance are allowed. The corresponding inputs for these should be "F", "cv", "db", and "var".If the rows are to be weighed equally, any other string will do.
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, default is FALSE. This is recommended if different distance types are used. More details on normalization in \code{Normalization}
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param gr A prespecified grouping of the samples to be used in calculating the F-statistic if stat="F".
#' @param bag Logical, indicating whether the columns should be bagged in each iteration. Defaults to TRUE.
#' @param numsim The number of iterations to be used in the ABC Algorithm. Defaults to 1000.
#' @param numvar The number of featurus to be used at each iteration to calculate the temporary clusters in the ABC Algorithm. Default is c(100,100) for two data sets.
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible"
#' @param NC Expected number of clusters in the data; passed to Wards Method in each iteration. 
#' @param NC2 Expected number of clusters in the data; passed to Wards Method in the final calculation of the clusters. By default set to NC. If NC2="syl", a silhouette will be used to determine the most likely number of clusters.
#' @param mds Logical, indicating whether the dissimilarities calculated in the ABC Algorithm should be plotted using Multi Dimensional Scaling. Defaults to FALSE
#' @references 
#' \insertRef{Amaratunga2008}{IntClust}
#' @return The returned value is a list of two elements:
#' \item{DistM}{The resulting distance matrix matrix}
#' \item{Clust}{The resulting clustering}
#' The value has class 'Ensemble'.
#' @examples
#' data(fingerprintMat)
#' data(targetMat)
#' L=list(fingerprintMat,targetMat)
#' 
#' MCF7_MABC=M_ABC(List=L,transpose=TRUE,distmeasure=c("tanimoto", "tanimoto"),
#' weighting=c(FALSE,FALSE),stat="var",normalize=c(FALSE,FALSE),method=c(NULL,NULL),
#' gr=c(),bag=TRUE, numsim=1000,numvar=c(100,100),linkage=c("flexible","flexible"),
#' alpha=0.625,NC=7, NC2=NULL, mds=FALSE)
M_ABC<- function(List,transpose=TRUE,distmeasure=c("tanimoto", "tanimoto"),weighting=c(FALSE,FALSE),stat="var",normalize=c(FALSE,FALSE),method=c(NULL,NULL), gr=c(),bag=TRUE, numsim=1000, numvar=c(100,100),linkage=c("flexible","flexible"),alpha=0.625,NC=NULL, NC2=NULL, mds=FALSE)
{
	Results=c()
	for(i in 1:length(List)){
		data=List[[i]]		
		res=ABC.SingleInMultiple(data,transpose=transpose,distmeasure=distmeasure[i],weighting=weighting[i],stat=stat,normalize=normalize[i],method=method[i],gr=gr, bag=bag, numsim=numsim, numvar=numvar[i],linkage=linkage[i],alpha=alpha,NC=NC, NC2=NC2, mds=mds)
		Results=rbind(Results,res)
		
	}
	
	f.clustABC.MultiSource(res=Results, numclust=NC2, distmeth=0,linkage="ward",alpha=alpha, mds=mds)
}

#' @title Single-source ABC clustering
#' 
#' @description The Aggregating Bundles of Clusters (ABC, \insertCite{Amaratunga2008}{IntClust}) was originally developed for a single 
#' gene expression data. ABC is an iterative algorithm in which for each iteration a random sample of objects and features is taken of each data set. A clustering algorithm is run 
#' on each subset and an incidence matrix $C$ is set up by dividing the resulting dendrogram in $k$ clusters. After $r$ iterations, 
#' all incidence matrices are summed and divided by number of times two objects were selected simultaneously. This similarity value is 
#' transformed into a dissimilarity measure expressing the number of times the objects are not clustered together when both are selected. 
#' The obtained matrix is used a input into a clustering algorithm.
#' @export
#' @param data A data matrix. It is assumed the rows are corresponding with the objects.
#' @param transpose Logical, whether the data should be transposed to have the ABC orginal format of rows being the variables and columns the samples. Defaults to TRUE.
#' @param distmeasure The distance measurs to be used for the data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to "euclidean".
#' @param weighting Logical value indicating whether the rows should be weighted in the resampling.
#' @param stat The statistic to be used in weighing the rows. Currently the Coefficient of Variation and Variance are allowed. The corresponding inputs for these should be, "cv" and "var". If the rows are to be weighed equally, any other string will do.
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, default is FALSE. This is recommended if different distance types are used. More details on normalization in \code{Normalization}
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is NULL.
#' @param gr A prespecified grouping of the samples to be used in calculating the F-statistic if stat="F".
#' @param bag Logical, indicating whether the columns should be bagged in each iteration. Defaults to TRUE.
#' @param numsim The number of iterations to be used in the ABC Algorithm. Default is 1000.
#' @param numvar The number of featurus to be used at each iteration to calculate the temporary clusters in the ABC Algorithm.
#' @param linkage Choice of inter group dissimilarity (character). Defaults to "ward".
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible"
#' @param NC Expected number of clusters in the data; passed to Wards Method in each iteration.  Default is NULL.
#' @param NC2 Expected number of clusters in the data; passed to Wards Method in the final calculation of the clusters. By default set to NULL such that NC2=NC. If NC2="syl", a silhouette will be used to determine the most likely number of clusters.
#' @param mds Logical, indicating whether the dissimilarities calculated in the ABC Algorithm should be plotted using Multi Dimensional Scaling. Defaults to FALSE.
#' @return The returned value is a list of two elements:
#' \item{DistM}{The resulting distance matrix matrix}
#' \item{Clust}{The resulting clustering}
#' The value has class 'Ensemble'.
#' @references 
#' \insertRef{Amaratunga2008}{IntClust}
#' @examples
#' data(fingerprintMat)
#' data(targetMat)
#' L=list(fingerprintMat,targetMat)
#' 
#' MCF7_ABC=ABC.SingleInMultiple(data=fingerprintMat,transpose=TRUE,distmeasure="tanimoto",
#' weighting=TRUE,stat="var", gr=c(),bag=TRUE, numsim=100,numvar=100,linkage="flexible",
#' alpha=0.625,NC=7, NC2=NULL, mds=FALSE)
ABC.SingleInMultiple<- function(data,transpose=TRUE,distmeasure="euclidean",weighting=FALSE,stat="var",normalize=FALSE,method=NULL,gr=c(), bag=TRUE, numsim=1000, numvar=100,linkage="ward",alpha=0.625,NC=NULL, NC2=NULL, mds=FALSE)
{
	if(transpose){
		#message("The data has been transposed. Make sure that the rows are the measured variables and the columns the samples.")
		data=t(data)
	}
	
	if(is.null(NC)){
		NC=length(table(colnames(data)))-1
	}
	
	if(is.null(NC2)){
		NC2=NC
	}
	
	dataORIG<-data
	
	nc<-ncol(data)
	
	res<-c()
	
	if(!bag){
		zf<- rep(1,nrow(data))		
		
		if(weighting){
		 	if(stat=="cv") zf<-sqrt(f.rmv(data)[,2])/f.rmv(data)[,1]
			else if(stat=="tt") zf<-apply(data,1,f.t)
			else if(stat=="var") zf<-f.rmv(data,varonly=TRUE)
			else zf<- rep(1,nrow(data))
			
		}
	}

	for(i in 1:numsim){
		if(bag){
			boot<-sample(1:nc,replace=TRUE)
			boot<-unique(boot)
			data2<-data[,boot]
			data3<-dataORIG[,boot]
			if(weighting){
				
				if(stat=="cv") z<-sqrt(f.rmv(data2)[,2])/f.rmv(data2)[,1]
				else if(stat=="tt") z<-apply(data2,1,f.t)
				else if(stat=="var")z<-f.rmv(data3,varonly=TRUE)
				else z<-rep(1,nrow(data))
			}
			else{
				z<-rep(1,nrow(data))
			}
		}
		
		else{
			boot<-1:nc
			z<-zf
			data2<-data
		}
		ggiv<-f.gsample(z, ng=numvar)
		
		Ok=FALSE
		while(Ok==FALSE){			
			if(class(as.vector(data2[ggiv,]))=="character"){
				if(distmeasure=="MCA_coord"){
					dat=t(data2[ggiv,])
					dat=as.data.frame(dat)
					dat=dat[, sapply(dat, nlevels) > 1,drop=FALSE]
					if(ncol(dat)<=1){
						ggiv<-f.gsample(z, ng=numvar)
					}
					else{
						Ok=TRUE
					}
				}
			}
			
			else if(any(colSums(data2[ggiv,])==0)){
				ggiv<-f.gsample(z, ng=numvar)
			}
			else{
				dat=t(data2[ggiv,])
				Ok=TRUE
			}				
		}			
		
		res1<-1:nc
		
		
		DistM=Distance(dat,distmeasure=distmeasure,normalize,method)  #altered metric
		Clust=cluster::agnes(DistM,diss=TRUE,method=linkage,par.method=alpha)
		
		HClust=as.hclust(Clust)
		res1[boot]<-cutree(HClust,NC)
		res1[-boot]<-0
		res<-rbind(res,res1)  #every row indicates an iteration: if 0: compound not selected ;; if a number: number of the cluster the compound was appointed too
		#cat(".")
		#if(i%%100==0){cat(paste(i,"\n" ))}
	}
	
	colnames(res)<-colnames(data)
	return(res)
}

#' @title ff
#' @param x data variable
#' @param his Logical. Whether a histogram should be plotted. Default is FALSE.
#' @description Internal function of M_ABC: determine weights with two sample t-test.
f.t <- function(x,his=FALSE) {
	n = length(x)
	tt= NULL
	x = sort(x)
	for(i in 4:(n-3)){ x1 = x[1:(i-1)]; x2 = x[i:n];
		m1 = mean(x1); m2 = mean(x2); 
		sp = sqrt(((i-1)*var(x1) + (n-i+1)*var(x2))/(n-2))*
				sqrt(1/(i-1) + 1/(n-i+1) )
		tt = c(tt,( m1 - m2)/sp)
	}
	if(his) hist(x)
	max(abs(tt))
}

#' @title f.rmv
#' @param x data variable
#' @param varonly Logical. Whether only the variance should be returned (TRUE) or the mean as well (FALSE). Default is FALSE. 
#' @description internal function of M_ABC: computes means and variances.
f.rmv = function(x,varonly=FALSE) {
	m <- ncol(x)
	n <- matrix(as.numeric(!is.na(x)),ncol=m)%*%rep(1,m)
	xx <- x
	xx[is.na(xx)]<- 0
	sx <- xx%*%rep(1,m)
	xx <- xx*xx
	sxx <- xx%*%rep(1,m)
	if(varonly) return((sxx-((sx*sx)/n))/(n-1))
	else return(cbind(sx/n,(sxx-((sx*sx)/n))/(n-1)))
}

#' @title f.gsample
#' @param zf Statistics to be used in the weighting of the rows.
#' @param ng The number to resample. Default is 100.
#' @description internal function of M_ABC: samples a predetermined number of rows according to specified weights.
f.gsample <- function(zf,ng=100) {
	w = 1/(rank(-zf)+100)
	w = w/sum(w)
	sample(length(zf),ng,replace=FALSE, prob=w)
}

#' @title f.clustABC.MultiSoucre
#' @param res Matrix object whose rows are the base clusters determined in each iteration of the ABC algorithm.
#' @param numclust The number of clusters. Default is NULL.
#' @param distmeth Binary, indicating whether the ABC dissimilarities should be scaled by the total number of simulations. By default distmeth=1 indicating the dissimilarities should be left as is.
#' @param linkage Choice of inter group dissimilarity (character) for the final clustering. Defaults to "ward.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible"
#' @param mds  Logical, indicating whether and MDS plot of the dissimilarities should be drawn. Default is FALSE
#' @description Internal function of M_ABC: performs the final clustering.
f.clustABC.MultiSource <- function(res, numclust=NULL, distmeth=1,linkage="ward",alpha=0.625, mds=FALSE){
	
	
	x<-matrix("NA",ncol=ncol(res),nrow=ncol(res))
	for(i in 1:ncol(x)){
		for(j in 1:ncol(x)){
			x[i,j]<-sum(res[((res[,i]!=0)*(res[,j]!=0))==1,i]!=
							res[((res[,i]!=0)*(res[,j]!=0))==1,j])  ### The number of times the objects are not clustered together when both are selected
		}
	}
	
	if(distmeth==1) x<-matrix(as.numeric(x),ncol=ncol(x))
	else{
		x<- (nrow(res)-matrix(as.numeric(x),ncol=ncol(x)))/nrow(res)
		x<- sqrt(1-x^2)
	}
	
	if(mds) plot(cmdscale(dist(x), 2))
	
#	cutree(hclust((as.dist(x)),method="ward"), numclust)
	
	Dist=x
	colnames(Dist)=colnames(res)
	rownames(Dist)=colnames(res)
	Clust=cluster::agnes(Dist,diss=TRUE,method=linkage,par.method=alpha)
	
	
	out=list(DistM=Dist,Clust=Clust)
	attr(out,'method')<-'Ensemble'
	return(out)
}

# Hierarchy-based procedures

#' @title Ensemble for hierarchical clustering
#' 
#' @description The Ensemble for Hierarchical Clustering (EHC, \insertCite{Hossain2012}{IntClust}) defines the strength of association between a 
#' pair of objects as a measure of how closely these are associated taking into account the levels of the dendrogram. Therefore, the sum of the 
#' normalized depths of the clusters in which both objects reside is taken as a measure of association. The depths are weighted by the intra-cluster 
#' proximity values. The resulting similarity matrix is seen as an adjacency matrix of a graph and the METIS algorithm is performed to cut the graph 
#' in $k$ clusters.
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or"clusters".	
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, default is FALSE. This is recommended if different distance types are used. More details on normalization in \code{Normalization}
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param clust Choice of clustering function (character). Defaults to "agnes". 
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible"
#' @param gap Logical. Whether the optimal number of clusters should be determined with the gap statistic. Defaults to FALSE.
#' @param maxK The maximal number of clusters to investigate in the gap statistic. Default is 15.
#' @param graphPartitioning The graph-partitioning method to be performed: "METIS" (implemented in MATLAB), "MST". 
#' @param optimalk An estimate of the final optimal number of clusters. Default is 7.
#' @param waitingtime The time in seconds to wait until the MATLAB results are generated. Defaults to 300.
#' @param file_number The specific file number to be placed as a tag in the file generated by MATLAB. Defaults to 00.
#' @param executable Logical. Whether the METIS MATLAB function is performed via an executable on the command line (TRUE, only possible for Linux systems) or by calling on MATLAB directly (FALSE). Defaults to FALSE. 
#' @return The returned value is a list of two elements:
#' \item{DistM}{The resulting distance matrix}
#' \item{Clust}{The resulting clusters}
#' The value has class 'Ensemble'.
#' @references \insertRef{Hossain2012}{IntClust}
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' L=list(fingerprintMat,targetMat)
#' 
#' MCF7_EHC=EHC(List=L,type="data",distmeasure=c("tanimoto", "tanimoto"),normalize=
#' c(FALSE,FALSE),method=c(NULL,NULL),clust="agnes",linkage = c("flexible","flexible"),
#' alpha=0.625,gap=FALSE,maxK=15,graphPartitioning="METIS",optimalk=7,
#' waitingtime=300,file_number=00,executable=FALSE)
#' 
#' }
EHC<-function(List,type=c("data","dist","clust"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust = "agnes", linkage = c("flexible","flexible"),alpha=0.625,gap = FALSE, maxK = 15,graphPartitioning=c("METIS","MST"),optimalk=7,waitingtime=300,file_number=00,executable=FALSE){
	
	## Step 1: perfom aggl clustering on each data set
	
	if(type=="data"){
		
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type="data",distmeasure[i],normalize[i],method[i],clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
	}
	
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize=FALSE,method=NULL,clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
	}
	else if(type=="clust"){
		
		Clusterings=List
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		
	}
	
	## Step 2: compute strengths of association
	theta_ab=matrix(0,nrow=nrow(Clusterings[[1]]$DistM),ncol=nrow(Clusterings[[1]]$DistM))
	rownames(theta_ab)=rownames(Clusterings[[1]]$DistM)
	colnames(theta_ab)=rownames(Clusterings[[1]]$DistM)
	
	for(i in 1:length(Clusterings)){
		Dend=Clusterings[[i]]$Clust
		nclusters=nrow(Dend$merge)
		MergingList=list()
		
		for(a in 1:nclusters){
			NewMerge=Dend$merge[a,]
			
			if(sign(NewMerge[1])==-1 & sign(NewMerge[2])==-1){
				MergingList[[length(MergingList)+1]]=abs(NewMerge)				
			}
			
			else if(sign(NewMerge[1])==1 & sign(NewMerge[2])==-1){
				MergingList[[length(MergingList)+1]]=c(MergingList[[NewMerge[1]]],abs(NewMerge[2]))					
			}
			
			else if(sign(NewMerge[1])==-1 & sign(NewMerge[2])==1){
				MergingList[[length(MergingList)+1]]=c(abs(NewMerge[1]),MergingList[[NewMerge[2]]])					
			}
			
			else if(sign(NewMerge[1])==1 & sign(NewMerge[2])==1){
				MergingList[[length(MergingList)+1]]=c(MergingList[[NewMerge[1]]],MergingList[[NewMerge[2]]])					
			}
		}
		
		#depths
		depths=c()
		for(b in 1:length(MergingList)){
			SubList=MergingList[[b]]
			depth=0
			for(c in 1:length(MergingList)){
				if(all(SubList%in%MergingList[[c]])){
					depth=depth+1
				}
			}	
			depths=c(depths,depth)
		}
		
		depths=depths-1 #root cluster has depth 0
		
		#max_depth
		max_depth=max(depths)
		
		#intra_cluster proximity values: based on cd ; compare with heights
		values=sort(Dend$height)
		
		#values=as.matrix(cophenetic(as.hclust(Clusterings[[i]]$Clust)))
		
		CheckDist<-function(Dist,StopRange){
			if(StopRange==FALSE & !(0<=min(Dist) & max(Dist)<=1)){
				#message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
				Dist=Normalization(Dist,method="Range")
			}
			else{
				Dist=Dist
			}
		}
		
		values=CheckDist(values,StopRange=FALSE)  #use 1-height as the intra-cluster proximity values after normalizing height between 0 and 1
		proximityvalues=1-values  #one for each cluster
		
		for(j in 1:length(MergingList)){
			Pairs=combn(MergingList[[j]],2)
			for(k in 1:ncol(Pairs)){
				ab=Pairs[,k]
				theta_ab[ab[1],ab[2]]=theta_ab[ab[1],ab[2]]+depths[j]*proximityvalues[j]/max_depth	
				theta_ab[ab[2],ab[1]]=theta_ab[ab[1],ab[2]]
			}			
		}
	}
	
	
	
	## Step 3: generate cluster association graph
	## theta_ab can be seen as an adjacency matrix of the graph: the higher the value, the closer the objects
	net=igraph::graph.adjacency(theta_ab,mode="undirected",weighted=TRUE,diag=FALSE)
	igraph::plot.igraph(net,vertex.label=igraph::V(net)$name,layout=igraph::layout.fruchterman.reingold, edge.color="black",edge.width=igraph::E(net)$weight)
	
	
	if(graphPartitioning=="METIS"){
		
		matlabdata=theta_ab
		utils::write.table(matlabdata,file=paste("matlabdata_",file_number,".csv",sep=""),sep=",",col.names=FALSE,row.names=FALSE)
		
		if(executable){
			system(paste("./MetisAlgorithm ",optimalk," ",file_number,sep=""),intern=TRUE)
		}
		else{
			matlabcode=c(paste("cls = csvread('matlabdata_",file_number,".csv')",sep=""),
				paste("Partition_",file_number,"=hmetis(cls,",optimalk,",",file_number,")",sep=""),
				paste("csvwrite('Partition_",file_number,".csv', Partition_",file_number,")",sep=""))
		
			writeLines(matlabcode, con=paste("MetisAlgorithm_",file_number,".m",sep=""))
		}
		Continue=FALSE
		time=0
		while(Continue==FALSE){
			Sys.sleep(15)
			time=time+15
			Continue=file.exists(paste("Partition_",file_number,".csv",sep=""))
			if(time>waitingtime & Continue==FALSE){
				stop(paste("Waited",waitingtime, "seconds for completion of the ensemble clustering procedure. Increase waiting time to continue.",sep=" "))
			}
		}	
		
		
		Partition <- utils::read.table(paste("Partition_",file_number,".csv",sep=""),sep=",")
		Partition=as.vector(as.matrix(Partition))
		names(Partition)=rownames(theta_ab)
		
		order=sort(Partition,index=TRUE)$ix
		order.lab=names(Partition[sort(Partition,index=TRUE)$ix])
		Clusters=Partition
		
		if(!executable){
			file.remove(paste("MetisAlgorithm_",file_number,".m",sep=""))
		}
		file.remove(paste("Partition_",file_number,".csv",sep=""))
		file.remove(paste("matlabdata_",file_number,".csv",sep=""))
		
	}
	
	else if(graphPartitioning=="MST"){
		
		#are assigning to the first encounter to break ties and get clusters: if changed  to join all if one is in common, all end up in 1 cluster
		
		Graph=igraph::graph.adjacency(adjmatrix=theta_ab, mode=c( "undirected"), weighted=TRUE, diag=TRUE,add.colnames=TRUE, add.rownames=NA)
		MST_Graph=igraph::as_data_frame(igraph::minimum.spanning.tree(Graph))
		
		Partition=list()
		Partition[[1]]=c(MST_Graph[1,1],MST_Graph[1,2])
		Placed=rep(FALSE,nrow(theta_ab))
		Placed[Partition[[1]]]=TRUE
		for(j in 2:nrow(MST_Graph)){
			Edge=as.numeric(MST_Graph[j,c(1:2)])		
			k=1
			
			while(any(Placed[Edge]==FALSE)){
				
				if(k>length(Partition)){
					Partition[[length(Partition)+1]]=Edge
					Placed[Edge]=TRUE
				}	
				
				else if(Edge[1]%in%Partition[[k]] | Edge[2]%in%Partition[[k]]){
					Partition[[k]]=c(unlist(Partition[[k]]),Edge)
					Partition[[k]]=unique(Partition[[k]])
					Placed[Edge]=TRUE
				}
				k=k+1
			}	
			
			
		}	
		
		order=unlist(Partition)
		order.lab=rownames(theta_ab)[order]
		t1<-sapply(1:length(Partition),function(i) rep(i,length(Partition[[i]])))
		t2<-sapply(1:length(Partition),function(i) rownames(theta_ab)[Partition[[i]]])
		t3<-cbind(unlist(t1),unlist(t2))
		clus=as.numeric(t3[,1])
		names(clus)=t3[,2]
		Clusters=clus[rownames(theta_ab)]
		
		
	}
	
	Out=list(DistM=theta_ab,Clust=list(order=order,order.lab=order.lab,Clusters=Clusters))
	attr(Out,"method")="Ensemble"
	return(Out)
	
}

#' @title Hierarchical ensemble clustering
#' 
#' @description \insertCite{Zheng2014}{IntClust} proposed the Hierarchical Ensemble Clustering (HEC) algorithm. For each dendrogram, the cophenetic 
#' distances between the object are calculated. The distances are aggregated across the data sets and an ultra-metric which is the closest to the 
#' distance matrix is determined. A final hierarchical clustering is based on the ultra-metric values.
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or"clusters".	
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, default is FALSE. This is recommended if different distance types are used. More details on normalization in \code{Normalization}
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param clust Choice of clustering function (character). Defaults to "agnes". 
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible".
#' @return The returned value is a list of two elements:
#' \item{DistM}{The resulting distance matrix}
#' \item{Clust}{The resulting hierarchical structure}
#' The value has class 'HEC'.
#' @references \insertRef{Zheng2014}{IntClust}
#' @examples
#' data(fingerprintMat)
#' data(targetMat)
#' L=list(fingerprintMat,targetMat)
#' 
#' MCF7_HEC=HierarchicalEnsembleClustering(List=L,type="data",distmeasure=
#' c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),
#' clust="agnes",linkage=c("flexible","flexible"),alpha=0.625)
HierarchicalEnsembleClustering<-function(List,type=c("data","dist","clust"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust = "agnes", linkage = c("flexible","flexible"),alpha=0.625){
	
	## Step 1: perfom aggl clustering on each data set
	
	if(type=="data"){
		
		if(is.null(rownames(List[[1]]))){
			for(i in 1:length(List)){
				rownames(List[[i]])=seq(1:nrow(List[[i]]))
			}
		}
		else{
			OrderNames=rownames(List[[1]])
			
			for(i in 1:length(List)){
				List[[i]]=List[[i]][OrderNames,,drop=FALSE]
			}
		}
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type="data",distmeasure[i],normalize[i],method[i],clust,linkage[i],alpha,gap=FALSE,maxK=5,StopRange=TRUE))
		#Clusterings=lapply(seq(length(List)),function(i) agnes(List[[i]],method=linkage[i],par.method=0.625))
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
	}
	
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize=FALSE,method=NULL,clust,linkage[i],alpha,gap=FALSE,maxK=5,StopRange=TRUE))
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
	}
	else if(type=="clust"){
		
		Clusterings=List
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		
		
	}
	
	## Step 2: cophenetic distance of each object in each dendrogram
	
	CD=list()
	for(j in 1:length(Clusterings)){
		CD[[j]]=as.matrix(cophenetic(as.hclust(Clusterings[[j]]$Clust)))
	}
	
	## Step 3: Aggregate the distance of the dendrograms and compute the matrix D
	## ? what is meant by aggregate?
	## Assumed: aggragete ==  added
	
	D_H=Reduce("+",CD)
	
	D=1/2*D_H
	
	## Step 4: Find an ultra-metric distance T closest to D
	## Apply the modified Floyd-Warshall algorithm on D
	
	Floyd_Warshall<-function(G){	
		M=G
		N=nrow(G)
		for(k in c(1:N)){
			for(i in c(1:N)){
				for(j in c(1:N)){
					M[i,j]=min(M[i,j],max(M[i,k],M[k,j]))
					
				}
			}
		}
		return(M)
	}
	
	T=Floyd_Warshall(D)  # Distance matrix
	
	## Step 5: Construct the final hierarchical clustering based on T
	## the alpha-cut method form Meyer et al 2004
	
	## ? defuzzify
	
#	## Should
#	test=cut(T, level = 1, type = c("alpha"), strict = FALSE)
#	
#	#Testing 
#	A=matrix(c(0,1,14,18,18,18,1,0,14,18,18,18,14,14,0,18,18,18,18,18,18,0,15,15,18,18,18,15,0,2,18,18,18,15,2,0),ncol=6,nrow=6)
#	
	
	
	#try alpha cut for each value: denote cluster
	#build tree from this?
	
	alpha=sort(unique(as.vector(T)),decreasing = TRUE)
	C=list()
#	for(i in alpha[-length(alpha)][c(1:11)]){
#		C[[length(C)+1]]=archi(as.dist(T),i)[[1]]
#		
#		if(length(C)>1){
#
#			if(max(C[[length(C)]])>max(C[[length(C)-1]])+1){
#				Changed=which(C[[length(C)]]+ C[[length(C)-1]]<2*C[[length(C)]]-1)
#				Groups=unique(C[[length(C)]][Changed])
#				
#				if(length(Groups)==1){
#					N=C[[length(C)]][-c(Changed)]
#					O=C[[length(C)-1]][-c(Changed)]
#					Joinedtemp=which(!as.vector(c(table(N)))==as.vector(c(table(O),0)))[c(2)]
#					Joined=as.numeric(names(table(O))[Joinedtemp])
#					New2=which(C[[length(C)]]==Joined[1])
#					Groups=c(Groups,unique(C[[length(C)]][New2]))
#					Changed=c(Changed,New2)
#				}
#				
#				Group1=Changed[which(C[[length(C)]][Changed]==Groups[1])]
#				Group2=Changed[which(C[[length(C)]][Changed]==Groups[2])]
#				
#				Clus=C[[length(C)-1]]
#				Clusnext=C[[length(C)]]
#				Clus1=Clusnext
#				ChangeAfter=unique(Clus1[Group2])
#				Clus1[Group2]=Clus[Group2]
#				Clus1[which(Clus1>ChangeAfter)]=Clus1[which(Clus1>ChangeAfter)]-1
#				
#				Clus2=Clusnext
#				
#				C=C[-c(length(C))]
#				C[[length(C)+1]]=Clus1
#				C[[length(C)+1]]=Clus2
#				
#			}
#		}
#	}
	
	for(k in c(1:nrow(List[[1]]))){
		Y=cutree(as.hclust(agnes(T,diss=TRUE)),k)
		names(Y)=NULL
		C[[length(C)+1]]=Y
	}
	#C[[length(C)+1]]=seq(1,nrow(T))
	
	Record=rep(0,length(C))
	names(Record)=seq(1,length(C))
	Merge=matrix(0,ncol=2,nrow=1)
	MergingList=list()
	
	for(i in length(C):1){
		if(i==length(C)){
			#All Leafs - Bottom of the tree
			Leafs=C[[length(C)]]	
			Record=Leafs
		}
		
		else{
			
			NewMerge=C[[i]]
			if(i==(length(C)-1)){
				NewOnes=which(duplicated(NewMerge))
				for(n in 1:length(NewOnes)){
					MergedTogether=which(NewMerge==NewMerge[NewOnes][n])
					if(n==1){
						Merge[1,]=MergedTogether*(-1)
						
					}
					else{
						Merge=rbind(Merge,MergedTogether*(-1))
						
					}
					MergingList[[length(MergingList)+1]]=MergedTogether
					Record=NewMerge
				}
			}
			
			else{
				
				Changed=which(NewMerge!=Record)
				if(length(Changed)!=0){
					
					PositionChanged=Position(function(x) identical(x,Changed), MergingList, nomatch = 0) 
					if(PositionChanged==0){
						
						formedafter=FALSE
						NewOnes=unique(c(which(duplicated(NewMerge)),which(duplicated(NewMerge,fromLast=TRUE))))[which(!unique(c(which(duplicated(NewMerge)),which(duplicated(NewMerge,fromLast=TRUE))))%in%unique(c(which(duplicated(Record)),which(duplicated(Record,fromLast=TRUE)))))]
						if(any(NewMerge[Changed]%in%NewMerge[NewOnes])){
							Temp=Changed[which(NewMerge[Changed]%in%NewMerge[NewOnes])]
							if(length(intersect(Temp,NewOnes))>0){
								NewOnes=Temp
								formedafter=TRUE
							}
						}
					}
					else{
						formedafter=FALSE
						NewOnes=which(NewMerge==unique(NewMerge[Changed]))[!which(NewMerge==unique(NewMerge[Changed]))%in%Changed]
					}
					if(length(NewOnes)==0){
						Joined=which(!as.vector(c(table(NewMerge),0))==as.vector(table(Record)))[c(1,2)]
						NewOnes=which(Record==Joined[1])
						
					}
					for(k in 1:length(unique(NewMerge[NewOnes]))){
						MergedTogether=which(NewMerge==unique(NewMerge[NewOnes])[k])
						
						Position1=Position(function(x) identical(x,(as.integer(MergedTogether[-which(MergedTogether%in%NewOnes)]))), MergingList, nomatch = 0) 
						if(Position1==0){
							Merge=rbind(Merge,MergedTogether*(-1))  #they are both leaves
							
						}
						else{
							#First group was formed before
							if(formedafter==TRUE){
								TrulyNew=unique(c(which(duplicated(NewMerge)),which(duplicated(NewMerge,fromLast=TRUE))))[which(!unique(c(which(duplicated(NewMerge)),which(duplicated(NewMerge,fromLast=TRUE))))%in%unique(c(which(duplicated(Record)),which(duplicated(Record,fromLast=TRUE)))))]
								Group=NewOnes[-c(which(NewOnes%in%TrulyNew))]
								Set=list(TrulyNew,Group)
								for(s in Set){
									if(length(s)>0){
										Position2=Position(function(x) identical(x,(as.integer(MergedTogether[which(MergedTogether%in%s)]))), MergingList, nomatch = 0) 
										
										if(Position2==0){  #first one was a group, newones is a leaf
											Merge=rbind(Merge,c(Position1,(intersect(s,MergedTogether)*(-1))))
											
										}
										else{  #both were grouped before
											Merge=rbind(Merge,c(Position1,Position2))
											
										}
									}
								}
							}
							else{
								Position2=Position(function(x) identical(x,(as.integer(MergedTogether[which(MergedTogether%in%NewOnes)]))), MergingList, nomatch = 0) 
								
								if(Position2==0){  #first one was a group, newones is a leaf
									Merge=rbind(Merge,c(Position1,(intersect(NewOnes,MergedTogether)*(-1))))
									
								}
								else{  #both were grouped before
									Merge=rbind(Merge,c(Position1,Position2))
									
								}
							}
						}
						
						MergingList[[length(MergingList)+1]]=MergedTogether
					}
					Record=NewMerge
				}
			}
		}
	}
	
	
	Start=which(Merge==-1,arr.ind=TRUE)[1]
	Checked=rep(0,(nrow(Merge)))
	out=c()
	Row=Start
	while(!all(Checked==1)){
		Joined=Merge[Row,]
		#print(Row)
		if(sign(Joined[1])==-1 & sign(Joined[2])==-1){
			
			Checked[Row]=1
			
			if(Row==Start){
				
				out=abs(Joined)
				RNext=which(Merge==Row,arr.ind=TRUE)[1]
				
			}
			else{
				
				if(any(Checked>1)){
					out=c(out,abs(Joined))
					RNext=which.max(Checked)
				}
				
				else if(any(Checked==0)){
					out=c(out,abs(Joined))
					RNext=which(Merge==Row,arr.ind=TRUE)[1]
				}
			}
			
		}
		
		else if(sign(Joined[1])==1 & sign(Joined[2])==1){
			
			if(all(Checked[Joined]==1)){
				Checked[Row]=1
				
				if(any(Checked>1)){
					RNext=which.max(Checked)
				}
				else if(!(all(Checked==1))){
					RNext=which(Merge==Row,arr.ind=TRUE)[1]
				}
				
			}
			
			else{
				Checked[Row]=max(Checked)+1
				RNext=Joined[which(!Joined%in%which(Checked==1))]
				
				if(length(RNext)==2){
					Checked[RNext[2]]=max(Checked)+1
					RNext=RNext[1]
				}
			}
			
		}
		
		else if(sign(Joined[1])==1 & sign(Joined[2])==-1){
			
			Checked[Row]=1
			
			if(Checked[Joined[1]]==0){
				out=c(out,abs(Joined[2]))
				RNext=Joined[1]
			}
			else{
				if(!(abs(Joined)[2]%in%out)){
					out=c(out,abs(Joined[2]))
				}
				RNext=which(Merge==Row,arr.ind=TRUE)[1]
				
			}			
		}
		
		else if(sign(Joined[1])==-1 & sign(Joined[2])==1){
			Checked[Row]=1
			
			if(Checked[Joined[2]]==0){
				out=c(out,abs(Joined[1])) #possibly to be shifted
				RNext=Joined[2]
				
			}
			else{
				if(!(abs(Joined)[1]%in%out)){
					out=c(out,abs(Joined[1]))
				}
#					print("301")
				RNext=which(Merge==Row,arr.ind=TRUE)[1]
				
			}
			
		}
		
		Row=RNext
	}
	
	
#	OrderObjects<-function(Row=Start,M=Merge,out=c(),Checked=rep(0,(nrow(T)-1))){
#		
#		if(!all(Checked==1)){
	##			print(out)
#			print(Row)
	##			print(Checked)
	##			print("212")
#			
#			Joined=M[Row,]
#
#			if(sign(Joined[1])==-1 & sign(Joined[2])==-1){
#				
#				Checked[Row]=1
#				
#				if(Row==Start){
#					
#					out=abs(Joined)
#					RNext=which(Merge==Row,arr.ind=TRUE)[1]
	##					print("218")
#					out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
#					
#				}
#				else{
	##					print("222")
#					
#					
#					if(any(Checked>1)){
#						out=c(out,abs(Joined))
#						#RNext=which(Merge==which.max(Checked),arr.ind=TRUE)[1]
#						RNext=which.max(Checked)
	##						print("227")
#						#Checked[which.max(Checked)]=1
#						out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
#					}
#					
#					else if(any(Checked==0)){
#						out=c(out,abs(Joined))
#						RNext=which(Merge==Row,arr.ind=TRUE)[1]
#						out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
#					}
	##					return(out)
#				}
#				
#				if((all(Checked==1))){
#					print(out)
#					out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
#				}
#				return(out)
#			}
#			
#			
#			else if(sign(Joined[1])==1 & sign(Joined[2])==1){
#				
#				if(all(Checked[Joined]==1)){
#					Checked[Row]=1
	##					print("239")
#					print(Checked)
#					
#					if(any(Checked>1)){
#						RNext=which.max(Checked)
	##						print("242")
#						out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
#					}
#					else if(!(all(Checked==1))){
#						RNext=which(Merge==Row,arr.ind=TRUE)[1]
	##						print("247")
#						out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
#					}
#					
#				}
#				
#				else{
#					Checked[Row]=max(Checked)+1
#					RNext=Joined[which(!Joined%in%which(Checked==1))]
	##					print("259")
#					
#					if(length(RNext)==2){
#						Checked[RNext[2]]=max(Checked)+1
#						RNext=RNext[1]
	##						print("263")
#						out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
#					}
#					else{
	##						print("269")
#						out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
#					}
#				}
#				
#				if((all(Checked==1))){
#					print(out)
#					out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
#				}
#				
#				return(out)
#			}
#			
#			else if(sign(Joined[1])==1 & sign(Joined[2])==-1){
#				
	##				print("276")
#				Checked[Row]=1
#
#				
#				if(Checked[Joined[1]]==0){
#					out=c(out,abs(Joined[2]))
	##					print("280")
#					RNext=Joined[1]
#					out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
#				}
#				else{
	##					print("288")
#					if(!(abs(Joined)[2]%in%out)){
#						out=c(out,abs(Joined[2]))
#					}
#					RNext=which(Merge==Row,arr.ind=TRUE)[1]
#					out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
#				}
#				if((all(Checked==1))){
	##					print(out)
	##					out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
#			}
#				return(out)
#			}
#			
#			else if(sign(Joined[1])==-1 & sign(Joined[2])==1){
#				Checked[Row]=1
	##				print("293")
#
#				if(Checked[Joined[2]]==0){
#					out=c(out,abs(Joined[1])) #possibly to be shifted
	##					print("296")
#					RNext=Joined[2]
#					out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
#				}
#				else{
#					if(!(abs(Joined)[1]%in%out)){
#						out=c(out,abs(Joined[1]))
#					}
	##					print("301")
#					RNext=which(Merge==Row,arr.ind=TRUE)[1]
#					out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
#				}
#				
#				if((all(Checked==1))){
#					out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
#				}
#				return(out)
#			}
#			
#		}
#		
#		else{
	##			print("here")
	##			print(out)
	##			print(Checked)
#			return(out)
#		}
#		
#		
#	}
#	Order=OrderObjects(Row=Start,M=Merge,out=c(),Checked=rep(0,(nrow(T)-1)))
#	
	
	Order=out
	Heights=c()	
	for(j in 1:nrow(Merge)){
		joined=Merge[j,]
		
		if(sign(joined[1])==-1){
			a=abs(joined[1])
		}
		else{
			a=MergingList[[joined[1]]]
			a=abs(a)
		}
		
		if(sign(joined[2])==-1){
			b=abs(joined[2])
		}
		else{
			b=MergingList[[joined[2]]]
			b=abs(b)
		}
		
		Heights=c(Heights,unique(as.vector(T[a,b]))) #should be put in the order of the names
		
	}
	
	
#	HeightsObjects<-function(Row=Start,M=Merge,H=Heights,heights=c(),Checked=rep(0,(nrow(T)-1))){
#		
#		if(!all(Checked==1)){
#			print(heights)
#			print(Row)
#			print(Checked)
#			print("212")
#			
#			Joined=M[Row,]
#			
#			if(sign(Joined[1])==-1 & sign(Joined[2])==-1){
#				
#				Checked[Row]=1
#				
#				if(Row==Start){
#					
#					heights=c(heights,H[Row])
#					RNext=which(Merge==Row,arr.ind=TRUE)[1]
#					print("218")
#					heights=HeightsObjects(Row=RNext,M=M,H=H,heights=heights,Checked=Checked)
#					
#				}
#				else{
#					print("222")
#					heights=c(heights,H[Row])
#					
#					if(any(Checked>1)){
#						#RNext=which(Merge==which.max(Checked),arr.ind=TRUE)[1]
#						RNext=which.max(Checked)
#						print("227")
#						#Checked[which.max(Checked)]=1
#						heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
#					}
	##					return(out)
#				}
#				
#				if((all(Checked==1))){
#					print(heights)
#					heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
#				}
#				return(heights)
#			}
#			
#			
#			else if(sign(Joined[1])==1 & sign(Joined[2])==1){
#				
#				if(all(Checked[Joined]==1)){
#					Checked[Row]=1
#					heights=c(heights,H[Row])
#					print("239")
#					
#					if(any(Checked>1)){
#						RNext=which.max(Checked)
#						print("242")
#						heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
#					}
#					else if(!(all(Checked==1))){
#						RNext=which(Merge==Row,arr.ind=TRUE)[1]
#						print("247")
#						heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
#					}
#					
#				}
#				
#				else{
#					Checked[Row]=max(Checked)+1
#					#Checked[Row]=1
#					#heights=c(heights,H[Row])
#					RNext=Joined[which(!Joined%in%which(Checked==1))]
#					print("259")
#					
#					if(length(RNext)==2){
#						Checked[RNext[2]]=max(Checked)+1
#						RNext=RNext[1]
#						print("263")
#						heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
#					}
#					else{
#						print("269")
#						heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
#					}
#				}
#				
#				if((all(Checked==1))){
#					print(heights)
#					heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
#				}
#				
#				return(heights)
#			}
#			
#			else if(sign(Joined[1])==1 & sign(Joined[2])==-1){
#				
#				print("276")
#				Checked[Row]=1
#				heights=c(heights,H[Row]) #possibly to be shifted
#								
#				if(Checked[Joined[1]]==0){
#					print("280")
#					RNext=Joined[1]
#					heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
#				}
#				else{
#					print("288")
#					RNext=which(Merge==Row,arr.ind=TRUE)[1]
#					heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
#				}
#				if((all(Checked==1))){
#					print(heights)
#					heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
#				}
#				return(heights)
#			}
#			
#			else if(sign(Joined[1])==-1 & sign(Joined[2])==1){
#				Checked[Row]=1
#				print("293")
#				heights=c(heights,H[Row]) #possibly to be shifted
#				if(Checked[Joined[2]]==0){
#					print("296")
#					RNext=Joined[2]
#					heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
#				}
#				else{
#					print("301")
#					RNext=which(Merge==Row,arr.ind=TRUE)[1]
#					heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
#				}
#				
#				if((all(Checked==1))){
#					print(heights)
#					heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
#				}
#				return(heights)
#			}
#			
#		}
#		
#		else{
#			print("here")
#			print(heights)
#			return(heights)
#		}
#		
#		
#	}
#	HeightsOrder<-HeightsObjects(Row=Start,M=Merge,H=Heights,heights=c(),Checked=rep(0,(nrow(T)-1)))
	
	Labels=rownames(T)
	
	out=list(height=Heights,merge=Merge,order=Order,labels=rownames(T))
	class(out)="hclust"
	plot(out)
	
	order.lab=rownames(T)[Order]
	Out=list(DistM=T,Clust=list(height=Heights,merge=Merge,order=Order,labels=rownames(T),order.lab=order.lab))
	attr(Out$Clust,"class")="hclust"
	attr(Out,"method")="HEC"
	return(Out)
	
}



## Aid Functions

#' @title Normalization of features
#' @description If data of different scales are being employed by the user, it is recommended to perform a normalization to make the data structures comparable.
#' @export
#' @param Data A data matrix. It is assumed the  rows are corresponding with the objects.
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates","standardize","Range" or any of the first letters of these names.
#' @details The method "Quantile" refers to the Quantile-Normalization widely used in 
#' omics data. The "Fisher-Yates" normalization has a similar approach as the Quantile-
#' Normalization but does not rely on the data, just on the number of rows present in
#' the data matrix. The "Standardize" method refers to the \code{stdize} function of the
#' \pkg{pls} package and centers and scales the data matrix. The method "Range" computes
#' the maximum and minimum value of the matrix and determines the range. Every value is then
#' reduced by the minimum and divided by the range of the data matrix. The latter normalization will
#' result in values between 0 and 1.
#' @return The returned value is a normalized matrix.
#' @examples 
#' x=matrix(rnorm(100),ncol=10,nrow=10)
#' Norm_x=Normalization(x,method="R")
Normalization<-function(Data,method=c("Quantile","Fisher-Yates","Standardize","Range","Q","q","F","f","S","s","R","r")){
	
	method=substring(method[1],1,1)  #Function also accepts begin letter of each method	
	method=match.arg(method)
	if(method=="S"|method=="s"){
		Data1<-pls::stdize(Data)#center and divide by standard deviation
		return(Data1)
	}
	
	else if(method=="R"|method=="r"){
#			rangenorm<-function(x){
#				minc=min(x)
#				maxc=max(x)
#				rx=(x-minc)/(maxc-minc)
#			}
#			
#			DataN=apply(Data,2,rangenorm)
#			Data=DataN
#			return(Data)
		
		minc=min(as.vector(Data))
		maxc=max(as.vector(Data))
		DataN=(Data-minc)/(maxc-minc)
		Data=DataN
		return(Data)
		
	}
	else if(method=="Q"|method=="q"){
		DataColSorted <- apply(Data, 2, sort)
		NValues=apply(DataColSorted,1,stats::median,na.rm=TRUE)
	}
	else{
		DataColSorted <- apply(Data, 2, sort)
		NValues=stats::qnorm(1:nrow(Data)/(nrow(Data)+1))
	}
	
	DataRanked=c(apply(Data,2,rank))
	
	Data=array(stats::approx(1:nrow(Data),NValues,DataRanked)$y,dim(Data),dimnames(Data))
	t=stats::approx(1:nrow(Data),NValues,DataRanked)
	return(Data)
	
}

#' @title Single source clustering
#' @description The function \code{Cluster} performs clustering on a single source of information, i.e one data matrix. The option is available to compute the gap statistic to determine the optimal number of clusters. 
#' @export 
#' @param Data A matrix containing the data. It is assumed the rows are corresponding with the objects.
#' @param type Type indicates whether the provided matrix in "Data" is either a data or a distance matrix obtained from the data. If type="dist" the calculation of the distance matrix is skipped. Type should be one of "data" or "dist".	
#' @param distmeasure Choice of metric for the dissimilarity matrix (character). Should be one of "tanimoto", "euclidean", "jaccard","hamming". Default is "tanimoto".
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, default is FALSE. This is recommended if different distance types are used. More details on normalization in \code{Normalization}
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is NULL.
#' @param clust Choice of clustering function (character). Defaults to "agnes". Note for now, the only option is to carry out agglomerative hierarchical clustering as it was implemented in the \code{agnes} function in the cluster package.
#' @param linkage Choice of inter group dissimilarity (character). Defaults to "flexible".
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible"
#' @param gap Logical. Whether the optimal number of clusters should be determined with the gap statistic. Default is TRUE.
#' @param maxK The maximal number of clusters to investigate in the gap statistic. Default is 15.
#' @param StopRange Logical. Indicates whether the distance matrices with values not between zero and one should be standardized to have so.
#' #' If FALSE the range normalization is performed. See \code{Normalization}. If TRUE, the distance matrices are not changed.
#' This is recommended if different types of data are used such that these are comparable. Default is TRUE.
#' @details The gap statistic is determined by the criteria described by the cluster package:
#' firstSEmax, globalSEmax, firstmax,globalmax, Tibs2001SEmax. The number of 
#' iterations is set to a default of 500. The implemented distances to be used for
#' the dissimilarity matrix are jaccard, tanimoto and euclidean. The jaccard distances
#' were computed with the \code{dist.binary(\ldots,method=1)} function in the ade4
#' package and the euclidean ones with the \code{daisy} function in again the cluster
#' package. The Tanimoto distances were implemented manually. 
#' @return The returned value is a list with two elements:
#'\item{DistM}{The distance matrix of the data matrix}
#'\item{Clust}{The resulting clustering}
#'If the gap option was indicated to be true, another 3 elements are joined to the list.
#'Clust\_gap contains the output from the function to compute the gap statistics and
#'gapdata is a subset of this output. Both can be used to make plots to visualize the 
#'gap statistic. The final component is k which is a matrix containing the optimal number
#'of clusters determined by each criterion mentioned earlier.
#' @examples
#' data(fingerprintMat)
#' data(targetMat)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' 		method=NULL,clust="agnes",linkage="flexible",alpha=0.625,gap=FALSE,maxK=55
#' 		,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' 		method=NULL,clust="agnes",linkage="flexible",alpha=0.625,gap=FALSE,maxK=55
#' 		,StopRange=FALSE)
Cluster<-function(Data,type=c("data","dist"),distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="flexible",alpha=0.625,gap=TRUE,maxK=15,StopRange=TRUE){	
	
	#STEP 1: Distance Matrices
	type<-match.arg(type)
	if(type=="data"){
		DistM=Distance(Data,distmeasure,normalize,method)
		if(StopRange==FALSE & !(0<=min(DistM) & max(DistM)<=1)){
			#message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
			DistM=Normalization(DistM,method="Range")
		}
	}
	else{
		DistM=Data
		if(StopRange==FALSE  & !(0<=min(DistM) & max(DistM)<=1)){
			#message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
			DistM=Normalization(DistM,method="Range")
		}
	}
	
	#STEP 2: Hierarchical Clustering
	
	Clust=cluster::agnes(DistM,diss=TRUE,method=linkage,par.method=alpha)
	
	func = function(x,k){
		return(list(cluster=stats::cutree(Clust,k=k)  )  )
	}
	
	#with optional gap statistic
	if(gap==TRUE){
		Clust_gap = cluster::clusGap(Data,FUNcluster=func,K.max=maxK,B=500)
		gapdata = as.data.frame(Clust_gap$Tab)
		
		k1 = cluster::maxSE(gapdata[-maxK,3],gapdata[-maxK,4],"firstSEmax")
		k2 = cluster::maxSE(gapdata[-maxK,3],gapdata[-maxK,4],"globalSEmax")
		k3 = cluster::maxSE(gapdata[-maxK,3],gapdata[-maxK,4],"firstmax")
		k4 = cluster::maxSE(gapdata[-maxK,3],gapdata[-maxK,4],"globalmax")
		k5 = cluster::maxSE(gapdata[-maxK,3],gapdata[-maxK,4],"Tibs2001SEmax")
		
		k = data.frame(firstSEmax=k1,globalSEmax=k2,firstmax=k3,globalmax=k4,Tibs2001SEmax=k5)
		
		out = list(DistM=DistM,Clust=Clust,Clust_gap=Clust_gap,gapdata=gapdata,k=k)
	}
	else{
		out=list(DistM=DistM,Clust=Clust)
	}
	attr(out,'method')<-'Single Clustering'
	return(out)
}


#' @title Distance calculation
#' @description The \code{Distance} function calculates the distances between the data objects. The included
#' distance measures are euclidean for continuous data and the tanimoto coefficient or jaccard index for binary data.
#' @export
#' @param Data A data matrix. It is assumed the  rows are corresponding with the objects.
#' @param distmeasure Choice of metric for the dissimilarity matrix (character). Should be one of "tanimoto", "euclidean", "jaccard","hamming","cont tanimoto","MCA_coord","gower","chi.squared" or "cosine"
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, default is FALSE. This is recommended if different distance types are used. More details on normalization in \code{Normalization}.
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is NULL.
#' @details The euclidean distance distance is included for continuous matrices while
#' for binary matrices, one has the choice of either the jaccard index, the
#' tanimoto coeffcient or the hamming distance. The hamming distance is obtained
#' by applying the \code{hamming.distance} function of the \pkg{e1071} package.
#' It will compute the hamming distance between the rows of the data matrix. The
#' hamming distance counts the number of times where two rows differ in their
#' zero and one values. The Jaccard index is calcaluted as determined by the
#' formula of the \code{dist.binary} function in the \pkg{a4} package and the
#' tanimoto coefficient as described by \cite{Li2011}. For both, first the
#' similarity is calculated as\deqn{s=frac{n11}{n11+n10+n01}}
#' with n11 the number of features the 2 objects have in common, n10
#' the number of features of the first compound and n01 the number of features
#' of the second compound. These similarities are converted to distances by:
#' \deqn{J=\sqrt{1-s}}
#' for the jaccard index and by:
#' \deqn{T=1-s}
#' for the tanimoto coefficient. The lower the similarity values s are, the more features are
#' shared between the two objects and the more alike they are. Since clustering is
#' based on dissimilarity, the conversion to distances is performed.
#' If normalize=TRUE and the distance meausure is euclidean, the data matrix is
#' normalized beforehand. Further, a version of the tanimoto coefficient is also available for continuous data. 
#' @return The returned value is a distance matrix.
#' @examples 
#' 	data(fingerprintMat)
#' 	Dist_F=Distance(fingerprintMat,distmeasure="tanimoto",normalize=FALSE,method=NULL)
Distance<-function(Data,distmeasure=c("tanimoto","jaccard","euclidean","hamming","cont tanimoto","MCA_coord","gower","chi.squared","cosine"), normalize=FALSE,method=NULL){
	if(distmeasure!="MCA_coord"){
		Data <- Data+0
	}
	distmeasure=match.arg(distmeasure)
	if((distmeasure=="euclidean") & normalize==TRUE){
		Data=Normalization(Data,method)				
	}
	
	tanimoto = function(m){
		S = matrix(0,nrow=dim(m)[1],ncol=dim(m)[1])
		
#		for(i in 1:dim(m)[1]){
#			for(j in 1:i){
#				N.A = sum(m[i,])
#				N.B = sum(m[j,])
#				N.C = sum(m[i,(m[i,]==m[j,])])
#				
#				if(N.A==0&N.B==0){
#					coef = 1				
#				}
#				else{
#					coef = N.C / (N.A+N.B-N.C)
#				}
#				S[i,j] = coef
#				S[j,i] = coef
#			}
#			
#		}
		#via matrix multiplication
		m=as.matrix(m)
		N.C=m %*% t(m)
		N.A=m %*% (1-t(m))
		N.B=(1-m) %*% t(m)
		S=N.C/(N.A+N.B+N.C)
		D = 1 - S
		return(D)
	}
	
	# Computing the distance matrices
	
	if(distmeasure=="jaccard"){
		dist = ade4::dist.binary(Data,method=1)
		dist = as.matrix(dist)
	}
	else if(distmeasure=="tanimoto"){
		dist = tanimoto(Data)
		dist = as.matrix(dist)
		rownames(dist) <- rownames(Data)
	}
	else if(distmeasure=="euclidean"|distmeasure=="gower"){
		dist = cluster::daisy(Data,metric=distmeasure)
		dist = as.matrix(dist)
	}
	else if(distmeasure=="hamming"){
		dist=e1071::hamming.distance(Data)
		dist=as.matrix(dist)
	}
	else if(distmeasure=="MCA_coord"){
		Data=as.data.frame(Data)
		Data=Data[, sapply(Data, nlevels) > 1]
		MCA_result<-FactoMineR::MCA(Data,graph=FALSE)
		dist=Distance(MCA_result$ind$coord,distmeasure="euclidean",normalize=FALSE,method=NULL)
		
	}
	else if(distmeasure=="chi.squared"){
		Temp=ade4::acm.disjonctif(Data)
		dist=analogue::distance(Temp,method="chi.distance")
	}
	else if(distmeasure=="cosine"){
		Temp=lsa::cosine(t(as.matrix(Data)))
		dist=1-(1-(2*acos(Temp))/pi)
	}
	else if(distmeasure=="cont tanimoto"){
		if(normalize==TRUE){
			Data=Normalization(Data,method)				
		}		
#		S=matrix(0,nrow=nrow(Data),ncol=nrow(Data))
#		for(i in 1:nrow(Data)){
#          for(j in 1:i){
#
#				N.AB= sum(Data[i,]*Data[j,])
#				N.A = sum(Data[i,]*Data[i,])
#				N.B = sum(Data[j,]*Data[j,])
#				
#			
#				coef = N.AB / (N.A+N.B-N.AB)
#				S[i,j] = coef
#				S[j,i] = coef
#			}
#			
#		}
		#matrix manipulation (generally faster)
		m=as.matrix(Data)
		X_A_X_B=m %*% t(m)
		temp=diag(X_A_X_B)
		X_A=matrix(rep(temp,nrow(Data)),nrow=nrow(Data),ncol=nrow(Data),byrow=FALSE)
		X_B=t(X_A)
		
		Denom=X_A+X_B-X_A_X_B
		
		S=X_A_X_B/Denom
		
		dist=1-S
		
		
	}
	
	else{
		stop("Incorrect choice of distmeasure. Must be one of: tanimoto, jaccard or euclidean.")
	}
	colnames(dist)=rownames(dist)
	
	return(dist)
}


## Secondary analysis

## Visualization tools

#' @title Visualization of characteristic binary features of a single data set
#' 
#' @description  A tool to visualize characteristic binary features of a set of objects in comparison with the remaining objects for a single data set. The result is a matrix with coloured cells. Columns represent 
#' objects and rows represent the specified features. A feature which is present is give a coloured cell while an absent feature is represented by a grey cell. The labels on the right indicate the names of the features while the labels on the bottom are the names of the objects.
#' @param leadCpds A character vector with the names of the objects in a first group, i.e., the group for which the specified features are characteristic. Default is NULL.
#' @param orderLab A character vector with the order of the objects. Default is NULL.
#' @param features A character vector with the names of the features to be visualized. Default is NULL.
#' @param data The data matrix. Default is NULL.
#' @param colorLab Optional. A clustering object if the objects are to be coloured accoring to their clustering order. Default is NULL.
#' @param nrclusters Optional. The number of clusters to divide the dendrogram of ColorLab. Default is NULL.
#' @param cols Optional. A character vector with the colours of the different clusters. Default is NULL.
#' @param name A character string with the name of the data. Default is "Data".
#' @param colors1 A character vector with the colours to indicate the presence (first element) or the absence of the features for the objects in LeadCpds. Default is c('gray90','blue').
#' @param colors2 A character vector with the colours to indicate the presence (first element) or the absence of the features for the objects in the remaining objects. Default is c('gray90','green').
#' @param highlightFeat Optional. A character vector with names of features to be highlighted. The names of the features are coloured purple. Default is NULL.
#' @param margins A vector with the margings of the plot. Default is c(5.5,3.5,0.5,5.5).
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location Optional. If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' Comps=FindCluster(list(MCF7_F),nrclusters=10,select=c(1,8)) 
#' 
#' MCF7_Char=CharacteristicFeatures(List=list(fingerprintMat),Selection=Comps,
#' binData=list(fingerprintMat),datanames=c("FP"),nrclusters=NULL,topC=NULL,
#' sign=0.05,fusionsLog=TRUE,weightclust=TRUE,names=c("FP"))Feat=MCF7_Char$
#' Selection$Characteristics$FP$TopFeat$Names[c(1:10)]
#' 
#' BinFeaturesPlot_SingleData(leadCpds=Comps,orderLab=MCF7_Char$Selection$
#' objects$OrderedCpds,features=Feat,data=fingerprintMat,colorLab=NULL,
#' nrclusters=NULL,cols=NULL,name=c("FP"),colors1=c('gray90','blue'),colors2=
#' c('gray90','green'),highlightFeat=NULL,margins=c(5.5,3.5,0.5,5.5),
#' plottype="new",location=NULL)
#' }
BinFeaturesPlot_SingleData<-function(leadCpds=c(),orderLab=c(),features=c(),data=NULL,colorLab=NULL,nrclusters=NULL,cols=NULL,name=c("Data"),colors1=c('gray90','blue'),
		colors2=c('gray90','green'),highlightFeat=NULL,margins=c(5.5,3.5,0.5,5.5),plottype="new",location=NULL){
	
	if(all(leadCpds%in%rownames(data))){
		data=t(data)
	}
	
	if(!is.null(orderLab)){
		if(class(orderLab)=="character"){
			orderlabs=orderLab
		}
		else{
			orderlabs=orderLab$Clust$order.lab
			data=data[,match(orderlabs,colnames(data))]
			
		}
	}
	else{
		orderlabs=colnames(data)
	}
	
	temp=orderlabs[which(!(orderlabs%in%leadCpds))]
	AllCpds=c(leadCpds,temp)
	
	
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	
	x<-c(1:length(AllCpds)) #x=comps
	y<-c(1:length(features)) #y=feat
	PlotData<-t(data[as.character(features),AllCpds,drop=FALSE])
	plottypein(plottype,location)
	graphics::par(mar=margins)
	graphics::image(x,y,PlotData,col=colors1,xlab="",axes=FALSE,ann=FALSE,xaxt='n')
	if(length(unique(as.vector(PlotData[1:length(leadCpds),])))==1){
		if(unique(as.vector(PlotData[1:length(leadCpds),]))==1){
			PlotData[1:length(leadCpds),]=2
			colors2=c("gray90","blue","green")
			graphics::image(x,y,PlotData,col=colors2,xlab="",axes=FALSE,ann=FALSE,xaxt='n')
		}
		else{
			colors2=c("gray90","blue")
			graphics::image(x,y,PlotData,col=colors2,xlab="",axes=FALSE,ann=FALSE,xaxt='n')
		}
	}
	else{
		graphics::image(x[1:length(leadCpds)],y,PlotData[1:length(leadCpds),,drop=FALSE],col=colors2,add=TRUE,xlab="",axes=FALSE,ann=FALSE,xaxt='n')		
	}
	
	if(!(is.null(colorLab)) & !is.null(nrclusters)){
		Data1 <- colorLab$Clust
		ClustData1=stats::cutree(Data1,nrclusters) 
		
		ordercolors=ClustData1[Data1$order]
		names(ordercolors)=Data1$order.lab
		
		ClustData1=ClustData1[Data1$order]	
		
		
		order=seq(1,nrclusters)
		
		for (k in 1:length(unique(ClustData1))){
			select=which(ClustData1==unique(ClustData1)[k])
			ordercolors[select]=order[k]
		}
		
		colors<- cols[ordercolors]
		names(colors) <-names(ordercolors)	
	}
	else{
		colors1<-rep("green",length(leadCpds))
		colors2<-rep("black",length(temp))
		colors=c(colors1,colors2)
		names(colors)=AllCpds
	}
	if(!is.null(highlightFeat)){
		colfeat=rep("black",ncol(PlotData))
		colfeat[which(colnames(PlotData)%in%highlightFeat)]="purple"
		graphics::mtext(colnames(PlotData), side = 4, at= c(1:ncol(PlotData)), line=0.2, las=2,cex=0.8,col=colfeat)
	}
	else{
		graphics::mtext(colnames(PlotData), side = 4, at= c(1:ncol(PlotData)), line=0.2, las=2,cex=0.8)
	}
	graphics::mtext(name, side = 2,  line=1, las=0, cex=1)
	graphics::mtext(rownames(PlotData), side = 1, at= c(1:nrow(PlotData)), line=0.2, las=2, cex=0.8,col=colors[AllCpds])
	plottypeout(plottype)
}

#' @title Visualization of characteristic binary features of multiple data sets
#' 
#' @description A tool to visualize characteristic binary features of a set of objects in comparison with the remaining objects for multiple data sets. The result is a matrix with coloured cells. Columns represent 
#' objects and rows represent the specified features. A feature which is present is give a coloured cell while an absent feature is represented by a grey cell. The labels on the right indicate the names of the features while the labels on the bottom are the names of the objects.
#' @param leadCpds A character vector with the names of the objects in a first group, i.e., the group for which the specified features are characteristic. Default is NULL.
#' @param orderLab A character vector with the order of the objects. Default is NULL.
#' @param features A list with as elements character vectors with the names of the features to be visualized for each data set. Default is NULL.
#' @param data A list with the different data sets. Default is NULL.
#' @param validate Optional. A list with validation data sets. If a feature has a validation reference, these are added in a red colour. Default is NULL.
#' @param colorLab Optional. A clustering object if the objects are to be coloured accoring to their clustering order. Default is NULL.
#' @param nrclusters Optional. The number of clusters to divide the dendrogram of ColorLab. Default is NULL.
#' @param cols Optional. A character vector with the colours of the different clusters. Default is NULL.
#' @param name A character string with the names of the data sets. Default is c("Data1" ,"Data2") for two data sets.
#' @param colors1 A character vector with the colours to indicate the presence (first element) or the absence of the features for the objects in LeadCpds. Default is c('gray90','blue').
#' @param colors2 A character vector with the colours to indicate the presence (first element) or the absence of the features for the objects in the remaining objects. Default is c('gray90','green').
#' @param margins A vector with the margings of the plot. Default is c(5.5,3.5,0.5,5.5).
#' @param cexB The font size of the labels on the bottom: the object labels. Default is 0.80.
#' @param cexL The font size of the labels on the left: the data labels. Default is 0.80.
#' @param cexR The font size of the labels on the right: the feature labels. Default is 0.80.
#' @param spaceNames A percentage of the height of the figure to be reserved for the names of the objects. Default is 0.20.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location Optional. If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' Comps=FindCluster(list(MCF7_F),nrclusters=10,select=c(1,8)) 
#' 
#' MCF7_Char=CharacteristicFeatures(List=NULL,Selection=Comps,binData=
#' list(fingerprintMat,targetMat),datanames=c("FP","TP"),nrclusters=NULL,
#' topC=NULL,sign=0.05,fusionsLog=TRUE,weightclust=TRUE,names=c("FP","TP"))
#' 
#' FeatFP=MCF7_Char$Selection$Characteristics$FP$TopFeat$Names[c(1:10)]
#' FeatTP=MCF7_Char$Selection$Characteristics$TP$TopFeat$Names[c(1:10)]
#' 
#' BinFeaturesPlot_MultipleData(leadCpds=Comps,orderLab=MCF7_Char$Selection$
#' objects$OrderedCpds,features=list(FeatFP,FeatTP),data=list(fingerprintMat,targetMat),
#' validate=NULL,colorLab=NULL,nrclusters=NULL,cols=NULL,name=c("FP","TP"),colors1=
#' c('gray90','blue'),colors2=c('gray90','green'),margins=c(5.5,3.5,0.5,5.5),cexB=0.80,
#' cexL=0.80,cexR=0.80,spaceNames=0.20,plottype="new",location=NULL)
#' }
BinFeaturesPlot_MultipleData<-function(leadCpds,orderLab,features=list(),data=list(),validate=NULL,colorLab=NULL,nrclusters=NULL,cols=NULL,name=c("Data1" ,"Data2"),colors1=c('gray90','blue'),colors2=c('gray90','green'),margins=c(5.5,3.5,0.5,5.5),cexB=0.80,cexL=0.80,cexR=0.80,spaceNames=0.20,plottype="new",location=NULL){
	
	if(all(leadCpds%in%rownames(Data))){
		Data=t(Data)
	}
	
	if(!is.null(orderLab)){
		if(class(orderLab)=="character"){
			orderlabs=orderLab
		}
		else{
			orderlabs=orderLab$Clust$order.lab
			data=data[,match(orderlabs,colnames(data))]
		}
	}
	else{
		orderlabs=rownames(data[[1]])
	}
	
	
	
	temp=orderlabs[which(!(orderlabs%in%leadCpds))]
	AllCpds=c(leadCpds,temp)
	
	
	if(!is.null(validate)){
		for(v in 1:length(validate)){
			validate[[v]]=validate[[v]][AllCpds,]
		}
	}
	
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			dev.off()
		}
	}
	
	x<-c(1:length(AllCpds)) #x=comps

	Data_new=list()
	for(j in 1:length(data)){
		tempD=data[[j]]
		Data_temp=tempD[AllCpds,as.character(features[[j]]),drop=FALSE]
		Data_new[[j]]=as.matrix(Data_temp)
#		 if(ncol(Data_new[[j]])==1){
#			 colnames(Data_new[[j]])=Features[[j]][which(as.character(Features[[j]])%in%colnames(tempD))]
#			 rownames(Data_new[[j]])=rownames(Data[[1]])
#		 }
	}
	names(Data_new)=name
	
	Draw=matrix(0,ncol=(length(name)+2),nrow=length(unique(unlist(features))))
	i=0
	for(f in unique(unlist(features))){
		Row=c(f,rep(0,length(name)+1))
		for(j in 1:length(Data_new)){
			if(f%in%colnames(Data_new[[j]])){
				Row[j+1]=name[j]
			}
		}
		Nr=length(which((Row[c(2:(ncol(Draw)-1))])!="0"))
		Row[ncol(Draw)]=Nr
		i=i+1
		Draw[i,]=Row
	} 
	Draw=Draw[order(as.numeric(Draw[,ncol(Draw)]),decreasing=TRUE),,drop=FALSE]
	
	
	Temp_Draw=Draw[,-c(1,ncol(Draw)),drop=FALSE]
	Fin_Draw=c()
	for(a in 1:nrow(Temp_Draw[!duplicated(Temp_Draw),,drop=FALSE])){
		example=Temp_Draw[!duplicated(Temp_Draw),,drop=FALSE][a,,drop=FALSE]
		join=c()
		for(b in 1:nrow(Temp_Draw)){
			if(all(Temp_Draw[b,]==example)){
				join=c(join,b)
			}				
		}
		
		t=paste(Draw[join,1],collapse=",")
		Fin_Draw=rbind(Fin_Draw,cbind(t,paste(example,collapse=",")))
	}
	
	plottypein(plottype,location)
	Mat <- matrix(c(1:(nrow(Fin_Draw)+1)),nrow = (nrow(Fin_Draw)+1),ncol = 1,byrow = TRUE)
	H=rep(0,nrow(Fin_Draw))
	Nrtotal=length(unlist(features))
	for(m in 1:nrow(Fin_Draw)){
		Methods=unlist(strsplit(Fin_Draw[m,2],","))
		NrofMethods=length(which(Methods!="0"))
		if(NrofMethods==1){
			Nrtotal=Nrtotal+1
		}
	}
	
	for(m in 1:nrow(Fin_Draw)){
		NrofFeat=length(unlist(strsplit(Fin_Draw[m,1],",")))
		Methods=unlist(strsplit(Fin_Draw[m,2],","))
		NrofMethods=length(which(Methods!="0"))
		if(!is.null(validate)){
			NrofMethods=NrofMethods+length(validate)
		}
		if(NrofFeat*NrofMethods==1){
			NrofMethods=1.5
		}
		H[m]=((NrofFeat*NrofMethods)/Nrtotal)*(1-spaceNames)
		
	}	
	layout(mat = Mat,heights=c(H,spaceNames))
	
	for(n in 1:nrow(Fin_Draw)){
		print(n)
		Shared=unlist(strsplit(Fin_Draw[n,1],","))
		Datasets=unlist(strsplit(Fin_Draw[n,2],","))
		V=0
		if(!is.null(validate)){
			V=length(validate)
			Datasets=c(Datasets,names(validate))
		}
		
		if(any("0"%in%Datasets)){
			Datasets=Datasets[-which(Datasets==0)]
		}
		
		
		if(V>0){
			Data_new=c(Data_new,validate)
		}
		ImageData=c()
		for(d in Datasets){	
			if(all(Shared%in%colnames(Data_new[[d]]))){
				ImageData=cbind(ImageData,Data_new[[d]][,Shared,drop=FALSE])
			}
			else{
				ImageData=cbind(ImageData,Data_new[[d]][,Shared[which(Shared%in%colnames(Data_new[[d]]))],drop=FALSE])
				Fill=Shared[which(!Shared%in%colnames(Data_new[[d]]))]
				for(f in 1:length(Fill)){
					Data_new[[d]]=cbind(Data_new[[d]],rep(0,nrow(Data_new[[d]])))
					colnames(Data_new[[d]])[ncol(Data_new[[d]])]=Fill[f]
				}
			}
		}
		
		
		if(V>0){
			for(i in 1:(length(Shared)*length(validate))){
				ImageData[,(ncol(ImageData)-(i-1))][which(ImageData[,(ncol(ImageData)-(i-1))]==1)]=3
			}
		}
		
		ImageData=as.matrix(ImageData)
		ImageData=ImageData[,order(colnames(ImageData)),drop=FALSE]
		colnames(ImageData)=paste(colnames(ImageData),rep(Datasets,length(Shared)),sep="_")
		rownames(ImageData)=rownames(Data_new[[1]])
		#ImageData=ImageData[,c(ncol(ImageData):1),drop=FALSE]
		Colors=matrix(0,nrow(ImageData),ncol(ImageData))
		for(nr in 1:nrow(ImageData)){
			for(nc in 1:ncol(ImageData)){
				if(rownames(ImageData)[nr]%in%leadCpds){
					if(ImageData[nr,nc]==1){
						Colors[nr,nc]="green"
					}
					else if(ImageData[nr,nc]==3){
						Colors[nr,nc]=adjustcolor("red", alpha.f = 0.3)  
					}
					else{
						Colors[nr,nc]="grey90"
					}
				}
				else{
					if(ImageData[nr,nc]==1){
						Colors[nr,nc]="blue"
					}
					else if(ImageData[nr,nc]==3){
						Colors[nr,nc]=adjustcolor("red", alpha.f = 0.3) 
					}
					else{
						Colors[nr,nc]="grey90"
					}
				}
			}
		}
		
		par(mar=margins)
		
		
		plotrix::color2D.matplot(t(ImageData),cellcolors=t(Colors),show.values=FALSE,axes=FALSE,xlab="",ylab="",border=NA)
		
		ColorsF=rep("black",ncol(ImageData))
		names(ColorsF)=colnames(ImageData)
		
		if(!(is.null(colorLab)) & !is.null(nrclusters)){
			Data1 <- colorLab$Clust
			ClustData1=cutree(Data1,nrclusters) 
			
			ordercolors=ClustData1[Data1$order]
			names(ordercolors)=Data1$order.lab
			
			ClustData1=ClustData1[Data1$order]	
			
			
			order=seq(1,nrclusters)
			
			for (k in 1:length(unique(ClustData1))){
				select=which(ClustData1==unique(ClustData1)[k])
				ordercolors[select]=order[k]
			}
			
			colors<- cols[ordercolors]
			names(colors) <-names(ordercolors)	
		}
		else{
			Colors1<-rep("green",length(leadCpds))
			Colors2<-rep("black",length(temp))
			Colors=c(Colors1,Colors2)
			names(Colors)=AllCpds
			
		}
		
		mtext(colnames(ImageData), side = 4, at= c(ncol(ImageData):1), line=0.2, las=2,cex=cexR,col=ColorsF)
		if(length(Datasets)==length(Data_new)){
			mtext("All", side = 2,  line=1, las=1, cex=2)
		}	 
		else{
			if(length(Datasets)>2){
				N=c()
				for(w in 1:length(Datasets)){
					if(w%%2==0&w!=length(Datasets)){
						N=c(N," & ",Datasets[w],"\n")
					}
					else if(w==1){
						
						N=c(N,Datasets[w])
					}
					else{
						N=c(N," & ",Datasets[w])
					}
				}
				mtext(paste(N,collapse=""), side = 2,  line=1, las=1, cex=cexL)
			}
			else{
				mtext(paste(Datasets,collapse=" & "), side = 2,  line=1, las=1, cex=cexL)
			}
		}
		
		if(n==nrow(Fin_Draw)){
			mtext(rownames(ImageData), side = 1, at= c(1:nrow(ImageData)), line=0.2, las=2, cex=cexB,col=Colors[AllCpds])
		}
		
	}
	
	plottypeout(plottype)
}

#' @title Comparison of clustering results over multiple results
#' 
#' @description A visual comparison of all methods is handy to see which objects will
#' always cluster together independent of the applied methods. To this aid the
#' function \code{ComparePlot} has been written. The function relies on methods
#' of the \code{circlize} package.
#' @export ComparePlot
#' @details This function makes use of the functions \code{ReorderToReference} and
#' \code{Colorsnames}.  Given a list with the outputs of several methods, the
#' first step is to call upon \code{ ReorderToReference} and to produce a
#' matrix of which the columns are ordered according to the ordering of the
#' objects of the first method in the list. Each cell represent the number of
#' the cluster the object belongs to for a specific method indicated by the
#' rows. The clusters are arranged in such a way that these correspond to that
#' one cluster of the referenced method that they have the most in common with.
#' The function \code{color2D.matplot} produces a plot of this matrix but needs
#' a vector indicating the names of the colors to be used. This is where
#' \code{ColorsNames} comes in. A vector of the color names of the output of
#' the \code{ReorderToReference} is created and handed to
#' \code{color2D.matplot}. It is optional to adjust the margins of the plot and
#' to give a vector with the names of the methods which will be used as labels
#' for the rows in the plot. The labels for the columns are the names of the
#' object in the order of clustering of the referenced method. Further, the
#' similarity measures of the methods compared to the reference will be
#' computed and shown on the right side of the plot.
#' @param List A list of the outputs from the methods to be compared. The first
#' element of the list will be used as the reference in
#' \code{ReorderToReference}.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param cols A character vector with the colours to be used. Default is NULL.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' @param names Optional. Names of the methods to be used as labels for the
#' columns. Default is NULL.
#' @param margins Optional. Margins to be used for the plot. Default is c(8.1,3.1,3.1,4.1).
#' @param circle Logical. Whether the figure should be circular (TRUE) or a rectangle (FALSE). Default is FALSE.
#' @param canvaslims The limits for the circular dendrogam. Default is c(-1.0,1.0,-1.0,1.0).
#' @param Highlight Optional. A list of character vectors of objects to be highlighted. The names of the elements in the list are the names to appear on the figure. The median similarities of the objects in each list elemented is computed. Default is NULL.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location Optional. If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return A plot which translates the matrix output of the function
#' \code{ReorderToReference} in which the columns represent the objects in the
#' ordering the referenced method and the rows the outputs of the given
#' methods. Each cluster is given a distinct color. This way it can be easily
#' observed which objects will cluster together. The labels on the right side
#' of the plot are the similarity measures computed by
#' \code{SimilarityMeasure}.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(Colors1)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' L=list(MCF7_F,MCF7_T)
#' N=c("FP","TP")
#' 
#' #rectangular
#' ComparePlot(List=L,nrclusters=7,cols=Colors1,fusionsLog=TRUE,weightclust=TRUE,
#' names=N,margins=c(9.1,4.1,4.1,4.1),plottype="new",location=NULL)
#' 
#' #circle
#' Comps_I=c("fluphenazine","trifluoperazine","prochlorperazine","chlorpromazine")  
#' Comps_II=c("butein","genistein","resveratrol")
#' 
#' ComparePlot(List=L,nrclusters=7,cols=c(Colors1), fusionsLog=TRUE,weightclust=FALSE,
#' names =N, margins = c(8.1, 3.1,3.1, 4.1),circle=TRUE,canvaslims=c(-1.1,1.1,-1.1,1.1),
#' Highlight=list("Comps I"=Comps_I,"Comps II"=Comps_II,"Cancer Treatments"=c("estradiol",
#' "fulvestrant")),plottype = "new")
#' }
ComparePlot<-function(List,nrclusters=NULL,cols=NULL,fusionsLog=FALSE,weightclust=FALSE,names=NULL,margins=c(8.1,3.1,3.1,4.1),circle=FALSE,canvaslims=c(-1.0,1.0,-1.0,1.0),Highlight=NULL,plottype="new",location=NULL){	
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	
	for(i in 1:length(List)){
		if(attributes(List[[i]])$method == "Weighted" & weightclust==TRUE){
			T=List[[i]]$Clust
			attr(T,"method")="Single Clustering"
			List[[i]]=T
		}
	}
	
	MatrixColors=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)
	#capture singletons
	for(i in 1:nrow(MatrixColors)){
		if(any(table(MatrixColors[i,])==1)){
			clusters=names(table(MatrixColors[i,]))[which(table(MatrixColors[i,])==1)]
			MatrixColors[i,][which(MatrixColors[i,]%in%as.numeric(clusters))]=5000
		}
	}
	
	Names=ColorsNames(MatrixColors,cols)
	colnames(Names)=colnames(MatrixColors)
	nobs=dim(MatrixColors)[2]
	nmethods=dim(MatrixColors)[1]
	
	if(is.null(names)){
		for(j in 1:nmethods){
			names[j]=paste("Method",j,sep=" ")	
		}
	}
	
	if(circle){
		plottypein(plottype, location)
		
		circlize::circos.initialize(factors =c(1:ncol(MatrixColors)) , xlim = c(0, ncol(MatrixColors)))
		for(i in 1:length(List)){
			circlize::circos.trackPlotRegion(factors = c(1:ncol(MatrixColors)), ylim = c(0,1),track.height = 0.05)
		}
		track=c(length(List):1)
		for(i in 1:nrow(MatrixColors)){
			for(j in 1:ncol(MatrixColors)){
				circlize::highlight.sector(sector.index=j, track.index = track[i],col=Names[i,j])
			}
		}
		hc=List[[1]]$Clust
		labels=substr(hc$order.lab,1,5)
		
		if(!is.null(Highlight)){
			for(h in 1:length(Highlight)){
				Name=names(Highlight)[h]
				HL=which(labels%in%substr(Highlight[[h]],1,5))
				
				Sims=c()
				for(i in 1:length(List)){
					Values=List[[i]]$DistM[lower.tri(List[[i]]$DistM)]
					Sims=c(Sims,as.numeric(1-List[[i]]$DistM[Highlight[[h]],Highlight[[h]]][lower.tri(List[[i]]$DistM[Highlight[[h]],Highlight[[h]]])]))
				}
				MedSim=round(median(Sims),2)
				
				circlize::draw.sector(get.cell.meta.data("cell.start.degree", sector.index = min(HL)),
						get.cell.meta.data("cell.end.degree", sector.index = max(HL)),
						rou1 = 1, col = "#00000020")
				
				circlize::highlight.sector(sector.index=c(min(HL):max(HL)), track.index = 1, text = paste(Name,": ",MedSim,sep=""),
						facing = "bending.inside", niceFacing = TRUE, text.vjust = -1.5)
				
			}
		}
		circlize::circos.clear()
		par(new = TRUE)
		hc=List[[1]]$Clust
		max_height=max(hc$height)
		dend=as.dendrogram(hc)
		labels=substr(hc$order.lab,1,5)
		#ct=cutree(dend,6)
		circlize::circos.par("canvas.xlim" = c(canvaslims[1], canvaslims[2]), "canvas.ylim" = c(canvaslims[3], canvaslims[4]))
		circlize::circos.initialize(factors =1 , xlim = c(0, ncol(MatrixColors)))
		circlize::circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.4,bg.border=NA,
				panel.fun = function(x, y) {
					for(i in seq_len(ncol(MatrixColors))) {
						circlize::circos.text(i-0.5, 0, labels[i], adj = c(0, 0.5),
								facing = "clockwise", niceFacing = TRUE,
								col = Names[1,colnames(MatrixColors)[i]], 
								cex = 0.8,font=2)
					}
				})
		circlize::circos.trackPlotRegion(ylim = c(0, max_height), bg.border = NA,track.height = 0.4, panel.fun = function(x, y) {
					circos.dendrogram(dend, max_height = max_height)})
		circlize::circos.clear()
		
		
		
		
		plottypeout(plottype)
		
		
	}
	
	else{
		#similar=round(SimilarityMeasure(MatrixColors),2)
		plottypein(plottype,location)
		graphics::par(mar=margins)
		plotrix::color2D.matplot(MatrixColors,cellcolors=Names,show.values=FALSE,axes=FALSE,xlab="",ylab="")
		graphics::axis(1,at=seq(0.5,(nobs-0.5)),labels=colnames(MatrixColors),las=2,cex.axis=1.5)
		graphics::axis(2,at=seq(0.5,(nmethods-0.5)),labels=rev(names),cex.axis=1.5,las=2)
		#axis(4,at=seq(0.5,(nmethods-0.5)),labels=rev(similar),cex.axis=0.65,las=2)
		plottypeout(plottype)
	}	
}

#' @title Function that annotates colors to their names
#' 
#' @description The \code{ColorsNames} function is used on the output of the
#' \code{ReorderToReference} and matches the cluster numbers indicated by the
#' cell with the names of the colors.  This is necessary to produce the plot of
#' the \code{ComparePlot} function and is therefore an internal function of
#' this function but can also be applied separately.
#' @export ColorsNames
#' 
#' @param matrixColors The output of the ReorderToReference function.
#' @param cols A character vector with the names of the colours to be used. Default is NULL.
#' @return A matrix containing the hex code of the color that corresponds to
#' each cell of the matrix to be colored. This function is called upon by the
#' \code{ComparePlot} function.
#' @examples
#' 
#' data(fingerprintMat)
#' data(targetMat)
#' data(Colors1)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' L=list(MCF7_F,MCF7_T)
#' names=c("FP","TP")
#' 
#' MatrixColors=ReorderToReference(List=L,nrclusters=7,fusionsLog=TRUE,weightclust=TRUE,
#' names=names)
#' 
#' Names=ColorsNames(matrixColors=MatrixColors,cols=Colors1)
ColorsNames<-function(matrixColors,cols=NULL){
	Names=matrix(0,nrow=dim(matrixColors)[1],ncol=dim(matrixColors)[2])
	for (i in 1:dim(matrixColors)[1]){
		for(j in 1:dim(matrixColors)[2]){
			if(is.na(matrixColors[i,j])){
				Names[i,j]="grey"
				
			}
			else if(length(unique(matrixColors[i,]))==1){
				Names[i,j]="grey"
				
			}
			else if(matrixColors[i,j]==5000){
				Names[i,j]="white"
				
			}
			else{
				temp=matrixColors[i,j]	
				Color=cols[temp]
				Names[i,j]=Color
			}
		}
	}
	return(Names)
}


#' @title Order the outputs of the clustering methods against a reference
#' 
#' @description When multiple methods are performed on a data set, it is interesting to
#' compare their results. However, a comparison is not easily done since
#' different methods leads to a different ordering of the objects. The
#' \code{ReorderToReference} rearranges the cluster to a reference method.
#' @export ReorderToReference
#' @details It is interesting to compare the results of the methods described in the
#' methodology. All methods result in a dendrogram which is cut into a specific
#' number of clusters with the \code{cutree} function. This results in an
#' numbering of cluster based on the ordering of the names in the data and not
#' on the order in which they are grouped into clusters. However, different
#' methods lead to different clusters and it is possible that cluster $1$ of
#' one method will not be the cluster that has the most in common with cluster
#' 1 of another method. This makes comparisons rather difficult. Therefore the
#' ReorderToReference function was written which takes one method as a
#' reference and rearranges the cluster numbers of the other methods to this
#' reference such that clusters are appointed to that cluster they have the
#' most in common with. The result of this function is a matrix of which the
#' columns are in the order of the clustering of the objects of the
#' referenced method and the rows represent the methods. Each cell contains the
#' number of the cluster the compound is in for that method compared to the
#' method used as a reference. This function is applied in the functions
#' \code{SimilarityMeasure}, \code{DiffGenes}, \code{Pathways} and
#' \code{ComparePlot}. It is a possibility that 2 or more clusters are fused
#' together compared to the reference method. If this is true, the function
#' will alert the user and will ask to put the parameter fusionsLog to true.
#' Since \code{ReorderToReference} is often used as an internal function, also
#' for visualization, it will print out how many more colors should be
#' specified for those clusters that did not find a suitable match. This can be
#' due to fusion or complete segregation of its objects into other clusters.
#' 
#' @param List A list of clustering outputs to be compared. The first element
#' of the list will be used as the reference.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param fusionsLog Logical. Indicator for the fusion of clusters. Default is FALSE.
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' @param names Optional. A character vector with the names of the methods.
#' @return A matrix of which the cells indicate to what cluster the objects
#' belong to according to the rearranged methods.
#' @note The \code{ReorderToReference} function was optimized for the
#' situations presented by the data sets at hand. It is noted that the function
#' might fail in a particular situation which results in a infinite loop.
#' @examples
#' 
#' data(fingerprintMat)
#' data(targetMat)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_ADC=ADC(list(fingerprintMat,targetMat),distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible")
#' 
#' L=list(MCF7_F,MCF7_ADC,MCF7_T)
#' names=c("FP","ADC","TP")
#' 
#' MCF7_Matrix=ReorderToReference(List=L,nrclusters = 7, fusionsLog = FALSE, weightclust = 
#' FALSE, names = names)
ReorderToReference<-function(List,nrclusters=NULL,fusionsLog=FALSE,weightclust=FALSE,names=NULL){
	
	matequal <- function(x, y)
		is.matrix(x) && is.matrix(y) && dim(x) == dim(y) && all(x == y)
	
	ListNew=list()
	element=0
	for(i in 1:length(List)){
		if(attributes(List[[i]])$method != "CEC" & attributes(List[[i]])$method != "Weighted" & attributes(List[[i]])$method!= "WeightedSim"){
			ResultsClust=list()
			ResultsClust[[1]]=list()
			ResultsClust[[1]][[1]]=List[[i]]
			names(ResultsClust[[1]])[1]="Clust"
			element=element+1					
			ListNew[[element]]=ResultsClust[[1]]
			#attr(ListNew[element],"method")="Weights"
		}
		else if(attributes(List[[i]])$method=="CEC" | attributes(List[[i]])$method=="Weighted" | attributes(List[[i]])$method == "WeightedSim"){
			ResultsClust=list()
			if(weightclust==TRUE){
				ResultsClust[[1]]=list()
				if(attributes(List[[i]])$method != "WeightedSim"){
					
					ResultsClust[[1]][[1]]=List[[i]]$Clust
					names(ResultsClust[[1]])[1]="Clust"
					attr(ResultsClust[[1]]$Clust,"method")="Weights"	
					element=element+1					
					ListNew[[element]]=ResultsClust[[1]]
					
				}
				else{
					ResultsClust[[1]]=list()
					ResultsClust[[1]][[1]]=List[[i]]
					names(ResultsClust[[1]])[1]="Clust"
					attr(ResultsClust[[1]]$Clust,"method")="Weights"	
					element=element+1					
					ListNew[[element]]=ResultsClust[[1]]
				}
			}
			else{
				for (j in 1:length(List[[i]]$Results)){
					ResultsClust[[j]]=list()
					ResultsClust[[j]][[1]]=List[[i]]$Results[[j]]
					names(ResultsClust[[j]])[1]="Clust"
					attr(ResultsClust[[j]]$Clust,"method")="Weights"	
					element=element+1					
					ListNew[[element]]=ResultsClust[[j]]
					#attr(ListNew[[element]],"method")="Weights"
				}		
			}		
		}	
	}
	
	if(is.null(names)){
		names=seq(1,length(ListNew),1)
		for(i in 1:length(ListNew)){
			names[i]=paste("Method",i,sep=" ")
		}
	}
	names(ListNew)=names
	List=ListNew
	
	Clusters=list()
	CutTree<-function(i,Data,nrclusters){
		print(i)
		if(attributes(Data$Clust)$method == "Ensemble"){
			Clusters=Data$Clust$Clust$Clusters
			names(Clusters)=NULL
		}
		else{
			Clusters=stats::cutree(Data$Clust$Clust,k=nrclusters)
		}
		return(Clusters)
	}
	Clusters=lapply(seq(1,length(List)),function(i) CutTree(i,Data=ListNew[[i]],nrclusters=nrclusters))
	
	
	xaxis=List[[1]]$Clust$Clust$order #order of the objects as for method 1.
	xaxis.names=List[[1]]$Clust$Clust$order.lab #might be that names of methods are not in the same order...
	
	ordercolors=Clusters[[1]][xaxis]
	order=seq(1,nrclusters)
	
	for (k in 1:length(unique(Clusters[[1]][xaxis]))){
		select=which(Clusters[[1]][xaxis]==unique(Clusters[[1]][xaxis])[k])
		ordercolors[select]=order[k]
	}
	
	cols=unique(ordercolors) #order of the colors as depicted by method 1
	
	Ordered=list()
	
	autograph=list()
	for(i in cols){
		autograph[[i]]=xaxis[which(ordercolors==i)]	
	}
	
	#for(j in 1:length(List)){		
	#		temp=Clusters[[j]][xaxis]  #put clusternumbers of the other method into the same order as those of method (1)
	#	clusternumbers=temp		   #problem:cutree is based on the ordering of the names as they are in the rownames not in the order of joined objects 
	#	for(k in 1:length(cols)){
	#		change=which(temp==unique(temp)[k])
	#		clusternumbers[change]=cols[which(cols==unique(temp)[k])]
	#	}
	#	Ordered[[j]]=clusternumbers
	#}
	
	for (j in 1:length(List)){
		message(j)
		#ordercolorsj=Clusters[[j]][xaxis]
		if(attributes(List[[j]]$Clust)$method=="Ensemble"){
			DistM=matrix(0,ncol=length(List[[j]]$Clust$Clust$order),nrow=length(List[[j]]$Clust$Clust$order))
			if(is.null(names(List[[j]]$Clust$Clust$Clusters))){
				names(List[[j]]$Clust$Clust$Clusters)=paste("Comp", c(1:nrow(DistM)),sep=" ")
			}
			colnames(DistM)=names(List[[j]]$Clust$Clust$Clusters)
			rownames(DistM)=names(List[[j]]$Clust$Clust$Clusters)
			List[[j]]$Clust$DistM=DistM
		}
		ordercolorsj=Clusters[[j]][match(xaxis.names,rownames(List[[j]]$Clust$DistM))]
		order=seq(1,nrclusters)
		
		for (k in 1:length(unique(Clusters[[j]][match(xaxis.names,rownames(List[[j]]$Clust$DistM))]))){
			select=which(Clusters[[j]][match(xaxis.names,rownames(List[[j]]$Clust$DistM))]==unique(Clusters[[j]][match(xaxis.names,rownames(List[[j]]$Clust$DistM))])[k])
			ordercolorsj[select]=order[k]
		}
		
		
		temp2=ordercolorsj
		#temp3=xaxis
		temp3=match(xaxis.names,rownames(List[[j]]$Clust$DistM))
		fan=list()
		for(i in cols){
			fan[[i]]=match(xaxis.names,rownames(List[[j]]$Clust$DistM))[which(temp2==i)]	
		}
		
		favors=matrix(0,length(autograph),length(fan))
		rownames(favors)=seq(1,length(autograph))
		colnames(favors)=seq(1,length(fan))
		
		for(a in 1:length(autograph)){
			for (b in 1:length(fan)){
				favorab=length(which(rownames(List[[j]]$Clust$DistM)[fan[[b]]] %in% rownames(List[[1]]$Clust$DistM)[autograph[[a]]]))/length(autograph[[a]])	
				favors[a,b]=favorab	
			}
		}
		
		#See function woman and men CB (put back what has value replaced)
		
		tempfavors=favors
		
		matched=c(rep("Free",nrclusters))
		proposed=c(rep("No",nrclusters))
		Switches=c(rep("Open",nrclusters))
		
		proposals=matrix(0,length(autograph),length(fan))
		
		#First match does "fans" that only have 1 element in their column: only one choice
		for(a in 1:dim(tempfavors)[1]){
			for (b in 1:dim(tempfavors)[2]){
				if(favors[a,b]==1){
					matched[a]=b
					proposed[b]="Yes"
					proposals[a,b]=1
					col=a
					
					change=which(xaxis.names %in% rownames(List[[j]]$Clust$DistM)[fan[[b]]])
					temp3[change]=col
					
					tempfavors[,b]=0
					tempfavors[a,]=0
					
					Switches[a]="Closed"
				}
				
			}
		}	
		
		
		#OneLeftC=FALSE
		#OneLeftR=FALSE
		for(b in 1:dim(tempfavors)[2]){
			if(length(which(tempfavors[,b]!=0))==1){
				match=which(tempfavors[,b]!=0)
				test=which(tempfavors[match,]==max(tempfavors[match,]))[1]
				if(length(which(tempfavors[,test]!=0))!=1 | b %in% which(tempfavors[match,]==max(tempfavors[match,])) ){
					matched[match]=b
					proposed[b]="Yes"
					proposals[match,b]=1
					col=match
					
					change=which(xaxis.names %in% rownames(List[[j]]$Clust$DistM)[fan[[b]]])
					temp3[change]=col
					
					tempfavors[,b]=0
					tempfavors[match,]=0
					
					Switches[match]="Closed"
				}
			}		
			#Unneccesary? 
			#if(length(which(tempfavors[,b]==1))>1){
			#	matches=which(tempfavors[,b]==1)
			#	matched[matches]=b
			#	proposed[b]="Yes"
			#	proposals[matches,b]=1
			#	col=matches[1]
			#	
			#	change=which(xaxis %in% fan[[b]])
			#	temp3[change]=col
			#	
			#	tempfavors[,b]=0
			#	tempfavors[matches,]=0
			#	
			#	Switches[matches]="Closed"
			#
			#	OneLeftC=TRUE
			#}
		}
		
		for(a in 1:dim(tempfavors)[1]){
			if(length(which(tempfavors[a,]!=0))==1){
				propose=which(tempfavors[a,]!=0)
				test=which(tempfavors[,propose]==max(tempfavors[,propose]))[1]
				if(length(which(tempfavors[test,]!=0))!=1 | a %in% which(tempfavors[,propose]==max(tempfavors[,propose]))){
					
					matched[a]=propose
					proposed[propose]="Yes"
					proposals[a,propose]=1
					col=a
					
					change=which(xaxis.names %in% rownames(List[[j]]$Clust$DistM)[fan[[propose]]])
					temp3[change]=col
					
					tempfavors[a,]=0
					tempfavors[,propose]=0
					
					Switches[a]="Closed"
				}
			}
			#Unnecessary?
			#if(length(which(tempfavors[a,]==1))>1){
			#	proposes=which(tempfavors[a,]==1)
			#	matched[a]="Left"
			#	proposed[proposes]="Yes"
			#	proposals[a,proposes]=1
			#	col=a
			#	
			#	change=which(xaxis %in% fan[[proposes]])
			#	temp3[change]=col
			#	
			#	tempfavors[a,]=0
			#	tempfavors[,proposes]=0
			#	
			#	Switches[a]="Closed"
			#	
			#	OneLeftR=TRUE
			#}
		}
		Continue=TRUE
		if(length(which(matched=="Free")) == 0){
			Continue=FALSE
		}	
		
		while(length(which(matched=="Free")) != 0 | !(matequal(proposals[which(matched=="Free"),], matrix(1, length(which(matched=="Free")), nrclusters))) | Continue!=FALSE){
			#for(a in which(matched=="Free")){
			#if(length(which(tempfavors[a,]!=0))==1){
			#	propose=which.max(tempfavors[a,])
			#
			#	matched[a]=propose
			#	proposed[propose]="Yes"
			#	proposals[a,propose]=1
			#	col=a
			#	
			#	change=which(xaxis %in% fan[[propose]])
			#	temp3[change]=col
			#	
			#	tempfavors[a,propose]=0
			#}
			
			#else{
			a=which(matched=="Free")[1]
			propose=which.max(tempfavors[a,])
			if(tempfavors[a,propose]==0){
				if(length(which(matched=="Free"))==1){
					Continue=FALSE
				}
				matched[a]="Left"
			}	
			else{
				if(proposed[propose]=="No"){
					proposed[propose]="Yes"
					matched[a]=propose
					proposals[a,propose]=1
					col=a
					
					change=which(xaxis.names %in% rownames(List[[j]]$Clust$DistM)[fan[[propose]]])
					temp3[change]=col
					
					tempfavors[a,propose]=0
					
					if(length(which(tempfavors[a,]==0))==dim(tempfavors)[2]){
						Switches[a]="Closed"
						tempfavors[,propose]=0
						
						c=1
						while(c < a){
							if(Switches[c] != "Closed" & length(which(tempfavors[c,]==0))==dim(tempfavors)[2]){
								Switches[c]="Closed"
								if(matched[c]=="Left"){
									tempfavors[c,]=0
								}
								else{
									tempfavors[,matched[c]]=0
								}
								c=1								
							}
							else{ 
								c=c+1
							}
						}
					}
				}
				else if(proposed[propose]=="Yes"){
					if(favors[a,propose] > max(favors[which(matched==propose),propose]) & Switches[which(matched==propose)]=="Open"){
						
						
						#first undo then replace
						#tempfavors[which(matched==propose),propose]=favors[which(matched==propose),propose]
						
						changeback=which(xaxis.names %in%  rownames(List[[j]]$Clust$DistM)[fan[[propose]]])
						temp3[changeback]=match(xaxis.names,rownames(List[[j]]$Clust$DistM))[changeback]
						matched[which(matched==propose)]="Free"
						
						matched[a]=propose
						proposals[a,propose]=1
						col=a
						change=which(xaxis.names %in% rownames(List[[j]]$Clust$DistM)[fan[[propose]]])
						temp3[change]=col
						
						tempfavors[a,propose]=0
					}
					else if(length(which(tempfavors[a,]!=0))==1){
						#if only 1 remains, these MUST BE matched
						changeback=which(xaxis.names %in% rownames(List[[j]]$Clust$DistM)[fan[[propose]]])
						temp3[changeback]=match(xaxis.names,rownames(List[[j]]$Clust$DistM))[changeback]
						matched[which(matched==propose)]="Free"
						
						matched[a]=propose
						proposals[a,propose]=1
						col=a
						
						change=which(xaxis.names %in% rownames(List[[j]]$Clust$DistM)[fan[[propose]]])
						temp3[change]=col
						
						tempfavors[a,propose]=0	
						
					}
					else{							
						proposals[a,propose]=1
						tempfavors[a,propose]=0	
					}
					
				}	
			}	
			if(length(which(matched=="Free"))==0){
				Continue=FALSE
			}
		}	
		fusions=0
		for( i in unique(matched)){
			if(length(which(!(seq(1,nrclusters) %in% matched)))>=1){
				fusions=length(which(!(seq(1,nrclusters) %in% matched)))
			}
		}
		
		if(fusions != 0 & fusionsLog==FALSE){
			message(paste("specify",fusions,"more color(s) and put fusionsLog equal to TRUE",sep=" "))
		}
		premiumcol=c()
		for (i in 1:(fusions)){
			premiumcol=c(premiumcol,length(matched)+i)
		}
		
		if((length(which(matched=="Left"))!=0) | (length(which(proposed=="No"))!=0)){						
			if(length(which(proposed=="No"))!=0){
				for(i in 1:length(which(proposed=="No"))){
					Left=which(proposed=="No")[1]
					maxLeft=which(favors[,Left]==max(favors[,Left]))
					
					proposed[Left]="Yes"
					proposals[maxLeft,Left]=1
					col=premiumcol[i]
					
					change=which(xaxis.names %in% rownames(List[[j]]$Clust$DistM)[fan[[Left]]])
					temp3[change]=col
					
					tempfavors[,Left]=0
					tempfavors[maxLeft,]=0
				}	
				
			}
			if(length(which(matched=="Left"))!=0){
				for (i in 1:length(which(matched=="Left")))
					Left=which(matched=="Left")[1]
				message(paste("Cluster",Left,"of the reference has found no suitable match.",sep=" "))
				#maxLeft=which(favors[Left,]==max(favors[Left,]))
				
			}					
		} 
		
		Ordered[[j]]=temp3
	}
	
	Matrix=c()
	for(j in 1:length(Ordered)){
		Matrix=rbind(Matrix,Ordered[[j]])		
	}
	colnames(Matrix)=List[[1]]$Clust$Clust$order.lab
	rownames(Matrix)=names
	return(Matrix)
	
}

#' @title Box plots of one distance matrix categorized against another distance
#' matrix.
#' 
#' @description Given two distance matrices, the function categorizes one distance matrix
#' and produces a box plot from the other distance matrix against the created
#' categories. The option is available to choose one of the plots or to have
#' both plots. The function also works on outputs from ADEC and CEC functions
#' which do not have distance matrices but incidence matrices.
#' 
#' 
#' @param Data1 The first data matrix, cluster outcome or distance matrix to be
#' plotted.
#' @param Data2 The second data matrix, cluster outcome or distance matrix to
#' be plotted.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clusters".	
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}. 
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param lab1 The label to plot for Data1.
#' @param lab2 The label to plot for Data2.
#' @param limits1 The limits for the categories of Data1.
#' @param limits2 The limits for the categories of Data2.
#' @param plot The type of plots: 1 - Plot the values of Data1 versus the
#' categories of Data2. 2 - Plot the values of Data2 versus the categories of
#' Data1. 3 - Plot both types 1 and 2.
#' @param StopRange Logical. Indicates whether the distance matrices with
#' values not between zero and one should be standardized to have so. If FALSE
#' the range normalization is performed. See \code{Normalization}. If TRUE, the
#' distance matrices are not changed. This is recommended if different types of
#' data are used such that these are comparable. Default is FALSE.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return One/multiple box plots.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' 
#' BoxPlotDistance(MCF7_F,MCF7_T,type="cluster",lab1="FP",lab2="TP",limits1=c(0.3,0.7),
#' limits2=c(0.3,0.7),plot=1,StopRange=FALSE,plottype="new", location=NULL)
#' }
#' @export BoxPlotDistance
BoxPlotDistance<-function(Data1,Data2,type=c('data','dist','clusters'),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),lab1,lab2,limits1=NULL,limits2=NULL,plot=1,StopRange=FALSE,plottype="new",location=NULL){
	type<-match.arg(type)
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	
	C1=D1=C2=D2=NULL
	if(type=='clusters'){
		
		Dist1<-Data1$DistM
		Dist2<-Data2$DistM
		
		
	}
	else if(type=='data'){
		Dist1<-Distance(Data1,distmeasure[1],normalize[1],method[1])
		Dist2<-Distance(Data2,distmeasure[2],normalize[2],method[2])
		DistL=list(Dist1,Dist2)
		for(i in 1:2){
			if(StopRange==FALSE & !(0<=min(DistL[[i]]) & max(DistL[[i]])<=1)){
				message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
				DistL[[i]]=Normalization(DistL[[i]],method="Range")
			}
		}
		
	}
	else if(type=='dist'){
		Dist1=Data1
		Dist2=Data2	
		DistL=list(Dist1,Dist2)
		for(i in 1:2){
			if(StopRange==FALSE &  !(0<=min(DistL[[i]]) & max(DistL[[i]])<=1)){
				message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
				DistL[[i]]=Normalization(DistL[[i]],method="Range")
			}
		}
	}
	
	OrderNames=rownames(Dist1)
	Dist2=Dist2[OrderNames,OrderNames]
	
	Dist1lower <- Dist1[lower.tri(Dist1)]
	Dist2lower <- Dist2[lower.tri(Dist2)]
	
	Categorize<-function(Distlower,limits){
		Cat=c(rep(0,length(Distlower)))
		for(j in 1:(length(limits)+1)){
			if(j==1){
				Cat[Distlower<=limits[j]]=j
			}
			else if(j<=length(limits)){
				Cat[Distlower>limits[j-1] & Distlower<=limits[j]]=j
			}	
			else{
				Cat[Distlower>limits[j-1]]=j
			}
		}
		Cat<-factor(Cat)
		return(Cat)
		
	}
	
	#plot2
	if(!(is.null(limits1))){
		Dist1cat<-Categorize(Dist1lower,limits1)
		
		dataBox2<-data.frame(D2=Dist2lower,C1=Dist1cat)
		p2<-ggplot2::ggplot(dataBox2,ggplot2::aes(factor(C1),D2)) #x,y
		p2<-p2+ggplot2::geom_boxplot(outlier.shape=NA)+ggplot2::geom_point(color="blue",size=2,shape=19,position="jitter",cex=1.5)+ggplot2::xlab(lab1)+ggplot2::ylab(lab2)
	}
	#plot1
	if(!(is.null(limits2))){
		Dist2cat<-Categorize(Dist2lower,limits2)
		dataBox1<-data.frame(D1=Dist1lower,C2=Dist2cat)
		p1<-ggplot2::ggplot(dataBox1,ggplot2::aes(factor(C2),D1)) #x,y
		p1<-p1+ggplot2::geom_boxplot(outlier.shape=NA)+ggplot2::geom_point(color="blue",size=2,shape=19,position="jitter",cex=1.5)+ggplot2::xlab(lab2)+ggplot2::ylab(lab1)
		
	}	
	if(plot==3){
		if(plottype=="pdf"){
			location=paste(location,'_type3.pdf',sep="")
		}
		plottypein(plottype,location)
		gridExtra::grid.arrange(p1, p2, ncol=2,nrow=1)		
		
	}
	else if(plot==1){
		if(plottype=="pdf"){
			location=paste(location,'_type1.pdf',sep="")
		}
		plottypein(plottype,location)
		print(p2)
		
	}
	else if(plot==2){
		if(plottype=="pdf"){
			location=paste(location,'_type2.pdf',sep="")
		}
		plottypein(plottype,location)
		print(p2)		
	}
	
}


#' @title Determining the characteristic features of a cluster
#' 
#' @description The function \code{CharacteristicFeatures} requires as input a list of one
#' or multiple clustering results. It is capable of selecting the binary
#' features which determine a cluster with the help of the fisher's exact test.
#' @export CharacteristicFeatures
#' @details The function rearranges the clusters of the methods to a reference method
#' such that a comparison is made easier.  Given a list of methods, it calls
#' upon \code{ReorderToReference} to rearrange the number of clusters according
#' to the first element of the list which will be used as the reference.
#' 
#' @param List A list of the clustering outputs to be compared. The first
#' element of the list will be used as the reference in
#' \code{ReorderToReference}.
#' @param Selection If differential gene expression should be investigated for
#' a specific selection of objects, this selection can be provided here.
#' Selection can be of the type "character" (names of the objects) or
#' "numeric" (the number of specific cluster). Default is NULL.
#' @param binData A list of the binary feature data matrices. These will be
#' evaluated with the fisher's extact test. Default is NULL.
#' @param contData A list of continuous data sets of the objects. These will
#' be evaluated with the t-test. Default is NULL.
#' @param datanames A vector with the names of the data matrices. Default is NULL.
#' @param nrclusters Optional. The number of clusters to cut the dendrogram in.
#' The number of clusters should not be specified if the interest lies only in
#' a specific selection of objects which is known by name.  Otherwise, it is
#' required. Default is NULL.
#' @param sign The significance level to be handled. Default is 0.05.
#' @param topChar Overrules sign. The number of features to display for each
#' cluster.  If not specified, only the significant genes are shown. Default is NULL.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' @param names Optional. Names of the methods. Default is NULL.
#' @return The returned value is a list with an element per method. Each
#' element contains a list per cluster with the following elements:
#' \item{objects}{A list with the elements LeadCpds (the objects of
#' interest) and OrderedCpds (all objects in the order of the clustering
#' result)} \item{Characteristics}{A list with an element per defined binary
#' data matrix in BinData and continuous data in ContData. Each element is
#' again a list with the elements TopFeat (a table with information on the top
#' features) and AllFeat (a table with information on all features)}
#' @examples
#' 
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(geneMat)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' L=list(MCF7_T ,MCF7_F)
#' 
#' MCF7_Char=CharacteristicFeatures(List=L,Selection=NULL,BinData=list(fingerprintMat,
#' targetMat),datanames=c("FP","TP"),nrclusters=7,topC=NULL,sign=0.05,fusionsLog=TRUE,
#' weightclust=TRUE,names=c("FP","TP"))
#' }
CharacteristicFeatures<-function(List,Selection=NULL,binData=NULL,contData=NULL,datanames=NULL,nrclusters=NULL,sign=0.05,topChar=NULL,fusionsLog=TRUE,weightclust=TRUE,names=NULL){
	if(is.null(datanames)){
		for(j in 1:(length(binData)+length(contData))){
			datanames[j]=paste("Data",j,sep=" ")	
		}
	}
	
	if(!(is.null(Selection))){
		ResultFeat=FeatSelection(List,Selection,binData,contData,datanames,nrclusters,topChar,sign,fusionsLog,weightclust)
	}
	else{
		ListNew=list()
		element=0
		for(i in 1:length(List)){
			if(attributes(List[[i]])$method != "CEC" & attributes(List[[i]])$method != "Weighted" & attributes(List[[i]])$method!= "WeightedSim"){
				ResultsClust=list()
				ResultsClust[[1]]=list()
				ResultsClust[[1]][[1]]=List[[i]]
				names(ResultsClust[[1]])[1]="Clust"
				element=element+1					
				ListNew[[element]]=ResultsClust[[1]]
				#attr(ListNew[element],"method")="Weights"
			}
			else if(attributes(List[[i]])$method=="CEC" | attributes(List[[i]])$method=="Weighted" | attributes(List[[i]])$method == "WeightedSim"){
				ResultsClust=list()
				if(weightclust==TRUE){
					ResultsClust[[1]]=list()
					if(attributes(List[[i]])$method != "WeightedSim"){
						ResultsClust[[1]][[1]]=List[[i]]$Clust
						names(ResultsClust[[1]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[1]]
						attr(ListNew[element],"method")="Weights"
					}
					else{
						ResultsClust[[1]]=list()
						ResultsClust[[1]][[1]]=List[[i]]
						names(ResultsClust[[1]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[1]]
					}
				}
				else{
					for (j in 1:length(List[[i]]$Results)){
						ResultsClust[[j]]=list()
						ResultsClust[[j]][[1]]=List[[i]]$Results[[j]]
						names(ResultsClust[[j]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[j]]
						attr(ListNew[element],"method")="Weights"
					}		
				}		
			}	
		}
		
		if(is.null(names)){
			names=seq(1,length(ListNew),1)
			for(i in 1:length(ListNew)){
				names[i]=paste("Method",i,sep=" ")
			}
		}
		names(ListNew)=names
		MatrixClusters=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)
		
		List=ListNew
		
		if(!is.null(binData)){
			cpdSet <- rownames(binData[[1]])
		}
		else if(!is.null(contData)){
			cpdSet <- rownames(contData[[1]])
		}
		else{
			stop("Specify a data set in binData and/or in contData")
		}
		
		ResultFeat=list()
		maxclus=0
		for (k in 1:dim(MatrixClusters)[1]){
			
			clusters=MatrixClusters[k,]
			if(max(clusters)>maxclus){
				maxclus=max(clusters)
			}
			Characteristics=list()
			clust=sort(unique(clusters)) #does not matter: Genes[i] puts right elements on right places
			hc<-stats::as.hclust(List[[k]]$Clust$Clust)
			OrderedCpds <- hc$labels[hc$order]
			for (i in clust){	
				
				temp=list()
				LeadCpds=names(clusters)[which(clusters==i)] 
				temp[[1]]=list(LeadCpds,OrderedCpds)
				names(temp[[1]])=c("LeadCpds","OrderedCpds") #names of the objects
				
				group <- factor(ifelse(cpdSet %in% LeadCpds, 1, 0)) #identify the group of interest
				
				#Determine characteristic features for the objects: fishers exact test
				result=list()
				if(!is.null(binData)){
					for(i in 1:length(binData)){
						binData[[i]]=binData[[i]]+0
						binData[[i]]<-binData[[i]][,which(colSums(binData[[i]]) != 0 & colSums(binData[[i]]) != nrow(binData[[i]]))]
					}
					
					
					for(j in 1: length(binData)){
						binMat=binData[[j]]
						
						pFish <- apply(binMat, 2, function(x) stats::fisher.test(table(x, group))$p.value)
						pFish <- sort(pFish)
						adjpFish<-stats::p.adjust(pFish, method = "fdr")
						
						AllFeat=data.frame(Names=names(pFish),P.Value=pFish,adj.P.Val=adjpFish)
						AllFeat$Names=as.character(AllFeat$Names)
						
						if(is.null(topChar)){
							topChar=length(which(pFish<sign))
						}
						
						TopFeat=AllFeat[0:topChar,]
						TopFeat$Names=as.character(TopFeat$Names)
						temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
						result[[j]]<-temp1
						names(result)[j]=datanames[length(binData)+j]
						
					}
				}
				
				resultC=list()
				if(!is.null(contData)){
					for(j in 1:length(contData)){
						contMat=contData[[j]]
						
						group1=which(group==1)
						group2=which(group==0)
						
						
						pTTest <- apply(contMat, 2, function(x) stats::t.test(x[group1],x[group2])$p.value)
						
						pTTest <- sort(pTTest)
						adjpTTest<-stats::p.adjust(pTTest, method = "fdr")
						
						AllFeat=data.frame(Names=as.character(names(pTTest)),P.Value=pTTest,adj.P.Val=adjpTTest)
						AllFeat$Names=as.character(AllFeat$Names)
						if(is.null(topChar)){
							topChar=length(which(pTTest<sign))
						}
						
						TopFeat=data.frame(Names=as.character(names(pTTest[0:topChar])),P.Value=pTTest[0:topChar],adj.P.Val=adjpTTest[0:topChar])
						TopFeat$Names=as.character(TopFeat$Names)
						temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
						resultC[[j]]<-temp1
						names(resultC)[j]=datanames[j]
						
					}
				}
				
				temp[[2]]=c(result,resultC)
				
				names(temp)=c("objects","Characteristics")
				
				Characteristics[[i]]=temp
				
				names(Characteristics)[i]=paste("Cluster",i,sep=" ")
			}
			ResultFeat[[k]]=Characteristics
			
		}
		names(ResultFeat)=names
		for(i in 1:length(ResultFeat)){
			for(k in 1:length(ResultFeat[[i]])){
				if(is.null(ResultFeat[[i]][[k]])[1]){
					ResultFeat[[i]][[k]]=NA
					names(ResultFeat[[i]])[k]=paste("Cluster",k,sep=" ")
				}			
			}
			if(length(ResultFeat[[i]]) != maxclus){
				extra=maxclus-length(ResultFeat[[i]])
				for(j in 1:extra){
					ResultFeat[[i]][[length(ResultFeat[[i]])+1]]=NA
					names(ResultFeat[[i]])[length(ResultFeat[[i]])]=paste("Cluster",length(ResultFeat[[i]]),sep=" ")
				}
			}
		} 	
		
	}
	
	return(ResultFeat)	
}


#' @title Interactive plot to determine DE Genes and DE features for a specific
#' cluster
#' 
#' @description If desired, the function produced a dendrogram of a clustering results. One
#' or multiple cluster can be indicated by a mouse click. From these clusters
#' DE genes and characteristic features are determined. It is also possible to
#' provide the objects of interest without producing the plot. Note, it is required to click on the dendrogram branches, not on the objects.
#' #' @export ChooseCluster
#' @details The DE genes are determined by testing for significance of the specified
#' cluster versus all other objects combined. This is performed by the limma
#' function. The binary features are evaluated with the fisher exact test while
#' the continuous features are tested with the t-test. Multiplicity correction
#' is included.
#' 
#' @param Interactive Logical. Whether an interactive plot should be made. Defaults to TRUE.
#' @param leadCpds A list of the objects of the clusters of interest. If
#' Interactive=TRUE, these are determined by the mouse-click and it defaults to
#' NULL.
#' @param clusterResult The output of one of the aggregated cluster functions,
#' The clustering result of interest. Default is NULL.
#' @param colorLab The clustering result the dendrogram should be colored after
#' as in \code{ClusterPlot}. It is the output of one of the clustering
#' functions.
#' @param binData A list of the binary feature data matrices. These will be
#' evaluated with the fisher's extact test. Default is NULL.
#' @param contData A list of continuous data sets of the objects. These will
#' be evaluated with the t-test. Default is NULL.
#' @param datanames A vector with the names of the data matrices. Default is NULL.
#' @param geneExpr A gene expression matrix, may also be an ExpressionSet. The
#' rows should correspond with the genes. Default is NULL.
#' @param topChar The number of top characteristics to return. If NULL, only
#' the significant characteristics are saved. Default is NULL.
#' @param topG The number of top genes to return. If NULL, only the significant
#' genes are saved. Default is NULL.
#' @param sign The significance level. Default is 0.05.
#' @param nrclusters Optional. The number of clusters to cut the dendrogram in.
#' If NULL, the dendrogram will be plotted without colors to discern the
#' different clusters. Default is NULL.
#' @param cols The colors to use in the dendrogram. Default is NULL.
#' @param n The number of clusters one wants to identify by a mouse click. Default is 1.
#' @return The returned value is a list with one element per cluster of
#' interest indicated by the prefix "Choice".  This element is again a list
#' with the following three elements: \item{objects}{A list with the elements
#' LeadCpds (the objects of interest) and OrderedCpds (all objects in the
#' order of the clustering result)} \item{Characteristics}{The found (top)
#' characteristics of the feature data} \item{Genes}{A list with the elements
#' TopDE (a table with information on the top genes) and AllDE (a table with
#' information on all genes)}
#' @examples
#' 
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(geneMat)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' MCF7_Interactive=ChooseCluster(Interactive=TRUE,leadCpds=NULL,clusterResult=MCF7_T,
#' colorLab=MCF7_F,binData=list(fingerprintMat),datanames=c("FP"),geneExpr=geneMat,
#' topChar = 20, topG = 20,nrclusters=7,n=1)
#' }
ChooseCluster=function(Interactive=TRUE,leadCpds=NULL,clusterResult=NULL,colorLab=NULL,binData=NULL,contData=NULL,datanames=c("FP"),geneExpr=NULL,topChar = 20, topG = 20,sign=0.05,nrclusters=NULL,cols=NULL,n=1){
	if(is.null(datanames)){
		for(j in 1:(length(binData)+length(contData))){
			datanames[j]=paste("Data",j,sep=" ")	
		}
	}
	OrInteractive=Interactive
	
	if(Interactive==TRUE){
		#windows()
		ClusterPlot(clusterResult,colorLab,nrclusters,cols)
		hc1<-stats::as.hclust(clusterResult$Clust)
		ClusterSpecs<-list()
		ClusterSpecs=graphics::identify(hc1, N=n, MAXCLUSTER = nrow(binData[[1]]), function(j) ChooseCluster(Interactive=FALSE,leadCpds=rownames(binData[[1]][j,]),clusterResult,colorLab=NULL,binData=binData,contData=contData,datanames=datanames,geneExpr=geneExpr,topChar=topChar,topG=topG,sign=sign,nrclusters=nrclusters,cols=cols))		
		
		names(ClusterSpecs)<-sapply(seq(1,n),FUN=function(x) paste("Choice",x,sep=" "))
		
	}
	else{

		DistM<-clusterResult$DistM
		Clust<-clusterResult$Clust
		if(is.null(Clust)){
			clusterResult$Clust=clusterResult
			Clust<-clusterResult$Clust
		}
		
		hc <- stats::as.hclust(Clust)
		OrderedCpds <- hc$labels[hc$order]
		
		if(class(leadCpds)=="character"){
			leadCpds=list(leadCpds)
		}
		
		Specs=list()
		
		if(!is.null(binData)){
			cpdSet <- rownames(binData[[1]])
		}
		else if(!is.null(contData)){
			cpdSet <- rownames(contData[[1]])
		}
		else if(!is.null(geneExpr)){
			cpdSet <- colnames(geneExpr)
		}
		else{
			stop("Specify a data set in binData, contData and/or geneExpr")
		}
		
		for(a in 1:length(leadCpds)){
			objects=list(leadCpds[[a]],OrderedCpds)
			names(objects)=c("LeadCpds","OrderedCpds")
			group <- factor(ifelse(cpdSet %in% leadCpds[[a]], 1, 0)) #identify the group of interest
			
			#Determine characteristic features for the objects: fishers exact test
			Characteristics=list()
			
			resultB=list()
			if(!is.null(binData)){
				
				if(class(binData)!="list"){
					stop("The binary data matrices must be put into a list")
				}
				for(i in 1:length(binData)){
					binData[[i]]=binData[[i]]+0
					binData[[i]]<-binData[[i]][,which(colSums(binData[[i]]) != 0 & colSums(binData[[i]]) != nrow(binData[[i]]))]
				}
				for(j in 1: length(binData)){
					
					binMat=binData[[j]]
					
					pFish <- apply(binMat, 2, function(x) stats::fisher.test(table(x, group))$p.value)
					
					pFish <- sort(pFish)
					adjpFish<-stats::p.adjust(pFish, method = "fdr")
					
					AllFeat=data.frame(Names=names(pFish),P.Value=pFish,adj.P.Val=adjpFish)
					AllFeat$Names=as.character(AllFeat$Names)
					
					if(is.null(topChar)){
						topChar=length(which(pFish<sign))
					}
					
					TopFeat=AllFeat[0:topChar,]
					TopFeat$Names=as.character(TopFeat$Names)
					temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
					resultB[[j]]<-temp1
					names(resultB)[j]=datanames[length(binData)+j]
					
				}
			}
			
			resultC=list()
			if(!is.null(contData)){
				
				if(class(contData)!="list"){
					stop("The continuous data matrices must be put into a list")
				}
				for(j in 1:length(contData)){
					contMat=contData[[j]]
					
					group1=which(group==1)
					group2=which(group==0)
					
					pTTest <- apply(contMat, 2, function(x) stats::t.test(x[group1],x[group2])$p.value)
					
					pTTest <- sort(pTTest)
					adjpTTest<-stats::p.adjust(pTTest, method = "fdr")
					
					AllFeat=data.frame(Names=as.character(names(pTTest)),P.Value=pTTest,adj.P.Val=adjpTTest)
					AllFeat$Names=as.character(AllFeat$Names)
					if(is.null(topChar)){
						topChar=length(which(pTTest<sign))
					}
					
					TopFeat=data.frame(Names=as.character(names(pTTest[0:topChar])),P.Value=pTTest[0:topChar],adj.P.Val=adjpTTest[0:topChar])
					TopFeat$Names=as.character(TopFeat$Names)
					temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
					resultC[[j]]<-temp1
					names(resultC)[j]=datanames[j]
					
				}
				
			}
			
			Characteristics=c(resultB,resultC)
			names(Characteristics)=datanames
			
			#Determine DE Genes with limma --> make difference between "regular" data matrix and "expression set"
			#GeneExpr.2=GeneExpr[,colnames(Matrix)]
			if(!is.null(geneExpr)){
				cpdSetG <-colnames(geneExpr)
				groupG <- factor(ifelse(cpdSetG %in% leadCpds[[a]], 1, 0)) #identify the group of interest
				if(class(geneExpr)[1]=="ExpressionSet"){
					geneExpr$LeadCmpds<-groupG		
					
					
					if (!requireNamespace("a4Base", quietly = TRUE)) {
						stop("a4Base needed for this function to work. Please install it.",
								call. = FALSE)
					}
					
					DElead <- a4Base::limmaTwoLevels(geneExpr,"LeadCpds")
					
					#allDE <- topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL), resort.by = "logFC",sort.by="p")
					allDE <- a4Core::topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL),sort.by="p")
					if(is.null(allDE$ID)){
						allDE$ID <- rownames(allDE)
					}
					else
					{
						allDE$ID=allDE$ID
					}
					if(is.null(topG)){
						topG=length(which(allDE$adj.P.Val<=sign))
					}
					TopDE <- allDE[0:topG, ]
					#TopAdjPval<-TopDE$adj.P.Val
					#TopRawPval<-TopDE$P.Value
					
					#RawpVal<-allDE$P.Value
					#AdjpVal <- allDE$adj.P.Val
					#genesEntrezId <- allDE$ENTREZID
					
					Genes<-list(TopDE,allDE)
					names(Genes)<-c("TopDE","AllDE")
					#Genes <- list(TopDE$SYMBOL,TopAdjPval,TopRawPval,genesEntrezId,RawpVal,AdjpVal)	
					#names(Genes)<-c("DE_Genes","DE_RawPvals","DE_AdjPvals", "All_Genes", "All_RawPvals","All_AdjPvals")
				}
				else{
					
					label.factor = factor(groupG)
					design = stats::model.matrix(~label.factor)
					fit = limma::lmFit(geneExpr,design=design)
					fit = limma::eBayes(fit)
					
					#allDE = topTable(fit,coef=2,adjust="fdr",n=nrow(GeneExpr),resort.by = "logFC", sort.by="p")
					allDE = limma::topTable(fit,coef=2,adjust="fdr",n=nrow(geneExpr), sort.by="p")
					if(is.null(allDE$ID)){
						allDE$ID <- rownames(allDE)
					}
					else
					{
						allDE$ID=allDE$ID
					}
					if(is.null(topG)){
						topG=length(which(allDE$adj.P.Val<=sign))
					}
					TopDE=allDE[0:topG,]
					#TopAdjPval<-TopDE$adj.P.Val
					#TopRawPval<-TopDE$P.Value
					
					#RawpVal<-allDE$P.Value
					#AdjpVal <- allDE$adj.P.Val
					
					Genes<-list(TopDE,allDE)
					names(Genes)<-c("TopDE","AllDE")
					#Genes <- list(TopDE[,1],TopAdjPval,TopRawPval,allDE[,1],RawpVal,AdjpVal)	
					#names(Genes)<-c("DE_Genes","DE_RawPvals","DE_AdjPvals", "All_Genes", "All_RawPvals","All_AdjPvals")
					
				}
				
				out=list(objects,Characteristics,Genes)
				names(out)=c("objects","Characteristics","Genes")
				Specs[[a]]=out
				names(Specs)[a]=paste("Choice",a,sep=" ")
			}	
			else{
				out=list(objects,Characteristics)
				names(out)=c("objects","Characteristics")
				Specs[[a]]=out
				names(Specs)[a]=paste("Choice",a,sep=" ")
			}
		}
		
		if(OrInteractive==TRUE|length(Specs)==1){
			return(out)
		}
		else{
			return(Specs)
		}
	}
	class(ClusterSpecs)="ChosenClusters"
	return(ClusterSpecs)
}


#' @title Matching clusters with colours
#' @param x The leaf of a dendrogram.
#' @param Data A clustering object.
#' @param nrclusters The number of clusters to divide the dendrogram in. Default is NULL.
#' @param cols A character vector with the colours to be used. Default is NULL.
#' @param colorComps A character vector of a specific set of objects to be coloured. Default is NULL.
#' @description Internal function of \code{ClusterPlot}.
ClusterCols <- function(x,Data,nrclusters=NULL,cols=NULL,colorComps=NULL) {
	
	if(is.null(nrclusters) & is.null(colorComps)){
		return(x)
	}
	else if(!is.null(nrclusters)){
		if(length(cols)<nrclusters){
			stop("Not for every cluster a color is specified")
		}
	}	
	
	if(!is.null(nrclusters)){
		Clustdata=stats::cutree(Data,nrclusters)
		Clustdata=Clustdata[Data$order]
		
		ordercolors=Clustdata
		order=seq(1,nrclusters)
		
		for (k in 1:length(unique(Clustdata))){
			select=which(Clustdata==unique(Clustdata)[k])
			ordercolors[select]=order[k]
		}
		names(ordercolors)=Data$order.lab
	}
	else{
		cols=rep("black",length(Data$order.lab))
		names(cols)=Data$order.lab
		cols[which(names(cols)%in%colorComps)]="red"
		ordercolors=cols
		
	}
	
	colfunc=function(x,cols,colorComps){
#		if(is.null(colorComps)){
#			color=cols[which(names(cols)==x)]
#			indextemp=which(attr(Data$diss,"Labels")==x)
#			index1=which(Data$order==indextemp)	
#			
#			index2=ordercolors[index1]
#			
#			color=cols[index2]
#			
#		}
#		else{
#			color=cols[which(names(cols)==x)]
#		}
		color=cols[which(names(cols)==x)]
		return(color)	
	}
	
	if (stats::is.leaf(x)) {
		## fetch label
		label <- attr(x, "label") 
		## set label color to clustercolor
		attr(x, "nodePar") <- list(pch=NA,lab.col=colfunc(label,ordercolors,colorComps),lab.cex=0.9,font=2)
		attr(x, "edgePar") <- list(lwd=2,col=colfunc(label,ordercolors,colorComps))
	}
	return(x)
}


#' @title Colouring clusters in a dendrogram
#' 
#' @description Plot a dendrogram with leaves colored by a result of choice.
#' @export ClusterPlot
#' @param Data1 The resulting clustering of a method which contains the
#' dendrogram to be colored.
#' @param Data2 Optional. The resulting clustering of another method , i.e. the
#' resulting clustering on which the colors should be based. Default is NULL.
#' @param nrclusters Optional. The number of clusters to cut the dendrogram in.
#' If not specified the dendrogram will be drawn without colours to discern the
#' different clusters. Default is NULL.
#' @param cols The colours for the clusters if nrclusters is specified. Default is NULL.
#' @param colorComps If only a specific set of objects needs to be
#' highlighted, this can be specified here. The objects should be given in a
#' character vector. If specified, all other compound labels will be colored
#' black. Default is NULL.
#' @param hangdend A specification for the length of the brances of the dendrogram. Default is 0.02.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @param \dots Other options which can be given to the plot function.
#' @return A plot of the dendrogram of the first clustering result with colored
#' leaves. If a second clustering result is given in Data2, the colors are
#' based on this clustering result.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(Colors1)
#' 
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' ClusterPlot(MCF7_T ,nrclusters=7,cols=Colors1,plottype="new",location=NULL,
#' main="Clustering on Target Predictions: Dendrogram",ylim=c(-0.1,1.8))
#' }
ClusterPlot<-function(Data1,Data2=NULL,nrclusters=NULL,cols=NULL,colorComps=NULL,hangdend=0.02,plottype="new",location=NULL,...){
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	
	cx=Data1$Clust
	if(is.null(Data2)){
		Data=Data1$Clust
	}
	else{
		Data=Data2$Clust
	}
	
	d_temp<- stats::dendrapply(stats::as.dendrogram(as.hclust(cx),hang=hangdend),ClusterCols,Data,nrclusters,cols,colorComps)
	plottypein(plottype,location)
	graphics::plot(d_temp,nodePar=list(pch=NA),edgePar=list(lwd=2),ylab="Height",font.axis=2,font.lab=2,font=2)
	graphics::axis(side = 2, lwd = 2)	
	plottypeout(plottype)
}


#' @title Create a color palette to be used in the plots
#' 
#' @description In order to facilitate the visualization of the influence of the different
#' methods on the clustering of the objects, colours can be used. The function
#' \code{ColorPalette} is able to pick out as many colours as there are
#' clusters. This is done with the help of the \code{ColorRampPalette} function
#' of the grDevices package
#' @export ColorPalette
#' @param colors A vector containing the colors of choice
#' @param ncols The number of colors to be specified. If higher than the number
#' of colors, it specifies colors in the region between the given colors.
#' @return A vector containing the hex codes of the chosen colors.
#' @examples
#' 
#' Colors1<-ColorPalette(c("cadetblue2","chocolate","firebrick2",
#' "darkgoldenrod2", "darkgreen","blue2","darkorchid3","deeppink2"), ncols=8)
ColorPalette<-function(colors=c("red","green"),ncols=5){
	my_palette=grDevices::colorRampPalette(colors)(ncols)
	
	return(my_palette)
	
}


#' @title Interactive comparison of clustering results for a specific cluster or
#' method.
#' 
#' @description A visual comparison of all methods is handy to see which objects will
#' always cluster together independent of the applied methods. The function
#' \code{CompareInteractive} plots the comparison over the specified methods. A
#' cluster or method can than be identified by clicking and is plotted
#' separately against the single source or other specified methods.
#' @export CompareInteractive
#' @param ListM A list of the multiple source clustering or other methods to be
#' compared and from which a cluster or method will be identified. The first
#' element of the list will be used as the reference in
#' \code{ReorderToReference}.
#' @param ListS A list of the single source clustering or other methods the
#' identified result will be compared to. The first element of the list will be
#' used as the reference in \code{ReorderToReference}.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param cols A character vector with the names of the colours.  Default is NULL.
#' @param fusionsLogM The fusionsLog parameter for the elements in ListM. To be
#' handed to \code{ReorderToReference}. Default is FALSE.
#' @param fusionsLogS The fusionslog parameter for the elements in ListS. To be
#' handed to \code{ReorderToReference}. Default is FALSE.
#' @param weightclustM The weightclust parameter for the elements in ListM. To
#' be handed to \code{ReorderToReference}. Default is FALSE.
#' @param weightclustS The weightclust parameter for the elements in ListS. To
#' be handed to \code{ReorderToReference}. Default is FALSE.
#' @param namesS Optional. Names of the single source clusterings to be used as
#' labels for the columns. Default is NULL.
#' @param namesM Optional. Names of the multiple source clusterings to be used
#' as labels for the columns. Default is NULL.
#' @param marginsM Optional. Margins to be used for the plot for the elements
#' is ListM after the identification. Default is c(2,2.5,2,2.5).
#' @param marginsS Optional. Margins to be used for the plot for the elements
#' is ListS after the identification. Default is c(8,2.5,2,2.5).
#' @param Interactive Optional. Do you want an interactive plot? Defaults to
#' TRUE, if not the function provides the same as \code{ComparePlot} for the
#' elements in ListM. Default is TRUE.
#' @param n The number of methods/clusters you want to identify. Default is 1.
#' @return The returned value is a plot of the comparison of the elements of
#' ListM. On this plot multiple clusters and/or methods can be identified.
#' Click on a cluster of a specific method to see how that cluster of that
#' method compares to the elements in ListS. Click left next to a row to
#' identify a all cluster of a specific method. A new plotting window will
#' appear for every identification.
#' @examples
#' 
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(Colors1)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' L=list(fingerprintMat,targetMat)
#' 
#' MCF7_W=WeightedClust(List=L,type="data",distmeasure=c("tanimoto","tanimoto"),
#' normalize=c(FALSE,FALSE),method=c(NULL,NULL),weight=seq(1,0,-0.1),weightclust=0.5,
#' clust="agnes",linkage="ward",StopRange=FALSE)
#' 
#' ListM=list(MCF7_W)
#' namesM=c(seq(1.0,0.0,-0.1))
#' 
#' ListS=list(MCF7_F,MCF7_T)
#' namesS=c("FP","TP")
#' 
#' CompareInteractive(ListM,ListS,nrclusters=7,cols=Colors1,fusionsLogM=FALSE,
#' fusionsLogS=FALSE,weightclustM=FALSE,weightclustS=TRUE,namesM,namesS,
#' marginsM=c(2,2.5,2,2.5),marginsS=c(8,2.5,2,2.5),Interactive=TRUE,n=1)
#' }
CompareInteractive<-function(ListM,ListS,nrclusters=NULL,cols=NULL,fusionsLogM=FALSE,fusionsLogS=FALSE,weightclustM=FALSE,weightclustS=FALSE,namesM=NULL,namesS=NULL,marginsM=c(2,2.5,2,2.5),marginsS=c(8,2.5,2,2.5),Interactive=TRUE,n=1){
	
	MatrixColorsM=ReorderToReference(ListM,nrclusters,fusionsLogM,weightclustM,namesM)
	
	NamesM=ColorsNames(MatrixColorsM,cols)
	
	nobsM=dim(MatrixColorsM)[2]
	nmethodsM=dim(MatrixColorsM)[1]
	
	if(is.null(namesM)){
		for(j in 1:nmethodsM){
			namesM[j]=paste("Method",j,sep=" ")	
		}
	}
	
	similarM=round(SimilarityMeasure(MatrixColorsM),2)
	
	grDevices::dev.new()
	graphics::par(mar=marginsM)
	plotrix::color2D.matplot(MatrixColorsM,cellcolors=NamesM,show.values=FALSE,axes=FALSE,xlab="",ylab="")	
	graphics::axis(1,at=seq(0.5,(nobsM-0.5)),labels=colnames(MatrixColorsM),las=2,cex.axis=0.70)
	graphics::axis(2,at=seq(0.5,(nmethodsM-0.5)),labels=rev(namesM),cex.axis=0.65,las=2)
	graphics::axis(4,at=seq(0.5,(nmethodsM-0.5)),labels=rev(similarM),cex.axis=0.65,las=2)
	
	if(Interactive==TRUE){
		yseq=c(seq(dim(MatrixColorsM)[1]-0.5,0.5,-1))
		for(i in seq(dim(MatrixColorsM)[1]-0.5,0.5,-1)){
			yseq=c(yseq,rep(i,dim(MatrixColorsM)[2]))
		}
		ids=graphics::identify(x=c(rep(-1,dim(MatrixColorsM)[1]),rep(seq(0.5,dim(MatrixColorsM)[2]-0.5),dim(MatrixColorsM)[1])),y=yseq,n=n,plot=FALSE)
		
		comparison<-function(id){
			if(id%in%seq(dim(MatrixColorsM)[1])){
				grDevices::dev.new()
				graphics::layout(matrix(c(1,2),nrow=2), heights=c(1,2))
				NamesMSel=NamesM[id,]
				namesMSel=namesM[id]
				
				graphics::par(mar=marginsS)
				plotrix::color2D.matplot(t(as.matrix(MatrixColorsM[id,])),cellcolors=NamesMSel,show.values=FALSE,axes=FALSE,xlab="",ylab="")
				graphics::axis(2,at=c(0.5),labels=rev(namesMSel),cex.axis=0.65,las=2)
				#axis(4,at=c(0.5),labels=rev(similarSel),cex.axis=0.65,las=2)
				
				#Find reference for MatrixColorsM
				if(weightclustM==FALSE){
					temp=FindElement("Results",ListM[1])
					if(!(is.null(temp))&length(temp)!=0){
						Ref=list(Clust=temp$Results_1[[1]])
						attr(Ref,'method')<-attributes(ListM[[1]])$method
					}
					else if(length(temp)==0){
						temp=FindElement("Clust",ListM[1])
						Ref=list(Clust=temp$Clust_1)
						attr(Ref,'method')<-attributes(ListM[[1]])$method
					}
					else{
						message('Cannot find a reference for the second plot, try: weightclust=TRUE')
					}
				}
				else{
					Ref=ListM[[1]]
					attr(Ref,'method')<-attributes(ListM[[1]])$method
				}
				
				L=c(Ref,ListS)
				for(i in 1:length(L)){
					if(i==1){
						attr(L[[1]],'method')<-"Ref"
					}
					else{
						attr(L[[i]],"method")<-attributes(ListS[[i-1]])$method
					}					
				}
				MatrixColorsS=ReorderToReference(L,nrclusters,fusionsLogS,weightclustS,names=c("Ref",namesS))
				MatrixColorsS=MatrixColorsS[-1,]
				NamesS=ColorsNames(MatrixColorsS,cols)
				
				nobs=dim(MatrixColorsS)[2]
				nmethodS=dim(MatrixColorsS)[1]
				
				if(is.null(namesS)){
					for(j in 1:nmethodS){
						namesS[j]=paste("Method",j,sep=" ")	
					}
				}
				
				similarS=round(SimilarityMeasure(MatrixColorsS),2)
				
				graphics::par(mar=marginsS)
				plotrix::color2D.matplot(MatrixColorsS,cellcolors=NamesS,show.values=FALSE,axes=FALSE,xlab="",ylab="")	
				graphics::axis(1,at=seq(0.5,(nobs-0.5)),labels=colnames(MatrixColorsM),las=2,cex.axis=0.70)
				graphics::axis(2,at=c(seq(0.5,nmethodS-0.5)),labels=rev(namesS),cex.axis=0.65,las=2)
				graphics::axis(4,at=c(seq(0.5,nmethodS-0.5)),labels=rev(similarS),cex.axis=0.65,las=2)			
			}
			else{
				SelCluster=t(MatrixColorsM)[id-nrow(MatrixColorsM)]
				Temp=sapply(seq(nrow(MatrixColorsM)),function(i) ncol(MatrixColorsM)*i)
				Row=which(Temp>(id-nrow(MatrixColorsM)))[1]
				Index=which(MatrixColorsM[Row,]!=SelCluster)
				
				grDevices::dev.new()
				graphics::layout(matrix(c(1,2),nrow=2), heights=c(1,2))
				NamesMSel=NamesM[Row,]
				NamesMSel[Index]="white"
				namesMSel=namesM[Row]
				
				graphics::par(mar=marginsM)
				plotrix::color2D.matplot(t(as.matrix(MatrixColorsM[Row,])),cellcolors=NamesMSel,show.values=FALSE,axes=FALSE,xlab="",ylab="")	
				graphics::axis(2,at=c(0.5),labels=rev(namesMSel),cex.axis=0.65,las=2)	
				
				#Find reference for MatrixColorsM
				if(weightclustM==FALSE){
					temp=FindElement("Results",ListM[1])
					if(!(is.null(temp))&length(temp)!=0){
						Ref=list(Clust=temp$Results_1[[1]])
						attr(Ref,'method')<-attributes(ListM[[1]])$method
					}
					else if(length(temp)==0){
						temp=FindElement("Clust",ListM[1])
						Ref=list(Clust=temp$Clust_1)
						attr(Ref,'method')<-attributes(ListM[[1]])$method
					}
					else{
						message('Cannot find a reference for the second plot, try: weightclust=TRUE')
					}
				}
				else{
					Ref=ListM[[1]]
					attr(Ref,'method')<-attributes(ListM[[1]])$method
				}
				
				L=c(Ref,ListS)
				for(i in 1:length(L)){
					if(i==1){
						attr(L[[1]],'method')<-"Ref"
					}
					else{
						attr(L[[i]],"method")<-attributes(ListS[[i-1]])$method
					}					
				}
				
				MatrixColorsS=ReorderToReference(L,nrclusters,fusionsLogS,weightclustS,names=c("Ref",namesS))
				MatrixColorsS=MatrixColorsS[-1,]
				
				IndexS=lapply(seq(nrow(MatrixColorsS)),function(i) which(MatrixColorsS[i,]!=SelCluster))
				
				NamesS=ColorsNames(MatrixColorsS,cols)
				for(i in 1:nrow(NamesS)){
					NamesS[i,IndexS[[i]]]="white"
				}
				
				nobs=dim(MatrixColorsS)[2]
				nmethodsS=dim(MatrixColorsS)[1]
				
				if(is.null(namesS)){
					for(j in 1:nmethodsS){
						namesS[j]=paste("Method",j,sep=" ")	
					}
				}
				
				similarS=round(SimilarityMeasure(MatrixColorsS),2)
				
				graphics::par(mar=marginsS)
				plotrix::color2D.matplot(MatrixColorsS,cellcolors=NamesS,show.values=FALSE,axes=FALSE,xlab="",ylab="")	
				graphics::axis(1,at=seq(0.5,(nobs-0.5)),labels=colnames(MatrixColorsM),las=2,cex.axis=0.70)
				graphics::axis(2,at=c(seq(0.5,nmethodsS-0.5)),labels=rev(namesS),cex.axis=0.65,las=2)
				graphics::axis(4,at=c(seq(0.5,nmethodsS-0.5)),labels=rev(similarS),cex.axis=0.65,las=2)		
				
				
			}
		}
		
		plots=sapply(seq(length(ids)),function(i) comparison(ids[i]))
		
	}
	
}

#' @title Compares medoid clustering results based on silhouette widths
#' 
#' @description The function \code{CompareSilCluster} compares the results of two medoid
#' clusterings. The null hypothesis is that the clustering is identical. A test
#' statistic is calcluated and a p-value obtained with bootstrapping. See
#' "Details" for a more elaborate description.
#' @export CompareSilCluster
#' @details For the data or distance matrices in List, medoid clustering with nrclusters
#' is set up by the \code{pam} function of the \pkg{cluster} and the silhouette
#' widths are retrieved. These widths indicate how well an object fits in its
#' current cluster. Values around one indicate an appropriate cluster while
#' values around zero indicate that the object might as well lie in its
#' neighbouring cluster. The silhouette widths are than regressed in function
#' of the cluster membership of the objects. First the widths are modelled
#' according to the cluster membership of object these were derived from. Next,
#' these are modeled in function of the membership determined by the other
#' object. The regression function is fit by the \code{lm} function and the
#' \code{r.squared} value is retrieved. The\code{r.squared} value indicates how
#' much of the variance of the silhouette widths is explained by the
#' membership. Optimally this value is high.
#' 
#' Next, a statistic is determined. Suppose that RXX is the \code{r.squared}
#' retrieved from regressing the silhouette widths of object X versus the
#' corresponding cluster membership of object X and RXY the \code{r.squared}
#' retrieved from regressing the silhouette widths of object X versus the
#' cluster membership determined by object Y and vice versa.  The statistic is
#' obtained as: \deqn{Stat=abs(\sum{RXX}-\sum{RXY})}
#' 
#' The lower the statistical value, the better the clustering is explained by
#' the sources. Via bootstrapping a p-value is obtained.
#' 
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clusters".	
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}. 
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param nrclusters The number of clusters to cut the dendrogram in. This is
#' necessary for the computation of the Jaccard coefficient. Default is NULL.
#' @param names The labels to give to the elements in List. Default is NULL.
#' @param nboot Number of bootstraps to be run. Default is 100.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return A plots are made of the density of the statistic under the null
#' hypotheses. The p-value is also indicated on this plot. Further, a list with
#' two elements is returned: \item{Observed Statistic}{The observed statistical
#' value} \item{P-Value}{The P-value of the obtained statistic retrieved after
#' bootstrapping}
#' @examples
#' 
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' 
#' List=list(fingerprintMat,targetMat)
#' 
#' Comparison=CompareSilCluster(List=List,type="data",
#' distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),
#' nrclusters=7,names=NULL,nboot=100,plottype="new",location=NULL)
#' 
#' Comparison
#' }
CompareSilCluster<-function(List,type=c("data","dist"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),nrclusters=NULL,names=NULL,nboot=100,plottype="new",location=NULL){
	
	type=match.arg(type)
	
	if(is.null(names)){
		names=c()
		for(i in 1:length(List)){
			names=c(names,paste("Method",i,sep=" "))
		}	
	}

	if(type=="data"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}
		Dist=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i],normalize[i],method[i]))
		silwidth=lapply(Dist,function(x) cluster::pam(x,nrclusters)$silinfo$widths)
		names(silwidth)=names
		
	}
	else{
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		Dist=List
		silwidth=lapply(Dist,function(x) cluster::pam(x,nrclusters)$silinfo$widths)
		names(silwidth)=names
	}
	
	
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	
	regressioncomb=gtools::permutations(n=length(List),r=2,repeats.allowed=T)
	
	StatRSq<-function(regressioncomb,silwidth,ordernames,names){
		
		regressRSq<-function(x,silwidth,ordernames,names){
			
			i1=x[1]
			i2=x[2]
			L1=silwidth[[i1]][,3][ordernames]
			L2=silwidth[[i2]][,1][ordernames]
			
			regress<-stats::lm(L1~L2)
			Rsq<-summary(regress)$r.squared
			#names(Rsq)=paste("RSquared_",names[i1],names[i2],sep="_")
			return(Rsq)
			
			#paste names on this object!!!
		}
		
		RSqs=apply(regressioncomb,1,function(x) regressRSq(x,silwidth,ordernames,names))
		
		#for (i in 1:nrow(regressioncomb)){
		#	names(RSqs)[i]=paste("RSquared_",names[regressioncomb[i,1]],names[regressioncomb[i,2]],sep="_")
		#}
		
		stat=0
		xx=0
		xy=0
		for(i in 1:nrow(regressioncomb)){
			if(regressioncomb[i,1]==regressioncomb[i,2]){
				xx=xx+RSqs[i]
			}
			else{
				xy=xy+RSqs[i]
			}
		}
		stat=abs(xx-xy)  #check this formula with Nolen
		names(stat)=NULL
		return(stat)
		
	}	
	
	StatRSqObs=StatRSq(regressioncomb,silwidth,ordernames=rownames(Dist[[1]]),names)
	
	
	#bootstrapping
	statNULL=c(1:nboot)
	perm.rowscols <- function (D, n) 
	{
		s <- sample(1:n)
		D=D[s, s]
		return(D)
	}
	
	for(i in 1:nboot){
		set.seed(i)
		DistNULL=Dist
		DistNULL[[1]] <- perm.rowscols(DistNULL[[1]],nrow(DistNULL[[1]]))
		
		silwidthNULL=lapply(DistNULL,function(x) cluster::pam(x,nrclusters)$silinfo$widths)
		
		statNULL[i]=StatRSq(regressioncomb,silwidthNULL,ordernames=rownames(DistNULL[[1]]),names)
		
	}
	
	pval=(sum(abs(statNULL)<=abs(StatRSqObs))+1)/(nboot+1)
	
	plottypein(plottype,location)
	graphics::plot(stats::density(statNULL),type="l",main="The Density of the Statistic under the H0")
	graphics::abline(v=StatRSqObs)
	
	out=list()
	out[[1]]=StatRSqObs
	#out[[2]]=statNULL
	out[[2]]=pval
	names(out)=c("Observed Statistic","P-Value")
	
	return(out)
}

#' @title Comparison of clustering results for the single and multiple source
#' clustering.
#' 
#' @description A visual comparison of all methods is handy to see which objects will
#' always cluster together independent of the applied methods. The function
#' \code{CompareSvsM} plots the \code{ComparePlot} of the single source
#' clustering results on the left and that of the multiple source clustering
#' results on the right such that a visual comparison is possible.
#' @export CompareSvsM
#' @param ListS A list of the outputs from the single source clusterings to be
#' compared. The first element of the list will be used as the reference in
#' \code{ReorderToReference}.
#' @param ListM A list of the outputs from the multiple source clusterings to
#' be compared. The first element of the list will be used as the reference.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param cols A character vector with the names of the colours. Default is NULL.
#' @param fusionsLogS The fusionslog parameter for the elements in ListS. To be
#' handed to \code{ReorderToReference}. Default is FALSE.
#' @param fusionsLogM The fusionsLog parameter for the elements in ListM. To be
#' handed to \code{ReorderToReference}. Default is FALSE.
#' @param weightclustS The weightclust parameter for the elements in ListS. To
#' be handed to \code{ReorderToReference}. Default is FALSE.
#' @param weightclustM The weightclust parameter for the elements in ListM. To
#' be handed to \code{ReorderToReference}. Default is FALSE.
#' @param namesS Optional. Names of the single source clusterings to be used as
#' labels for the columns. Default is NULL.
#' @param namesM Optional. Names of the multiple source clusterings to be used
#' as labels for the columns. Default is NULL.
#' @param margins Optional. Margins to be used for the plot. Default is c(8.1,3.1,3.1,4.1).
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return The returned value is a plot with on the left the comparison over
#' the objects in ListS and on the right a comparison over the objects in
#' ListM.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(Colors1)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' L=list(fingerprintMat,targetMat)
#' 
#' MCF7_W=WeightedClust(List=L,type="data", distmeasure=c("tanimoto","tanimoto"),
#' normalize=c(FALSE,FALSE),method=c(NULL,NULL),weight=seq(1,0,-0.1),weightclust=0.5
#' ,clust="agnes",linkage="ward",StopRange=FALSE)
#' 
#' ListM=list(MCF7_W)
#' namesM=seq(1.0,0.0,-0.1)
#' 
#' ListS=list(MCF7_F,MCF7_T)
#' namesS=c("FP","TP")
#' 
#' CompareSvsM(ListS,ListM,nrclusters=7,cols=Colors1,fusionsLogS=FALSE,
#' fusionsLogM=FALSE,weightclustS=FALSE,weightclustM=FALSE,namesS,
#' namesM,plottype="new",location=NULL)
#' }
CompareSvsM<-function(ListS,ListM,nrclusters=NULL,cols=NULL,fusionsLogS=FALSE,fusionsLogM=FALSE,weightclustS=FALSE,weightclustM=FALSE,namesS=NULL,namesM=NULL,margins=c(8.1,3.1,3.1,4.1),plottype="new",location=NULL){
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new(wdith=14,height=7)
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	nmethodsS=0
	nmethodsM=0
	
	MatrixColorsS=ReorderToReference(ListS,nrclusters,fusionsLogS,weightclustS,namesS)
	MatrixColorsM=ReorderToReference(c(ListS[1],ListM),nrclusters,fusionsLogM,weightclustM,c("ref",namesM))
	
	similarS=round(SimilarityMeasure(MatrixColorsS),2)
	similarM=round(SimilarityMeasure(MatrixColorsM),2)
	
	MatrixColorsM=MatrixColorsM[-c(1),]
	
	NamesM=ColorsNames(MatrixColorsM,cols)
	NamesS=ColorsNames(MatrixColorsS,cols)
	
	nobsM=dim(MatrixColorsM)[2]
	nmethodsM=dim(MatrixColorsM)[1]
	
	nobsS=dim(MatrixColorsS)[2]
	nmethodsS=dim(MatrixColorsS)[1]
	
	if(is.null(namesS)){
		for(j in 1:nmethodsS){
			namesS[j]=paste("Method",j,sep=" ")	
		}
	}
	
	if(is.null(namesM)){
		for(j in 1:nmethodsM){
			namesM[j]=paste("Method",j,sep=" ")	
		}
	}
	
	
	plottypein(plottype,location)
	graphics::par(mfrow=c(1,2),mar=margins)
	plotrix::color2D.matplot(MatrixColorsS,cellcolors=NamesS,show.values=FALSE,axes=FALSE,xlab="",ylab="")
	graphics::axis(1,at=seq(0.5,(nobsS-0.5)),labels=colnames(MatrixColorsS),las=2,cex.axis=0.70)
	graphics::axis(2,at=seq(0.5,(nmethodsS-0.5)),labels=rev(namesS),cex.axis=0.65,las=2)
	graphics::axis(4,at=seq(0.5,(nmethodsS-0.5)),labels=rev(similarS),cex.axis=0.65,las=2)
	
	plotrix::color2D.matplot(MatrixColorsM,cellcolors=NamesM,show.values=FALSE,axes=FALSE,xlab="",ylab="")
	graphics::axis(1,at=seq(0.5,(nobsM-0.5)),labels=colnames(MatrixColorsM),las=2,cex.axis=0.70)
	graphics::axis(2,at=seq(0.5,(nmethodsM-0.5)),labels=rev(namesM),cex.axis=0.65,las=2)
	graphics::axis(4,at=seq(0.5,(nmethodsM-0.5)),labels=rev(similarM[-1]),cex.axis=0.65,las=2)
	plottypeout(plottype)
	
}

#' @title Plot of continuous features
#' 
#' @description The function \code{ContFeaturesPlot} plots the values of continuous features. 
#' It is possible to separate between objects of interest and the
#' other objects. 
#' @export ContFeaturesPlot
#' @param leadCpds A character vector containing the objects one wants to
#' separate from the others.
#' @param data The data matrix.
#' @param nrclusters Optional. The number of clusters to consider if colorLab
#' is specified. Default is NULL.
#' @param orderLab Optional. If the objects are to set in a specific order of
#' a specific method. Default is NULL.
#' @param colorLab The clustering result that determines the color of the
#' labels of the objects in the plot. If NULL, the labels are black. Default is NULL.
#' @param cols The colors for the labels of the objects. Default is NULL.
#' @param ylab The lable of the y-axis. Default is "features".
#' @param addLegend Logical. Indicates whether a legend should be added to the
#' plot. Default is TRUE.
#' @param margins Optional. Margins to be used for the plot. Default is c(5.5,3.5,0.5,8.7).
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return A plot in which the values of the features of the leadCpds are
#' separeted from the others.
#' @examples
#' \dontrun{
#' data(Colors1)
#' Comps=c("Cpd1", "Cpd2", "Cpd3", "Cpd4", "Cpd5")
#' 
#' Data=matrix(sample(15, size = 50*5, replace = TRUE), nrow = 50, ncol = 5)
#' colnames(Data)=colnames(Data, do.NULL = FALSE, prefix = "col")
#' rownames(Data)=rownames(Data, do.NULL = FALSE, prefix = "row")
#' for(i in 1:50){
#' 	rownames(Data)[i]=paste("Cpd",i,sep="")
#' }
#' 
#' ContFeaturesPlot(leadCpds=Comps,orderLab=rownames(Data),colorLab=NULL,data=Data,
#' nrclusters=7,cols=Colors1,ylab="features",addLegend=TRUE,margins=c(5.5,3.5,0.5,8.7),
#' plottype="new",location=NULL)
#' }
ContFeaturesPlot<-function(leadCpds,data,nrclusters=NULL,orderLab=NULL,colorLab=NULL,cols=NULL,ylab="features",addLegend=TRUE,margins=c(5.5,3.5,0.5,8.7),plottype="new",location=NULL){
	
	if(all(leadCpds%in%rownames(data))){
		data=t(data)
	}
	
	if(!is.null(orderLab)){
		if(class(orderLab)=="character"){
			orderlabs=orderLab
		}
		else{
			orderlabs=orderLab$Clust$order.lab
			data=data[,match(orderlabs,colnames(data))]
		}
	}
	else{
		orderlabs=colnames(data)
	}
	
	temp=orderlabs[which(!(orderlabs%in%leadCpds))]
	AllCpds=c(leadCpds,temp)
	
	if(is.null(dim(data))){
		data=t(as.matrix(data))
		rownames(data)="Feature"
	}
	data=data[,AllCpds,drop=FALSE]
	
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	
	
	plottypein(plottype,location)
	graphics::par(mar=margins)
	graphics::plot(x=0,y=0,xlim=c(0,(ncol(data)+3)),ylim=c(min(data)-0.5,max(data)+0.5),type="n",ylab=ylab,xlab='',xaxt='n')
	for(i in c(1:nrow(data))){	
		graphics::lines(x=seq(1,length(leadCpds)),y=data[i,which(colnames(data)%in%leadCpds)],col=i)
	}
	for(i in c(1:nrow(data))){	
		graphics::lines(x=seq(length(leadCpds)+4,(ncol(data)+3)),y=data[i,which(!(colnames(data)%in%leadCpds))],col=i)
	}
	
	
	
	if(!(is.null(colorLab)) | is.null(nrclusters)){
		Data1 <- colorLab$Clust
		ClustData1=stats::cutree(Data1,nrclusters) 
		
		ordercolors=ClustData1[Data1$order]
		names(ordercolors)=Data1$order.lab
		
		ClustData1=ClustData1[Data1$order]	
		
		
		order=seq(1,nrclusters)
		
		for (k in 1:length(unique(ClustData1))){
			select=which(ClustData1==unique(ClustData1)[k])
			ordercolors[select]=order[k]
		}
		
		colors<- cols[ordercolors]
		names(colors) <-names(ordercolors)	
	}
	else{
		colors1<-rep("green",length(leadCpds))
		colors2<-rep("black",length(temp))
		colors=c(colors1,colors2)
		names(colors)=AllCpds
	}
	
	graphics::mtext(leadCpds,side=1,at=seq(1,length(leadCpds)),line=0.6,las=2,cex=0.70,col=colors[leadCpds])
	graphics::mtext(temp, side = 1, at=c(seq(length(leadCpds)+4,(ncol(data)+3))), line=0.5, las=2, cex=0.70,col=colors[temp])
	if(addLegend==TRUE){
		
		labels=rownames(data)
		colslegend=seq(1,length(rownames(data)))
		
		graphics::par(xpd=T,mar=margins)
		graphics::legend(ncol(data)+5,mean(c(min(c(min(data)-0.5,max(data)+0.5)),max(c(min(data)-0.5,max(data)+0.5)))),legend=c(labels),col=c(colslegend),lty=1,lwd=3,cex=0.8)
		
	}
	plottypeout(plottype)
}

#' @title Determines an optimal weight for weighted clustering by silhouettes widths.
#' 
#' @description The function \code{DetermineWeight_SilClust} determines an optimal weight
#' for weighted similarity clustering by calculating silhouettes widths. See
#' "Details" for a more elaborate description.
#' 
#' @details For each given weight, a linear combination of the distance matrices of the
#' single data sources is obtained. For these distance matrices, medoid
#' clustering with nrclusters is set up by the \code{pam} function of the
#' \pkg{cluster} and the silhouette widths are retrieved. These widths
#' indicates how well an object fits in its current cluster. Values around one
#' indicate an appropriate cluster. The silhouette widths are regressed in
#' function of the cluster membership determined by the objects. First, in
#' function of the cluster membership determined by the weighted combination.
#' Then, also in function of the cluster membership determined by the single
#' source clustering. The regression function is fit by the \code{lm} function
#' and the \code{r.squared} value is retrieved. The\code{r.squared} value
#' indicates how much of the variance of the silhouette widths is explained by
#' the membership. Optimally this value is high.
#' 
#' Next, a statistic is determined. Suppose that RWW is the \code{r.squared}
#' retrieved from regressing the weighted silhouette widths versus the weighted
#' cluster membership and RWX the \code{r.squared} retrieved from regressing
#' the weighted silhouette widths versus the cluster membership determined by
#' data X.  If M is total number of data sources, than statistic is obtained
#' as: \deqn{Stat=abs(M*RWW-\sum{RWX})}
#' 
#' The lower the statistical value, the better the weighted clustering is
#' explained by the single data sources. The goal is to obtain the weights for
#' which this value is minimized. Via bootstrapping a p-value is obtained for
#' every statistic.
#' 
#' The weight combinations should be provided as elements in a list. For three
#' data matrices an example could be:
#' weights=list(c(0.5,0.2,0.3),c(0.1,0.5,0.4)). To provide a fixed weight for
#' some data sets and let it vary randomly for others, the element "x"
#' indicates a free parameter. An example is weights=list(c(0.7,"x","x")). The
#' weight 0.7 is now fixed for the first data matrix while the remaining 0.3
#' weight will be divided over the other two data sets. This implies that every
#' combination of the sequence from 0 to 0.3 with steps of 0.1 will be reported
#' and clustering will be performed for each.
#' 
#' @param List A list of matrices of the same type. It is assumed the rows are
#' corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clusters".	
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}. 
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param weight Optional. A list of different weight combinations for the data sets in List.
#' If NULL, the weights are determined to be equal for each data set.
#' It is further possible to fix weights for some data matrices and to
#' let it vary randomly for the remaining data sets. Defaults to seq(1,0,-0.1). An example is provided in the details.
#' @param nrclusters The number of clusters to cut the dendrogram in. This is
#' necessary for the computation of the Jaccard coefficient. Default is NULL.
#' @param names The labels to give to the elements in List. Default is NULL.
#' @param nboot Number of bootstraps to be run. Default is 10.
#' @param StopRange Logical. Indicates whether the distance matrices with
#' values not between zero and one should be standardized to have so. If FALSE
#' the range normalization is performed. See \code{Normalization}. If TRUE, the
#' distance matrices are not changed. This is recommended if different types of
#' data are used such that these are comparable. Default is FALSE.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is FALSE.
#' @return Two plots are made: one of the statistical values versus the weights
#' and one of the p-values versus the weights. Further, a list with two
#' elements is returned: \item{Result}{A data frame with the statistic for each
#' weight combination} \item{Weight}{The optimal weight}
#' @examples
#' 
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' L=list(MCF7_F,MCF7_T)
#' 
#' MC7_Weight=DetermineWeight_SilClust(List=L,type="clusters",distmeasure=
#' c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),
#' weight=seq(0,1,by=0.01),nrclusters=c(7,7),names=c("FP","TP"),nboot=10,
#' StopRange=FALSE,plottype="new",location=NULL)
#' }
#' 
#' @export DetermineWeight_SilClust
DetermineWeight_SilClust<-function(List,type=c("data","dist","clusters"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),weight=seq(0,1,by=0.01),nrclusters=NULL,names=NULL,nboot=10,StopRange=FALSE,plottype="new",location=NULL){
	
	type=match.arg(type)
	
	if(is.null(names)){
		names=c()
		for(i in 1:length(List)){
			names=c(names,paste("Method",i,sep=" "))
		}	
	}
	
	CheckDist<-function(Dist,StopRange){
		if(StopRange==FALSE  & !(0<=min(Dist) & max(Dist)<=1)){
			message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
			Dist=Normalization(Dist,method="Range")
		}
		else{
			Dist=Dist
		}
	}
	
	
	if(type=="data"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}
		Dist=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i],normalize[i],method[i]))
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		
		silwidth=lapply(Dist,function(x) cluster::pam(x,nrclusters[i])$silinfo$widths)
		names(silwidth)=names
		
	}
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		Dist=List
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		silwidth=lapply(Dist,function(x) cluster::pam(x,nrclusters[i])$silinfo$widths)
		names(silwidth)=names
	}
	else{
		Dist=lapply(seq(length(List)),function(i) return(List[[i]]$DistM))
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		OrderNames=rownames(Dist[[1]])
		for(i in 1:length(Dist)){
			Dist[[i]]=Dist[[i]][OrderNames,OrderNames]
		}
		silwidth=lapply(Dist,function(x) cluster::pam(x,nrclusters[i])$silinfo$widths)
		names(silwidth)=names
	}
	
	
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	
	
	namesw=c()
	for(i in 1:length(names)){
		namesw=c(namesw,paste("w_",names[i],sep=""))
	}
	
	labels<-c(namesw,"Observed Statistic","P-Value")
	
	ResultsWeight<-matrix(0,ncol=length(labels),nrow=length(weight))
	colnames(ResultsWeight)=labels
	
	if(is.null(weight)){
		equalweights=1/length(List)
		weight=list(rep(equalweights,length(List)))		
	}
	else if(class(weight)=='list' & length(weight[[1]])!=length(List)){
		stop("Give a weight for each data matrix or specify a sequence of weights")
	}
	else{
		message('The weights are considered to be a sequence, each situation is investigated')
	}
	
	if(class(weight)!="list"){
		condition<-function(l){		
			l=as.numeric(l)
			if( sum(l)==1 ){  #working wit characters since with the numeric values of comb or permutations something goes not the way is should: 0.999999999<0.7+0.3<1??
				#return(row.match(l,t1))
				return(l)
			}
			else(return(0))
		}
		t1=gtools::permutations(n=length(weight),r=length(List),v=as.character(weight),repeats.allowed = TRUE)
		t2=lapply(seq_len(nrow(t1)), function(i) if(sum(as.numeric(t1[i,]))==1) return(as.numeric(t1[i,])) else return(0)) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
		t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
		weight=t2[which(t3!=0)]
	}
	
	if(class(weight)=="list" & "x" %in% weight[[1]]){ #x indicates a free weight
		for(i in 1:length(weight)){
			w=weight[[i]]
			weightsfordata=which(w!="x") #position of the provided weight = position of the data to which the weight is given
			givenweights=as.numeric(w[weightsfordata])
			
			stilltodistribute=1-sum(givenweights)
			
			newweights=seq(stilltodistribute,0,-0.1)
			
			t1=gtools::permutations(n=length(newweights),r=length(List)-length(weightsfordata),v=as.character(newweights),repeats.allowed = TRUE)
			Input1=as.list(seq_len(nrow(t1)))
			Input2=lapply(seq(length(Input1)),function(i) {Input1[[i]][length(Input1[[i]])+1]=stilltodistribute
						return(Input1[[i]])})
			t2=lapply(seq(length(Input2)), FUN=function(i){if(sum(as.numeric(t1[Input2[[i]][1],])+0.00000000000000002775)==Input2[[i]][2]) return(as.numeric(t1[i,])) else return(0)}) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
			t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
			weightsforotherdata=t2[which(t3!=0)]
			
			new=list()
			for(i in 1:length(weightsforotherdata)){
				w1=weightsforotherdata[[i]]
				new[[i]]=rep(0,length(List))
				new[[i]][weightsfordata]=givenweights
				new[[i]][which(new[[i]]==0)]=w1
			}
			
			weight=new
		}
	}
	
	weightedcomb<-function(w,Dist){
		temp=lapply(seq_len(length(Dist)),function(i) w[i]*Dist[[i]])
		temp=Reduce("+",temp)		
	}
	
	DistW=lapply(weight,weightedcomb,Dist)
	
	nrclus=ceiling(mean(nrclusters))
	
	silwidthW=lapply(DistW,function(x) cluster::pam(x,nrclus)$silinfo$widths)
	
	
	StatRSq<-function(silwidthW,silwidth,ordernames,names){
		n=length(silwidth)
		
		regressRSq<-function(silwidthweight,silwidth,ordernames,names){
			
			L1=silwidthweight[,3][ordernames]		
			L2W=silwidthweight[,1][ordernames]
			
			regressWW<-stats::lm(L1~L2W)
			RsqWW<-summary(regressWW)$r.squared
			
			regressWX<-function(silw,L1,ordernames){
				L2<-silw[,1][ordernames]
				
				regresWX<-stats::lm(L1~L2)
				RsqWX<-summary(regresWX)$r.squared
				return(RsqWX)
			}
			
			RsqWX<-sapply(silwidth,regressWX,L1=L1,ordernames=ordernames)
			
			Rsq=c(RsqWW,RsqWX)
			return(Rsq)
		}
		
		
		RSqs=lapply(c(1:length(silwidthW)),function(x) regressRSq(silwidthW[[x]],silwidth,ordernames,names))
		
		statfunction<-function(RS){
			stat=0
			xx=0
			xy=0
			for(i in 1:length(RS)){
				if(i==1){
					xx=xx+RS[i]
				}
				else{
					xy=xy+RS[i]
				}
			}
			stat=abs(n*xx-xy)  #check this formula with Nolen
			names(stat)=NULL
			return(stat)
		}
		
		Stats=sapply(RSqs,statfunction)	
		
	}	
	
	StatRSqObs=StatRSq(silwidthW,silwidth,ordernames=rownames(Dist[[1]]),names)
	
	#bootstrapping
	statNULL=matrix(0,nrow=length(weight),ncol=nboot)
	perm.rowscols <- function (D, n) 
	{
		s <- sample(1:n)
		D=D[s, s]
		return(D)
	}
	
	for(i in 1:nboot){
		set.seed(i)
		DistNULL=Dist
		DistNULL[[1]] <- perm.rowscols(DistNULL[[1]],nrow(DistNULL[[1]]))
		silwidthNULL=lapply(1:length(DistNULL),function(x) cluster::pam(DistNULL[[x]],nrclusters[x])$silinfo$widths)
		
		
		DistWNULL=lapply(weight,weightedcomb,DistNULL)
		silwidthWNULL=lapply(DistWNULL,function(x) cluster::pam(x,nrclus)$silinfo$widths)
		
		statNULL[,i]=StatRSq(silwidthWNULL,silwidthNULL,ordernames=rownames(DistNULL[[1]]),names)
	}
	
	PVals=lapply(c(1:nrow(statNULL)),function(x) (1+sum(abs(statNULL[x,])<=abs(StatRSqObs[x])))/(nboot+1))
	
	ResultsWeight=t(mapply(c,weight,StatRSqObs,PVals))
	colnames(ResultsWeight)=labels
	
	#Choose weight with smallest observed test statistic
	
	Weight=ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),c(1:length(List))]
	
	plottypein(plottype,location)
	graphics::plot(x=ResultsWeight[,1],y=ResultsWeight[,"Observed Statistic"],xlim=c(0,max(ResultsWeight[,1])),ylim=c(min(ResultsWeight[,"Observed Statistic"]),max(ResultsWeight[,"Observed Statistic"])),xlab="",ylab="Observed Statistic",pch=19,col="black")
	graphics::points(ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),1],ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),"Observed Statistic"],pch=19,col="red")
	graphics::mtext("Weight Combinations", side=1, line=4)
	graphics::axis(1,labels=paste("Optimal weights:", paste(Weight,collapse=", "),sep=" "), at=ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),1],line=2)
	plottypeout(plottype)
	
	plottypein(plottype,location)
	graphics::plot(x=ResultsWeight[,1],y=ResultsWeight[,"P-Value"],xlim=c(0,max(ResultsWeight[,1])),ylim=c(min(ResultsWeight[,"P-Value"]),max(ResultsWeight[,"P-Value"])),xlab="",ylab="P-Value",pch=19,col="black")
	graphics::points(ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),1],ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),"P-Value"],pch=19,col="red")
	graphics::mtext("Weight Combinations", side=1, line=4)
	graphics::axis(1,labels=paste("Optimal weights:", paste(Weight,collapse=", "),sep=" "), at=ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),1],line=2)
	plottypeout(plottype)
	
	
	out=list()
	out[[1]]=ResultsWeight
	out[[2]]=Weight
	names(out)=c("Result","Weight")
	
	return(out)
	
}


#' @title Determines an optimal weight for weighted clustering by similarity weighted
#' clustering.
#' 
#' @description The function \code{DetermineWeight_SimClust} determines an optimal weight
#' for performing weighted similarity clustering on by applying similarity
#' clustering. For each given weight, is each separate clustering compared to
#' the clustering on a weighted dissimilarity matrix and a Jaccard coefficient
#' is calculated. The ratio of the Jaccard coefficients closets to one
#' indicates an optimal weight.
#' 
#' @details If the type of List is data, an hierarchical clustering is performed on each
#' data matrix separately. After obtaining clustering results for the two data
#' matrices, the distance matrices are extracted. If these are not calculated
#' with the same distance measure, they are normalized to be in the same range.
#' For each weight, a weighted linear combination of the distance matrices is
#' taken and hierarchical clustering is performed once again. The resulting
#' clustering is compared to each of the separate clustering results and a
#' Jaccard coefficient is computed. The ratio of the Jaccard coefficients
#' closets to one, indicates an optimal weight. A plot of all the ratios is
#' produced with an extra indication for the optimal weight.
#' 
#' The weight combinations should be provided as elements in a list. For three
#' data matrices an example could be:
#' weights=list(c(0.5,0.2,0.3),c(0.1,0.5,0.4)). To provide a fixed weight for
#' some data sets and let it vary randomly for others, the element "x"
#' indicates a free parameter. An example is weights=list(c(0.7,"x","x")). The
#' weight 0.7 is now fixed for the first data matrix while the remaining 0.3
#' weight will be divided over the other two data sets. This implies that every
#' combination of the sequence from 0 to 0.3 with steps of 0.1 will be reported
#' and clustering will be performed for each.
#' 
#' @param List A list of matrices of the same type. It is assumed the rows are
#' corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clusters".	
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}. 
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param weight Optional. A list of different weight combinations for the data sets in List.
#' If NULL, the weights are determined to be equal for each data set.
#' It is further possible to fix weights for some data matrices and to
#' let it vary randomly for the remaining data sets. Defaults to seq(1,0,-0.1). An example is provided in the details.
#' @param nrclusters The number of clusters to cut the dendrogram in. This is
#' necessary for the computation of the Jaccard coefficient. Default is NULL.
#' @param clust Choice of clustering function (character). Defaults to "agnes". 
#' @param linkage Choice of inter group dissimilarity (character) for the individual clusterings. Defaults to c("flexible","flexible").
#' @param linkageF Choice of inter group dissimilarity (character) for the final clustering. Defaults to "ward".
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible".
#' @param gap Logical. Whether or not to calculate the gap statistic in the
#' clustering on each data matrix separately. Only if type="data". Default is FALSE.
#' @param maxK The maximal number of clusters to consider in calculating the
#' gap statistic. Only if type="data". Default is 15.
#' @param names The labels to give to the elements in List. Default is NULL.
#' @param StopRange Logical. Indicates whether the distance matrices with
#' values not between zero and one should be standardized to have so. If FALSE
#' the range normalization is performed. See \code{Normalization}. If TRUE, the
#' distance matrices are not changed. This is recommended if different types of
#' data are used such that these are comparable. Default is FALSE.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is FALSE.
#' @return The returned value is a list with three elements:
#' \item{ClustSep}{The result of \code{Cluster} for each single element of
#' List} \item{Result}{A data frame with the Jaccard coefficients and their
#' ratios for each weight} \item{Weight}{The optimal weight}
#' @references 
#' \insertRef{PerualilaTan2016}{IntClust}
#' @examples
#' 
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",alpha=0.625,gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",alpha=0.625,gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' L=list(MCF7_F,MCF7_T)
#' 
#' MCF7_Weight=DetermineWeight_SimClust(List=L,type="clusters",weight=seq(0,1,by=0.01),
#' nrclusters=c(7,7),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),
#' method=c(NULL,NULL),clust="agnes",linkage=c("flexible","flexible"),linkageF="ward",
#' alpha=0.625,gap=FALSE,maxK=50,names=c("FP","TP"),StopRange=FALSE,plottype="new",location=NULL)
#' }
#' 
#' @export DetermineWeight_SimClust
DetermineWeight_SimClust<-function(List,type=c("data","dist","clusters"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),weight=seq(0,1,by=0.01),nrclusters=NULL,clust="agnes",linkage=c("flexible","flexible"),linkageF="ward",alpha=0.625,gap=FALSE,maxK=15,names=NULL,StopRange=FALSE,plottype="new",location=NULL){
	
	CheckDist<-function(Dist,StopRange){
		if(StopRange==FALSE & !(0<=min(Dist) & max(Dist)<=1)){
			message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
			Dist=Normalization(Dist,method="Range")
		}
		else{
			Dist=Dist
		}
	}
	
	
	type<-match.arg(type)
	if(type=="data"){

		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize[i],method[i],clust,linkage[i],alpha,gap,maxK,StopRange))
		
		Dist=lapply(seq(length(List)),function(i) Clusterings[[i]]$DistM)
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		
		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=ceiling(mean(clusters))
			}
		}
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		out<-list(Clusterings)
		names(out)="ClusterSep"		
	}
	
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		
		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize[i],method[i],clust,linkage[i],alpha,gap,maxK,StopRange))
		
		Dist=List
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		
		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters or put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=ceiling(mean(clusters))
			}
		}
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		out<-list(Clusterings)
		names(out)="ClusterSep"		
		
	}
	else{
		
		Dist=lapply(seq(length(List)),function(i) return(List[[i]]$DistM))
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		
		OrderNames=rownames(Dist[[1]])
		for(i in 1:length(Dist)){
			Dist[[i]]=Dist[[i]][OrderNames,OrderNames]
		}
		
		Clusterings=List
		
		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')			
		}
		out<-list(Clusterings)
		names(out)="ClusterSep"		
		
	}
	
	
	namesw=c()
	for(i in 1:length(names)){
		namesw=c(namesw,paste("w_",names[i],sep=""))
	}
	namesJ=c()
	for(i in 1:length(names)){
		namesJ=c(namesJ,paste("J(sim",names[i],",simW)",sep=""))
	}
	namesR=c()
	combs=utils::combn(seq(length(List)),m=2,simplify=FALSE)
	for(i in 1:length(combs)){
		namesR=c(namesR,paste("J_",names[combs[[i]][1]],"/J_",names[combs[[i]]][2],sep=""))
	}
	
	labels<-c(namesw,namesJ,namesR)
	
	ResultsWeight<-matrix(0,ncol=length(labels),nrow=length(weight))
	#data.frame(col1=numeric(),col2=numeric(),col3=numeric(),col4=numeric())
	colnames(ResultsWeight)=labels
	
	if(is.null(weight)){
		equalweights=1/length(List)
		weight=list(rep(equalweights,length(List)))		
	}
	else if(class(weight)=='list' & length(weight[[1]])!=length(List)){
		stop("Give a weight for each data matrix or specify a sequence of weights")
	}
	else{
		message('The weights are considered to be a sequence, each situation is investigated')
	}
	
	if(class(weight)!="list"){
		condition<-function(l){		
			l=as.numeric(l)
			if( sum(l)==1 ){  #working wit characters since with the numeric values of comb or permutations something goes not the way is should: 0.999999999<0.7+0.3<1??
				#return(row.match(l,t1))
				return(l)
			}
			else(return(0))
		}
		t1=gtools::permutations(n=length(weight),r=length(List),v=as.character(weight),repeats.allowed = TRUE)
		t2=lapply(seq_len(nrow(t1)), function(i) if(sum(as.numeric(t1[i,]))==1) return(as.numeric(t1[i,])) else return(0)) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
		t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
		weight=t2[which(t3!=0)]
	}
	
	if(class(weight)=="list" & "x" %in% weight[[1]]){ #x indicates a free weight
		for(i in 1:length(weight)){
			w=weight[[i]]
			weightsfordata=which(w!="x") #position of the provided weight = position of the data to which the weight is given
			givenweights=as.numeric(w[weightsfordata])
			
			stilltodistribute=1-sum(givenweights)
			
			newweights=seq(stilltodistribute,0,-0.1)
			
			t1=gtools::permutations(n=length(newweights),r=length(List)-length(weightsfordata),v=as.character(newweights),repeats.allowed = TRUE)
			Input1=as.list(seq_len(nrow(t1)))
			Input2=lapply(seq(length(Input1)),function(i) {Input1[[i]][length(Input1[[i]])+1]=stilltodistribute
						return(Input1[[i]])})
			t2=lapply(seq(length(Input2)), FUN=function(i){if(sum(as.numeric(t1[Input2[[i]][1],])+0.00000000000000002775)==Input2[[i]][2]) return(as.numeric(t1[i,])) else return(0)}) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
			t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
			weightsforotherdata=t2[which(t3!=0)]
			
			new=list()
			for(i in 1:length(weightsforotherdata)){
				w1=weightsforotherdata[[i]]
				new[[i]]=rep(0,length(List))
				new[[i]][weightsfordata]=givenweights
				new[[i]][which(new[[i]]==0)]=w1
			}
			
			weight=new
		}
	}
	weightedcomb<-function(w,Dist){
		temp=lapply(seq_len(length(Dist)),function(i) w[i]*Dist[[i]])
		temp=Reduce("+",temp)		
	}
	DistM=lapply(weight,weightedcomb,Dist)
	
	hclustOr=lapply(seq(length(List)),function(i) stats::cutree(Clusterings[[i]]$Clust,nrclusters[i]))
	nrclus=ceiling(mean(nrclusters))
	hclustW=lapply(seq(length(weight)),function(i) stats::cutree(cluster::agnes(DistM[[i]],diss=TRUE,method=linkageF),nrclus))
	
	Counts=function(clusterlabs1,clusterlabs2){
		index=c(1:length(clusterlabs1))
		allpairs=utils::combn(index,2,simplify=FALSE)  #all pairs of indices: now check clutserlabels for every pair==> only 1 for loop
		n11=n10=n01=n00=0
		
		counts<-function(pair){
			if(clusterlabs1[pair[1]]==clusterlabs1[pair[2]]){
				if(clusterlabs2[pair[1]]==clusterlabs2[pair[2]]){
					n11=n11+1
				}
				else{
					n10=n10+1						
				}
			}
			else{
				if(clusterlabs2[pair[1]]==clusterlabs2[pair[2]]){
					n01=n01+1
				}
				else{
					n00=n00+1
				}
				
			}
			return(c(n11,n10,n01,n00))		
		}
		
		n=lapply(seq(length(allpairs)),function(i) counts(allpairs[[i]]))
		nn=Reduce("+",n)
		#2: compute jaccard coefficient	
		Jac=nn[1]/(nn[1]+nn[2]+nn[3])
		return(Jac)
	}
	
	Jaccards<-function(hclust){
		jacs=lapply(seq(length(hclustOr)),function(i) Counts(clusterlabs1=hclustOr[[i]],clusterlabs2=hclust))
		return(unlist(jacs))
	}
	
	AllJacs=lapply(hclustW,Jaccards)  #make this faster:lapply + transfrom to data frame with plyr package
	
	Ratios<-function(Jacs){	
		combs=utils::combn(seq(length(List)),m=2,simplify=FALSE)
		ratio<-function(v,Jacs){
			return(Jacs[v[1]]/Jacs[v[2]])
		}
		
		ratios=lapply(seq(length(combs)),function(i) ratio(v=combs[[i]],Jacs=Jacs))
		
	}
	
	AllRatios=lapply(seq(length(AllJacs)),function(i) unlist(Ratios(AllJacs[[i]])))
	
	
	ResultsWeight=t(mapply(c,weight,AllJacs,AllRatios))
	colnames(ResultsWeight)=labels
	
	#Choose weight with ratio closest to one==> smallest where this happens: ##### START HERE WITH OPTIMIZATION #####
	ResultsWeight=cbind(ResultsWeight,rep(0,nrow(ResultsWeight)))
	colnames(ResultsWeight)[ncol(ResultsWeight)]="trick"
	Weight=ResultsWeight[which.min(rowSums(abs(ResultsWeight[,c(namesR,"trick")]-1))),c(1:length(List))]
	
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	
	plottypein(plottype,location)
	graphics::plot(x=0,y=0,type="n",xlim=c(0,dim(ResultsWeight)[1]),ylim=c(min(ResultsWeight[,namesR]),max(ResultsWeight[,namesR])),xlab="",ylab="Ratios")
	if(is.null(ncol(ResultsWeight[,namesR]))){
		L=1
	}
	else{
		L=ncol(ResultsWeight[,namesR])
	}
	for(i in 1:L){
		graphics::lines(1:dim(ResultsWeight)[1],y=ResultsWeight[,namesR[i]],col=i)
	}
	graphics::abline(h=0,v=which.min(rowSums(abs(ResultsWeight[,c(namesR,"trick")]-1))),col="black",lwd=2)
	graphics::mtext("Weight Combinations", side=1, line=3)
	graphics::axis(1,labels=paste("Optimal weights:", paste(Weight,collapse=", "),sep=" "), at=which.min(rowSums(abs(ResultsWeight[,c(namesR,"trick")]-1))),line=1,tck=1,lwd=2)
	plottypeout(plottype)
	
	ResultsWeight=ResultsWeight[,-ncol(ResultsWeight)]
	out[[2]]=ResultsWeight
	out[[3]]=Weight
	names(out)=c("ClusterSep","Result","Weight")
	
	
	
	return(out)
	
}

#' @title Differential gene expressions for multiple results
#' 
#' @description The function \code{DiffGenes} will, given the output of a certain method,
#' look for genes that are differentially expressed for each cluster by
#' applying the limma function to that cluster and compare it to all other
#' clusters simultaneously. If a list of outputs of several methods is
#' provided, DiffGenes will perform the limma function for each method.
#' @export DiffGenes
#' @details The function rearranges the clusters of the methods to a reference method
#' such that a comparison is made easier.  Given a list of methods, it calls
#' upon \code{ReorderToReference} to rearrange the number of clusters according
#' to the first element of the list which will be used as the reference.
#' 
#' @param List A list of the clustering outputs to be compared. The first
#' element of the list will be used as the reference in
#' \code{ReorderToReference}.
#' @param Selection If differential gene expression should be investigated for
#' a specific selection of objects, this selection can be provided here.
#' Selection can be of the type "character" (names of the objects) or
#' "numeric" (the number of specific cluster). Default is NULL.
#' @param geneExpr The gene expression matrix or ExpressionSet of the objects.
#' The rows should correspond with the genes.
#' @param nrclusters Optional. The number of clusters to cut the dendrogram in.
#' The number of clusters should not be specified if the interest lies only in
#' a specific selection of objects which is known by name.  Otherwise, it is
#' required. Default is NULL.
#' @param method The method to applied to look for DE genes. For now, only the
#' limma method is available. Default is "limma".
#' @param sign The significance level to be handled. Default is 0.05.
#' @param topG Overrules sign. The number of top genes to be shown. Default is NULL.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' @param names Optional. Names of the methods. Default is NULL.
#' @return The returned value is a list with an element per method. Each
#' element contains a list per cluster with the following elements:
#' \item{objects}{A list with the elements LeadCpds (the objects of
#' interest) and OrderedCpds (all objects in the order of the clustering
#' result)} \item{Genes}{A list with the elements TopDE (a table with
#' information on the top genes) and AllDE (a table with information on all
#' genes)}
#' @references SMYTH, G. K. (2004). Linear models and empirical Bayes methods
#' for assessing differential expression in microarray experiments. Statistical
#' Applications in Genetics and Molecular Biology. 3(1).
#' @examples
#' 
#' data(fingerprintMat)
#' data(targetMat)
#' data(geneMat)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' L=list(MCF7_T ,MCF7_F)
#' 
#' MCF7_FT_DE = DiffGenes(List=L,geneExpr=geneMat,nrclusters=7,method="limma",
#' sign=0.05,topG=10,fusionsLog=TRUE,weightclust=TRUE)
DiffGenes=function(List,Selection=NULL,geneExpr=NULL,nrclusters=NULL,method="limma",sign=0.05,topG=NULL,fusionsLog=TRUE,weightclust=TRUE,names=NULL){
	if(method != "limma"){
		stop("Only the limma method is implemented to find differentially expressed genes")
	} 	
	if(!is.null(Selection)){
		ResultLimma=DiffGenesSelection(List,Selection,geneExpr,nrclusters,method,sign,topG,fusionsLog,weightclust,names)
	}
	else{
		
		ListNew=list()
		element=0
		for(i in 1:length(List)){
			if(attributes(List[[i]])$method != "CEC" & attributes(List[[i]])$method != "Weighted" & attributes(List[[i]])$method!= "WeightedSim"){
				ResultsClust=list()
				ResultsClust[[1]]=list()
				ResultsClust[[1]][[1]]=List[[i]]
				names(ResultsClust[[1]])[1]="Clust"
				element=element+1					
				ListNew[[element]]=ResultsClust[[1]]
				#attr(ListNew[element],"method")="Weights"
			}
			else if(attributes(List[[i]])$method=="CEC" | attributes(List[[i]])$method=="Weighted" | attributes(List[[i]])$method == "WeightedSim"){
				ResultsClust=list()
				if(weightclust==TRUE){
					ResultsClust[[1]]=list()
					if(attributes(List[[i]])$method != "WeightedSim"){
						ResultsClust[[1]][[1]]=List[[i]]$Clust
						names(ResultsClust[[1]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[1]]
						attr(ListNew[element],"method")="Weights"
					}
					else{
						ResultsClust[[1]]=list()
						ResultsClust[[1]][[1]]=List[[i]]
						names(ResultsClust[[1]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[1]]
					}
				}
				else{
					for (j in 1:length(List[[i]]$Results)){
						ResultsClust[[j]]=list()
						ResultsClust[[j]][[1]]=List[[i]]$Results[[j]]
						names(ResultsClust[[j]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[j]]
						attr(ListNew[element],"method")="Weights"
					}		
				}		
			}	
		}
		
		if(is.null(names)){
			names=seq(1,length(ListNew),1)
			for(i in 1:length(ListNew)){
				names[i]=paste("Method",i,sep=" ")
			}
		}
		
		if(is.null(names)){
			for(j in 1:length(List)){
				names[j]=paste("Method",j,sep=" ")	
			}
		}
		
		names(ListNew)=names
		
		if(is.null(topG)){
			top1=FALSE
		}
		else{
			top1=TRUE
		}	
		
		
		MatrixClusters=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)
		List=ListNew
		ResultLimma=list()
		maxclus=0
		for (k in 1:dim(MatrixClusters)[1]){
			clusters=MatrixClusters[k,]
			if(max(clusters)>maxclus){
				maxclus=max(clusters)
			}
			Genes=list()
			clust=sort(unique(clusters)) #does not matter: Genes[i] puts right elements on right places
			hc<-stats::as.hclust(List[[k]]$Clust$Clust)
			OrderedCpds <- hc$labels[hc$order]
			for (i in clust){
				
				temp=list()
				LeadCpds=names(clusters)[which(clusters==i)] 
				temp[[1]]=list(LeadCpds,OrderedCpds)
				names(temp[[1]])=c("LeadCpds","OrderedCpds") #names of the objects
				
				label = rep(0,length(names(clusters)))
				label[which(clusters==i)] = 1
				
				label.factor = factor(label)
				GeneExpr.2=geneExpr[,names(clusters)]
				
				if(class(GeneExpr.2)[1]=="ExpressionSet"){
					
					if (!requireNamespace("a4Base", quietly = TRUE)) {
						stop("a4Base needed for this function to work. Please install it.",
								call. = FALSE)
					}
					
					
					GeneExpr.2$LeadCmpds<-label.factor	
					DElead <- a4Base::limmaTwoLevels(GeneExpr.2,"LeadCpds")
					
					allDE <-a4Core::topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL),sort.by="p")
					
					if(is.null(allDE$ID)){
						allDE$ID<- rownames(allDE)
					}
					else
					{
						allDE$ID=allDE$ID
					}
					
					if(top1==TRUE){
						result = list(allDE[1:topG,],allDE)
						names(result)=c("TopDE","AllDE")
						
					}
					else if(top1==FALSE){
						topG=length(which(allDE$adj.P.Val<=sign))
						result = list(allDE[1:topG,],allDE)
						names(result)=c("TopDE","AllDE")
						
					}
					
				}
				else{
					
					design = stats::model.matrix(~label.factor)
					fit = limma::lmFit(GeneExpr.2,design=design)
					fit = limma::eBayes(fit)
					allDE=limma::topTable(fit,n=dim(geneExpr)[1],coef=2,adjust="fdr",sort.by="P")
					
					if(is.null(allDE$ID)){
						allDE$ID <- rownames(allDE)
					}
					else
					{
						allDE$ID=allDE$ID
					}
					
					if(top1==TRUE){
						result = list(allDE[1:topG,],allDE)
						names(result)=c("TopDE","AllDE")
					}
					else if(top1==FALSE){
						topG=length(which(allDE$adj.P.Val<=sign))
						result = list(allDE[1:topG,],allDE)
						names(result)=c("TopDE","AllDE")
					}
				}	
				
				temp[[2]]=result
				
				names(temp)=c("objects","Genes")
				
				Genes[[i]]=temp
				
				names(Genes)[i]=paste("Cluster",i,sep=" ")
			}
			ResultLimma[[k]]=Genes
			
		}
		names(ResultLimma)=names
		for(i in 1:length(ResultLimma)){
			for(k in 1:length(ResultLimma[[i]])){
				if(is.null(ResultLimma[[i]][[k]])[1]){
					ResultLimma[[i]][[k]]=NA
					names(ResultLimma[[i]])[k]=paste("Cluster",k,sep=" ")
				}			
			}
			if(length(ResultLimma[[i]]) != maxclus){
				extra=maxclus-length(ResultLimma[[i]])
				#temp=length(ResultLimma[[i]])
				for(j in 1:extra){
					ResultLimma[[i]][[length(ResultLimma[[i]])+1]]=NA
					names(ResultLimma[[i]])[length(ResultLimma[[i]])]=paste("Cluster",length(ResultLimma[[i]]),sep=" ")
				}
			}
		} 	
		
	}
	return(ResultLimma)
}

#' @title Differential expression for a selection of objects
#' @param List A list of the clustering outputs to be compared. The first
#' element of the list will be used as the reference in
#' \code{ReorderToReference}.
#' @param Selection If differential gene expression should be investigated for
#' a specific selection of objects, this selection can be provided here.
#' Selection can be of the type "character" (names of the objects) or
#' "numeric" (the number of specific cluster). Default is NULL.
#' @param geneExpr The gene expression matrix or ExpressionSet of the objects.
#' The rows should correspond with the genes.
#' @param nrclusters Optional. The number of clusters to cut the dendrogram in.
#' The number of clusters should not be specified if the interest lies only in
#' a specific selection of objects which is known by name.  Otherwise, it is
#' required. Default is NULL.
#' @param method The method to applied to look for DE genes. For now, only the
#' limma method is available. Default is "limma".
#' @param sign The significance level to be handled. Default is 0.05.
#' @param topG Overrules sign. The number of top genes to be shown. Default is NULL.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' @param names Optional. Names of the methods. Default is NULL.
#' @description Internal function of \code{DiffGenes}.
DiffGenesSelection=function(List,Selection,geneExpr=NULL,nrclusters=NULL,method="limma",sign=0.05,topG=NULL,fusionsLog=TRUE,weightclust=TRUE,names=NULL){
	if(method != "limma"){
		stop("Only the limma method is implemented to find differentially expressed genes")
	} 
	
	if(is.null(topG)){
		top1=FALSE
	}
	else{
		top1=TRUE
	}	
	
	if(class(Selection)=="character"){
		ResultLimma=list()
		Genes=list()
		temp=list()
		
		LeadCpds=Selection #names of the objects
		OrderedCpds=colnames(geneExpr)
		temp[[1]]=list(LeadCpds,OrderedCpds)
		names(temp[[1]])=c("LeadCpds","OrderedCpds")
		
		label = rep(0,dim(geneExpr)[2])
		label[which(colnames(geneExpr)%in%Selection)] = 1
		label.factor = factor(label)
		
		
		if(class(geneExpr)[1]=="ExpressionSet"){
			
			if (!requireNamespace("a4Base", quietly = TRUE)) {
				stop("a4Base needed for this function to work. Please install it.",
						call. = FALSE)
			}
			
			geneExpr$LeadCmpds<-label.factor 
			DElead <- a4Base::limmaTwoLevels(geneExpr,"LeadCpds")
			
			allDE <- a4Core::topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL),sort.by="p")
			if(is.null(allDE$ID)){
				allDE$ID <- rownames(allDE)
			}
			else
			{
				allDE$ID=allDE$ID
			}
			if(top1==TRUE){
				result = list(allDE[1:topG,],allDE)
				names(result)=c("TopDE","AllDE")
				
			}
			else if(top1==FALSE){
				topG=length(which(allDE$adj.P.Val<=sign))
				result = list(allDE[0:topG,],allDE)
				names(result)=c("TopDE","AllDE")
				
			}
			
		}
		else{
			
			design = stats::model.matrix(~label.factor)
			fit = limma::lmFit(geneExpr,design=design)
			fit = limma::eBayes(fit)
			
			allDE=limma::topTable(fit,coef=2,n=dim(geneExpr)[1],adjust="fdr",sort.by="P")
			if(is.null(allDE$ID)){
				allDE$ID <- rownames(allDE)
			}
			else
			{
				allDE$ID=allDE$ID
			}
			if(top1==TRUE){
				result = list(allDE[0:topG,],allDE)
				names(result)=c("TopDE","AllDE")
				
			}
			else if(top1==FALSE){
				topG=length(which(allDE$adj.P.Val<=sign))
				result = list(allDE[0:topG,],allDE)
				names(result)=c("TopDE","AllDE")
				
			}	
			
		}
		temp[[2]]=result
		
		names(temp)=c("objects","Genes")
		ResultLimma[[1]]=temp
		names(ResultLimma)="Selection"
		
	}
	else if(class(Selection)=="numeric" & !(is.null(List))){
		
		ListNew=list()
		element=0
		
		ListNew=list()
		element=0
		for(i in 1:length(List)){
			if(attributes(List[[i]])$method != "CEC" & attributes(List[[i]])$method != "Weighted" & attributes(List[[i]])$method!= "WeightedSim"){
				ResultsClust=list()
				ResultsClust[[1]]=list()
				ResultsClust[[1]][[1]]=List[[i]]
				names(ResultsClust[[1]])[1]="Clust"
				element=element+1					
				ListNew[[element]]=ResultsClust[[1]]
				#attr(ListNew[element],"method")="Weights"
			}
			else if(attributes(List[[i]])$method=="CEC" | attributes(List[[i]])$method=="Weighted" | attributes(List[[i]])$method == "WeightedSim"){
				ResultsClust=list()
				if(weightclust==TRUE){
					ResultsClust[[1]]=list()
					if(attributes(List[[i]])$method != "WeightedSim"){
						ResultsClust[[1]][[1]]=List[[i]]$Clust
						names(ResultsClust[[1]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[1]]
						attr(ListNew[element],"method")="Weights"
					}
					else{
						ResultsClust[[1]]=list()
						ResultsClust[[1]][[1]]=List[[i]]
						names(ResultsClust[[1]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[1]]
					}
				}
				else{
					for (j in 1:length(List[[i]]$Results)){
						ResultsClust[[j]]=list()
						ResultsClust[[j]][[1]]=List[[i]]$Results[[j]]
						names(ResultsClust[[j]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[j]]
						attr(ListNew[element],"method")="Weights"
					}		
				}		
			}	
		}
		
		if(is.null(names)){
			names=seq(1,length(ListNew),1)
			for(i in 1:length(ListNew)){
				names[i]=paste("Method",i,sep=" ")
			}
		}
		names(ListNew)=names
		
		Matrix=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)
		
		List=ListNew
		ResultLimma=list()
		for(k in 1:dim(Matrix)[1]){	
			cluster=Selection
			
			hc<-stats::as.hclust(List[[k]]$Clust$Clust)
			OrderedCpds <- hc$labels[hc$order]
			
			Genes=list()
			temp=list()
			LeadCpds=colnames(Matrix)[which(Matrix[k,]==cluster)] #names of the objects
			temp[[1]]=list(LeadCpds,OrderedCpds)
			names(temp[[1]])=c("LeadCpds","OrderedCpds")
			
			label = rep(0,dim(Matrix)[2])
			label[which(Matrix[k,]==cluster)] = 1
			label.factor = factor(label)
			
			GeneExpr.2=geneExpr[,colnames(Matrix)]
			
			if(class(GeneExpr.2)[1]=="ExpressionSet"){
				
				if (!requireNamespace("a4Base", quietly = TRUE)) {
					stop("a4Base needed for this function to work. Please install it.",
							call. = FALSE)
				}
				
				GeneExpr.2$LeadCmpds<-label.factor 
				DElead <- a4Base::limmaTwoLevels(GeneExpr.2,"LeadCpds")
				
				allDE <- a4Core::topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL),sort.by="p")
				if(is.null(allDE$ID)){
					allDE$ID <- rownames(allDE)
				}
				else
				{
					allDE$ID=allDE$ID
				}
				if(top1==TRUE){
					result = list(allDE[0:topG,],allDE)
					names(result)=c("TopDE","AllDE")
					
				}
				else if(top1==FALSE){
					topG=length(which(allDE$adj.P.Val<=sign))
					result = list(allDE[0:topG,],allDE)
					names(result)=c("TopDE","AllDE")
					
				}
				
			}
			else{
				
				
				design = stats::model.matrix(~label.factor)
				fit = limma::lmFit(GeneExpr.2,design=design)
				fit = limma::eBayes(fit)
				
				allDE=limma::topTable(fit,coef=2,n=dim(geneExpr)[1],adjust="fdr",sort.by="P")
				if(is.null(allDE$ID)){
					allDE$ID <- rownames(allDE)
				}
				else
				{
					allDE$ID=allDE$ID
				}
				if(top1==TRUE){
					result = list(allDE[1:topG,],allDE)
					names(result)=c("TopDE","AllDE")
					
				}
				else if(top1==FALSE){
					topG=length(which(allDE$adj.P.Val<=sign))
					
					result = list(allDE[0:topG,],allDE)
					names(result)=c("TopDE","AllDE")
					
				}	
				
			}
			temp[[2]]=result
			
			names(temp)=c("objects","Genes")
			ResultLimma[[k]]=temp
			names(ResultLimma)[k]=paste(names[k],": Cluster", cluster, sep="")
		}		
	}
	
	else{
		message("If a specific cluster is specified, clustering results must be provided in List")
	}
	return(ResultLimma)
	
}

#' @title Determine the distance in a heatmap
#' @param Data1 The resulting clustering of method 1.
#' @param Data2 The resulting clustering of method 2. 
#' @param names The names of the objects in the data sets. Default is NULL.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @description Internal function of \code{HeatmapPlot}
distanceheatmaps<-function(Data1,Data2,names=NULL,nrclusters=7){
	ClustData1=stats::cutree(Data1,nrclusters) #original clusters (aggregated data clustering)
	ClustData2=stats::cutree(Data2,nrclusters) #clusters of changed method
	
	ClustData1=ClustData1[Data1$order]
	ClustData2=ClustData2[Data2$order]
	
	trueorder1=sort(Data1$order,index.return = TRUE)
	trueorder2=sort(Data2$order,index.return = TRUE)
	
	ordercolors=ClustData1
	order=seq(1,nrclusters)
	
	for (k in 1:length(unique(ClustData1))){
		select=which(ClustData1==unique(ClustData1)[k])
		ordercolors[select]=order[k]
	}
	
	ordercolors2=ClustData2
	
	
	for (k in 1:length(unique(ClustData2))){
		select=which(ClustData2==unique(ClustData2)[k])
		ordercolors2[select]=order[k]
	}
	
	ClustData1=ordercolors[trueorder1$ix]
	ClustData2=ordercolors2[trueorder2$ix]
	
	out=matrix(0,length(ClustData1),length(ClustData2)) #want the rows to be the other method and the columns to be the aggregated data clustering
	
	#names=names[Data2$order]
	
	rownames(out)=names
	colnames(out)=names
	
	for(i in 1:length(names)){
		focus=names[i] #defines the column	
		
		label=ClustData2[i] #color of the cluster is defined by the original cluster that contains focus (1 to 7)		
		
		for(j in 1:length(names)){ #go over the rows
			other=names[j]		
			found=FALSE  #find cluster of other
			k=1
			while(found==FALSE & k<=nrclusters){
				label2=k
				if(other %in% names[ClustData1==label2]){
					othercluster=names[ClustData1==label2]
					found=TRUE
				}
				k=k+1
			}	
			
			if(focus %in% othercluster){ #if other and focus still together: give it color of cluster defined by focus
				out[j,i]=label
			}				
		}				
	}
	return(out)	
}

#' @title feature selection for a selection of objects
#' @param List A list of the clustering outputs to be compared. The first
#' element of the list will be used as the reference in
#' \code{ReorderToReference}.
#' @param Selection If differential gene expression should be investigated for
#' a specific selection of objects, this selection can be provided here.
#' Selection can be of the type "character" (names of the objects) or
#' "numeric" (the number of specific cluster). Default is NULL.
#' @param binData A list of the binary feature data matrices. These will be
#' evaluated with the fisher's extact test. Default is NULL.
#' @param contData A list of continuous data sets of the objects. These will
#' be evaluated with the t-test. Default is NULL.
#' @param datanames A vector with the names of the data matrices. Default is NULL.
#' @param nrclusters Optional. The number of clusters to cut the dendrogram in.
#' The number of clusters should not be specified if the interest lies only in
#' a specific selection of objects which is known by name.  Otherwise, it is
#' required. Default is NULL.
#' @param sign The significance level to be handled. Default is 0.05.
#' @param topChar Overrules sign. The number of features to display for each
#' cluster.  If not specified, only the significant genes are shown. Default is NULL.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' @description Internal function of \code{CharacteristicFeatures}.
FeatSelection<-function(List,Selection=NULL,binData=NULL,contData=NULL,datanames=NULL,nrclusters=NULL,topChar=NULL,sign=0.05,fusionsLog=TRUE,weightclust=TRUE){
	
	
	if(is.null(datanames)){
		for(j in 1:(length(binData)+length(contData))){
			datanames[j]=paste("Data",j,sep=" ")	
		}
	}
	
	if(class(Selection)=="character"){
		
		
		if(!is.null(binData)){
			cpdSet <- rownames(binData[[1]])
		}
		else if(!is.null(contData)){
			cpdSet <- rownames(contData[[1]])
		}
		else{
			stop("Specify a data set in binData and/or in contData")
		}
		ResultFeat=list()
		Characteristics=list()
		temp=list()
		
		LeadCpds=Selection #names of the objects
		OrderedCpds=cpdSet
		temp[[1]]=list(LeadCpds,OrderedCpds)
		names(temp[[1]])=c("LeadCpds","OrderedCpds")
		
		group <- factor(ifelse(cpdSet %in% LeadCpds, 1, 0)) #identify the group of interest
		
		#Determine characteristic features for the objects: fishers exact test
		result=list()
		
		if(!is.null(binData)){
			for(i in 1:length(binData)){
				binData[[i]]=binData[[i]]+0
				binData[[i]]<-binData[[i]][,which(colSums(binData[[i]]) != 0 & colSums(binData[[i]]) != nrow(binData[[i]]))]
			}
			for(j in 1: length(binData)){
				binMat=binData[[j]]
				
				
				if(length(LeadCpds)==1){
					FP=which(binMat[LeadCpds,]==1)
					Ranks=sort(colSums(binMat[,FP]),)
					
					N=c(names(Ranks),colnames(binMat)[which(!colnames(binMat)%in%names(Ranks))])
					
					AllFeat=data.frame(Names=as.character(N))
					AllFeat$Names=as.character(AllFeat$Names)
					if(is.null(topC)){
						topC=length(Ranks)
					}
					else{
						Ranks=Ranks[c(1:topC)]
					}
					TopFeat=data.frame(Names=as.character(names(Ranks)))
					TopFeat$Names=as.character(TopFeat$Names)
					
					
					temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
					result[[j]]<-temp1
					names(result)[j]=datanames[j]
					
					
				}
				else{
					pFish <- apply(binMat, 2, function(x) stats::fisher.test(table(x, group))$p.value)
					
					pFish <- sort(pFish)
					adjpFish<-stats::p.adjust(pFish, method = "fdr")
					
					AllFeat=data.frame(Names=as.character(names(pFish)),P.Value=pFish,adj.P.Val=adjpFish)
					AllFeat$Names=as.character(AllFeat$Names)
					if(is.null(topC)){
						topC=length(which(pFish<sign))
					}
					
					TopFeat=data.frame(Names=as.character(names(pFish[0:topC])),P.Value=pFish[0:topC],adj.P.Val=adjpFish[0:topC])
					TopFeat$Names=as.character(TopFeat$Names)
					temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
					result[[j]]<-temp1
					names(result)[j]=datanames[j]
				}
			}
			
		}
		
		resultC=list()
		if(!is.null(contData)){
			for(j in 1:length(contData)){
				contMat=contData[[j]]
				
				group1=which(group==1)
				group2=which(group==0)
				
				
				pTTest <- apply(contMat, 2, function(x) stats::t.test(x[group1],x[group2])$p.value)
				
				pTTest <- sort(pTTest)
				adjpTTest<-stats::p.adjust(pTTest, method = "fdr")
				
				AllFeat=data.frame(Names=as.character(names(pTTest)),P.Value=pTTest,adj.P.Val=adjpTTest)
				AllFeat$Names=as.character(AllFeat$Names)
				if(is.null(topC)){
					topC=length(which(pTTest<sign))
				}
				
				TopFeat=data.frame(Names=as.character(names(pTTest[0:topC])),P.Value=pTTest[0:topC],adj.P.Val=adjpTTest[0:topC])
				TopFeat$Names=as.character(TopFeat$Names)
				temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
				resultC[[j]]<-temp1
				names(resultC)[j]=datanames[length(binData)+j]
				
			}
		}
		
		temp[[2]]=c(result,resultC)
		
		
		names(temp)=c("objects","Characteristics")
		
		ResultFeat[[1]]=temp
		names(ResultFeat)="Selection"
		
		
		
	}
	
	else if(class(Selection)=="numeric" & !(is.null(List))){
		
		ListNew=list()
		element=0
		for(i in 1:length(List)){
			if(attributes(List[[i]])$method != "CEC" & attributes(List[[i]])$method != "Weighted" & attributes(List[[i]])$method!= "WeightedSim"){
				ResultsClust=list()
				ResultsClust[[1]]=list()
				ResultsClust[[1]][[1]]=List[[i]]
				names(ResultsClust[[1]])[1]="Clust"
				element=element+1					
				ListNew[[element]]=ResultsClust[[1]]
				#attr(ListNew[element],"method")="Weights"
			}
			else if(attributes(List[[i]])$method=="CEC" | attributes(List[[i]])$method=="Weighted" | attributes(List[[i]])$method == "WeightedSim"){
				ResultsClust=list()
				if(weightclust==TRUE){
					ResultsClust[[1]]=list()
					if(attributes(List[[i]])$method != "WeightedSim"){
						ResultsClust[[1]][[1]]=List[[i]]$Clust
						names(ResultsClust[[1]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[1]]
						attr(ListNew[element],"method")="Weights"
					}
					else{
						ResultsClust[[1]]=list()
						ResultsClust[[1]][[1]]=List[[i]]
						names(ResultsClust[[1]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[1]]
					}
				}
				else{
					for (j in 1:length(List[[i]]$Results)){
						ResultsClust[[j]]=list()
						ResultsClust[[j]][[1]]=List[[i]]$Results[[j]]
						names(ResultsClust[[j]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[j]]
						attr(ListNew[element],"method")="Weights"
					}		
				}		
			}	
		}
		
		if(is.null(names)){
			names=seq(1,length(ListNew),1)
			for(i in 1:length(ListNew)){
				names[i]=paste("Method",i,sep=" ")
			}
		}
		names(ListNew)=names
		
		Matrix=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)
		List=ListNew
		
		
		if(!is.null(binData)){
			cpdSet <- rownames(binData[[1]])
		}
		else{
			cpdSet <- rownames(contData[[1]])
		}
		
		ResultFeat=list()
		for(k in 1:dim(Matrix)[1]){	
			cluster=Selection
			
			hc<-stats::as.hclust(List[[k]]$Clust$Clust)
			OrderedCpds <- hc$labels[hc$order]
			
			Genes=list()
			temp=list()
			
			LeadCpds=colnames(Matrix)[which(Matrix[k,]==cluster)] #names of the objects
			temp[[1]]=list(LeadCpds,OrderedCpds)
			names(temp[[1]])=c("LeadCpds","OrderedCpds")
			
			group <- factor(ifelse(cpdSet %in% LeadCpds, 1, 0)) #identify the group of interest
			
			#Determine characteristic features for the objects: fishers exact test
			result=list()
			if(!is.null(binData)){
				for(i in 1:length(binData)){
					binData[[i]]=binData[[i]]+0
					binData[[i]]<-binData[[i]][,which(colSums(binData[[i]]) != 0 & colSums(binData[[i]]) != nrow(binData[[i]]))]
				}	
				
				for(j in 1: length(binData)){
					binMat=binData[[j]]
					pFish <- apply(binMat, 2, function(x) stats::fisher.test(table(x, group))$p.value)
					
					pFish <- sort(pFish)
					adjpFish<-stats::p.adjust(pFish, method = "fdr")
					
					AllFeat=data.frame(Names=as.character(names(pFish)),P.Value=pFish,adj.P.Val=adjpFish)
					AllFeat$Names=as.character(AllFeat$Names)
					if(is.null(topC)){
						topC=length(which(pFish<0.05))
					}
					
					TopFeat=AllFeat[0:topC,]
					TopFeat$Names=as.character(TopFeat$Names)
					temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
					result[[j]]<-temp1
					names(resultC)[j]=datanames[length(binData)+j]
					
				}
			}
			
			resultC=list()
			if(!is.null(contData)){
				for(j in 1:length(contData)){
					contMat=contData[[j]]
					
					group1=which(group==1)
					group2=which(group==0)
					
					
					pTTest <- apply(contMat, 2, function(x) stats::t.test(x[group1],x[group2])$p.value)
					
					pTTest <- sort(pTTest)
					adjpTTest<-stats::p.adjust(pTTest, method = "fdr")
					
					AllFeat=data.frame(Names=as.character(names(pTTest)),P.Value=pTTest,adj.P.Val=adjpTTest)
					AllFeat$Names=as.character(AllFeat$Names)
					if(is.null(topC)){
						topC=length(which(pTTest<sign))
					}
					
					TopFeat=data.frame(Names=as.character(names(pTTest[0:topC])),P.Value=pTTest[0:topC],adj.P.Val=adjpTTest[0:topC])
					TopFeat$Names=as.character(TopFeat$Names)
					temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
					resultC[[j]]<-temp1
					names(resultC)[j]=datanames[j]
					
				}
			}
			
			temp[[2]]=c(result,resultC)
			
			names(temp)=c("objects","Characteristics")
			ResultFeat[[k]]=temp
			names(ResultFeat)[k]=paste(names[k],": Cluster",cluster, sep=" ")
		}		
	}
	
	else{
		message("If a specific cluster is specified, clustering results must be provided in List")
	}
	return(ResultFeat)
	
}

#' @title List all features present in a selected cluster of objects
#' 
#' @description The function \code{FeaturesOfCluster} lists the number of features objects
#' of the cluster have in common. A threshold can be set selecting among how
#' many objects of the cluster the features should be shared. An optional
#' plot of the features is available.
#' @export FeaturesOfCluster
#' @param leadCpds A character vector containing the objects one wants to
#' investigate in terms of features.
#' @param data The data matrix.
#' @param threshold The number of objects the features at least should be
#' shared amongst. Default is set to 1 implying that the features should be
#' present in at least one of the objects specified in leadCpds.
#' @param plot Logical. Indicates whether or not a BinFeaturesPlot should be
#' set up for the selectcion of objects and discovered features. Default is FALSE.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return A plot indicating the values of the features of the LeadCpds in
#' green and those of the others in blue. It list all features which are
#' present in at least the threshold number of objects. By including all
#' other objects as well, one can see whether features are common in the
#' objects or rather specific for the cluster.
#' 
#' Further, it returns a list with 2 items. The first indicates the number of
#' shared features among the objects. This provides an overview of which
#' objects are more similar than others. The second item is a character
#' vector of the plotted features such that these can be retrieved for further
#' investigation.
#' @examples
#' 
#' \dontrun{
#' data(fingerprintMat)
#' 
#' Lead=rownames(fingerprintMat)[1:5]
#' 
#' FeaturesOfCluster(leadCpds=Lead,data=fingerprintMat,
#' threshold=1,plot=TRUE,plottype="new",location=NULL)
#' }
FeaturesOfCluster<-function(leadCpds,data,threshold=1,plot=TRUE,plottype="new",location=NULL){
	SubsetData=as.matrix(data[which(rownames(data)%in%leadCpds),])
	
	Common=SubsetData%*%t(SubsetData)
	
	if(threshold>length(leadCpds)){
		stop("threshold is larger than the number of LeadCpds. This number can maximally be the number of LeadCpds.")
	}
	
	
	Features=colnames(SubsetData[,which(apply(SubsetData,2,sum)>=threshold)])
	
	if(is.null(Features)){
		Features=""
		Plot=FALSE
	}
	
	#SharedAmongLeadCpds=SubsetData[,Features]
	
	if(plot==TRUE){
		BinFeaturesPlot_SingleData(leadCpds=leadCpds,orderLab=rownames(data),
				features=as.character(Features),data=data,colorLab=NULL,nrclusters=NULL,cols=NULL,name=c("Shared Features Among Selected objects"),
				margins=c(8.5,2.0,0.5,9.5),plottype=plottype,location=location)	
	}
	
	Out=list("Number_of_Common_Features"=Common,"SharedFeatures"=Features)
	
	return(Out)
}


#' @title Find a selection of objects in the output of \code{ReorderToReference}
#' 
#' @description \code{FindCluster} selects the objects belonging to a cluster after the
#' results of the methods have been rearranged by the
#' \code{ReorderToReference}.
#' @export FindCluster
#' 
#' @param List A list of the clustering outputs to be compared. The first
#' element of the list will be used as the reference in
#' \code{ReorderToReference}.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param select The row (the method) and the number of the cluster to select. Default is c(1,1).
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' @param names Optional. Names of the methods. Default is NULL.
#' @return A character vector containing the names of the objects in the
#' selected cluster.
#' @examples
#' 
#' data(fingerprintMat)
#' data(targetMat)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' L=list(MCF7_F,MCF7_T)
#' names=c("FP","TP")
#' 
#' Comps=FindCluster(List=L,nrclusters=7,select=c(1,4))
#' Comps
FindCluster<-function(List,nrclusters=NULL,select=c(1,1),fusionsLog=TRUE, weightclust=TRUE,names=NULL){
	if(length(List)==1 & attributes(List[[1]])$method == "Weighted" & weightclust==TRUE){
		T=List[[1]]$Clust
		attr(T,"method")="Single Clustering"
		List=list(T)
	}
	
	if(length(List)==1){
		Matrix=stats::cutree(List[[1]]$Clust,nrclusters)
		names(Matrix)=rownames(List[[1]]$DistM)
		clusternr=select[2]
		Comps=names(which(Matrix==clusternr))
	}
	
	else{
		Matrix=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)
		methodnr=select[1]
		clusternr=select[2]
		Comps=names(which(Matrix[methodnr,]==clusternr))
	}
	return(Comps)
}	

#' @title Find an element in a data structure
#' 
#' @description The function \code{FindElement} is used internally in the
#' \code{PreparePathway} function but might come in handy for other uses as
#' well. Given the name of an object, the function searches for that object in
#' the data structure and extracts it. When multiple objects have the same
#' name, all are extracted.
#' @export FindElement
#' @param what A character string indicating which object to look for. Default is NULL.
#' @param object The data structure to look into. Only the classes data frame
#' and list are supported. Default is NULL.
#' @param element Not to be specified by the user.
#' @return The returned value is a list with an element for each object found.
#' The element contains everything the object contained in the original data
#' structure.
#' @examples
#' 
#' data(fingerprintMat)
#' data(targetMat)
#' data(geneMat)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' MCF7_DiffGenes_FandT10=DiffGenes(list(MCF7_F,MCF7_T),Selection=NULL,geneExpr=geneMat,
#' nrclusters=7,method="limma",sign=0.05,top=10,fusionsLog = TRUE, weightclust = TRUE, 
#' names = NULL)
#' 
#' Find=FindElement(what='TopDE',object=MCF7_DiffGenes_FandT10)
FindElement<-function(what=NULL,object=NULL,element=list()){
	#str(Object)
	if(class(object)=="data.frame"){
		#search in columns
		if(what %in% colnames(object)){			
			element[[length(element)+1]]<-object[,what]
			names(element)[length(element)]=paste(what,"_",length(element),sep="")
		}
		else if(what %in% rownames(object)){
			element[[length(element)+1]]<-object[what,]
			names(element)[length(element)]=paste(what,"_",length(element),sep="")
		}	
	}
	if(class(object)=="list"){
		#Element=list()
		
		for(i in 0:length(object)){
			if(i==0){
				Names=names(object)
				if(what%in%Names){
					for(j in which(what==Names)){
						element[length(element)+1]=object[j]
						names(element)[length(element)]=paste(what,"_",length(element),sep="")
						return(element)
					}
				}
			}
			else if(class(object[[i]])[1]=="list"){
				#Names=names(Object[[i]])
				#if(What%in%Names){
				#	for(j in which(What==Names)){
				#			Element[length(Element)+1]=Object[[i]][j]
				#		names(Element)[length(Element)]=paste(What,"_",length(Element),sep="")
				#
				#	}
				#}
				element=FindElement(what,object[[i]],element=element)
				#for(j in 1:length(temp)){
				#	Element[length(Element)+1]=temp[j]
				#	names(Element)[length(Element)]=paste(What,"_",length(Element),sep="")
				#}
				
			}
			else if(class(object[[i]])[1]=="data.frame"){
				element=FindElement(what,object[[i]],element=element)
				
			}
		}
	}	
	return(element)
}

#' @title Investigates whether genes are differential expressed in multiple clusters
#' 
#' @description Due to the shifting of objects over the clusters for the different
#' methods, it is possible that the same gene is found significant for a
#' different cluster in another method. These can be tracked with the
#' \code{FindGenes} function. Per method and per cluster, it will take note of
#' the genes found significant and investigate if these were also find for
#' another cluster in another method.
#' @export FindGenes
#' 
#' @param dataLimma Preferably an output of the \code{DiffGenes} function. If
#' not, an ID element of the top genes must be present for each cluster of each
#' method specified in the data structure.
#' @param names Optional. Names of the methods. Default is NULL.
#' @return The returned value is a list with an element per cluster and per
#' cluster one for every gene.  Per gene, a vector is given which contain the
#' methods for which the gene was found. If the cluster is changed compared to
#' the reference method of DataLimma, this is indicated with an underscore.
#' @author Marijke Van Moerbeke
#' @examples
#' 
#' data(fingerprintMat)
#' data(targetMat)
#' data(geneMat)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' MCF7_DiffGenes_FandT10=DiffGenes(list(MCF7_F,MCF7_T),Selection=NULL,geneExpr=geneMat,
#' nrclusters=7,method="limma",sign=0.05,top=10,fusionsLog = TRUE, weightclust = TRUE, 
#' names = NULL)
#' 
#' MCF7_SharedGenes=FindGenes(dataLimma=MCF7_DiffGenes_FandT10,names=c("FP","TP"))
FindGenes<-function(dataLimma,names=NULL){
	
	FoundGenes=list()
	
	if(is.null(names)){
		for(j in 1:length(dataLimma)){
			names[j]=paste("Method",j,sep=" ")	
		}
	}
	
	nrclusters=length(dataLimma[[1]])
	
	
	for(j in 1:nrclusters){
		FoundGenes[[j]]=list()
	}
	
	for(i in 1:length(dataLimma)){ #i == method
		
		for(j in 1:nrclusters){ #j == cluster
			if(!(is.na(dataLimma[[i]][[j]]))[1]){
				
				tempgenes=dataLimma[[i]][[j]]$Genes$TopDE$ID		
				
				
				if(!(is.null(tempgenes))){
					tempgenes=tempgenes[!(is.na(tempgenes))]
					for(k in 1:length(tempgenes)){
						if(!(tempgenes[k] %in% names(FoundGenes[[j]]))){	
							#FoundGenes[[j]][[length(FoundGenes[[j]])+1]]=c()
							FoundGenes[[j]][[length(FoundGenes[[j]])+1]]=names[i]
							names(FoundGenes[[j]])[length(FoundGenes[[j]])]=tempgenes[k]
							
						}
						else if (tempgenes[k] %in% names(FoundGenes[[j]])){
							found=which(names(FoundGenes[[j]])==tempgenes[k])
							FoundGenes[[j]][[found]]=c(FoundGenes[[j]][[found]],names[i])
						}	
					}					
				}	
			}
		}
	}
	
	for(l in 1:nrclusters){
		namesl=names(FoundGenes[[l]])
		if(!(is.null(namesl))){
			
			for(m in l:nrclusters){
				namesm=names(FoundGenes[[m]])	
				Templist=FoundGenes[[m]]
				if(length(namesm)!=0){
					if(l != m){
						
						for(k in 1:length(namesm)){
							
							if (namesm[k] %in% namesl){
								found=which(namesl==namesm[k])
								methods=Templist[[k]]
								del=which(names(FoundGenes[[m]])==namesm[k])
								FoundGenes[[m]][[del]]=c()
								for(a in 1:length(methods)){
									methods[a]=paste(methods[a],"_",m,sep="")
								}
								FoundGenes[[l]][[found]]=c(FoundGenes[[l]][[found]],methods)
							}
							
						}
					}					
				}
			}
		}			
	}
	
	
	for(i in 1:length(FoundGenes)){
		names(FoundGenes)[i]=paste("Cluster",i,sep=" ")
	}
	
	return(FoundGenes)	
}

#' @title Intersection over resulting gene sets of \code{PathwaysIter} function
#' 
#' @description The function \code{Geneset.intersect} collects the results of the
#' \code{PathwaysIter} function per method for each cluster and takes the
#' intersection over the iterations per cluster per method. This is to see if
#' over the different resamplings of the data, similar pathways were
#' discovered.
#' @export Geneset.intersect
#' 
#' @param PathwaysOutput The output of the \code{PathwaysIter} function.
#' @param Selection Logical. Indicates whether or not the output of the
#' pathways function were concentrated on a specific selection of objects. If
#' this was the case, Selection should be put to TRUE. Otherwise, it should be
#' put to FALSE. Default is TRUE.
#' @param sign The significance level to be handled for cutting of the
#' pathways. Default is 0.05.
#' @param names Optional. Names of the methods. Default is NULL.
#' @param seperatetables Logical. If TRUE, a separate element is created per
#' cluster containing the pathways for each iteration. Default is FALSE.
#' @param separatepvals Logical. If TRUE, the p-values of the each iteration of
#' each pathway in the intersection is given. If FALSE, only the mean p-value
#' is provided. Default is FALSE.
#' @return The output is a list with an element per method. For each method, it
#' is portrayed per cluster which pathways belong to the intersection over all
#' iterations and their corresponding mean p-values.
#' @examples
#' 
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(geneMat)
#' data(GeneInfo)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' L=list(MCF7_F,MCF7_T)
#' 
#' MCF7_Paths_FandT=PathwaysIter(List=L, geneExpr = geneMat, nrclusters = 7, method = 
#' c("limma", "MLP"), geneInfo = GeneInfo, geneSetSource = "GOBP", topP = NULL, 
#' topG = NULL, GENESET = NULL, sign = 0.05,niter=2,fusionsLog = TRUE, 
#' weightclust = TRUE, names =names)
#' 
#' MCF7_Paths_intersection=Geneset.intersect(PathwaysOutput=MCF7_Paths_FandT,
#' sign=0.05,names=c("FP","TP"),seperatetables=FALSE,separatepvals=FALSE)
#' 
#' str(MCF7_Paths_intersection)
#' }
Geneset.intersect<-function(PathwaysOutput,Selection=FALSE,sign=0.05,names=NULL,seperatetables=FALSE,separatepvals=FALSE){
	
	if(Selection==TRUE){
		if(length(PathwaysOutput$'Iteration 1')==1){
			names="Selection"
		}
		Intersect=Geneset.intersectSelection(PathwaysOutput,sign,names,seperatetables,separatepvals)	
	}
	
	else{
		
		if(is.null(names)){
			for(j in 1:length(PathwaysOutput$"Iteration 1")){
				names[j]=paste("Method",j,sep=" ")	
			}
		}
		
		
		#put all of same method together:preparation of lists
		subsets=list()
		nmethods=length(PathwaysOutput$"Iteration 1") 
		for(i in 1:nmethods){
			subsets[[i]]=list()
			
		}
		names(subsets)=names
		
		#put all of same method together: go through PathwaysOutput
		for(j in 1:length(PathwaysOutput)){
			name1=names(PathwaysOutput)[j]
			for(k in 1:nmethods){
				name2=names[k]
				subsets[[name2]][[name1]]=PathwaysOutput[[name1]][[name2]]
				
			}
			
		}	
		
		#for every subset (= every method) take intersection over the interations per cluster
		Intersect=list()
		
		for(i in 1:length(subsets)){
			Method=subsets[[i]]
			Clusters=list()
			nclus=length(Method[[1]])
			for(j in 1:length(Method)){
				name3=paste("Iteration",j,sep=" ")
				for(k in 1:nclus){
					name4=paste("Cluster",k,sep=" ")
					if(!(is.na(Method[[name3]][[name4]])[1])){
						Clusters[[name4]][[name3]]=Method[[name3]][[name4]]
					}
					else{
						Clusters[[name4]][[name3]]=NA
					}
				}
				
			}
			
			IntersectM=list()
			
			for(a in 1:length(Clusters)){ #per cluster
				if(!(is.na(Clusters[[a]])[1])){
					result.out=list()
					result.name = c()
					for(b in 1:length(Clusters[[a]])){#per iteration
						if(b==1){
							objects=Clusters[[a]][[1]]$objects
							Genes=Clusters[[a]][[1]]$Genes	
							Names=data.frame("description"=Clusters[[a]][[1]]$Pathways$AllPaths$geneSetDescription,"genesetcode"=rownames(Clusters[[a]][[1]]$Pathways$AllPaths))
							Names$description=as.character(Names$description)	
							Names$genesetcode=as.character(Names$genesetcode)	
						}				
						cut = Clusters[[a]][[b]]$Pathways$AllPaths[Clusters[[a]][[b]]$Pathways$AllPaths$geneSetPValue<=sign,]
						colnames(cut)[4] = paste("pvalues.",b,sep="")
						colnames(cut)[2] = paste("testedgenesetsize.",b,sep="")
						colnames(cut)[3] = paste("genesetstatistic.",b,sep="")
						cut=cut[,c(1,5,2,3,4)]
						result.out[[b]] = cut
						result.name = c(result.name,paste("genesettable",b,sep=""))
						
					}
					
					
					
					names(result.out) = result.name
					
					genesets.table.intersect = plyr::join_all(result.out,by=c("totalGeneSetSize","geneSetDescription"),type="inner")
					genesets.table.intersect$mean_testedGeneSetSize=round(apply(genesets.table.intersect[,which(substring(colnames(genesets.table.intersect),1,nchar(colnames(genesets.table.intersect))-nchar(".1"))=='testedgenesetsize')],1,mean),1)
					genesets.table.intersect$mean_geneSetStatistic=apply(genesets.table.intersect[,which(substring(colnames(genesets.table.intersect),1,nchar(colnames(genesets.table.intersect))-nchar(".1"))=='genesetstatistic')],1,mean)
					genesets.table.intersect$mean_geneSetPValue=apply(genesets.table.intersect[,which(substring(colnames(genesets.table.intersect),1,nchar(colnames(genesets.table.intersect))-nchar(".1"))=='pvalues')],1,mean)
					
					rownames(genesets.table.intersect)=as.character(Names[which(genesets.table.intersect$geneSetDescription%in%Names[,1]),2])
					
					class(genesets.table.intersect)=c("MLP","data.frame")
					attr(genesets.table.intersect,'geneSetSource')=attributes(Clusters[[1]][[1]]$Pathways$AllPaths)$geneSetSource
					
					
					result.out$genesets.table.intersect = genesets.table.intersect
					
					
					
					if(separatepvals==FALSE){
						result.out$genesets.table.intersect=genesets.table.intersect[,c(1,2,(ncol(genesets.table.intersect)-2):ncol(genesets.table.intersect))]
						class(result.out$genesets.table.intersect)=c("MLP","data.frame")
						attr(result.out$genesets.table.intersect,'geneSetSource')=attributes(Clusters[[1]][[1]]$Pathways$AllPaths)$geneSetSource
					}
					
					
					if(seperatetables==FALSE){
						result.out=result.out$genesets.table.intersect
						class(result.out)=c("MLP","data.frame")
						attr(result.out,'geneSetSource')=attributes(Clusters[[1]][[1]]$Pathways$AllPaths)$geneSetSource
					}
					
					newresult=list(objects=objects,Genes=Genes,Pathways=result.out)
					
					
					IntersectM[[a]]=newresult	
					names(IntersectM)[a]=names(Clusters)[[a]]
				}
				else{
					IntersectM[[a]]=NA	
					names(IntersectM)[a]=names(Clusters)[[a]]
				}
			}
			
			Intersect[[i]]=IntersectM
			
		}
	}
	names(Intersect)=names
	return(Intersect)
}

#' @title Intersection over resulting gene sets of \code{PathwaysIter} function for a selection of objects
#' @param list.output The output of the \code{PathwaysIter} function.
#' @param sign The significance level to be handled for cutting of the
#' pathways. Default is 0.05.
#' @param names Optional. Names of the methods. Default is NULL.
#' @param seperatetables Logical. If TRUE, a separate element is created per
#' cluster containing the pathways for each iteration. Default is FALSE.
#' @param separatepvals Logical. If TRUE, the p-values of the each iteration of
#' each pathway in the intersection is given. If FALSE, only the mean p-value
#' is provided. Default is FALSE.
#' @description Internal function of \code{Geneset.intersect}.
Geneset.intersectSelection<-function(list.output,sign=0.05,names=NULL,seperatetables=FALSE,separatepvals=FALSE){
	if(is.null(names)){
		for(j in 1:length(list.output$"Iteration 1")){
			names[j]=paste("Method",j,sep=" ")	
		}
	}
	
	#put all of same method together:preparation of lists
	subsets=list()
	nmethods=length(list.output$"Iteration 1") 
	for(i in 1:nmethods){
		subsets[[i]]=list()		
	}
	names(subsets)=names
	
	#put all of same method together: go through list.output
	for(j in 1:length(list.output)){
		name1=names(list.output)[j]
		for(k in 1:nmethods){
			name2=k
			subsets[[name2]][[name1]]=list.output[[name1]][[name2]]
			
		}
		
	}	
	
	#for every subset (= every method) take intersection over the interations per cluster
	Intersect=list()
	
	for(i in 1:length(subsets)){
		Method=subsets[[i]]
		Clusters=list()
		nclus=1
		for(j in 1:length(Method)){
			name3=paste("Iteration",j,sep=" ")
			Clusters[[name3]]=Method[[name3]]			
		}
		
		IntersectM=list()
		
		result.out=list()
		result.name = c()
		for(a in 1:length(Clusters)){ #per cluster
			if(a==1){
				objects=Clusters[[a]]$objects
				Genes=Clusters[[a]]$Genes
				Names=data.frame("description"=Clusters[[a]]$Pathways$AllPaths$geneSetDescription,"genesetcode"=rownames(Clusters[[a]]$Pathways$AllPaths))
				Names$description=as.character(Names$description)	
				Names$genesetcode=as.character(Names$genesetcode)	
			}				
			cut = Clusters[[a]]$Pathways$AllPaths[Clusters[[a]]$Pathways$AllPaths$geneSetPValue<=sign,]
			colnames(cut)[4] = paste("pvalues.",a,sep="")
			colnames(cut)[2] = paste("testedgenesetsize.",a,sep="")
			colnames(cut)[3] = paste("genesetstatistic.",a,sep="")
			cut=cut[,c(1,5,2,3,4)]
			result.out[[a]] = cut
			result.name = c(result.name,paste("genesettable",a,sep=""))
		}
		
		names(result.out) = result.name
		
		genesets.table.intersect = plyr::join_all(result.out,by=c("totalGeneSetSize","geneSetDescription"),type="inner")
		genesets.table.intersect$mean_testedGeneSetSize=round(apply(genesets.table.intersect[,which(substring(colnames(genesets.table.intersect),1,nchar(colnames(genesets.table.intersect))-nchar(".1"))=='testedgenesetsize')],1,mean),1)
		genesets.table.intersect$mean_geneSetStatistic=apply(genesets.table.intersect[,which(substring(colnames(genesets.table.intersect),1,nchar(colnames(genesets.table.intersect))-nchar(".1"))=='genesetstatistic')],1,mean)
		genesets.table.intersect$mean_geneSetPValue=apply(genesets.table.intersect[,which(substring(colnames(genesets.table.intersect),1,nchar(colnames(genesets.table.intersect))-nchar(".1"))=='pvalues')],1,mean)
		
		rownames(genesets.table.intersect)=as.character(Names[which(genesets.table.intersect$geneSetDescription%in%Names[,1]),2])
		
		class(genesets.table.intersect)=c("MLP","data.frame")
		attr(genesets.table.intersect,'geneSetSource')=attributes(Clusters[[1]]$Pathways$AllPaths)$geneSetSource
		
		
		result.out$genesets.table.intersect = genesets.table.intersect
		
		if(separatepvals==FALSE){
			result.out$genesets.table.intersect=genesets.table.intersect[,c(1,2,(ncol(genesets.table.intersect)-2):ncol(genesets.table.intersect))]
			class(result.out$genesets.table.intersect)=c("MLP","data.frame")
			attr(result.out$genesets.table.intersect,'geneSetSource')=attributes(Clusters[[1]]$Pathways$AllPaths)$geneSetSource
		}
		
		
		if(seperatetables==FALSE){
			result.out=result.out$genesets.table.intersect
			class(result.out)=c("MLP","data.frame")
			attr(result.out,'geneSetSource')=attributes(Clusters[[1]]$Pathways$AllPaths)$geneSetSource
		}
		
		
		newresult=list(objects=objects,Genes=Genes,Pathways=result.out)
		
		
		#IntersectM[[a]]=
		#names(IntersectM)[a]=names(Clusters)[[a]]
		Intersect[[i]]=newresult	
		
	}
	names(Intersect)=names
	return(Intersect)
}

#' @title Comparing two clustering results with a heatmap
#' 
#' @description The \code{HeatmapCols} function calculates the distance between two outputs
#' of clustering methods and plots the resulting heatmap. The function
#' heatmap.2 is called upon to make the actual plot of the heatmap. It is noted
#' that for thi s function the number of colors should be one more than the
#' number of clusters to color the so called zero cells in the distance matrix.
#' 
#' Another way to compare to methods is via an adaptation of heatmaps. The
#' input of this function is the resulting clustering (the Clust element of the
#' list) of two methods and can be seen as: method 1 versus method 2. The
#' dendrograms are cut into a specific number of clusters. Each cluster of
#' method 2 and its members are given a distinct color represented by a number.
#' These are the clusters to which a comparison is made. A matrix is set up of
#' which the columns are determined by the ordering of clustering of method 2
#' and the rows by the ordering of method 1. Every column represent one object
#' just as every row and every column represent the color of its cluster. A
#' function visits every cell of the matrix. If the objects represented by the
#' cell are still together in a cluster, the color of the column is passed to
#' the cell. This creates the distance matrix which can be given to the
#' HeatmapCols function to create the heatmap.
#' 
#' @param Data1 The resulting clustering of method 1.
#' @param Data2 The resulting clustering of method 2.
#' @param names The names of the objects in the data sets. Default is NULL.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param cols A character vector with the colours for the clusters. Default is NULL.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return A heatmap based on the distance matrix created by the function with
#' the dendrogram of method 2 on top of the plot and the one from method 1 on
#' the left. The names of the objects are depicted on the bottom in the order
#' of clustering of method 2 and on the right by the ordering of method 1.
#' Vertically the cluster of method 2 can be seen while horizontally those of
#' method 1 are portrayed.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(Colors1)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' clust="agnes",linkage="flexible",gap=FALSE,maxK=15)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' clust="agnes",linkage="flexible",gap=FALSE,maxK=15)
#' 
#' 
#' L=list(MCF7_F,MCF7_T)
#' names=c("FP","TP")
#' 
#' HeatmapPlot(Data1=MCF7_T,Data2=MCF7_F,names=rownames(fingerprintMat)
#' ,nrclusters=7,cols=Colors1,plottype="new", location=NULL)
#' }
#' 
#' @export HeatmapPlot
HeatmapPlot<-function(Data1,Data2,names=NULL,nrclusters=NULL,cols=NULL,plottype="new",location=NULL){
	data1=Data1$Clust
	data2=Data2$Clust
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	DistM=distanceheatmaps(data1,data2,names,nrclusters)
	plottypein(plottype,location)
	gplots::heatmap.2(DistM,Rowv =stats::as.dendrogram(data1), Colv=stats::as.dendrogram(data2),trace="none",col=cols,key=FALSE)
	plottypeout(plottype)
}

#' @title A function to select a group of objects via the similarity heatmap.
#' 
#' @description The function \code{HeatmapSelection} plots the similarity values between
#' objects. The plot is similar to the one produced by
#' \code{SimilarityHeatmap} but without the dendrograms on the sides. The
#' function is rather explorative and experimental and is to be used with some
#' caution. By clicking in the plot, the user can select a group of objects
#' of interest. See more in \code{Details}.
#' 
#' A similarity heatmap is created in the same way as in
#' \code{SimilarityHeatmap}. The user is now free to select two points on the
#' heatmap. It is advised that these two points are in opposite corners of a
#' square that indicates a high similarity among the objects. The points do
#' not have to be the exact corners of the group of interest, a little
#' deviation is allowed as rows and columns of the selected subset of the
#' matrix with sum equal to 1 are filtered out. A sum equal to one, implies
#' that the compound is only similar to itself.
#' 
#' The function is meant to be explorative but is experimental. The goal was to
#' make the selection of interesting objects easier as sometimes the labels
#' of the dendrograms are too distorted to be read. If the figure is exported
#' to a pdf file with an appropriate width and height, the labels can be become
#' readable again.
#' 
#' @param Data The data of which a heatmap should be drawn.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clusters".
#' @param distmeasure The distance measure. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to "tanimoto".
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}. 
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is NULL.
#' @param linkage Choice of inter group dissimilarity (character). Defaults to "flexible".
#' @param cutoff Optional. If a cutoff value is specified, all values lower are
#' put to zero while all other values are kept. This helps to highlight the
#' most similar objects. Default is NULL.
#' @param percentile Logical. The cutoff value can be a percentile. If one want
#' the cutoff value to be the 90th percentile of the data, one should specify
#' cutoff = 0.90 and percentile = TRUE. Default is FALSE.
#' @param dendrogram Optional. If the clustering results of the data is already
#' available and should not be recalculated, this results can be provided here.
#' Otherwise, it will be calculated given the data. This is necessary to have
#' the objects in their order of clustering on the plot. Default is NULL.
#' @param width The width of the plot to be made. This can be adjusted since
#' the default size might not show a clear picture. Default is 7.
#' @param height The height of the plot to be made. This can be adjusted since
#' the default size might not show a clear picture. Default is 7.
#' @return A heatmap with the names of the objects on the right and bottom.
#' Once points are selected, it will return the names of the objects that are
#' in the selected square provided that these show similarity among each other.
#' @examples
#' 
#' \dontrun{
#' data(fingerprintMat)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55)
#' 
#' HeatmapSelection(Data=MCF7_F$DistM,type="dist",cutoff=0.90,percentile=TRUE,
#' dendrogram=MCF7_F,width=7,height=7)
#' }
#' 
#' @export HeatmapSelection
HeatmapSelection<-function(Data,type=c("data","dist","clust","sim"),distmeasure="tanimoto",normalize=FALSE,method=NULL,linkage="flexible",cutoff=NULL,percentile=FALSE,dendrogram=NULL,width=7,height=7){
	
	#create binary similarity heatmap first
	if(type=="data"){
		ClustData<-Cluster(Data=Data,distmeasure=distmeasure,normalize=normalize,method=method,clust="agnes",linkage=linkage,gap=FALSE,maxK=55,StopRange=FALSE)
		Data=ClustData$DistM
		type="dist"
	}
	
	
	if(type=="clust"){
		Dist=Data$DistM
		if(0<=min(Dist) & max(Dist)<=1){
			SimData=1-Dist
		}
		else{
			NormData=Normalization(Dist,method="Range")
			SimData=1-NormData
		}
		if(is.null(dendrogram)){
			dendrogram=Data
		}
	}
	
	else if(type=="dist"){
		if(0<=min(Data) & max(Data)<=1){
			SimData=1-Data
			if(is.null(dendrogram)){
				dendrogram=Cluster(Data=Data,type="dist",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
			}
		}
		else{
			NormData=Normalization(Data,method="Range")
			SimData=1-NormData
			if(is.null(dendrogram)){
				dendrogram=Cluster(Data=Data,type="dist",distmeasure="tanimoto",normalize=TRUE,method="Q",clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
			}
		}
		
		
	}
	else if(type=="sim"){
		SimData=Data
		if(0<=min(SimData) & max(SimData)<=1){
			if(is.null(dendrogram)){
				DistData=1-Data
				ClustData=Cluster(Data=DistData,type="dist",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
			}
		}
		else{
			if(is.null(dendrogram)){
				NormData=Normalization(Dist,method="Range")
				DistData=1-Data
				ClustData=Cluster(Data=DistData,type="dist",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)		
			}
		}
	}	
	
	
	if(!is.null(cutoff)){
		if(percentile==TRUE){
			cutoff=stats::quantile(SimData[lower.tri(SimData)], cutoff)
		}
		
		SimData_bin <- ifelse(SimData<=cutoff,0,SimData) # Every value higher than the 90ieth percentile is kept, all other are put to zero
	}
	
	else{
		SimData_bin=SimData
	}
	
	
	
	dend <- stats::as.dendrogram(dendrogram$Clust)
	Ind <- stats::order.dendrogram(dend)
	
	SimData_bin=SimData_bin[Ind,Ind]
	
	#Layout<-rbind(4:3, 2:1)
	#lhei <- c(0.4, 4)	
	#lwid <- c(0.4, 4)
	#layout(Layout, widths = lwid, heights = lhei, respect = FALSE)
	grDevices::dev.new(width=width,height=height)
	graphics::par(mar = c(9,7, 7, 9))
	graphics::image(x=1:nrow(SimData_bin),y=1:ncol(SimData_bin),z=t(SimData_bin),col=(grDevices::gray(seq(0.9,0,len=1000))),axes=FALSE,xlab="",ylab="")
	graphics::axis(1, 1:ncol(SimData_bin), labels = colnames(SimData_bin), las = 2, line =0, tick = 0, cex.axis = 0.6)
	graphics::axis(4, 1:nrow(SimData_bin), labels = rownames(SimData_bin), las = 2, line = 0, tick = 0, cex.axis = 0.6)
	
	points=graphics::locator(n=2,type="l")	
	cols=c(floor(points$x[1]),ceiling(points$x[2]))
	rows=c(floor(points$y[1]),ceiling(points$y[2]))
	
	if(cols[1]>cols[2]){
		colseq=seq(cols[2],cols[1],1)
	}
	else{
		colseq=seq(cols[1],cols[2],1)
	}
	
	if(rows[1]>rows[2]){
		rowseq=seq(rows[2],rows[1],1)
	}
	else{
		rowseq=seq(rows[1],rows[2],1)
	}
	
	print(rowseq)
	print(colseq)
	SubsetData=SimData_bin[rowseq,colseq]
#	DelRows=rownames(SubsetData)[which(rowSums(SubsetData)==1)]
#	DelCols=colnames(SubsetData)[which(colSums(SubsetData)==1)]
#	
#	if(length(DelRows)!=0 & length(DelCols)!=0){
#		Subset=SubsetData[-which(rownames(SubsetData)%in%c(DelRows,DelCols)),-which(colnames(SubsetData)%in%c(DelRows,DelCols))]
#	}
#	else if(length(DelRows)!=0 & length(DelCols)==0){
#		Subset=SubsetData[-which(rownames(SubsetData)%in%c(DelRows)),]
#		
#	}
#	else if(length(DelRows)==0 & length(DelCols)!=0){
#		Subset=SubsetData[,-which(colnames(SubsetData)%in%c(DelCols))]
#		
#	}	
#	else if(length(DelRows)==0 & length(DelCols)==0){
#		Subset=SubsetData
#		
#	}
#	SelComps=colnames(Subset)
	SelComps=colnames(SubsetData)
	
	return(SelComps)
}

#' @title Colouring labels
#' @param x The leaf of a dendrogram.
#' @param Sel1 The selection of objects to be colored. Default is NULL.
#' @param Sel2 An optional second selection to be colored. Default is NULL.
#' @param col1 The color for the first selection. Default is NULL.
#' @param col2 The color for the optional second selection. Default is NULL.
#' @description Internal function of \code{LabelPlot}.
LabelCols <- function(x,Sel1,Sel2=NULL,col1=NULL,col2=NULL) {
	colfunc=function(x,Sel1,Sel2,col1,col2){
		if (x %in% Sel1){
			return(col1)
		}
		else if(x %in% Sel2){
			return(col2)
		}
		else{
			return("black")
		}
	}
	
	if (stats::is.leaf(x)) {
		## fetch label
		label <- attr(x, "label") 
		## set label color to red for SelF, to black otherwise
		attr(x, "nodePar") <- list(pch=NA,lab.col=colfunc(label,Sel1,Sel2,col1,col2),lab.cex=0.9,font=2)
		attr(x, "edgePar") <- list(lwd=2,col=colfunc(label,Sel1,Sel2,col1,col2))
	}
	return(x)
}

#' @title Coloring specific leaves of a dendrogram
#' 
#' @description The function plots a dendrogrmam of which specific leaves are coloured.
#' 
#' @param Data The result of a method which contains the dendrogram to be
#' colored.
#' @param sel1 The selection of objects to be colored. Default is NULL.
#' @param sel2 An optional second selection to be colored. Default is NULL.
#' @param col1 The color for the first selection. Default is NULL.
#' @param col2 The color for the optional second selection. Default is NULL.
#' @return A plot of the dendrogram of which the leaves of the selection(s) are
#' colored.
#' @examples
#' 
#' data(fingerprintMat)
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' ClustF_6=cutree(MCF7_F$Clust,6)
#' 
#' SelF=rownames(fingerprintMat)[ClustF_6==6]
#' SelF
#' 
#' LabelPlot(Data=MCF7_F,sel1=SelF,sel2=NULL,col1='darkorchid')
#' 
#' 
#' @export LabelPlot
LabelPlot<-function(Data,sel1,sel2=NULL,col1=NULL,col2=NULL){
	x=Data$Clust
	
	d_temp <- stats::dendrapply(stats::as.dendrogram(x,hang=0.02),LabelCols,sel1,sel2,col1,col2)
	
	graphics::plot(d_temp,nodePar=list(pch=NA),edgePar=list(lwd=2),ylab="Height",font.axis=2,font.lab=2,font=2)
	graphics::axis(side = 2, lwd = 2)
}

#' @title Pathway Analysis
#' 
#' @description The \code{PathwayAnalysis} function combines the functions
#' \code{PathwaysIter} and \code{Geneset.intersect} such that only one function
#' should be called.
#' 
#' 
#' @param List A list of clustering outputs or output of the\code{DiffGenes}
#' function. The first element of the list will be used as the reference in
#' \code{ReorderToReference}. The output of \code{ChooseFeatures} is also
#' accepted.
#' @param Selection If pathway analysis should be conducted for a specific
#' selection of objects, this selection can be provided here. Selection can
#' be of the type "character" (names of the objects) or "numeric" (the number
#' of specific cluster). Default is NULL.
#' @param geneExpr The gene expression matrix of the objects. The rows should
#' correspond with the genes.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param method The method to applied to look for differentially expressed genes and related pathways. For now, only the
#' limma method is available for gene analysis and the MLP method for pathway analysis. Default is c("limma","MLP").
#' @param geneInfo A data frame with at least the columns ENTREZID and SYMBOL.
#' This is necessary to connect the symbolic names of the genes with their
#' EntrezID in the correct order. The order of the gene is here not in the
#' order of the rownames of the gene expression matrix but in the order of
#' their significance. Default is NULL.
#' @param geneSetSource The source for the getGeneSets function, defaults to
#' "GOBP".
#' @param topP Overrules sign. The number of pathways to display for each
#' cluster. If not specified, only the significant genes are shown. Default is NULL.
#' @param topG Overrules sign. The number of top genes to be returned in the
#' result. If not specified, only the significant genes are shown. Default is NULL.
#' @param GENESET Optional. Can provide own candidate gene sets. Default is NULL.
#' @param sign The significance level to be handled. Default is 0.05.
#' @param niter The number of times to perform pathway analysis. Default is 10.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' @param names Optional. Names of the methods. Default is NULL.
#' @param seperatetables Logical. If TRUE, a separate element is created per
#' cluster. containing the pathways for each iteration. Default is FALSE.
#' @param separatepvals Logical. If TRUE, the p-values of the each iteration of
#' each pathway in the intersection is given. If FALSE, only the mean p-value
#' is provided. Default is FALSE.
#' @return The output is a list with an element per method. For each method, it
#' is portrayed per cluster which pathways belong to the intersection over all
#' iterations and their corresponding mean p-values.
#' @examples
#' 
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(geneMat)
#' data(GeneInfo)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' L=list(MCF7_F,MCF7_T)
#' names=c('FP','TP')
#' 
#' MCF7_PathsFandT=PathwayAnalysis(List=L, geneExpr = geneMat, nrclusters = 7, method = c("limma", 
#' "MLP"), geneInfo = GeneInfo, geneSetSource = "GOBP", topP = NULL, 
#' topG = NULL, GENESET = NULL, sign = 0.05,niter=2,fusionsLog = TRUE, weightclust = TRUE, 
#'  names =names,seperatetables=FALSE,separatepvals=FALSE)
#' }
#' 
#' @export PathwayAnalysis
PathwayAnalysis<-function(List,Selection=NULL,geneExpr=NULL,nrclusters=NULL,method=c("limma", "MLP"),geneInfo=NULL,geneSetSource = "GOBP",topP=NULL,topG=NULL,GENESET=NULL,sign=0.05,niter=10,fusionsLog=TRUE,weightclust=TRUE,names=NULL,seperatetables=FALSE,separatepvals=FALSE){
	if (!requireNamespace("MLP", quietly = TRUE)) {
		stop("MLP needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	if (!requireNamespace("biomaRt", quietly = TRUE)) {
		stop("biomaRt needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
		stop("org.Hs.eg.db needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	Pathways=PathwaysIter(List,Selection,geneExpr,nrclusters,method,geneInfo,geneSetSource,topP,topG,GENESET,sign,niter,fusionsLog,weightclust,names)
	
	if(is.null(Selection)){
		Selection=FALSE
	}
	else{
		Selection=TRUE
	}
	
	Intersection=Geneset.intersect(PathwaysOutput=Pathways,Selection,sign,names,seperatetables,separatepvals)
	
	return(Intersection)
	
}

#' @title Pathway analysis for multiple clustering results
#' 
#' @description A pathway analysis per cluster per method is conducted.
#' @export Pathways
#' @details After finding differently expressed genes, it can be investigated whether
#' pathways are related to those genes. This can be done with the help of the
#' function \code{Pathways} which makes use of the \code{MLP} function of the
#' MLP package. Given the output of a method, the cutree function is performed
#' which results into a specific number of clusters. For each cluster, the
#' limma method is performed comparing this cluster to the other clusters. This
#' to obtain the necessary p-values of the genes. These are used as the input
#' for the \code{MLP} function to find interesting pathways. By default the
#' candidate gene sets are determined by the \code{AnnotateEntrezIDtoGO}
#' function. The default source will be GOBP, but this can be altered.
#' Further, it is also possible to provide own candidate gene sets in the form
#' of a list of pathway categories in which each component contains a vector of
#' Entrez Gene identifiers related to that particular pathway. The default
#' values for the minimum and maximum number of genes in a gene set for it to
#' be considered were used. For MLP this is respectively 5 and 100. If a list
#' of outputs of several methods is provided as data input, the cluster numbers
#' are rearranged according to a reference method. The first method is taken as
#' the reference and ReorderToReference is applied to get the correct ordering.
#' When the clusters haven been re-appointed, the pathway analysis as described
#' above is performed for each cluster of each method.
#' @param List A list of clustering outputs or output of the\code{DiffGenes}
#' function. The first element of the list will be used as the reference in
#' \code{ReorderToReference}. The output of \code{ChooseFeatures} is also
#' accepted.
#' @param Selection If pathway analysis should be conducted for a specific
#' selection of objects, this selection can be provided here. Selection can
#' be of the type "character" (names of the objects) or "numeric" (the number
#' of specific cluster). Default is NULL.
#' @param geneExpr The gene expression matrix or ExpressionSet of the objects.
#' The rows should correspond with the genes.
#' @param nrclusters Optional. The number of clusters to cut the dendrogram in.
#' The number of clusters should not be specified if the interest lies only in
#' a specific selection of objects which is known by name.  Otherwise, it is
#' required. Default is NULL.
#' @param method The method to applied to look for differentially expressed genes and related pathways. For now, only the
#' limma method is available for gene analysis and the MLP method for pathway analysis. Default is c("limma","MLP").
#' @param geneInfo A data frame with at least the columns ENTREZID and SYMBOL.
#' This is necessary to connect the symbolic names of the genes with their
#' EntrezID in the correct order. The order of the gene is here not in the
#' order of the rownames of the gene expression matrix but in the order of
#' their significance. Default is NULL.
#' @param geneSetSource The source for the getGeneSets function, defaults to
#' "GOBP".
#' @param topP Overrules sign. The number of pathways to display for each
#' cluster. If not specified, only the significant genes are shown. Default is NULL.
#' @param topG Overrules sign. The number of top genes to be returned in the
#' result. If not specified, only the significant genes are shown. Defaults is NULL.
#' @param GENESET Optional. Can provide own candidate gene sets. Default is NULL.
#' @param sign The significance level to be handled. Default is 0.05.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' @param names Optional. Names of the methods. Default is NULL.
#' @return The returned value is a list with an element per cluster per method.
#' This element is again a list with the following four elements:
#' \item{objects}{A list with the elements LeadCpds (the objects of
#' interest) and OrderedCpds (all objects in the order of the clustering
#' result)} \item{Characteristics}{The found (top) characteristics of the
#' feauture data} \item{Genes}{A list with the elements TopDE (a table with
#' information on the top genes) and AllDE (a table with information on all
#' genes)} \item{Pathways}{A list with the element ranked.genesets.table which
#' is a data frame containing the genesets, their p-values and their
#' descriptions. The second element is nr.genesets and contains the used and
#' total number of genesets.}
#' @examples
#' 
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(geneMat)
#' data(GeneInfo)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' L=list(MCF7_F,MCF7_T)
#' names=c('FP','TP')
#' 
#' MCF7_PathsFandT=Pathways(List=L, geneExpr = geneMat, nrclusters = 7, method = c("limma", 
#' "MLP"), geneInfo = GeneInfo, geneSetSource = "GOBP", topP = NULL, 
#' topG = NULL, GENESET = NULL, sign = 0.05,fusionsLog = TRUE, weightclust = TRUE, 
#'  names =names)
#'  }
Pathways<-function(List,Selection=NULL,geneExpr=NULL,nrclusters=NULL,method=c("limma", "MLP"),geneInfo=NULL,geneSetSource = "GOBP",topP=NULL,topG=NULL,GENESET=NULL,sign=0.05,fusionsLog=TRUE,weightclust=TRUE,names=NULL){
	
	if (!requireNamespace("MLP", quietly = TRUE)) {
		stop("MLP needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	if (!requireNamespace("biomaRt", quietly = TRUE)) {
		stop("biomaRt needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
		stop("org.Hs.eg.db needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	
	
	if(!(is.null(Selection))){
		ResultMLP=PathwaysSelection(List,Selection,geneExpr,nrclusters,method,geneInfo,geneSetSource,topP,topG,GENESET,sign,fusionsLog,weightclust,names)
		
	}
	else if(class(List)=="ChosenClusters"){
		ResultMLP=list()
		for(i in 1:length(List)){	
			Selection=List[[i]]$objects$LeadCpds
			L=List[i]
			ResultMLP[[i]]=PathwaysSelection(List=L,Selection,geneExpr,nrclusters,method,geneInfo,geneSetSource,topP,topG,GENESET,sign,fusionsLog,weightclust,names)	
			names(ResultMLP)=paste("Choice",i,sep=' ')
		}	
	}
	else{
		
		#Check for gene expression data: if not, reordering to ListNew not necessary
		DataPrepared<-plyr::try_default(PreparePathway(List[[1]],geneExpr,topG,sign),NULL,quiet=TRUE)
		if(is.null(DataPrepared)){
			
			ListNew=list()
			element=0
			for(i in 1:length(List)){
				if(attributes(List[[i]])$method != "CEC" & attributes(List[[i]])$method != "Weighted" & attributes(List[[i]])$method!= "WeightedSim"){
					ResultsClust=list()
					ResultsClust[[1]]=list()
					ResultsClust[[1]][[1]]=List[[i]]
					names(ResultsClust[[1]])[1]="Clust"
					element=element+1					
					ListNew[[element]]=ResultsClust[[1]]
					#attr(ListNew[element],"method")="Weights"
				}
				else if(attributes(List[[i]])$method=="CEC" | attributes(List[[i]])$method=="Weighted" | attributes(List[[i]])$method == "WeightedSim"){
					ResultsClust=list()
					if(weightclust==TRUE){
						ResultsClust[[1]]=list()
						if(attributes(List[[i]])$method != "WeightedSim"){
							ResultsClust[[1]][[1]]=List[[i]]$Clust
							names(ResultsClust[[1]])[1]="Clust"
							element=element+1					
							ListNew[[element]]=ResultsClust[[1]]
							attr(ListNew[element],"method")="Weights"
						}
						else{
							ResultsClust[[1]]=list()
							ResultsClust[[1]][[1]]=List[[i]]
							names(ResultsClust[[1]])[1]="Clust"
							element=element+1					
							ListNew[[element]]=ResultsClust[[1]]
						}
					}
					else{
						for (j in 1:length(List[[i]]$Results)){
							ResultsClust[[j]]=list()
							ResultsClust[[j]][[1]]=List[[i]]$Results[[j]]
							names(ResultsClust[[j]])[1]="Clust"
							element=element+1					
							ListNew[[element]]=ResultsClust[[j]]
							attr(ListNew[element],"method")="Weights"
						}		
					}		
				}	
			}
			
			if(is.null(names)){
				names=seq(1,length(ListNew),1)
				for(i in 1:length(ListNew)){
					names[i]=paste("Method",i,sep=" ")
				}
			}
			names(ListNew)=names
			MatrixClusters=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)
			List=ListNew	
			
			maxclus=0
			DataPrepared=list()
			for (k in 1:dim(MatrixClusters)[1]){
				message(k)
				clusters=MatrixClusters[k,]
				
				if(max(clusters)>maxclus){
					maxclus=max(clusters)
				}
				
				check<-plyr::try_default(PreparePathway(List[[k]],geneExpr,topG,sign),NULL,quiet=TRUE)
				if(is.null(check)){
					Temp=List[[k]]
					
					for(i in unique(clusters)){
						objects=list()
						objects$LeadCpds=names(clusters)[which(clusters==i)] 
						objects$OrderedCpds=stats::as.hclust(List[[k]]$Clust$Clust)$labels[stats::as.hclust(List[[k]]$Clust$Clust)$order]
						
						Temp[[i+1]]=list(objects=objects)
						names(Temp)[i+1]=paste("Cluster",i,sep=" ")
						
					}
					DataPrepared[[k]]<-PreparePathway(Temp,geneExpr,topG,sign)
					
				}	
			}
		}
		
		else{
			for(k in 1:length(List)){
				DataPrepared[[k]]<-plyr::try_default(PreparePathway(List[[k]],geneExpr,topG,sign),NULL,quiet=TRUE)
				if(is.null(DataPrepared[[k]])){
					Temp=List[[k]]
					
					for(i in unique(clusters)){
						objects=list()
						objects$LeadCpds=List[[k]]$objects$LeadCpds
						objects$OrderedCpds=List[[k]]$objects$OrderedCpds
						
						Temp[[i+1]]=list(objects=objects)
						names(Temp)[i+1]=paste("Cluster",i,sep=" ")
						
					}
					DataPrepared[[k]]<-PreparePathway(Temp,geneExpr,topG,sign)
					
				}	
			}
			
			
		}
		
		
		
		method.test = function(sign.method,path.method){
			method.choice = FALSE
			
			if( sign.method=="limma"  & path.method=="MLP"  ){
				method.choice = TRUE
			}
			if(method.choice==TRUE){
				return(list(sign.method=sign.method,path.method=path.method))
			}	
			else{
				stop("Incorrect choice of method.")
			}
			
		}
		
		method.out = method.test(method[1],method[2])
		
		sign.method = method.out$sign.method
		path.method = method.out$path.method
		
		if(length(geneInfo$ENTREZID)==1){
			geneInfo$ENTREZID = colnames(geneExpr)
		}
		
		# Determining the genesets if they were not given with the function input
		if((class(GENESET)=="geneSetMLP")[1] ){
			geneSet <- GENESET
		}
		else{
			geneSet <- MLP::getGeneSets(species = "Human",geneSetSource = geneSetSource,entrezIdentifiers = geneInfo$ENTREZID)
		}
		
		if(is.null(topP)){
			top1=FALSE
		}
		else{
			top1=TRUE
		}
		
		
		ResultMLP=list()
		
		for (k in 1:length(DataPrepared)){
			message(k)
			
			PathwaysResults=list()
			
			for (i in 1:length(DataPrepared[[k]]$pvalsgenes)){
				message(paste(k,i,sep='.'))
				temp=list()
				temp[[1]]=DataPrepared[[k]]$objects[[i]] # the objects
				temp[[2]]=DataPrepared[[k]]$Genes[[i]]		# the genes		
				pvalscluster=DataPrepared[[k]]$pvalsgenes[[i]]
				
				Entrezs=sapply(names(pvalscluster),function(x) return(geneInfo$ENTREZID[which(geneInfo$SYMBOL==x)]))
				
				if(path.method=="MLP"){
					## WE WILL USE THE RAW P-VALUES TO PUT IN MLP -> LESS GRANULAR
					
					names(pvalscluster) = Entrezs
					
					out.mlp <- MLP::MLP(
							geneSet = geneSet,
							geneStatistic = pvalscluster,
							minGenes = 5,
							maxGenes = 100,
							rowPermutations = TRUE,
							nPermutations = 100,
							smoothPValues = TRUE,
							probabilityVector = c(0.5, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999),df = 9,addGeneSetDescription=TRUE)
					
					output = list()
					#output$gene.p.values = p.adjust(p.values,method="fdr")
					
					#ranked.genesets.table = data.frame(genesets = (rownames(out.mlp)),p.values = as.numeric(out.mlp$geneSetPValue),descriptions = out.mlp$geneSetDescription)
					#ranked.genesets.table$genesets = as.character(ranked.genesets.table$genesets)
					#ranked.genesets.table$descriptions = as.character(ranked.genesets.table$descriptions)
					
					#if(is.null(topP)){
					#		topP=length(ranked.genesets.table$p.values<=sign)
					#}
					
					if(is.null(topP)){
						topP=length(which(out.mlp$geneSetPValue<=sign))
					}
					
					#TopPaths=ranked.genesets.table[1:topP,]
					#AllPaths=ranked.genesets.table
					
					TopPaths=out.mlp[1:topP,]
					AllPaths=out.mlp
					
					output$TopPaths=TopPaths
					attr(output$TopPaths,'geneSetSource')=geneSetSource
					output$AllPaths=AllPaths
					attr(output$AllPaths,'geneSetSource')=geneSetSource
					#output$ranked.genesets.table = ranked.genesets.table[ranked.genesets.table$p.values<=sign,]
					
					
					#nr.genesets = c( dim(ranked.genesets.table)[1]  ,  length(geneSet) 	)
					#names(nr.genesets) = c("used.nr.genesets","total.nr.genesets")
					#output$nr.genesets = nr.genesets
					
					#output$object = out.mlp
					#output$method = "MLP"
					
					temp[[3]]=output				
				}
				names(temp)=c("objects","Genes","Pathways")
				PathwaysResults[[i]]=temp
				names(PathwaysResults)[i]=paste("Cluster",i,sep=" ")
				
			}
			
			ResultMLP[[k]]=PathwaysResults	
		}
		names(ResultMLP)=names
		for(i in 1:length(ResultMLP)){
			for(k in 1:length(ResultMLP[[i]])){
				if(is.null(ResultMLP[[i]][[k]])[1]){
					ResultMLP[[i]][[k]]=NA
					names(ResultMLP[[i]])[k]=paste("Cluster",k,sep=" ")
				}			
			}
			if(length(ResultMLP[[i]]) != maxclus){
				extra=maxclus-length(ResultMLP[[i]])
				for(j in 1:extra){
					ResultMLP[[i]][[length(ResultMLP[[i]])+j]]=NA
					names(ResultMLP[[i]])[length(ResultMLP[[i]])]=paste("Cluster",length(ResultMLP[[i]]),sep=" ")
				}
			}
		} 	
	}
	return(ResultMLP)	
}

#' @title Iterations of the pathway analysis
#' 
#' @description The MLP method to perform pathway analysis is based on resampling of the
#' data. Therefore it is recommended to perform the pathway analysis multiple
#' times to observe how much the results are influenced by a different
#' resample. The function \code{PathwaysIter} performs the pathway analysis as
#' described in \code{Pathways} a specified number of times. The input can be
#' one data set or a list as in \code{Pathway.2} and \code{Pathways}.
#' @export PathwaysIter
#' 
#' @param List A list of clustering outputs or output of the\code{DiffGenes}
#' function. The first element of the list will be used as the reference in
#' \code{ReorderToReference}. The output of \code{ChooseFeatures} is also
#' accepted.
#' @param Selection If pathway analysis should be conducted for a specific
#' selection of objects, this selection can be provided here. Selection can
#' be of the type "character" (names of the objects) or "numeric" (the number
#' of specific cluster). Default is NULL.
#' @param geneExpr The gene expression matrix of the objects. The rows should
#' correspond with the genes.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param method The method to applied to look for differentially expressed genes and related pathways. For now, only the
#' limma method is available for gene analysis and the MLP method for pathway analysis. Default is c("limma","MLP").
#' @param geneInfo A data frame with at least the columns ENTREZID and SYMBOL.
#' This is necessary to connect the symbolic names of the genes with their
#' EntrezID in the correct order. The order of the gene is here not in the
#' order of the rownames of the gene expression matrix but in the order of
#' their significance. Default is NULL.
#' @param geneSetSource The source for the getGeneSets function ("GOBP",
#' "GOMF","GOCC", "KEGG" or "REACTOME"). Default is "GOBP".
#' @param topP Overrules sign. The number of pathways to display for each
#' cluster. If not specified, only the significant genes are shown. Default is NULL.
#' @param topG Overrules sign. The number of top genes to be returned in the
#' result. If not specified, only the significant genes are shown. Default is NULL.
#' @param GENESET Optional. Can provide own candidate gene sets. Default is NULL.
#' @param sign The significance level to be handled. Default is 0.05.
#' @param niter The number of times to perform pathway analysis. Default is 10.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' @param names Optional. Names of the methods. Default is NULL.
#' @return This element is again a list with the following four elements:
#' \item{objects}{A list with the elements LeadCpds (the objects of
#' interest) and OrderedCpds (all objects in the order of the clustering
#' result)} \item{Characteristics}{The found (top) characteristics of the
#' feauture data} \item{Genes}{A list with the elements TopDE (a table with
#' information on the top genes) and AllDE (a table with information on all
#' genes)} \item{Pathways}{A list with the element ranked.genesets.table which
#' is a data frame containing the genesets, their p-values and their
#' descriptions. The second element is nr.genesets and contains the used and
#' total number of genesets.}
#' @examples
#' 
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(geneMat)
#' data(GeneInfo)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' L=list(MCF7_F,MCF7_T)
#' names=c('FP','TP')
#' 
#' MCF7_Paths_FandT=PathwaysIter(List=L, geneExpr = geneMat, nrclusters = 7, method = 
#' c("limma", "MLP"), geneInfo = GeneInfo, geneSetSource = "GOBP", topP = NULL, 
#' topG = NULL, GENESET = NULL, sign = 0.05,niter=2,fusionsLog = TRUE, 
#' weightclust = TRUE, names =names)
#' }
PathwaysIter<-function(List,Selection=NULL,geneExpr=NULL,nrclusters=NULL,method=c("limma", "MLP"),geneInfo=NULL,geneSetSource = "GOBP",topP=NULL,topG=NULL,GENESET=NULL,sign=0.05,niter=10,fusionsLog=TRUE,weightclust=TRUE,names=NULL){
	if (!requireNamespace("MLP", quietly = TRUE)) {
		stop("MLP needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	if (!requireNamespace("biomaRt", quietly = TRUE)) {
		stop("biomaRt needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
		stop("org.Hs.eg.db needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	PathwaysOutput = list() 
	for (i in 1:niter){
		message(paste("Iteration",i,sep=" "))
		mlp = Pathways(List,Selection,geneExpr,nrclusters,method,geneInfo,geneSetSource,topP,topG,GENESET,sign=sign,fusionsLog,weightclust,names)
		PathwaysOutput [[length(PathwaysOutput )+1]] = mlp
		names(PathwaysOutput )[i]=paste("Iteration",i,sep=" ")
	}
	
	
	return(PathwaysOutput)
}

#' @title Pathway analysis for a selection of objects
#' @param List A list of clustering outputs or output of the\code{DiffGenes}
#' function. The first element of the list will be used as the reference in
#' \code{ReorderToReference}. The output of \code{ChooseFeatures} is also
#' accepted.
#' @param Selection If pathway analysis should be conducted for a specific
#' selection of objects, this selection can be provided here. Selection can
#' be of the type "character" (names of the objects) or "numeric" (the number
#' of specific cluster). Default is NULL.
#' @param geneExpr The gene expression matrix or ExpressionSet of the objects.
#' The rows should correspond with the genes.
#' @param nrclusters Optional. The number of clusters to cut the dendrogram in.
#' The number of clusters should not be specified if the interest lies only in
#' a specific selection of objects which is known by name.  Otherwise, it is
#' required. Default is NULL.
#' @param method The method to applied to look for differentially expressed genes and related pathways. For now, only the
#' limma method is available for gene analysis and the MLP method for pathway analysis. Default is c("limma","MLP").
#' @param geneInfo A data frame with at least the columns ENTREZID and SYMBOL.
#' This is necessary to connect the symbolic names of the genes with their
#' EntrezID in the correct order. The order of the gene is here not in the
#' order of the rownames of the gene expression matrix but in the order of
#' their significance. Default is NULL.
#' @param geneSetSource The source for the getGeneSets function, defaults to
#' "GOBP".
#' @param topP Overrules sign. The number of pathways to display for each
#' cluster. If not specified, only the significant genes are shown. Default is NULL.
#' @param topG Overrules sign. The number of top genes to be returned in the
#' result. If not specified, only the significant genes are shown. Defaults is NULL.
#' @param GENESET Optional. Can provide own candidate gene sets. Default is NULL.
#' @param sign The significance level to be handled. Default is 0.05.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' @param names Optional. Names of the methods. Default is NULL.
#' @description Internal function of \code{Pathways}.
PathwaysSelection<-function(List=NULL,Selection,geneExpr=NULL,nrclusters=NULL,method=c("limma", "MLP"),geneInfo=NULL,geneSetSource = "GOBP",topP=NULL,topG=NULL,GENESET=NULL,sign=0.05,fusionsLog=TRUE,weightclust=TRUE,names=NULL){
	if (!requireNamespace("MLP", quietly = TRUE)) {
		stop("MLP needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	if (!requireNamespace("biomaRt", quietly = TRUE)) {
		stop("biomaRt needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
		stop("org.Hs.eg.db needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	
	method.test = function(sign.method,path.method){
		method.choice = FALSE
		
		if( sign.method=="limma"  & path.method=="MLP"  ){
			method.choice = TRUE
		}
		if(method.choice==TRUE){
			return(list(sign.method=sign.method,path.method=path.method))
		}	
		else{
			stop("Incorrect choice of method.")
		}
		
	}
	
	method.out = method.test(method[1],method[2])
	
	sign.method = method.out$sign.method
	path.method = method.out$path.method
	
	if(length(geneInfo$ENTREZID)==1){
		geneInfo$ENTREZID = colnames(geneExpr)
	}
	
	# Determining the genesets if they were not given with the function input
	if((class(GENESET)=="geneSetMLP")[1] ){
		geneSet <- GENESET
	}
	else{
		geneSet <- MLP::getGeneSets(species = "Human",geneSetSource = geneSetSource,entrezIdentifiers = geneInfo$ENTREZID)
	}
	
	
	if(class(Selection)=="character"){
		ResultMLP=list()
		
		DataPrepared<-plyr::try_default(PreparePathway(List[[1]],geneExpr,topG,sign),NULL,quiet=TRUE)
		if(is.null(DataPrepared)){
			Temp=List[[1]]
			objects=list()
			objects$LeadCpds=Selection
			objects$OrderedCpds=colnames(geneExpr)
			Temp[[length(Temp)+1]]=list(objects=objects)
			names(Temp)[length(Temp)]=paste("Cluster")
			
			DataPrepared<-PreparePathway(Temp,geneExpr,topG,sign)
		}
		
		
		temp=list()
		temp[[1]]=DataPrepared$objects[[1]] #names of the objects
		temp[[2]]=DataPrepared$Genes[[1]]
		
		
		
		
		if(path.method=="MLP"){
			## WE WILL USE THE RAW P-VALUES TO PUT IN MLP -> LESS GRANULAR
			p.values=DataPrepared$pvalsgenes[[1]]
			
			Entrezs=sapply(names(p.values),function(x) return(geneInfo$ENTREZID[which(geneInfo$SYMBOL==x)]))
			
			names(p.values) = Entrezs
			out.mlp <- MLP::MLP(
					geneSet = geneSet,
					geneStatistic = p.values,
					minGenes = 5,
					maxGenes = 100,
					rowPermutations = TRUE,
					nPermutations = 100,
					smoothPValues = TRUE,
					probabilityVector = c(0.5, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999),df = 9)
			
			output = list()
			#output$gene.p.values = p.adjust(p.values,method="fdr")
			
			#ranked.genesets.table = data.frame(genesets = (rownames(out.mlp)),p.values = as.numeric(out.mlp$geneSetPValue),descriptions = out.mlp$geneSetDescription)
			#ranked.genesets.table$genesets = as.character(ranked.genesets.table$genesets)
			#ranked.genesets.table$descriptions = as.character(ranked.genesets.table$descriptions)
			
			#if(is.null(topP)){
			#		topP=length(ranked.genesets.table$p.values<=sign)
			#}
			
			if(is.null(topP)){
				topP=length(which(out.mlp$geneSetPValue<=sign))
			}
			
			#TopPaths=ranked.genesets.table[1:topP,]
			#AllPaths=ranked.genesets.table
			
			TopPaths=out.mlp[1:topP,]
			AllPaths=out.mlp
			
			output$TopPaths=TopPaths
			attr(output$TopPaths,'geneSetSource')=geneSetSource
			output$AllPaths=AllPaths
			attr(output$AllPaths,'geneSetSource')=geneSetSource
			#output$ranked.genesets.table = ranked.genesets.table[ranked.genesets.table$p.values<=sign,]
			
			#nr.genesets = c( dim(ranked.genesets.table)[1]  ,  length(geneSet) )
			#names(nr.genesets) = c("used.nr.genesets","total.nr.genesets")
			#output$nr.genesets = nr.genesets
			
			#output$object = out.mlp
			#output$method = "MLP"
			
			temp[[3]]=output				
			
			
		}
		names(temp)=c("objects","Genes","Pathways")			
		ResultMLP[[1]]=temp
		names(ResultMLP)="Selection"
	}
	
	
	
	else if(class(Selection)=="numeric" & !(is.null(List))){
		check<-plyr::try_default(PreparePathway(List[[1]],geneExpr,topG,sign),NULL,quiet=TRUE)
		if(is.null(check)){
			
			ListNew=list()
			element=0
			for(i in 1:length(List)){
				if(attributes(List[[i]])$method != "CEC" & attributes(List[[i]])$method != "Weighted" & attributes(List[[i]])$method!= "WeightedSim"){
					ResultsClust=list()
					ResultsClust[[1]]=list()
					ResultsClust[[1]][[1]]=List[[i]]
					names(ResultsClust[[1]])[1]="Clust"
					element=element+1					
					ListNew[[element]]=ResultsClust[[1]]
					#attr(ListNew[element],"method")="Weights"
				}
				else if(attributes(List[[i]])$method=="CEC" | attributes(List[[i]])$method=="Weighted" | attributes(List[[i]])$method == "WeightedSim"){
					ResultsClust=list()
					if(weightclust==TRUE){
						ResultsClust[[1]]=list()
						if(attributes(List[[i]])$method != "WeightedSim"){
							ResultsClust[[1]][[1]]=List[[i]]$Clust
							names(ResultsClust[[1]])[1]="Clust"
							element=element+1					
							ListNew[[element]]=ResultsClust[[1]]
							attr(ListNew[element],"method")="Weights"
						}
						else{
							ResultsClust[[1]]=list()
							ResultsClust[[1]][[1]]=List[[i]]
							names(ResultsClust[[1]])[1]="Clust"
							element=element+1					
							ListNew[[element]]=ResultsClust[[1]]
						}
					}
					else{
						for (j in 1:length(List[[i]]$Results)){
							ResultsClust[[j]]=list()
							ResultsClust[[j]][[1]]=List[[i]]$Results[[j]]
							names(ResultsClust[[j]])[1]="Clust"
							element=element+1					
							ListNew[[element]]=ResultsClust[[j]]
							attr(ListNew[element],"method")="Weights"
						}		
					}		
				}	
			}
			
			if(is.null(names)){
				names=seq(1,length(ListNew),1)
				for(i in 1:length(ListNew)){
					names[i]=paste("Method",i,sep=" ")
				}
			}
			names(ListNew)=names
			Matrix=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)
			List=ListNew	
			
			DataPrepared=list()
			for (k in 1:dim(Matrix)[1]){
				
				cluster=Selection
				
				check<-plyr::try_default(PreparePathway(List[[k]],geneExpr,topG,sign),NULL,quiet=TRUE)
				if(is.null(check)){
					Temp=List[[k]]
					objects=list()
					objects$LeadCpds=colnames(Matrix)[which(Matrix[k,]==cluster)]
					objects$OrderedCpds=stats::as.hclust(List[[k]]$Clust$Clust)$labels[stats::as.hclust(List[[k]]$Clust$Clust)$order]
					Temp[[length(Temp)+1]]=list(objects=objects)
					names(Temp)[length(Temp)]=paste("Cluster")
					
					DataPrepared[[k]]<-PreparePathway(Temp,geneExpr,topG,sign)
				}
			}
			names(DataPrepared)=names
		}
		
		else{
			for(k in 1:length(List)){
				DataPrepared[[k]]<-plyr::try_default(PreparePathway(List[[k]],geneExpr,topG,sign),NULL,quiet=TRUE)
				if(is.null(DataPrepared[[k]])){
					Temp=List[[k]]
					
					
					objects=list()
					objects$LeadCpds=List[[k]]$objects$LeadCpds
					objects$OrderedCpds=List[[k]]$objects$OrderedCpds
					
					Temp[[length(Temp)+1]]=list(objects=objects)
					names(Temp)[length(Temp)]=paste("Cluster")
					
					DataPrepared[[k]]<-PreparePathway(Temp,geneExpr,topG,sign)
					
				}	
			}
			
			
		}
		
		ResultMLP=list()
		for (k in 1:length(DataPrepared)){
			message(k)
			cluster=Selection
			
			
			temp=list()
			temp[[1]]=DataPrepared[[k]]$objects[[1]] #names of the objects
			temp[[2]]=DataPrepared[[k]]$Genes[[1]]
			
			
			if(path.method=="MLP"){
				## WE WILL USE THE RAW P-VALUES TO PUT IN MLP -> LESS GRANULAR
				p.values=DataPrepared[[k]]$pvalsgenes[[1]]
				Entrezs=sapply(names(p.values),function(x) return(geneInfo$ENTREZID[which(geneInfo$SYMBOL==x)]))
				
				names(p.values) = Entrezs
				
				out.mlp <- MLP::MLP(
						geneSet = geneSet,
						geneStatistic = p.values,
						minGenes = 5,
						maxGenes = 100,
						rowPermutations = TRUE,
						nPermutations = 100,
						smoothPValues = TRUE,
						probabilityVector = c(0.5, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999),df = 9)
				
				output = list()
				#output$gene.p.values = p.adjust(p.values,method="fdr")
				
				#ranked.genesets.table = data.frame(genesets = (rownames(out.mlp)),p.values = as.numeric(out.mlp$geneSetPValue),descriptions = out.mlp$geneSetDescription)
				#ranked.genesets.table$genesets = as.character(ranked.genesets.table$genesets)
				#ranked.genesets.table$descriptions = as.character(ranked.genesets.table$descriptions)
				
				#if(is.null(topP)){
				#		topP=length(ranked.genesets.table$p.values<=sign)
				#}
				
				if(is.null(topP)){
					topP=length(which(out.mlp$geneSetPValue<=sign))
				}
				
				#TopPaths=ranked.genesets.table[1:topP,]
				#AllPaths=ranked.genesets.table
				
				TopPaths=out.mlp[1:topP,]
				AllPaths=out.mlp
				
				output$TopPaths=TopPaths
				attr(output$TopPaths,'geneSetSource')=geneSetSource
				output$AllPaths=AllPaths
				attr(output$AllPaths,'geneSetSource')=geneSetSource
				#output$ranked.genesets.table = ranked.genesets.table[ranked.genesets.table$p.values<=sign,]
				
				#nr.genesets = c( dim(ranked.genesets.table)[1]  ,  length(geneSet) )
				#names(nr.genesets) = c("used.nr.genesets","total.nr.genesets")
				#output$nr.genesets = nr.genesets
				
				#output$object = out.mlp
				#output$method = "MLP"
				
				temp[[3]]=output				
			}
			names(temp)=c("objects","Genes","Pathways")			
			ResultMLP[[k]]=temp
			
		}
		names(ResultMLP)=names
		
	}
	else{
		message("If a specific cluster is specified, clustering results must be provided in List")
	}
	
	return(ResultMLP)	
}

#' @title A GO plot of a pathway analysis output.
#' 
#' @description The \code{PlotPathways} function takes an output of the
#' \code{PathwayAnalysis} function and plots a GO graph with the help of the
#' \code{plotGOgraph} function of the MLP package.
#' @param Pathways One element of the output list returned by
#' \code{PathwayAnalysis} or \code{Geneset.intersect}.
#' @param nRow Number of GO IDs for which to produce the plot. Default is 5.
#' @param main Title of the plot. Default is NULL.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return The output is a GO graph.
#' @examples
#' 
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(geneMat)
#' data(GeneInfo)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' L=list(MCF7_F,MCF7_T)
#' names=c('FP','TP')
#' 
#' MCF7_PathsFandT=PathwayAnalysis(List=L, geneExpr = geneMat, nrclusters = 7, method = c("limma", 
#' "MLP"), geneInfo = GeneInfo, geneSetSource = "GOBP", topP = NULL, 
#' topG = NULL, GENESET = NULL, sign = 0.05,niter=2,fusionsLog = TRUE, weightclust = TRUE, 
#'  names =names,seperatetables=FALSE,separatepvals=FALSE)
#'  
#' PlotPathways(MCF7_PathsFandT$FP$"Cluster 1"$Pathways,nRow=5,main=NULL)
#' }
#' 
#' @export PlotPathways
PlotPathways<-function(Pathways,nRow=5,main=NULL,plottype="new",location=NULL){	
	
	if (!requireNamespace("MLP", quietly = TRUE)) {
		stop("MLP needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	if (!requireNamespace("biomaRt", quietly = TRUE)) {
		stop("biomaRt needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
		stop("org.Hs.eg.db needed for this function to work. Please install it.",
				call. = FALSE)
	}
	
	
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	#preparing data structure for the plotGOGraph
	colnames(Pathways)[3:ncol(Pathways)]=sub("mean_","",colnames(Pathways)[3:ncol(Pathways)])	
	#plot GOgraph
	plottypein(plottype,location)
	MLP::plotGOgraph(Pathways,nRow=nRow,main=main)
	plottypeout(plottype)
	
}

#' @title Preparing a data set for pathway analysis
#' 
#' @description The functions for pathway analysis in this package can also work on results
#' of the integrated data functions. However, a differential gene expression
#' needs to be conducted to perform pathway analysis. The function
#' \code{PreparePathway} checks if the necessary elements are present in the
#' data structures and if not, the elements such as p-values are created. It is
#' an internal function to all pathway analysis functions but can be used
#' separately as well.
#' 
#' 
#' @param Object A list with at least an element with the name "objects" such
#' that the function knows which objects to test for differential gene
#' expression. If the elements "Genes" and "pvalsgenes" are present as well,
#' these will be collected and the gene expression is not analyzed.
#' @param geneExpr The gene expression matrix or ExpressionSet of the objects.
#' The rows should correspond with the genes.
#' @param topG Overrules sign. The number of top genes to be returned in the
#' result. If not specified, only the significant genes are shown. Default is NULL.
#' @param sign The significance level to be handled. Default is 0.05.
#' @return The returned value is a list with three elements: \item{pvalsgenses
#' }{This is a list with that contains a vector of raw p-values for every group
#' of tested objects.} \item{objects}{This is a list with that contains
#' another list per group of tested objects. Every list contains the lead
#' objects and the ordered objects.} \item{Genes }{This is a list with that
#' contains contains another list per group of tested objects. Every list
#' contains two data frames, one with information on the top genes and one with
#' information on all genes.}
#' @examples
#' 
#' data(fingerprintMat)
#' data(geneMat)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' 
#' L1=list(MCF7_F)
#' 
#' Comps1=FindCluster(L1,nrclusters=7,select=c(1,1))
#' Comps2=FindCluster(L1,nrclusters=7,select=c(1,2))
#' Comps3=FindCluster(L1,nrclusters=7,select=c(1,3))
#' 
#' 
#' L2=list()
#' 
#' L2$'Cluster 1'$objects$LeadCpds=Comps1
#' L2$'Cluster 2'$objects$LeadCpds=Comps2
#' L2$'Cluster 3'$objects$LeadCpds=Comps2
#' 
#' MCF7_PreparePaths=PreparePathway(Object=L2,geneExpr=geneMat,topG=NULL,sign=0.05)
#' str(MCF7_PreparePaths)
#' 
#' @export PreparePathway
PreparePathway<-function(Object,geneExpr,topG,sign){
	FoundGenes=NULL
	FoundComps=NULL
	
	FoundGenes=FindElement("Genes",Object)
	
	if(is.null(FoundGenes)|(is.list(FoundGenes) & length(FoundGenes) == 0)){
		FoundComps=FindElement("objects",Object)
		if(is.null(FoundComps)|(is.list(FoundComps) & length(FoundComps) == 0)){
			stop("Specify either the p-values of the genes or a selection of objects to test for DE genes.")
		}
		
		pvalsgenes=list()
		FoundGenes=list()
		CompsP=list()
		TopDEP=list()
		for(i in 1:length(FoundComps)){
			LeadCpds=FoundComps[[i]]$LeadCpds
			CompsP[[i]]=FoundComps[[i]]
			names(CompsP)[[i]]=paste("objects_",i,sep="")
			if(is.null(LeadCpds)){
				stop("In the objects element, specify an element LeadCpds")
			} 
			
			group <- factor(ifelse(colnames(geneExpr) %in% LeadCpds, 1, 0))
			
			if(class(geneExpr)[1]=="ExpressionSet"){
				geneExpr$LeadCmpds<-group		
				if (!requireNamespace("a4Base", quietly = TRUE)) {
					stop("a4Base needed for this function to work. Please install it.",
							call. = FALSE)
				}
				
				DElead <- a4Base::limmaTwoLevels(geneExpr,"LeadCmpds")
				
				allDE <- a4Core::topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL),sort.by="p")
				
				if(is.null(allDE$ID)){
					allDE$Genes <- rownames(allDE)
				}
				else
				{
					allDE$Genes=allDE$ID
				}
				
				if(is.null(topG)){
					topG=length(which(allDE$adj.P.Val<=sign))
				}
				TopDE <- allDE[1:topG, ]
				
				Genes <- list(TopDE,allDE)	
				names(Genes)<-c("TopDE","AllDE") 
			}
			else{				
				label.factor = factor(group)
				design = stats::model.matrix(~label.factor)
				fit = limma::lmFit(geneExpr,design=design)
				fit = limma::eBayes(fit)
				allDE = limma::topTable(fit,coef=2,adjust="fdr",n=dim(geneExpr)[1], sort.by="p")
				
				if(is.null(allDE$ID)){
					allDE$ID <- rownames(allDE)
				}
				
				if(is.null(topG)){
					topG=length(which(allDE$adj.P.Val<=sign))
				}
				
				
				TopDE<-allDE[1:topG,]
				
				Genes <- list(TopDE,allDE)	
				names(Genes)<-c("TopDE","AllDE") 
				
			}
			FoundGenes[[i]]=Genes
			TopDEP[[i]]=FoundGenes[[i]]
			names(TopDEP)[i]=paste("genes_",i,sep="")
			names(FoundGenes)=paste("Genes_",i,sep="")	
			pvalsgenes[[i]]=Genes$AllDE$P.Value
			names(pvalsgenes[[i]])=Genes$AllDE$ID
			names(pvalsgenes)[i]=paste("pvals_",i,sep="")
		}
	}	
	
	else{
		pvalsgenes=list()
		TopDEP=list()
		for(i in 1:length(FoundGenes)){
			#names(FoundGenes)[i]=paste("Genes_",i,sep="")
			TopDEP[[i]]=FoundGenes[[i]]
			names(TopDEP)[i]=paste("genes_",i,sep="")
			pvalsgenes[[i]]=FoundGenes[[i]]$AllDE$P.Value
			names(pvalsgenes[[i]])=FoundGenes[[i]]$AllDE$ID
			names(pvalsgenes)[i]=paste("pvals_",i,sep="")
		}
		
		FoundComps=FindElement("objects",Object)
		CompsP=list()
		for(i in 1:length(FoundComps)){
			if(is.null(FoundComps[[i]]$LeadCpds)){
				CompsP[[i]]="No LeadCpds specified"
			}
			else{
				CompsP[[i]]=FoundComps[[i]]
				names(CompsP)[[i]]=paste("objects_",i,sep="")
			}
		}
	}
	
	
	return(list(pvalsgenes=pvalsgenes,objects=CompsP,Genes=TopDEP))
}

#' @title Plotting gene profiles
#' 
#' @description In \code{ProfilePlot}, the gene profiles of the significant genes for a
#' specific cluster are shown on 1 plot. Therefore, each gene is normalized by
#' subtracting its the mean.
#' 
#' @param Genes The genes to be plotted.
#' @param Comps The objects to be plotted or to be separated from the other
#' objects.
#' @param geneExpr The gene expression matrix or ExpressionSet of the objects.
#' @param raw Logical. Should raw p-values be plotted? Default is FALSE.
#' @param orderLab Optional. If the objects are to set in a specific order of
#' a specific method. Default is NULL.
#' @param colorLab The clustering result that determines the color of the
#' labels of the objects in the plot. Default is NULL.
#' @param nrclusters Optional. The number of clusters to cut the dendrogram in.
#' @param cols Optional. The color to use for the objects in Clusters for each
#' method.
#' @param addLegend Optional. Whether a legend of the colors should be added to
#' the plot.
#' @param margins Optional. Margins to be used for the plot. Default is margins=c(8.1,4.1,1.1,6.5).
#' @param extra The space between the plot and the legend. Default is 5.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return A plot which contains multiple gene profiles. A distinction is made
#' between the values for the objects in Comps and the others.
#' @examples
#' 
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(geneMat)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#'		method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#'		method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#'
#' L=list(MCF7_F,MCF7_T)
#' names=c('FP','TP')
#'
#' MCF7_FT_DE = DiffGenes(List=L,geneExpr=geneMat,nrclusters=7,method="limma",sign=0.05,topG=10,
#' fusionsLog=TRUE,weightclust=TRUE)
#'
#' Comps=SharedComps(list(MCF7_FT_DE$`Method 1`$"Cluster 1",MCF7_FT_DE$`Method 2`$"Cluster 1"))[[1]]
#'
#' MCF7_SharedGenes=FindGenes(dataLimma=MCF7_FT_DE,names=c("FP","TP"))
#'
#' Genes=names(MCF7_SharedGenes[[1]])[-c(2,4,5)]
#'
#' colscl=ColorPalette(colors=c("red","green","purple","brown","blue","orange"),ncols=9)
#'
#' ProfilePlot(Genes=Genes,Comps=Comps,geneExpr=geneMat,raw=FALSE,orderLab=MCF7_F,
#' colorLab=NULL,nrclusters=7,cols=colscl,addLegend=TRUE,margins=c(16.1,6.1,1.1,13.5),
#' extra=4,plottype="sweave",location=NULL)
#' }
#' 
#' @export ProfilePlot
ProfilePlot<-function(Genes,Comps,geneExpr=NULL,raw=FALSE,orderLab=NULL,colorLab=NULL,nrclusters=NULL,cols=NULL,addLegend=TRUE,margins=c(8.1,4.1,1.1,6.5),extra=5,plottype="new",location=NULL){
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	
	
	
	if(class(geneExpr)[1]=="ExpressionSet"){
		geneExpr <- Biobase::exprs(geneExpr)
		
	}
	
	if(!is.null(orderLab)){
		if(class(orderLab)=="character"){
			orderlabs=orderLab
		}
		else{
			orderlabs=orderLab$Clust$order.lab
			geneExpr=geneExpr[,match(orderlabs,colnames(geneExpr))]
		}
	}
	else{
		orderlabs=colnames(geneExpr)
	}
	
	if(!is.null(colorLab)){
		Data1 <- colorLab$Clust
		ClustData1=stats::cutree(Data1,nrclusters) 
		
		ordercolors=ClustData1[Data1$order]
		names(ordercolors)=Data1$order.lab
		
		ClustData1=ClustData1[Data1$order]	
		
		order=seq(1,nrclusters)
		
		for (k in 1:length(unique(ClustData1))){
			select=which(ClustData1==unique(ClustData1)[k])
			ordercolors[select]=order[k]
		}
		
		if(!is.null(orderLab)){
			if(class(orderLab)=="character"){
				ordernames=orderLab
			}
			else{
				ordernames=orderLab$Clust$order.lab
			}	
			ordercolors=ordercolors[ordernames]
		}
		
		colors<- cols[ordercolors]
		names(colors) <-names(ordercolors)	
		
	}
	else{
		colors1<-rep("green",length(Comps))
		colors2<-rep("black",length(orderlabs[which(!(orderlabs%in%Comps))]))
		colors=c(colors1,colors2)
		AllCpds=c(Comps,orderlabs[which(!(orderlabs%in%Comps))])
		names(colors)=AllCpds
	}
	
	#yvalues=c()
	#allvalues=c()
	#for(i in 1:length(Genes)){
	#	yvalues=as.vector(GeneExpr[which(rownames(GeneExpr)==Genes[i]),])
	#	allvalues=c(allvalues,yvalues-mean(yvalues))
	#}	
	yvalues=geneExpr[Genes,]
	if(raw==FALSE & class(yvalues) != "numeric"){
		allvalues=as.vector(apply(yvalues,1,function(c) c-mean(c)))
	}
	else if(raw==FALSE & class(yvalues) == "numeric"){
		allvalues=as.vector(sapply(yvalues,function(c) c-mean(yvalues)))
	}
	else{
		allvalues=as.vector(yvalues)
	}	
	ylims=c(min(allvalues)-0.1,max(allvalues)+0.1)
	
	plottypein(plottype,location)
	graphics::par(mar=margins,xpd=TRUE)
	graphics::plot(type="n",x=0,y=0,xlim=c(0,ncol(geneExpr)),ylim=ylims,ylab=expression(log[2] ~ paste("fold ", "change")),xlab=" ",xaxt="n",cex.axis=1.5,cex.lab=2)
	#ylims=c()
	Indices=c(colnames(geneExpr)[which(colnames(geneExpr)%in%Comps)],colnames(geneExpr)[which(!colnames(geneExpr)%in%Comps)])
	
	for(i in 1:length(Genes)){
		GenesComps=as.numeric(geneExpr[which(rownames(geneExpr)==Genes[i]),colnames(geneExpr)%in%Comps])
		Others=as.numeric(geneExpr[which(rownames(geneExpr)==Genes[i]),!(colnames(geneExpr)%in%Comps)])
		if(length(Others)==0){
			Continue=FALSE
		}else{Continue=TRUE}
		
		#ylims=c(ylims,c(GenesComps,Others)-mean(c(GenesComps,Others)))	
		if(raw==FALSE){
			yvalues1=GenesComps-mean(c(GenesComps,Others))	
		}
		else{
			yvalues1=GenesComps
		}
		
		graphics::lines(x=seq(1,length(GenesComps)),y=yvalues1,lty=1,col=i,lwd=1.6)
		#points(x=seq(1,length(GenesComps)),y=yvalues1,pch=19,col=i)
		graphics::segments(x0=1,y0=mean(yvalues1[1:length(GenesComps)]),x1=length(GenesComps),y1=mean(yvalues1[1:length(GenesComps)]),lwd=1.5,col=i)
		
		
		if(Continue==TRUE){
			if(raw==FALSE){
				yvalues2=Others-mean(c(GenesComps,Others))	
			}
			else{
				yvalues2=Others
			}
			
			
			graphics::lines(x=seq(length(GenesComps)+1,ncol(geneExpr)),y=yvalues2,lty=1,col=i,lwd=1.6)
			graphics::segments(x0=length(GenesComps)+1,y0=mean(yvalues2[1:length(Others)]),x1=ncol(geneExpr),y1=mean(yvalues2[1:length(Others)]),lwd=1.5,col=i)
			
		}
		
		
	}	
	#Indices=c(colnames(GeneExpr)[which(colnames(GeneExpr)%in%Comps)],colnames(GeneExpr)[which(!colnames(GeneExpr)%in%Comps)])
	if(!is.null(colorLab)){
		graphics::axis(1, labels=FALSE)
		#box("outer")
		graphics::mtext(substr(Indices,1,15), side = 1,  at=seq(0.5,(ncol(geneExpr)-0.5)), line=0.2, las=2, cex=1.5,col=colors[Indices])
		#mtext(substr(Indices,1,15), side = 1,  at=seq(0.5,(ncol(GeneExpr)-0.5)), line=0.2, las=2, cex=0.70,col=c(rep("blue",7),rep("black",(56-7))))
	}
	else{
		#axis(1,at=seq(0.5,(ncol(GeneExpr)-0.5)),labels=Indices,las=2,cex.axis=0.70,xlab=" ",col=colors[Indices])
		#axis(1,at=seq(0.5,(ncol(GeneExpr)-0.5)), labels=FALSE)
		graphics::mtext(Indices, side = 1,  at=seq(0.5,(ncol(geneExpr)-0.5)),line=0.2, las=2, cex=1.5,col=colors[Indices])
	}
	graphics::axis(2,ylab=expression(log[2] ~ paste("fold ", "change")),cex=2,cex.axis=1.5)
	
	if(addLegend==TRUE){
		
		labels=Genes
		colslegend=seq(1,length(Genes))
		
		graphics::par(xpd=T,mar=margins)
		graphics::legend(ncol(geneExpr)+extra,max(ylims),legend=c(labels),col=c(colslegend),lty=1,lwd=3,cex=1.5)
		
	}
	plottypeout(plottype)
}

#' @title Determines an optimal number of clusters based on silhouette widths
#' 
#' @description The function \code{SelectnrClusters} determines an optimal optimal number of
#' clusters based by calculating silhouettes widths for a sequence of clusters.
#' See "Details" for a more elaborate description.
#' 
#' If the object provided in List are data or distance matrices clustering
#' around medoids is performed with the \code{pam} function of the
#' \pkg{cluster} package. Of the obtained pam objects, average silhouette
#' widths are retrieved. A silhouette width represents how well an object lies
#' in its current cluster. Values around one are an indication of an
#' appropriate clustering while values around zero show that the object might
#' as well lie in the neighbouring cluster. The average silhouette width is a
#' measure of how tightly grouped the data is.  This is performed for every
#' number of cluster for every object provided in List. Then the average is
#' taken for every number of clusters over the provided objects. This results
#' in one average value per number of clusters. The number width the maximal
#' average silhouette width is chosen as the optimal number of clusters.
#' 
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clusters".	
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}. 
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param nrclusters A sequence of numbers of clusters to cut the dendrogram in. Default is a sequence of 5 to 25.
#' @param names The labels to give to the elements in List. Default is NULL.
#' @param StopRange Logical. Indicates whether the distance matrices with
#' values not between zero and one should be standardized to have so. If FALSE
#' the range normalization is performed. See \code{Normalization}. If TRUE, the
#' distance matrices are not changed. This is recommended if different types of
#' data are used such that these are comparable. Default is FALSE.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return A plots are made showing the average silhouette widths of the
#' provided objects for each number of clusters. Further, a list with two
#' elements is returned: \item{Silhouette_Widths}{A data frame with the
#' silhouette widths for each object and the average silhouette widths per
#' number of clusters} \item{Optimal_Nr_of_CLusters}{The determined optimal
#' number of cluster }
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' 
#' L=list(fingerprintMat,targetMat)
#' 
#' NrClusters=SelectnrClusters(List=L,type="data",distmeasure=c("tanimoto",
#' "tanimoto"),nrclusters=seq(5,10),normalize=c(FALSE,FALSE),method=c(NULL,NULL),
#' names=c("FP","TP"),StopRange=FALSE,plottype="new",location=NULL)
#' 
#' NrClusters
#' }
#' @export SelectnrClusters
SelectnrClusters<-function(List,type=c("data","dist","pam"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),nrclusters = seq(5, 25, 1),names=NULL,StopRange=FALSE,plottype="new",location=NULL){
	
	type=match.arg(type)
	avsilwidth<-matrix(0,ncol=length(List),nrow=length(nrclusters))
	pamfunction<-function(DistM,nrclusters){
		asw=sapply(nrclusters, function(x) cluster::pam(DistM,x)$silinfo$avg.width)
		return(asw)
	}
	
	CheckDist<-function(Dist,StopRange){
		if(StopRange==FALSE & !(0<=min(Dist) & max(Dist)<=1)){
			message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
			Dist=Normalization(Dist,method="Range")
		}
		else{
			Dist=Dist
		}
	}
	
	
	if(type=="data"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}
		Dist=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i],normalize[i],method[i]))
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		
		avsilwidth=sapply(Dist,function(x) pamfunction(x,nrclusters=nrclusters))
		rownames(avsilwidth)=nrclusters
	}
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		Dist=List
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		
		avsilwidth=sapply(Dist,function(x) pamfunction(x,stats::start,stats::end))
		rownames(avsilwidth)=nrclusters
	}
	else{
		avsilwidth=sapply(List,function(x) return(x$silinfo$avg.width))
	}
	
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	
	
	if(is.null(names)){
		names1=c()
		names2=c()
		for(i in 1:length(List)){
			names1=c(names1,paste("Silhouette widths for Data",i,sep=" "))
			names2=c(names2,paste("Nr Clusters for Data",i,sep=' '))
		}
		names1=c(names1,"Average Silhoutte Widths")
		names2=c(names2,"Optimal nr of clusters")
		
	}
	else{
		names1=c()
		names2=c()
		for(i in 1:length(List)){
			names1=c(names1,paste("Silhouette widths for",names[i],sep=" "))
			names2=c(names2,paste("Nr Clusters for",names[i],sep=' '))
		}
		names1=c(names1,"Average Silhoutte Widths")
		names2=c(names2,"Optimal nr of clusters")
	}
	
	
	rownames(avsilwidth)=nrclusters
	
	
	avsil=apply(avsilwidth,1,mean)
	avsilwidth=cbind(avsilwidth,avsil)
	colnames(avsilwidth)=names1
	
	
	plotsil<-function(sils,plottype,location,name){
		k.best=as.numeric(names(sils)[which.max(sils)])
		cat("silhouette-optimal number of clusters:", k.best, "\n")
		plottypein(plottype,location)
		graphics::plot(nrclusters, sils, type= "h", main = name,
				xlab= "k  (# clusters)", ylab = "average silhouette width")
		graphics::axis(1, k.best, paste("best",k.best,sep="\n"), col = "red", col.axis = "red")
		plottypeout(plottype)
	}
	
	sapply(c(1:ncol(avsilwidth)),function(x) plotsil(avsilwidth[,x],plottype,location,names1[x]))
	
	Output=list()
	Output[[1]]=avsilwidth
	nrclusters=apply(avsilwidth,2,function(x) return(as.numeric(names(x)[which.max(x)])))
	nrclusters=as.data.frame(t(nrclusters))
	colnames(nrclusters)=names2
	rownames(nrclusters)="NrClusters"
	
	
	Output[[2]]=nrclusters
	
	names(Output)=c("Silhoutte_Widths","Optimal_Nr_of_CLusters")
	return(Output)
}

#' @title Intersection of clusters across multiple methods
#' 
#' @description The \code{SharedComps} function is an easy way to select the objects that
#' are shared over clusters of different methods.
#' @param List A list of clustering outputs or the output of the
#' \code{DiffGenes} function. The first element of the list will be used as a
#' reference in \code{ReorderToReference}.
#' @param nrclusters If List is the output several clustering methods, it has
#' to be provided in how many clusters to cut the dendrograms in. Default is NULL.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is FALSE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is FALSE.
#' @param names Names of the methods or clusters. Default is NULL.
#' @return A vector containing the shared objects of all listed elements.
#' @examples
#' 
#' data(fingerprintMat)
#' data(targetMat)
#' data(geneMat)
#' data(GeneInfo)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' L=list(MCF7_F,MCF7_T)
#' names=c('FP','TP')
#' 
#' Comps=SharedComps(List=L,nrclusters=7,fusionsLog=FALSE,weightclust=FALSE,names=names)
#' 
#' 
#' @export SharedComps
SharedComps<-function(List,nrclusters=NULL,fusionsLog=FALSE,weightclust=FALSE,names=NULL){
	FoundComps=NULL
	
	FoundComps=FindElement("objects",List)
	
	if(is.null(FoundComps)|(is.list(FoundComps) & length(FoundComps) == 0)){
		ListNew=list()
		element=0
		for(i in 1:length(List)){
			if(attributes(List[[i]])$method != "CEC" & attributes(List[[i]])$method != "Weighted" & attributes(List[[i]])$method!= "WeightedSim"){
				ResultsClust=list()
				ResultsClust[[1]]=list()
				ResultsClust[[1]][[1]]=List[[i]]
				names(ResultsClust[[1]])[1]="Clust"
				element=element+1					
				ListNew[[element]]=ResultsClust[[1]]
				#attr(ListNew[element],"method")="Weights"
			}
			else if(attributes(List[[i]])$method=="CEC" | attributes(List[[i]])$method=="Weighted" | attributes(List[[i]])$method == "WeightedSim"){
				ResultsClust=list()
				if(weightclust==TRUE){
					ResultsClust[[1]]=list()
					if(attributes(List[[i]])$method != "WeightedSim"){
						ResultsClust[[1]][[1]]=List[[i]]$Clust
						names(ResultsClust[[1]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[1]]
						attr(ListNew[element],"method")="Weights"
					}
					else{
						ResultsClust[[1]]=list()
						ResultsClust[[1]][[1]]=List[[i]]
						names(ResultsClust[[1]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[1]]
					}
				}
				else{
					for (j in 1:length(List[[i]]$Results)){
						ResultsClust[[j]]=list()
						ResultsClust[[j]][[1]]=List[[i]]$Results[[j]]
						names(ResultsClust[[j]])[1]="Clust"
						element=element+1					
						ListNew[[element]]=ResultsClust[[j]]
						attr(ListNew[element],"method")="Weights"
					}		
				}		
			}	
		}
		
		if(is.null(names)){
			names=seq(1,length(ListNew),1)
			for(i in 1:length(ListNew)){
				names[i]=paste("Method",i,sep=" ")
			}
		}
		names(ListNew)=names
		
		
		MatrixClusters=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)
		List=ListNew
		Comps=list()
		for (k in 1:dim(MatrixClusters)[1]){
			clusters=MatrixClusters[k,]
			
			clust=sort(unique(clusters)) #does not matter: Genes[i] puts right elements on right places
			hc<-stats::as.hclust(List[[k]]$Clust$Clust)
			OrderedCpds <- hc$labels[hc$order]
			clusters=MatrixClusters[k,]
			Method=list()
			for (i in 1:nrclusters){
				
				temp=list()
				LeadCpds=names(clusters)[which(clusters==i)] 
				temp[[1]]=list(LeadCpds,OrderedCpds)
				names(temp[[1]])=c("LeadCpds","OrderedCpds")
				names(temp)="objects"
				Method[[i]]=temp
				names(Method)[i]=paste("Cluster",i,sep=" ")
			}
			Comps[[k]]=Method
			
		}
		names(Comps)=names
		List=Comps
	}
	
	for(a in 1:length(List)){
		if(a==1){
			Comps1=FindElement("LeadCpds",List[[a]])		
		}
		else if(is.list(Comps1) & length(Comps1) != 0){
			Comps2=FindElement("LeadCpds",List[[a]])
			for(b in 1:length(Comps1)){
				Comps1[[b]]=intersect(Comps1[[b]],Comps2[[b]])
			}			
		}
		
	}
	
	namesCl=c()
	for(c in 1:length(Comps1)){
		namesCl=c(namesCl,paste("Cluster",c,sep=" "))
	}
	
	names(Comps1)=namesCl
	
	
	return(Comps1)
}

#' @title Intersection of genes and pathways over multiple methods
#' 
#' @description It is interesting to investigate exactly which and how many differently
#' expressed genes, pathways and characteristics are shared by the clusters
#' over the different methods. The function \code{SharedGenesPathsFeat} will
#' provide this information. Given the outputs of the \code{DiffGenes}, the
#' \code{Geneset.intersect} function and/or \code{CharacteristicFeatures}, it
#' investigates how many genes, pathways and/or characteristics are expressed
#' by each cluster per method, how many of these are shared over the methods
#' and which ones are shared including their respective p-values of each method
#' and a mean p-value. This is very handy to look into the shared genes and
#' pathways of clusters that share many objects but also of those that only
#' share only a few. Further, the result also includes the number of objects
#' per cluster per method and how many of these are shared over the methods.
#' The input can also be focused for a specific selection of objects or a
#' specific cluster.
#' 
#' 
#' @param DataLimma Optional. The output of a \code{DiffGenes} function. Default is NULL.
#' @param DataMLP Optional. The output of \code{Geneset.intersect} function. Default is NULL.
#' @param DataFeat Optional. The output of \code{CharacteristicFeatures}
#' function. Default is NULL.
#' @param names Optional. Names of the methods or "Selection" if it only
#' considers a selection of objects. Default is NULL.
#' @param Selection Logical. Do the results concern only a selection of
#' objects or a specific cluster?  If yes, then Selection should be put to
#' TRUE. Otherwise all objects and clusters are considered.Default is FALSE.
#' @return The result of the \code{SharedGenesPathsFeat} function is a list
#' with two elements. The first element Table is a table indicating how many
#' genes, pathways and/or characteristics were found to be differentially
#' expressed and how many of these are shared. The table also contains the
#' number of objects shared between the clusters of the different methods.
#' The second element Which is another list with a component per cluster. Each
#' component consists of four vectors: SharedComps indicating which objects
#' were shared across the methods, SharedGenes represents the shared genes,
#' SharedPaths shows the shared pathways and SharedFeat the shared features.
#' @examples
#' 
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(geneMat)
#' data(GeneInfo)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' L=list(MCF7_F,MCF7_T)
#' names=c('FP','TP')
#' 
#' MCF7_Paths_FandT=PathwaysIter(List=L, geneExpr=geneMat, nrclusters=7, method=
#' c("limma", "MLP"), geneInfo=GeneInfo, geneSetSource="GOBP", topP=NULL,
#' topG=NULL, GENESET=NULL, sign=0.05,niter=2,fusionsLog=TRUE, 
#' weightclust=TRUE, names =names)
#'
#' MCF7_Paths_intersection=Geneset.intersect(MCF7_Paths_FandT,0.05,names=names,
#' seperatetables=FALSE,separatepvals=FALSE)
#'
#' MCF7_DiffGenes_FandT10=DiffGenes(list(MCF7_F,MCF7_T),Selection=NULL,geneExpr=geneMat,
#' nrclusters=7,method="limma",sign=0.05,top=10,fusionsLog=TRUE,weightclust=TRUE,names=NULL)
#'
#' MCF7_Char=CharacteristicFeatures(list(MCF7_F,MCF7_T),Selection=NULL,binData=
#' list(fingerprintMat,targetMat),datanames=c("FP","TP"),nrclusters=7,top=NULL,
#' sign=0.05,fusionsLog=TRUE,weightclust=TRUE,names=c("FP","TP")) 
#'
#' MCF7_Shared=SharedGenesPathsFeat(DataLimma=MCF7_DiffGenes_FandT10,
#' DataMLP=MCF7_Paths_intersection,DataFeat=MCF7_Char)
#' 
#' str(MCF7_Shared)
#' }
#' 
#' @export SharedGenesPathsFeat
SharedGenesPathsFeat<-function(DataLimma=NULL,DataMLP=NULL,DataFeat=NULL,names=NULL,Selection=FALSE){  #Input=result of DiffGenes and Geneset.intersect
	#Include sharedLimma and SharedMLP inside the function
	if(Selection==TRUE){
		ResultShared=SharedSelection(DataLimma,DataMLP,DataFeat,names)
	}	
	
	else if(is.null(DataLimma) & is.null(DataMLP) & is.null(DataFeat)){	
		stop("At least one Data set should be specified")
	}
	
	else{
		
		List=list(DataLimma,DataMLP,DataFeat)
		AvailableData=sapply(seq(length(List)),function(i) if(!(is.null(List[[i]]))) return(i))
		AvailableData=unlist(AvailableData)
		len=c()
		for(i in AvailableData){
			len=c(len,length(List[[i]]))
		}
		if(length(unique(len))!=1){
			stop("Unequal number of methods for limma and MLP")
		}
		else{
			DataSets=lapply(AvailableData,function(i)  return(List[[i]]))
			nmethods=length(DataSets[[1]])
			nclusters=length(DataSets[[1]][[1]])
		}
		
		if(is.null(names)){
			for(j in 1:length(DataSets[[1]])){
				names[j]=paste("Method",j,sep=" ")	
			}
		}
		
		which=list()	
		table=c()
		
		
		for (i in 1:nclusters){
			
			name=paste("Cluster",i,sep=" ")
			
			comps=c()
			
			temp1g=c()
			temp1p=c()
			temp1f=list()
			
			
			for(j in 1:nmethods){
				
				if(!(is.na(DataSets[[1]][[j]][[i]])[1])){
					comps=c(comps,length(DataSets[[1]][[j]][[i]]$objects$LeadCpds))
					names(comps)[j]=paste("Ncomps", names[j],sep=" ")
				}
				else{
					comps=c(comps,"-")
				}
				
				if(!(is.null(DataLimma))){
					if(!(is.na(DataLimma[[j]][[i]])[1])){
						temp1g=c(temp1g,length(DataLimma[[j]][[i]]$Genes$TopDE$ID))
					}
					else{
						temp1g=c(temp1g,"-")
					}	
					names(temp1g)[j]=names[j]
				}
				else{
					temp1g=NULL
				}
				
				if(!(is.null(DataMLP))){
					if(!(is.na(DataMLP[[j]][[i]])[1])){
						temp1p=c(temp1p,length(DataMLP[[j]][[i]][[3]]$geneSetDescription))
					}
					else{
						temp1p=c(temp1g,"-")
					}	
					names(temp1p)[j]=names[j]
				}
				else{
					temp1p=NULL
				}
				
				if(!(is.null(DataFeat))){
					temp=c()
					for(f in 1:length(DataFeat[[j]][[i]]$Characteristics)){			
						if(!(is.na(DataFeat[[j]][[i]])[1])){
							temp=c(temp,length(DataFeat[[j]][[i]]$Characteristics[[f]]$TopFeat$Names))
						}
						else{
							temp=c(temp,"-")
						}
						
						names(temp)[f]=names(DataFeat[[j]][[i]]$Characteristics)[f]
					}	
					temp1f[[j]]=temp
					names(temp1f)[j]=names[j]	
				}
				else{
					temp1f=NULL
				}
			}	
			
			j=1		
			Continue=TRUE
			while (Continue==TRUE){
				cont=c()
				for(d in 1:length(DataSets)){
					cont=c(cont,!(is.na(DataSets[[d]][[j]][[i]])[1]))
				}	
				
				if(any(cont)){
					
					sharedcomps=DataSets[[1]][[j]][[i]]$objects$LeadCpds
					nsharedcomps=length(sharedcomps)
					names(nsharedcomps)="nsharedcomps"
					
					
					if(!(is.null(DataLimma))){
						sharedgenes=DataLimma[[j]][[i]]$Genes$TopDE$ID
						nsharedgenes=length(sharedgenes)
						names(nsharedgenes)="Nshared"
						#pvalsg=DataLimma[[j]][[i]]$Genes$TopDE$adj.P.Val
					}
					else{
						sharedgenes=NULL
						nsharedgenes=NULL
					}
					
					
					if(!(is.null(DataMLP))){
						sharedpaths=DataMLP[[j]][[i]][[3]]$geneSetDescription
						nsharedpaths=length(sharedpaths)
						names(nsharedpaths)="Nshared"
						#pvalsp=DataMLP[[j]][[i]][[3]]$mean_geneSetPValue
						
					}
					else{
						sharedpaths=NULL
						nsharedpaths=NULL
					}
					if(!(is.null(DataFeat))){
						sharedfeat=list()
						nsharedfeat=list()
						#pvalsf=list()
						for(f in 1:length(DataFeat[[j]][[i]]$Characteristics)){	
							sharedfeat[[f]]=DataFeat[[j]][[i]]$Characteristics[[f]]$TopFeat$Names	
							names(sharedfeat)[f]=paste("shared: ",names(DataFeat[[j]][[i]]$Characteristics)[f],sep="")
							
							nsharedfeat[[f]]=length(sharedfeat[[f]])
							names(nsharedfeat)[f]=paste("Nshared: ",names(DataFeat[[j]][[i]]$Characteristics)[f],sep="")
							
							#pvalsf[[f]]=DataFeat[[j]][[i]]$Characteristics[[f]]$TopFeat$adj.P.Val	
							#names(pvalsf)[f]=paste("shared: ",names(DataFeat[[j]][[i]]$Characteristics)[f],sep="")
							
						}
						
					}
					else{
						sharedfeat=NULL
						nsharedfeat=NULL
					}
					
					Continue=FALSE
				}
				j=j+1
			}
			
			if(nmethods>=2){
				for (j in 2:nmethods){
					cont=c()
					for(d in 1:length(DataSets)){
						cont=c(cont,!(is.na(DataSets[[d]][[j]][[i]])[1]))
					}	
					if(any(cont)){
						
						sharedcomps=intersect(sharedcomps,DataSets[[1]][[j]][[i]]$objects$LeadCpds)
						nsharedcomps=length(sharedcomps)
						names(nsharedcomps)="Nsharedcomps"
						
						
						if(!(is.null(DataLimma))){
							sharedgenes=intersect(sharedgenes,DataLimma[[j]][[i]]$Genes$TopDE$ID)
							nsharedgenes=length(sharedgenes)
							names(nsharedgenes)="Nshared"
							
						}
						if(!(is.null(DataMLP))){
							sharedpaths=intersect(sharedpaths,DataMLP[[j]][[i]][[3]]$geneSetDescription)
							nsharedpaths=length(sharedpaths)
							names(nsharedpaths)="Nshared"
							
						}
						if(!(is.null(DataFeat))){
							
							for(f in 1:length(DataFeat[[j]][[i]]$Characteristics)){	
								sharedfeat[[f]]=intersect(sharedfeat[[f]],DataFeat[[j]][[i]]$Characteristics[[f]]$TopFeat$Names)	
								names(sharedfeat)[f]=paste("shared: ",names(DataFeat[[j]][[i]]$Characteristics)[f],sep="")
								
								nsharedfeat[[f]]=length(sharedfeat[[f]])
								names(nsharedfeat)[f]=paste("Nshared: ",names(DataFeat[[j]][[i]]$Characteristics)[f],sep="")
								
							}
							
						}
						
					}
				}			
			}
			
			pvalsgenes=list()
			meanpvalsgenes=c()
			
			pvalspaths=list()
			meanpvalspaths=c()
			
			pvalsfeat=list()
			meanpvalsfeat=c()
			
			
			if(!(is.null(sharedgenes))&length(sharedgenes)!=0){
				for(c in 1:nmethods){
					pvalsg=c()
					for(g in sharedgenes){
						if(!(is.na(DataLimma[[c]][[i]])[1])){
							pvalsg=c(pvalsg,DataLimma[[c]][[i]]$Genes$TopDE$adj.P.Val[DataLimma[[c]][[i]]$Genes$TopDE$ID==g])	
						}	
					}
					
					pvalsgenes[[c]]=pvalsg
					names(pvalsgenes)[c]=paste("P.Val.",names[c],sep="")
				}	
				
				for(g1 in 1:length(sharedgenes)){
					pvalstemp=c()			
					for(c in 1:nmethods){
						if(!(is.na(DataLimma[[c]][[i]])[1])){
							pvalstemp=c(pvalstemp,pvalsgenes[[c]][[g1]])
						}
					}			
					meanpvalsgenes=c(meanpvalsgenes,mean(pvalstemp))			
				}
				pvalsgenes[[nmethods+1]]=meanpvalsgenes	
				names(pvalsgenes)[nmethods+1]="Mean pvals genes"
			}
			else{pvalsgenes=NULL}
			
			if(!(is.null(sharedpaths))&length(sharedpaths)!=0){
				for(c in 1:nmethods){
					pvalsp=c()
					if(!(is.na(DataMLP[[c]][[i]])[1])){
						for(p in sharedpaths){
							pvalsp=c(pvalsp,DataMLP[[c]][[i]][[3]][DataMLP[[c]][[i]][[3]]$geneSetDescriptions==p,5][1])
						}
					}
					
					pvalspaths[[c]]=pvalsp
					names(pvalspaths)[c]=paste("P.Val.",names[c],sep="")
				}
				
				
				for(p1 in 1:length(sharedpaths)){
					pvalstemp1=c()
					for(c in 1:nmethods){
						if(!(is.na(DataMLP[[c]][[i]])[1])){
							pvalstemp1=c(pvalstemp1,pvalspaths[[c]][[p1]])
							
						}
						
					}			
					
					meanpvalspaths=c(meanpvalspaths,mean(pvalstemp1))
				}
				pvalspaths[[nmethods+1]]=meanpvalspaths	
				names(pvalspaths)[nmethods+1]="Mean pvals paths"
			}
			else{pvalpaths=NULL}
			
			
			if(!(is.null(sharedfeat))){
				for(f in 1:length(DataFeat[[j]][[i]]$Characteristics)){
					
					if(length(sharedfeat[[f]])!=0){
						pvalschar=list()
						for(c in 1:nmethods){
							pvalsf=c()
							if(!(is.na(DataFeat[[c]][[i]])[1])){
								for(s in sharedfeat[[f]]){
									pvalsf=c(pvalsf,DataFeat[[c]][[i]]$Characteristics[[f]]$TopFeat$adj.P.Val[DataFeat[[c]][[i]]$Characteristics[[f]]$TopFeat$Names==s])
								}
							}
							pvalschar[[c]]=pvalsf
							names(pvalschar)[c]=paste("P.Val.",names[c],sep="")
						}
					}
					else{
						pvalschar=list()
						pvalschar[1:nmethods]=0
						for(c in 1:nmethods){
							names(pvalschar)[c]=paste("P.Val.",names[c],sep="")
						}
					}
					pvalsfeat[[f]]=pvalschar
					names(pvalsfeat)[f]=names(DataFeat[[j]][[i]]$Characteristics)[f]
					
				}
				
				for(f in 1:length(DataFeat[[j]][[i]]$Characteristics)){
					meanpvalsfeat=c()
					
					if((length(sharedfeat[[f]])!=0)){				
						for(f1 in 1:length(sharedfeat[[f]])) {
							
							pvalstemp=c()			
							for(c in 1:nmethods){
								if(!(is.na(DataFeat[[c]][[i]])[1])){
									pvalstemp=c(pvalstemp,pvalsfeat[[f]][[c]][[f1]])
								}
							}			
							meanpvalsfeat=c(meanpvalsfeat,mean(pvalstemp))			
						}
					}
					else{
						meanpvalsfeat=0
					}
					
					pvalsfeat[[f]][[nmethods+1]]=meanpvalsfeat
					names(pvalsfeat[[f]])[nmethods+1]="Mean pvals feat"
				}
				lenchar=length(DataFeat[[j]][[i]]$Characteristics)
				
			}
			else{
				pvalsfeat=NULL
				lenchar=0
			}
			
			
			if(!(is.null(temp1f))){
				temp1f=do.call(rbind.data.frame, temp1f)	
				colnames(temp1f)=names(DataFeat[[j]][[i]]$Characteristics)
				temp1f=as.matrix(temp1f)
				nsharedfeat=do.call(cbind.data.frame, nsharedfeat)
				nsharedfeat=as.matrix(nsharedfeat)
			}
			part1=cbind(cbind(temp1g,temp1p),temp1f)
			part1=as.matrix(part1)
			if(is.null(nsharedgenes) & is.null(nsharedpaths) &is.null(nsharedfeat)){
				if(!(is.null(temp1g))){
					nsharedgenes=0
				}
				if(!(is.null(temp1p))){
					nsharedpaths=0
				}
				if(!(is.null(temp1f))){
					nsharedfeat=rep(0,length(temp1f))
				}
			}
			part2=cbind(cbind(nsharedgenes,nsharedpaths),nsharedfeat)
			part2=as.matrix(part2)
			colnames(part1)=NULL
			colnames(part2)=NULL
			rownames(part2)="NShared"
			part3=c()
			for(r in 1:length(comps)){
				part3=rbind(part3,rep(comps[r],dim(part1)[2]))
			}
			colnames(part3)=NULL
			rownames(part3)=names(comps)
			part4=rep(nsharedcomps,dim(part1)[2])
			names(part4)=NULL
			temp=rbind(part1,part2,part3,part4)	
			rownames(temp)[nrow(temp)]="Nsharedcomps"
			
			
			table=cbind(table,temp)
			colnames(table)=seq(1,ncol(table))
			
			if(!(is.null(pvalsgenes))){
				SharedGenes=cbind(sharedgenes,do.call(cbind.data.frame, pvalsgenes))
			}
			else{
				SharedGenes=NULL
			}
			if(!(is.null(pvalspaths))){
				SharedPaths=cbind(sharedpaths,do.call(cbind.data.frame, pvalspaths))
			}
			else{
				SharedPaths=NULL
			}
			if(!(is.null(pvalsfeat))){
				SharedFeat=list()
				for(f in 1:lenchar){
					if(length(sharedfeat[[f]])==0){
						SharedFeat[[f]]=NULL
					}
					else{
						SharedFeat[[f]]=cbind(sharedfeat[[f]],do.call(cbind.data.frame, pvalsfeat[[f]]))
						names(SharedFeat)[f]=names(pvalsfeat[[1]])[f]
					}
					
				}
			}
			else{
				SharedFeat=NULL
			}
			which[[i]]=list(SharedComps=sharedcomps,SharedGenes=SharedGenes,SharedPaths=SharedPaths,SharedFeat=SharedFeat)
			names(which)[i]=paste("Cluster",i,sep=" ")
			
		}
		#Sep for all situations?
		if(all(!(is.null(DataLimma)),!(is.null(DataMLP)),!(is.null(DataFeat)))){
			for (i in 1:length(seq(1,dim(table)[2],lenchar+2))){
				number=seq(1,dim(table)[2],2+lenchar)[i]
				colnames(table)[number]=paste("G.Cluster",i,sep=" ")	
				colnames(table)[number+1]=paste("P.Cluster",i,sep=" ")	
				for(u in seq(1:lenchar)){
					colnames(table)[number+1+u]=paste(paste('Feat.',names(DataFeat[[1]][[1]]$Characteristics)[u],sep=""),paste(".Cluster",i,sep=" "),sep="")
				}
				
			}
		}	
		
		else if(all(!(is.null(DataLimma)),!(is.null(DataMLP)),is.null(DataFeat))){
			for (i in 1:length(seq(1,dim(table)[2],2))){
				number=seq(1,dim(table)[2],2)[i]
				colnames(table)[number]=paste("G.Cluster",i,sep=" ")	
				colnames(table)[number+1]=paste("P.Cluster",i,sep=" ")	
				
			}
		}	
		
		else if(all(!(is.null(DataLimma)),is.null(DataMLP),!(is.null(DataFeat)))){
			for (i in 1:length(seq(1,dim(table)[2],lenchar+1))){
				number=seq(1,dim(table)[2],1+lenchar)[i]
				colnames(table)[number]=paste("G.Cluster",i,sep=" ")	
				for(u in seq(1:lenchar)){
					colnames(table)[number+u]=paste(paste('Feat.',names(DataFeat[[1]][[1]]$Characteristics)[u],sep=""),paste(".Cluster",i,sep=" "),sep="")
				}
				
			}
		}
		
		else if(all((is.null(DataLimma)),!(is.null(DataMLP)),!(is.null(DataFeat)))){
			for (i in 1:length(seq(1,dim(table)[2],lenchar+1))){
				number=seq(1,dim(table)[2],lenchar+1)[i]	
				colnames(table)[number]=paste("P.Cluster",i,sep=" ")	
				for(u in seq(1:lenchar)){
					colnames(table)[number+u]=paste(paste('Feat.',names(DataFeat[[1]][[1]]$Characteristics)[u],sep=""),paste(".Cluster",i,sep=" "),sep="")
				}
				
			}
		}	
		
		else if(all(!(is.null(DataLimma)),(is.null(DataMLP)),(is.null(DataFeat)))){
			for (i in 1:length(seq(1,dim(table)[2],1))){
				colnames(table)[i]=paste("G.Cluster",i,sep=" ")				
			}
		}	
		
		else if(all((is.null(DataLimma)),!(is.null(DataMLP)),(is.null(DataFeat)))){
			for (i in 1:length(seq(1,dim(table)[2],1))){
				colnames(table)[i]=paste("P.Cluster",i,sep=" ")				
			}
		}	
		
		else if(all((is.null(DataLimma)),(is.null(DataMLP)),!(is.null(DataFeat)))){
			for (i in 1:length(seq(1,dim(table)[2],lenchar))){
				number=seq(1,dim(table)[2],2)[i]
				for(u in c(0,1)){
					colnames(table)[number+u]=paste(paste('Feat.',names(DataFeat[[1]][[1]]$Characteristics)[u+1],sep=""),paste(".Cluster",i,sep=" "),sep="")
				}
				
			}
		}	
		
		ResultShared=list(Table=table,Which=which)	
		
	}
	return(ResultShared)
	
}

#' @title Intersection of genes and pathways over multiple methods for a selection of objects.
#' @param DataLimma Optional. The output of a \code{DiffGenes} function. Default is NULL.
#' @param DataMLP Optional. The output of \code{Geneset.intersect} function. Default is NULL.
#' @param DataFeat Optional. The output of \code{CharacteristicFeatures}
#' function. Default is NULL.
#' @param names Optional. Names of the methods or "Selection" if it only
#' considers a selection of objects. Default is NULL.
#' @description Internal function of \code{SharedGenesPathsFeat}.
SharedSelection<-function(DataLimma=NULL,DataMLP=NULL,DataFeat=NULL,names=NULL){  #Input=result of DiffGenes.2 and Geneset.intersect
	if(is.null(DataLimma) & is.null(DataMLP) & is.null(DataFeat)){	
		stop("At least one Data set should be specified")
	}
	
	List=list(DataLimma,DataMLP,DataFeat)
	AvailableData=sapply(seq(length(List)),function(i) if(!(is.null(List[[i]]))) return(i))
	AvailableData=unlist(AvailableData)
	len=c()
	for(i in AvailableData){
		len=c(len,length(List[[i]]))
	}
	if(length(unique(len))!=1){
		stop("Unequal number of methods for limma and MLP")
	}
	else{
		DataSets=lapply(AvailableData,function(i)  return(List[[i]]))
		nmethods=length(DataSets[[1]])
		#nclusters=length(DataSets[[1]][[1]])
	}
	
	if(is.null(names)){
		for(j in 1:nmethods){
			names[j]=paste("Method",j,sep=" ")	
		}
	}
	
	which=list()	
	table=c()
	
	comps=c()
	
	temp1g=c()
	temp1p=c()
	temp1f=list()
	
	for (i in 1:nmethods){			
		
		if(!(is.na(DataSets[[1]][[i]])[1])){
			comps=c(comps,length(DataSets[[1]][[i]]$objects$LeadCpds))
			names(comps)[i]=paste("Ncomps", names[i],sep=" ")
		}
		else{
			comps=c(comps,"-")
		}
		
		if(!(is.null(DataLimma))){
			if(!(is.na(DataLimma[[i]])[1])){
				temp1g=c(temp1g,length(DataLimma[[i]]$Genes$TopDE$ID))
			}
			else{
				temp1g=c(temp1g,"-")
			}	
			names(temp1g)[i]=names[i]
		}
		else{
			temp1g=NULL
		}
		
		if(!(is.null(DataMLP))){
			if(!(is.na(DataMLP[[i]])[1])){
				temp1p=c(temp1p,length(DataMLP[[i]][[3]]$geneSetDescription))
			}
			else{
				temp1p=c(temp1g,"-")
			}	
			names(temp1p)[i]=names[i]
		}
		else{
			temp1p=NULL
		}
		
		if(!(is.null(DataFeat))){
			temp=c()
			for(f in 1:length(DataFeat[[i]]$Characteristics)){			
				if(!(is.na(DataFeat[[i]])[1])){
					temp=c(temp,length(DataFeat[[i]]$Characteristics[[f]]$TopFeat$Names))
				}
				else{
					temp=c(temp,"-")
				}
				
				names(temp)[f]=names(DataFeat[[i]]$Characteristics)[f]
			}	
			temp1f[[i]]=temp
			names(temp1f)[i]=names[i]	
		}
		else{
			temp1f=NULL
		}
		
		if (i==1){
			cont=c()
			for(d in 1:length(DataSets)){
				cont=c(cont,!(is.na(DataSets[[d]][[i]])[1]))
			}	
			
			if(any(cont)){
				
				sharedcomps=DataSets[[1]][[i]]$objects$LeadCpds
				nsharedcomps=length(sharedcomps)
				names(nsharedcomps)="Nsharedcomps"
				
				
				if(!(is.null(DataLimma))){
					sharedgenes=DataLimma[[i]]$Genes$TopDE$ID
					nsharedgenes=length(sharedgenes)
					names(nsharedgenes)="Nshared"
				}
				else{
					sharedgenes=NULL
					nsharedgenes=0
				}
				
				
				if(!(is.null(DataMLP))){
					sharedpaths=DataMLP[[i]][[3]]$geneSetDescription
					nsharedpaths=length(sharedpaths)
					names(nsharedpaths)="Nshared"
					
				}
				else{
					sharedpaths=NULL
					nsharedpaths=0
				}
				if(!(is.null(DataFeat))){
					sharedfeat=list()
					nsharedfeat=list()
					for(f in 1:length(DataFeat[[i]]$Characteristics)){	
						sharedfeat[[f]]=DataFeat[[i]]$Characteristics[[f]]$TopFeat$Names	
						names(sharedfeat)[f]=paste("Shared: ",names(DataFeat[[1]]$Characteristics)[f],sep="")
						
						nsharedfeat[[f]]=length(sharedfeat[[f]])
						names(nsharedfeat)[f]=paste("Nshared: ",names(DataFeat[[1]]$Characteristics)[f],sep="")
						
					}
					
				}
				else{
					sharedfeat=NULL
					nsharedfeat=0
				}
			}		
		}
		else{	
			
			sharedcomps=intersect(sharedcomps,DataSets[[1]][[i]]$objects$LeadCpds)
			nsharedcomps=length(sharedcomps)
			names(nsharedcomps)="Nsharedcomps"
			
			
			if(!(is.null(DataLimma))){
				sharedgenes=intersect(sharedgenes,DataLimma[[i]]$Genes$TopDE$ID)
				nsharedgenes=length(sharedgenes)
				names(nsharedgenes)="Nshared"
				
			}
			if(!(is.null(DataMLP))){
				sharedpaths=intersect(sharedpaths,DataMLP[[i]][[3]]$geneSetDescription)
				nsharedpaths=length(sharedpaths)
				names(nsharedpaths)="Nshared"
				
			}
			if(!(is.null(DataFeat))){
				
				for(f in 1:length(DataFeat[[i]]$Characteristics)){	
					sharedfeat[[f]]=intersect(sharedfeat[[f]],DataFeat[[i]]$Characteristics[[f]]$TopFeat$Names)	
					names(sharedfeat)[f]=paste("Shared: ",names(DataFeat[[1]]$Characteristics)[f],sep="")
					
					nsharedfeat[[f]]=length(sharedfeat[[f]])
					names(nsharedfeat)[f]=paste("Nshared: ",names(DataFeat[[i]]$Characteristics)[f],sep="")
					
					
				}
				
			}
		}
	}
	
	pvalsgenes=list()
	meanpvalsgenes=c()
	
	pvalspaths=list()
	meanpvalspaths=c()
	
	pvalsfeat=list()
	meanpvalsfeat=c()
	
	
	if(!(is.null(sharedgenes)) & nsharedgenes != 0 ){
		for(c in 1:nmethods){
			pvalsg=c()
			for(g in sharedgenes){
				if(!(is.na(DataLimma[[c]])[1])){
					pvalsg=c(pvalsg,DataLimma[[c]]$Genes$TopDE$adj.P.Val[DataLimma[[c]]$Genes$TopDE$ID==g])	
				}	
			}
			
			pvalsgenes[[c]]=pvalsg
			names(pvalsgenes)[c]=paste("P.Val.",names[c],sep="")
		}	
		
		for(g1 in 1:length(sharedgenes)){
			pvalstemp=c()			
			for(c in 1:nmethods){
				if(!(is.na(DataLimma[[c]])[1])){
					pvalstemp=c(pvalstemp,pvalsgenes[[c]][[g1]])
				}
			}			
			meanpvalsgenes=c(meanpvalsgenes,mean(pvalstemp))			
		}
		pvalsgenes[[nmethods+1]]=meanpvalsgenes	
		names(pvalsgenes)[nmethods+1]="Mean pvals genes"
	}
	else{
		pvalsgenes=NULL
		nsharedgenes=NULL
	}
	
	if(!(is.null(sharedpaths)) & nsharedpaths != 0){
		for(c in 1:nmethods){
			pvalsp=c()
			if(!(is.na(DataMLP[[c]])[1])){
				for(p in sharedpaths){
					pvalsp=c(pvalsp,DataMLP[[c]][[3]][DataMLP[[c]][[3]]$geneSetDescription==p,5][1])
				}
			}
			
			pvalspaths[[c]]=pvalsp
			names(pvalspaths)[c]=paste("P.Val.",names[c],sep="")
		}
		
		
		for(p1 in 1:length(sharedpaths)){
			pvalstemp1=c()
			for(c in 1:nmethods){
				if(!(is.na(DataMLP[[c]])[1])){
					pvalstemp1=c(pvalstemp1,pvalspaths[[c]][[p1]])
					
				}
				
			}			
			
			meanpvalspaths=c(meanpvalspaths,mean(pvalstemp1))
		}
		pvalspaths[[nmethods+1]]=meanpvalspaths	
		names(pvalspaths)[nmethods+1]="Mean pvals paths"
	}
	else{
		pvalpaths=NULL
		nsharedpaths=NULL
	}
	
	if(!(is.null(sharedfeat)) & !(is.null(Reduce("+",nsharedfeat)))){
		for(f in 1:length(DataFeat[[1]]$Characteristics)){
			pvalschar=list()
			for(c in 1:nmethods){
				pvalsf=c()
				if(!(is.na(DataFeat[[c]])[1])){
					for(s in sharedfeat[[f]]){
						pvalsf=c(pvalsf,DataFeat[[c]]$Characteristics[[f]]$TopFeat$adj.P.Val[DataFeat[[c]]$Characteristics[[f]]$TopFeat$Names==s])
					}
				}
				
				pvalschar[[c]]=pvalsf
				names(pvalschar)[c]=paste("P.Val.",names[c],sep="")				
			}
			pvalsfeat[[f]]=pvalschar
			names(pvalsfeat)[f]=names(DataFeat[[1]]$Characteristics)[f]
			
		}
		
		for(f in 1:length(DataFeat[[1]]$Characteristics)){
			meanpvalsfeat=c()
			for(f1 in 1:length(sharedfeat[[f]])){
				pvalstemp=c()			
				for(c in 1:nmethods){
					if(!(is.na(DataFeat[[c]])[1]) & pvalsfeat[[f]]){
						pvalstemp=c(pvalstemp,pvalsfeat[[f]][[c]][[f1]])
					}
				}			
				meanpvalsfeat=c(meanpvalsfeat,mean(pvalstemp))			
			}
			
			pvalsfeat[[f]][[nmethods+1]]=meanpvalsfeat
			names(pvalsfeat[[f]])[nmethods+1]="Mean pvals feat"
		}
		lenchar=length(DataFeat[[1]]$Characteristics)
		
	}
	else{
		if(is.null(sharedfeat)){
			pvalsfeat=NULL
			nsharedfeat=NULL
			lenchar=0
		}
		else if(!(is.null(sharedfeat))){
			pvalsfeat=NULL
			lenchar=lenchar=length(DataFeat[[1]]$Characteristics)
		}
	}
	
	if(!(is.null(temp1f))){
		temp1f=do.call(rbind.data.frame, temp1f)	
		colnames(temp1f)=names(DataFeat[[1]]$Characteristics)
		temp1f=as.matrix(temp1f)
		nsharedfeat=do.call(cbind.data.frame, nsharedfeat)
		nsharedfeat=as.matrix(nsharedfeat)
	}
	part1=cbind(cbind(temp1g,temp1p),temp1f)
	part1=as.matrix(part1)
	if(is.null(nsharedgenes) | is.null(nsharedpaths)  | is.null(nsharedfeat)){
		if(!(is.null(temp1g))){
			nsharedgenes=0
		}
		if(!(is.null(temp1p))){
			nsharedpaths=0
		}
		if(!(is.null(temp1f))){
			nsharedfeat=rep(0,length(temp1f))
		}
	}
	part2=cbind(cbind(nsharedgenes,nsharedpaths),nsharedfeat)
	part2=as.matrix(part2)
	rownames(part2)="NShared"
	#print(str(part2))
	colnames(part1)=NULL
	colnames(part2)=NULL
	
	part3=c()
	for(r in 1:length(comps)){
		part3=rbind(part3,rep(comps[r],dim(part1)[2]))
	}
	colnames(part3)=NULL
	rownames(part3)=names(comps)
	part4=rep(nsharedcomps,dim(part1)[2])
	names(part4)=NULL
	temp=rbind(part1,part2,part3,part4)	
	rownames(temp)[nrow(temp)]="Nsharedcomps"
	
	
	table=cbind(table,temp)
	colnames(table)=seq(1,ncol(table))
	
	if(!(is.null(pvalsgenes))){
		SharedGenes=cbind(SharedGenes=sharedgenes,do.call(cbind.data.frame, pvalsgenes))
	}
	else{
		SharedGenes=NULL
	}
	if(!(is.null(pvalspaths))){
		SharedPaths=cbind(SharedPaths=sharedpaths,do.call(cbind.data.frame, pvalspaths))
	}
	else{
		SharedPaths=NULL
	}
	if(!(is.null(pvalsfeat))){
		SharedFeat=list()
		for(f in 1:lenchar){
			SharedFeat[[f]]=cbind(SharedFeat=sharedfeat[[f]],do.call(cbind.data.frame, pvalsfeat[[f]]))
			names(SharedFeat)[f]=names(pvalsfeat)[f]
		}
	}
	else{
		SharedFeat=NULL
	}
	
	which[[1]]=list(SharedComps=sharedcomps,SharedGenes=SharedGenes,SharedPaths=SharedPaths,SharedFeat=SharedFeat)
	names(which)[1]="Selection"
	
	#Sep for all situations?
	if(all(!(is.null(DataLimma)),!(is.null(DataMLP)),!(is.null(DataFeat)))){
		for (i in 1:length(seq(1,dim(table)[2],lenchar+2))){
			number=seq(1,dim(table)[2],2+lenchar)[i]
			colnames(table)[number]=paste("G.Cluster",sep=" ")	
			colnames(table)[number+1]=paste("P.Cluster",sep=" ")	
			for(u in seq(1:lenchar)){
				colnames(table)[number+1+u]=paste(paste('Feat.',names(DataFeat[[1]]$Characteristics)[u],sep=""),paste(".Cluster",sep=" "),sep="")
			}
			
		}
	}	
	
	else if(all(!(is.null(DataLimma)),!(is.null(DataMLP)),is.null(DataFeat))){
		for (i in 1:length(seq(1,dim(table)[2],2))){
			number=seq(1,dim(table)[2],2)[i]
			colnames(table)[number]=paste("G.Cluster",i,sep=" ")	
			colnames(table)[number+1]=paste("P.Cluster",i,sep=" ")	
			
		}
	}	
	
	else if(all(!(is.null(DataLimma)),is.null(DataMLP),!(is.null(DataFeat)))){
		for (i in 1:length(seq(1,dim(table)[2],lenchar+1))){
			number=seq(1,dim(table)[2],1+lenchar)[i]
			colnames(table)[number]=paste("G.Cluster",i,sep=" ")	
			for(u in seq(1:lenchar)){
				colnames(table)[number+u]=paste(paste('Feat.',names(DataFeat[[1]]$Characteristics)[u],sep=""),paste(".Cluster",i,sep=" "),sep="")
			}
			
		}
	}
	
	else if(all((is.null(DataLimma)),!(is.null(DataMLP)),!(is.null(DataFeat)))){
		for (i in 1:length(seq(1,dim(table)[2],lenchar+1))){
			number=seq(1,dim(table)[2],lenchar+1)[i]	
			colnames(table)[number]=paste("P.Cluster",i,sep=" ")	
			for(u in seq(1:lenchar)){
				colnames(table)[number+u]=paste(paste('Feat.',names(DataFeat[[1]]$Characteristics)[u],sep=""),paste(".Cluster",i,sep=" "),sep="")
			}
			
		}
	}	
	
	else if(all(!(is.null(DataLimma)),(is.null(DataMLP)),(is.null(DataFeat)))){
		for (i in 1:length(seq(1,dim(table)[2],1))){
			colnames(table)[i]=paste("G.Cluster",sep=" ")				
		}
	}	
	
	else if(all((is.null(DataLimma)),!(is.null(DataMLP)),(is.null(DataFeat)))){
		for (i in 1:length(seq(1,dim(table)[2],1))){
			colnames(table)[i]=paste("P.Cluster",sep=" ")				
		}
	}	
	
	else if(all((is.null(DataLimma)),(is.null(DataMLP)),!(is.null(DataFeat)))){
		for (i in 1:length(seq(1,dim(table)[2],lenchar))){
			number=seq(1,dim(table)[2],2)[i]
			for(u in c(0,1)){
				colnames(table)[number+u]=paste(paste('Feat.',names(DataFeat[[1]]$Characteristics)[u+1],sep=""),paste(".Cluster",sep=" "),sep="")
			}
			
		}
	}	
	
	
	
	ResultShared=list(Table=table,Which=which)
	return(ResultShared)
	
	
}

#' @title Intersection of genes over multiple methods for a selection of objects.
#' @param DataLimma Optional. The output of a \code{DiffGenes} function. Default is NULL.
#' @param names Optional. Names of the methods or "Selection" if it only
#' considers a selection of objects. Default is NULL.
#' @description Internal function of \code{SharedGenesPathsFeat}.
SharedSelectionLimma<-function(DataLimma=NULL,names=NULL){  #Input=result of DiffGenes.2 and Geneset.intersect
	
	which=list()	
	table=c()
	
	if(is.null(names)){
		for(j in 1:length(DataLimma)){
			names[j]=paste("Method",j,sep=" ")	
		}
	}
	
	nmethods=length(DataLimma)
	
	temp1g=c()
	comps=c()
	
	pvalsg=c()
	for (i in 1:nmethods){			
		
		temp1g=c(temp1g,length(DataLimma[[i]]$Genes$TopDE$Genes))
		comps=c(comps,length(DataLimma[[i]]$objects$LeadCpds))
		
		
		names(temp1g)[i]=names[i]
		names(comps)[i]=paste("Ncomps",names[i],i,sep=" ")
		
		
		
		
		if (i==1){
			if(!(is.na(DataLimma[[i]])[1])){
				sharedcomps=DataLimma[[i]]$objects$LeadCpds
				sharedgenes=DataLimma[[i]]$Genes$TopDE$Genes
				
				
				pvalsg=c(pvalsg,DataLimma[[i]]$Genes$TopDE$adj.P.Val)
				
				
				nsharedcomps=length(DataLimma[[i]]$objects$LeadCpds)
				nsharedgenes=length(DataLimma[[i]]$Genes$TopDE$Genes)
				
				names(nsharedgenes)="nshared"
				
				names(nsharedcomps)="nsharedcomps"
			}
			
		}
		else{			
			sharedcomps=intersect(sharedcomps,DataLimma[[i]]$objects$LeadCpds)
			sharedgenes=intersect(sharedgenes,DataLimma[[i]]$Genes$TopDE$Genes)
			
			
			nsharedcomps=length(intersect(sharedcomps,DataLimma[[i]]$objects$LeadCpds))
			nsharedgenes=length(intersect(sharedgenes,DataLimma[[i]]$Genes$TopDE$Genes))
			
			names(nsharedgenes)="nshared"
			
			names(nsharedcomps)="nsharedcomps"
		}
		
	}	
	pvalsgenes=list()
	meanpvalsgenes=c()
	
	if(nsharedgenes != 0){
		for(c in 1:nmethods){
			pvalsg=c()
			for(g in sharedgenes){
				if(!(is.na(DataLimma[[c]])[1])){
					pvalsg=c(pvalsg,DataLimma[[c]]$Genes$TopDE$adj.P.Val[DataLimma[[c]]$Genes$TopDE$Genes==g])	
				}	
			}
			
			pvalsgenes[[c]]=pvalsg
			names(pvalsgenes)[c]=paste("Method",c,sep=" ")
		}	
		
		for(g1 in 1:length(sharedgenes)){
			pvalstemp=c()			
			for(c in 1:nmethods){
				if(!(is.na(DataLimma[[c]])[1])){
					pvalstemp=c(pvalstemp,pvalsgenes[[c]][[g1]])
				}
			}			
			meanpvalsgenes=c(meanpvalsgenes,mean(pvalstemp))			
		}
		pvalsgenes[[nmethods+1]]=meanpvalsgenes	
		names(pvalsgenes)[nmethods+1]="Mean pvals genes"
	}
	else{pvalsgenes=0}
	
	
	temp=rbind(temp1g,nsharedgenes,comps,nsharedcomps)	
	
	table=cbind(table,temp)
	
	which[[1]]=list(sharedcomps=sharedcomps,sharedgenes=sharedgenes,pvalsgenes=pvalsgenes)
	#names(which)[1]=paste("Cluster",i,sep=" ")
	
	
	ResultShared=list(Table=table,Which=which)
	return(ResultShared)	
}

#' @title Intersection of pathways over multiple methods for a selection of objects.
#' @param DataMLP Optional. The output of \code{Geneset.intersect} function. Default is NULL.
#' @param names Optional. Names of the methods or "Selection" if it only
#' considers a selection of objects. Default is NULL.
#' @description Internal function of \code{SharedGenesPathsFeat}.
SharedSelectionMLP<-function(DataMLP=NULL,names=NULL){  #Input=result of DiffGenes.2 and Geneset.intersect
	
	which=list()	
	table=c()
	
	
	nmethods=length(DataMLP)
	
	if(is.null(names)){
		for(j in 1:length(DataMLP)){
			names[j]=paste("Method",j,sep=" ")	
		}
	}
	
	temp1g=c()
	temp1p=c()
	comps=c()
	
	pvalsg=c()
	pvalsp=c()	
	for (i in 1:nmethods){			
		
		temp1g=c(temp1g,length(DataMLP[[i]]$Genes$TopDE$Genes))
		temp1p=c(temp1p,length(DataMLP[[i]][[3]]$geneSetDescription))
		comps=c(comps,length(DataMLP[[i]]$objects$LeadCpds))
		
		
		names(temp1g)[i]=names[i]
		names(temp1p)[i]=names[i]
		names(comps)[i]=paste("Ncomps",names[i],i,sep=" ")
		
		
		if (i==1){
			if(!(is.na(DataMLP[[i]])[1]) | !(is.na(DataMLP[[i]])[1])){
				sharedcomps=DataMLP[[i]]$objects$LeadCpds
				sharedgenes=DataMLP[[i]]$Genes$TopDE$Genes
				sharedpaths=DataMLP[[i]][[3]]$geneSetDescription
				
				pvalsg=c(pvalsg,DataMLP[[i]]$Genes$TopDE$adj.P.Val)
				pvalsp=c(pvalsp,DataMLP[[i]]$mean_geneSetPValue)
				
				nsharedcomps=length(DataMLP[[i]]$objects$LeadCpds)
				nsharedgenes=length(DataMLP[[i]]$Genes$TopDE$Genes)
				nsharedpaths=length(DataMLP[[i]][[3]]$geneSetDescription)
				names(nsharedgenes)="nshared"
				names(nsharedpaths)="nshared"
				names(nsharedcomps)="nsharedcomps"
			}
			
		}
		else{			
			sharedcomps=intersect(sharedcomps,DataMLP[[i]]$objects$LeadCpds)
			sharedgenes=intersect(sharedgenes,DataMLP[[i]]$Genes$TopDE$Genes)
			sharedpaths=intersect(sharedpaths,DataMLP[[i]][[3]]$geneSetDescription)
			
			nsharedcomps=length(intersect(sharedcomps,DataMLP[[i]]$objects$LeadCpds))
			nsharedgenes=length(intersect(sharedgenes,DataMLP[[i]]$Genes$TopDE$Genes))
			nsharedpaths=length(intersect(sharedpaths,DataMLP[[i]][[3]]$geneSetDescription))
			names(nsharedgenes)="nshared"
			names(nsharedpaths)="nshared"
			names(nsharedcomps)="nsharedcomps"
		}
		
	}	
	pvalsgenes=list()
	meanpvalsgenes=c()
	meanpvalspaths=c()
	pvalspaths=list()
	
	if(nsharedgenes != 0){
		for(c in 1:nmethods){
			pvalsg=c()
			for(g in sharedgenes){
				if(!(is.na(DataMLP[[c]])[1])){
					pvalsg=c(pvalsg,DataMLP[[c]]$Genes$TopDE$adj.P.Val[DataMLP[[c]]$Genes$TopDE$Genes==g])	
				}	
			}
			
			pvalsgenes[[c]]=pvalsg
			names(pvalsgenes)[c]=paste("Method",c,sep=" ")
		}	
		
		for(g1 in 1:length(sharedgenes)){
			pvalstemp=c()			
			for(c in 1:nmethods){
				if(!(is.na(DataMLP[[c]])[1])){
					pvalstemp=c(pvalstemp,pvalsgenes[[c]][[g1]])
				}
			}			
			meanpvalsgenes=c(meanpvalsgenes,mean(pvalstemp))			
		}
		pvalsgenes[[nmethods+1]]=meanpvalsgenes	
		names(pvalsgenes)[nmethods+1]="Mean pvals genes"
	}
	else{pvalsgenes=0}
	
	if(nsharedpaths!=0){
		for(c in 1:nmethods){
			pvalsp=c()
			if(!(is.na(DataMLP[[c]])[1])){
				for(p in sharedpaths){
					pvalsp=c(pvalsp,DataMLP[[c]][[3]][DataMLP[[c]][[3]]$geneSetDescription==p,5])
				}
			}
			
			pvalspaths[[c]]=pvalsp
			
			
			names(pvalspaths)[c]=paste("Method",c,sep=" ")
		}
		
		
		for(p1 in 1:length(sharedpaths)){
			pvalstemp1=c()
			for(c in 1:nmethods){
				if(!(is.na(DataMLP[[c]])[1])){
					pvalstemp1=c(pvalstemp1,pvalspaths[[c]][[p1]])
					
				}
				
			}			
			
			meanpvalspaths=c(meanpvalspaths,mean(pvalstemp1))
			
		}
		pvalspaths[[nmethods+1]]=meanpvalspaths	
		names(pvalspaths)[nmethods+1]="Mean pvals paths"
	}
	else{pvalpaths=0}
	
	temp=rbind(cbind(temp1g,temp1p),cbind(nsharedgenes,nsharedpaths),cbind(comps,comps),cbind(nsharedcomps,nsharedcomps))	
	
	table=cbind(table,temp)
	
	which[[1]]=list(sharedcomps=sharedcomps,sharedgenes=sharedgenes,pvalsgenes=pvalsgenes,sharedpaths=sharedpaths,pvalspaths=pvalspaths)
	#names(which)[1]=paste("Cluster",i,sep=" ")
	
	
	for (i in 1:length(seq(1,dim(table)[2],2))){
		number=seq(1,dim(table)[2],2)[i]
		colnames(table)[number]=c("G.Cluster")
		colnames(table)[number+1]=c("P.Cluster")	
		
	}
	ResultShared=list(Table=table,Which=which)
	return(ResultShared)
	
	
}

#' @title A heatmap of similarity values between objects
#' 
#' @description The function \code{SimilarityHeatmap} plots the similarity values between
#' objects. The darker the shade, the more similar objects are. The option
#' is available to set a cutoff value to highlight the most similar objects.
#' 
#' @details If data is of type "clust", the distance matrix is extracted from the result
#' and transformed to a similarity matrix. Possibly a range normalization is
#' performed. If data is of type "dist", it is also transformed to a similarity
#' matrix and cluster is performed on the distances. If data is of type "sim",
#' the data is tranformed to a distance matrix on which clustering is
#' performed. Once the similarity mattrix is obtained, the cutoff value is
#' applied and a heatmap is drawn. If no cutoff value is desired, one can leave
#' the default NULL specification.
#' 
#' @param Data The data of which a heatmap should be drawn.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" ,"sim" or "clusters".
#' @param distmeasure The distance measure. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to "tanimoto".
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}. 
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is NULL.
#' @param linkage Choice of inter group dissimilarity (character). Defaults to "flexible".
#' @param cutoff Optional. If a cutoff value is specified, all values lower are
#' put to zero while all other values are kept. This helps to highlight the
#' most similar objects. Default is NULL.
#' @param percentile Logical. The cutoff value can be a percentile. If one want
#' the cutoff value to be the 90th percentile of the data, one should specify
#' cutoff = 0.90 and percentile = TRUE. Default is FALSE.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return A heatmap with the names of the objects on the right and bottom
#' and a dendrogram of the clustering at the left and top.
#' @examples
#' 
#' \dontrun{
#' data(fingerprintMat)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55)
#' 
#' SimilarityHeatmap(Data=MCF7_F,type="clust",cutoff=0.90,percentile=TRUE)
#' SimilarityHeatmap(Data=MCF7_F,type="clust",cutoff=0.75,percentile=FALSE)
#' 
#' }
#' 
#' @export SimilarityHeatmap
SimilarityHeatmap<-function(Data,type=c("data","clust","sim","dist"),distmeasure="tanimoto",normalize=FALSE,method=NULL,linkage="flexible",cutoff=NULL,percentile=FALSE,plottype="new",location=NULL){
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	
	if(type=="data"){
		ClustData<-Cluster(Data=Data,distmeasure=distmeasure,normalize=normalize,method=method,clust="agnes",linkage=linkage,gap=FALSE,maxK=55,StopRange=FALSE)
		Data=ClustData$DistM
		type="dist"
	}
	
	
	if(type=="clust"){
		Dist=Data$DistM
		if(0<=min(Dist) & max(Dist)<=1){
			SimData=1-Dist
		}
		else{
			NormData=Normalization(Dist,method="Range")
			SimData=1-NormData
		}
		
		ClustData=Data
	}
	
	else if(type=="dist"){
		if(0<=min(Data) & max(Data)<=1){
			SimData=1-Data
			ClustData=Cluster(Data=Data,type="dist",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
			
		}
		else{
			NormData=Normalization(Data,method="Range")
			SimData=1-NormData
			ClustData=Cluster(Data=Data,type="dist",distmeasure="tanimoto",normalize=TRUE,method="Q",clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
			
		}
		
		
	}
	else if(type=="sim"){
		SimData=Data
		if(0<=min(SimData) & max(SimData)<=1){
			DistData=1-Data
			ClustData=Cluster(Data=DistData,type="dist",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
			
		}
		else{
			NormData=Normalization(Dist,method="Range")
			DistData=1-Data
			ClustData=Cluster(Data=DistData,type="dist",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)		
		}
	}	
	
	
	if(!is.null(cutoff)){
		if(percentile==TRUE){
			cutoff=stats::quantile(SimData[lower.tri(SimData)], cutoff)
		}
		
		SimData_bin <- ifelse(SimData<=cutoff,0,SimData) # Every value higher than the 90ieth percentile is kept, all other are put to zero
	}
	
	else{
		SimData_bin=SimData
	}
	
	plottypein(plottype,location)
	gplots::heatmap.2(SimData_bin, 
			Rowv = stats::as.dendrogram(stats::as.hclust(ClustData$Clust)), Colv=stats::as.dendrogram(stats::as.hclust(ClustData$Clust)),trace="none",
			col=(grDevices::gray(seq(0.9,0,len=1000))),
			cexRow=0.6, cexCol=0.6,  
			margins=c(9,9),
			key=FALSE,
			keysize=0.4,
			symkey=FALSE,
			sepwidth=c(0.01,0.01),
			sepcolor="black",
			colsep=c(0,ncol(SimData_bin)),
			rowsep=c(0,nrow(SimData_bin))
	)
	plottypeout(plottype)
	
}

#' @title A measure of similarity for the outputs of the different methods
#' 
#' @description The function \code{SimilarityMeasure} computes the similarity of the
#' methods.  Given a list of outputs as input, the first element will be seen
#' as the reference.  Function \code{MatrixFunction} is called upon and the
#' cluster numbers are rearranged according to the reference. Per method,
#' \code{SimilarityMeasure} investigates which objects have the same cluster
#' number in reference and said method. This number is divided by the total
#' number of objects and used as a similarity measure.
#' 
#' @param List A list of clustering outputs to be compared. The first element
#' of the list will be used as the reference in \code{ReorderToReference}.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' 
#' @param names Optional. Names of the methods.
#' @return A vector of similarity measures, one for each method given as input.
#' @examples
#' 
#' data(fingerprintMat)
#' data(targetMat)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' L=list(MCF7_F,MCF7_T)
#' names=c("FP","TP")
#' 
#' MCF7_SimFandT=SimilarityMeasure(List=L,nrclusters=7,fusionsLog=TRUE,weightclust=TRUE,
#' names=names)
#' 
#' 
#' @export SimilarityMeasure
SimilarityMeasure<-function(List,nrclusters=NULL,fusionsLog=TRUE,weightclust=TRUE,names=NULL){
	
	if(class(List)!="list"){
		MatrixColors=List
	}	
	else{
		MatrixColors=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)
	}
	
	#Compare every row to the first row
	Similarity=c()
	for(i in 1:dim(MatrixColors)[1]){
		Shared=0
		for(j in 1:dim(MatrixColors)[2])
			if(MatrixColors[i,j]==MatrixColors[1,j]){
				Shared=Shared+1
			}
		Similarity=c(Similarity,Shared/ncol(MatrixColors))
		
		
	}
	
	return(Similarity)
}

#' @title Follow a cluster over multiple methods
#' 
#' @description It is often desired to track a specific selection of object over the
#' different methods and/or weights. This can be done with the
#' \code{ClusterDistribution}. For every method, it is tracked where the
#' objects of the selections are situated.
#' 
#' @details The result is provided with extra information as which objects of the
#' original selection can be found in this cluster and which are extra.
#' Further, plots of the distribution of the objects can be produced.  One
#' plot follows the complete distribution of the cluster while another one
#' focuses on either the maximum number of objects or a specific cluster,
#' whatever is specified. It are the number of objects that are plotted and
#' the first element indicated the number of objects in the selection. A
#' table can be produced as well, that separates the objects that are shared
#' over all methods from those extra in the original selection and extra for
#' the other methods. The \code{ReorderToReference} is applied to make sure
#' that the clusters are comparable over the methods.
#' 
#' The function is experimental and might not work in specific cases. Please
#' let us know such that we can improve its functionality.
#' 
#' @param List A list of the clustering outputs.The first element of the list
#' will be used as the reference in \code{ReorderToReference}.
#' @param Selection The selection of objects to follow or a specific cluster
#' number. Default is NULL.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param followMaxComps Logical for plot. Whether to follow the maximum of
#' objects. Default is FALSE.
#' @param followClust Logical for plot. Whether to follow the specific cluster. Default is TRUE.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' @param names Optional. Names of the methods. Default is NULL.
#' @param selectionPlot Logical. Should a plot be produced. Depending on
#' followMaxComps and followClust it focuses on the maximum of objects or a
#' cluster. It will not be indicated to which cluster objects moved. Default is FALSE.
#' @param table Logical. Should a table with the objects per method and the
#' shared objects be produced? Default is FALSE.
#' @param completeSelectionPlot Logical. Should the complete distribution of
#' the selection be plotted? This implies that it will be indicated to which
#' cluster objects will move. Default is FALSE.
#' @param ClusterPlot Logical. Plot of specific cluster. Default is FALSE.
#' @param cols The colors used for the different clusters. Default is NULL.
#' @param legendposx The x-coordinate of the legend on all plots. Default is 0.5.
#' @param legendposy The y-coordinate of the legend on all plots. Default is 2.4.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return The returned value is a list with an element for every method. This
#' element is another list with the following elements: \item{Selection}{The
#' selection of objects to follow} \item{nr.clusters}{the number of clusters
#' the selection is divided over} \item{nr.min.max.together }{the minimum and
#' maximum number of objects found together}
#' \item{perc.min.max.together}{minimum and maximum percentage of objects
#' found together} \item{AllClusters}{A list with an element per cluster that
#' contains at least one of the objects in Selection. The list contains the
#' cluster number, the complete cluster, the objects that originally could be
#' found in this cluster and which object were joined extra to it.} Depending
#' on whether followMaxComps or followClust is specified, the cluster of
#' interest is mentioned separately as well for easy access. If the option was
#' specified to create a table, this can be found under the Table element. Each
#' plot that was specified to be created is plotted in a new window in the
#' graphics console.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(Colors1)
#' 
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' 
#' L=list(MCF7_F,MCF7_T)
#' names=c("FP","TP")
#' 
#' Comps=FindCluster(List=L,nrclusters=7,select=c(1,4))
#' Comps
#' 
#' CompsFPAll=TrackCluster(List=L,Selection=Comps,nrclusters=7,followMaxComps=TRUE,
#' followClust=FALSE,fusionsLog=TRUE,weightclust=TRUE,names=names,selectionPlot=TRUE,
#' table=TRUE,completeSelectionPlot=TRUE,cols=Colors1,plottype="new",location=NULL)
#' }
#' 
#' @export TrackCluster
TrackCluster <- function(List,Selection,nrclusters=NULL,followMaxComps=FALSE,followClust=TRUE,fusionsLog=TRUE,weightclust=TRUE,names=NULL,selectionPlot=FALSE,table=FALSE,completeSelectionPlot=FALSE,ClusterPlot=FALSE,cols=NULL,legendposx=0.5,legendposy=2.4,plottype="sweave",location=NULL){
	
	ClusterDistribution.2<-function(List,Selection,nrclusters,followMaxComps,followClust,fusionsLog,weightclust,names){
		FoundClusters=list()	
		FoundCl=list()
		Matrix=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)
		
		ClusterNumber=NULL
		if(class(Selection)=="numeric"){
			ClusterNumber=Selection
			Selection=colnames(Matrix)[which(Matrix[1,]==Selection)]
		}
		
		
		if(followClust==TRUE){			
			cluster.interest=NULL
			m=1
			while(m<=dim(Matrix)[1] & is.null(cluster.interest)){
				clustersfound=unique(Matrix[m,which(colnames(Matrix)%in%Selection)])
				if(length(clustersfound)==1){
					cluster.interest=clustersfound
				}
				else{
					m=m+1
				}
			}
			if(is.null(cluster.interest)){
				message("This selection is not found to be part of a cluster. FollowMaxComps will be put to true and the function will proceed")
				followMaxComps=TRUE
				followClust=FALSE
			}
			
		}	
		
		for(i in 1:dim(Matrix)[1]){
			
			
			temp=list()
			temp[[1]]=Selection
			names(temp)[1]="Selection"
			
			clusternumbers=unique(Matrix[i,which(colnames(Matrix)%in%Selection)])
			
			nr.clusters=length(clusternumbers)
			temp[[2]]=nr.clusters
			names(temp)[2]="nr.clusters"
			
			min.together=min(table(Matrix[i,which(colnames(Matrix)%in%Selection)]))
			max.together=max(table(Matrix[i,which(colnames(Matrix)%in%Selection)]))
			
			nr.min.max.together=c(min.together,max.together)
			temp[[3]]=nr.min.max.together
			names(temp)[3]="nr.min.max.together"
			
			min.perc.together <- min.together/length(Selection) *100
			max.perc.together <- max.together/length(Selection) *100
			
			perc.min.max.together =c(min.perc.together,max.perc.together) 
			temp[[4]]=perc.min.max.together	
			names(temp)[4]="perc.min.max.together"
			
			temp[[5]]=list()
			names(temp)[5]="AllClusters"
			
			for(a in 1:length(clusternumbers)){
				temp[[5]][[a]]=list()
				names(temp[[5]])[a]=paste("Cluster",clusternumbers[a],sep=" ")
				
				temp[[5]][[a]][[1]]=clusternumbers[a]
				temp[[5]][[a]][[2]]=names(which(Matrix[i,]==clusternumbers[a])) #complete cluster
				temp[[5]][[a]][[3]]=intersect(Selection,temp[[5]][[a]][[2]]) #Objects from original selection in this cluster
				temp[[5]][[a]][[4]]=temp[[5]][[a]][[2]][which(!(temp[[5]][[a]][[2]] %in% Selection))] #Objects extra to this cluster
				names(temp[[5]][[a]])=c("clusternumber","Complete cluster","Objects from original selection in this cluster","Objects extra to this cluster")					
			}
			
			if(followMaxComps==TRUE){
				
				maxcluster=names(which(table(Matrix[i,which(colnames(Matrix)%in%Selection)])==max(table(Matrix[i,which(colnames(Matrix)%in%Selection)]))))
				temp[[6]]=maxcluster
				names(temp)[6]="Cluster with max Objects"
				complabels=rownames(Matrix[i,which(colnames(Matrix)%in%Selection & Matrix[i,]==as.numeric(maxcluster)),drop=FALSE])
				temp[[7]]=complabels
				names(temp)[7]="Complabels"
				complete.new.cluster=names(Matrix[i,which(Matrix[i,]==as.numeric(maxcluster))])
				temp[[8]]=complete.new.cluster
				names(temp)[8]="Complete.new.cluster"
				extra.new.cluster=complete.new.cluster[which(!(complete.new.cluster %in% Selection))]
				temp[[9]]=extra.new.cluster
				names(temp)[9]="Extra.new.cluster"
				
			}
			
			if(followClust==TRUE){
				
				temp[[6]]=cluster.interest
				names(temp)[6]="Cluster"
				complabels=rownames(Matrix[i,which(colnames(Matrix)%in%Selection & Matrix[i,]==as.numeric(cluster.interest)),drop=FALSE])
				temp[[7]]=complabels
				names(temp)[7]="Complabels"
				complete.new.cluster=names(Matrix[i,which(Matrix[i,]==as.numeric(cluster.interest))])
				temp[[8]]=complete.new.cluster
				names(temp)[8]="Complete.new.cluster"
				extra.new.cluster=complete.new.cluster[which(!(complete.new.cluster %in% Selection))]
				temp[[9]]=extra.new.cluster
				names(temp)[9]="Extra.new.cluster"
				
			}
			
			FoundClusters[[i]]=temp	
		}
		
		if(ClusterPlot==TRUE & !(is.null(ClusterNumber))){
			
			for(i in 1:dim(Matrix)[1]){
				temp=list()
				if(i==1){
					temp[[1]]=names(Matrix[i,which(Matrix[i,]==ClusterNumber)])
					names(temp)[1]=paste("Cluster ", ClusterNumber,sep="")
					PrevCluster=temp[[1]]
				}
				else{
					temp[[1]]=names(Matrix[i,which(Matrix[i,]==ClusterNumber)])
					names(temp)[1]=paste("Cluster ",  ClusterNumber,sep="")
					
					Diss=list()
					DissComps=NULL
					if(length((which(!(PrevCluster%in%temp[[1]]))!=0))){
						DissComps=PrevCluster[(which(!(PrevCluster%in%temp[[1]])))]
					}				
					if(!(is.null(DissComps))){ #Dissapeared comps: what is missing form PrevClust, where to did they move?
						cl=i
						for(t in 1:nrclusters){
							TempComps=which(DissComps%in%names(which(Matrix[cl,]==t)))
							disscl=NULL
							if(length(TempComps)!=0){
								disscl=list()
								disscl[[1]]=DissComps[TempComps]
								disscl[[2]]=t
								names(disscl)=c("Cpds","Cluster")
							}
							Diss[[t]]=disscl
						}
					}
					if(length(Diss)!=0){
						r=c()
						for(l in 1:length(Diss)){
							if(is.null(Diss[[l]])){
								r=c(r,l)
							}
						}
						if(!is.null(r)){
							Diss=Diss[-r]
						}
						for(l in 1:length(Diss)){
							names(Diss)[l]=paste("Comps Diss To Cluster " ,Diss[[l]][[2]],sep="")
							
						}
					}
					temp[[2]]=Diss
					names(temp)[2]="Dissapeared"
					
					Joined=list()
					JoinComps=NULL
					if(length((which(!(temp[[1]]%in%PrevCluster))!=0))){
						JoinComps=temp[[1]][(which(!(temp[[1]]%in%PrevCluster)))]
					}				
					if(!(is.null(JoinComps))){ #Dissapeared comps: what is missing form PrevClust, where to did they move?
						prevcl=i-1
						for(t in 1:nrclusters){
							TempComps=which(JoinComps%in%names(which(Matrix[prevcl,]==t)))
							joincl=NULL
							if(length(TempComps)!=0){
								joincl=list()
								joincl[[1]]=JoinComps[TempComps]
								joincl[[2]]=t
								names(joincl)=c("Cpds","Cluster")
								
							}
							Joined[[t]]=joincl
						}
					}
					if(length(Joined)!=0){
						r=c()
						for(l in 1:length(Joined)){
							if(is.null(Joined[[l]])){
								r=c(r,l)
							}
						}
						if(!is.null(r)){
							Joined=Joined[-r]
						}
						
						for(l in 1:length(Joined)){
							names(Joined)[l]=paste("Comps Joined From Cluster " ,Joined[[l]][[2]],sep="")
							
						}
					}
					temp[[3]]=Joined
					names(temp)[3]="Joined"	
					PrevCluster=temp[[1]]
				}
				FoundCl[[i]]=temp
			}
			
		}
		
		FoundClusters[[length(FoundClusters)+1]]=FoundCl
		
		FoundClusters[[length(FoundClusters)+1]]=Matrix
		
		return(FoundClusters)
	}	
	
	Found=ClusterDistribution.2(List,Selection,nrclusters,followMaxComps,followClust,fusionsLog,weightclust,names=names)
	
	Matrix=Found[[length(Found)]]
	Found=Found[-length(Found)]
	
	FoundCl=Found[[length(Found)]]
	Found=Found[-length(Found)]
	
	if(is.null(names)){
		for(j in 1:dim(Matrix)[1]){
			names[j]=paste("Method",j,1)
		}
	}
	
	ClusterNumber=NULL
	if(class(Selection)=="numeric"){
		ClusterNumber=Selection
		Selection=colnames(Matrix)[which(Matrix[1,]==Selection)]
		
	}
	
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	
	
	if(selectionPlot==TRUE){
		
		
		if(followMaxComps==TRUE){
			lab1=c("Maximum of objects of original cluster together")
			labelcluster=c()
			for(z in 1:length(Found)){
				labelcluster=c(labelcluster,Found[[z]]$Cluster)				
			}
		}
		else{lab1=c("Number of objects still in original cluster")}
		
		nrcluster=c()
		nrcomps=c()
		for(j in 1:length(Found)){
			nrcluster=c(nrcluster,Found[[j]]$nr.clusters)
			nrcomps=c(nrcomps,length(Found[[j]]$Complabels))
		}
		
		if(is.null(ClusterNumber)){
			xl=c(0,length(Found)+0.5)
		}
		else{
			xl=c(1,length(Found)+0.5)
		}
		
		if(!(is.null(location))){
			location=paste(location,"_SelectionPlot.pdf",sep="")
			
		}
		plottypein(plottype,location)
		
		graphics::plot(type="n",x=0,y=0,xlim=xl,ylim=c(0,max(nrcluster,nrcomps)+legendposy+0),xlab="",ylab="",xaxt="n",yaxt="n",cex.lab=1.25)
		graphics::lines(x=seq(1,length(Found)),y=nrcluster,lty=1,col="red",lwd=1.5)
		graphics::points(x=seq(1,length(Found)),y=nrcluster,pch=19,col="red",cex=1.5)
		
		if(is.null(ClusterNumber)){
			graphics::lines(x=c(0,1),y=c(length(Selection),nrcomps[1]),lty=1,col="black",lwd=1.5)
			graphics::points(x=0,y=length(Selection),pch=19,col="black",cex=1.5)
		}
		
		graphics::lines(x=seq(1,length(Found)),y=nrcomps,lty=1,col="blue",lwd=1.5)
		graphics::points(x=seq(1,length(Found)),y=nrcomps,pch=19,col="blue",cex=1.5)
		
		if(is.null(ClusterNumber)){
			graphics::text(0,length(Selection), "S",cex=1.5,pos=1,col="black",font=2)	
		}
		if(selectionPlot==TRUE & followMaxComps==TRUE)
			graphics::text(seq(1,length(Found)),nrcomps, labelcluster,cex=1.5,pos=1,col="black",font=2)
		
		if(is.null(ClusterNumber)){
			graphics::axis(1,at=seq(0,length(Found)),labels=c("Selection",names),las=2,cex=1)
		}
		else{
			graphics::axis(1,at=seq(1,length(Found)),labels=c(names),las=2,cex=1)
		}
		graphics::axis(2,at=seq(0,max(nrcluster,nrcomps)+2.0),labels=seq(0,max(nrcluster,nrcomps)+2.0),cex=1)
		
		
		graphics::legend(legendposx,max(nrcluster,nrcomps)+legendposy,lty=c(1,1,0),pch=c(19,19,0),col=c("blue","red","black"),legend=c(lab1,"Number of clusters original cluster divided amongst","Cluster number"),bty="n",cex=1.2)
		
		plottypeout(plottype)
	}
	if(completeSelectionPlot==TRUE){
		if(!(is.null(location))){
			location=paste(location,"_CompleteSelectionPlot.pdf",sep="")
			
		}
		plottypein(plottype,location)
		nrcluster=c(1)
		nrcomps=c(length(Selection))
		for(j in 1:length(Found)){
			nrcluster=c(nrcluster,Found[[j]]$nr.clusters)
			nrcomps=c(nrcomps,length(Found[[j]]$Complabels))
		}
		
		
		
		graphics::plot(type="n",x=0,y=0,xlim=c(0,length(Found)+0.5),ylim=c(0,max(nrcluster,nrcomps)+legendposy+0),xlab="",ylab="Number of objects",xaxt="n",yaxt="n")
		
		xnext=c()
		ynext=c()
		colorsp=c()
		
		if(is.null(ClusterNumber)){
			p=seq(0,length(Found))
		}
		else{
			p=seq(0,length(Found))
		}
		
		for(m in p){
			if(m==0){
				if(!(is.null(ClusterNumber))){
					howmany=length(Found[[1]]$AllClusters)
				}
				else{
					howmany=1
				}
			}
			else{
				howmany=length(Found[[m]]$AllClusters)
			}
			
			if(m==0){
				xnext=0.1
				ynext=c(length(Selection))
				if(!(is.null(ClusterNumber))){
					colorsp=cols[ClusterNumber]
				}
				else{
					colorsp=c("black")
				}
				
				
				graphics::points(x=xnext,y=ynext,col=colorsp,pch=19,cex=1.25)
				if(!(is.null(ClusterNumber))){
					labelcluster=ClusterNumber
				}
				else{
					labelcluster="S"
				}
				position=3
				if(!(is.integer(ynext[length(ynext)]))){
					position=1
				}
				graphics::text(xnext,ynext,labelcluster,cex=1.5,pos=position,col="black",font=2)					
				
				L1=list()
				for(n in 1:howmany){
					if(is.null(ClusterNumber)){				
						L1[[n]]=Selection
					}
					else{
						L1[[n]]=Found[[1]]$AllClusters[[n]][[3]]
					}
				}
				
			}
			else{
				xprev=xnext
				yprev=ynext
				ynext=c()
				colorsp=c()
				xnext=rep(seq(1,length(Found))[m],howmany)
				for(n in 1:howmany){										
					ynext=c(ynext,length(Found[[m]]$AllClusters[[n]][[3]]))
					
					colorsp=c(colorsp,cols[Found[[m]]$AllClusters[[n]]$clusternumber])
					
					if(length(ynext)>1){
						for(t in 1:(length(ynext)-1)){
							if(ynext[t]==ynext[length(ynext)]){
								ynext[length(ynext)]=ynext[length(ynext)]-0.3
							}
						}	
					}
					
					graphics::points(x=xnext[n],y=ynext[length(ynext)],col=colorsp[length(colorsp)],pch=19,cex=1.25)
					labelcluster=Found[[m]]$AllClusters[[n]]$clusternumber
					position=3
					if(!(is.integer(ynext[length(ynext)]))){
						position=1
					}
					graphics::text(xnext[n],ynext[length(ynext)],labelcluster,cex=1.5,pos=position,col="black",font=2)					
				}		
				
				
				L2=L1				
				L1=list()
				for(n in 1:howmany){
					L1[[n]]=Found[[m]]$AllClusters[[n]][[3]]					
				}
				
				for(q in 1:length(L1)){
					for(p in 1:length(L2)){
						if(length(which(L2[[p]] %in% L1[[q]])) != 0){
							graphics::segments(x0=xprev[p],y0=yprev[p],x1=xnext[q],y1=ynext[q],col=colorsp[q],lwd=2)
						}							
					}
				}
			}
		}
		if(is.null(ClusterNumber)){
			graphics::axis(1,at=seq(0,length(Found)),labels=c("Selection",names),las=2,cex=1.5)
		}
		else{
			graphics::axis(1,at=seq(0,length(Found)-1),labels=c(names),las=2,cex=1.5)
		}
		graphics::axis(2,at=seq(0,max(nrcluster,nrcomps)),labels=seq(0,max(nrcluster,nrcomps)),cex=1.5)
		graphics::legend(legendposx,max(nrcluster,nrcomps)+legendposy,pch=c(0),col=c("black"),legend=c("Cluster number"),bty="n",cex=1.2)
		
		
		plottypeout(plottype)
	}
	
	if(ClusterPlot==TRUE & !(is.null(ClusterNumber))){
		
		#FoundCl=Found[[length(Found)]]
		
		if(!(is.null(location))){
			location=paste(location,"_SelectionPlot.pdf",sep="")
			
		}
		plottypein(plottype,location)
		nrcluster=c(1)
		nrcomps=c(length(Selection))
		for(j in 1:length(FoundCl)){
			nrcomps=c(nrcomps,length(FoundCl[[j]][[1]]))
		}
		
		
		graphics::plot(type="n",x=0,y=0,xlim=c(0,length(FoundCl)-0.5),ylim=c(0,max(nrcluster,nrcomps)+legendposy+0),xlab="",ylab="Number of objects",xaxt="n",yaxt="n")
		
		xnext=c()
		ynext=c()
		colorsp=c()
		
		
		p=seq(1,length(FoundCl))
		
		for(m in p){
			
			if(m==1){
				xnext=0.1
				ynext=length(FoundCl[[1]][[1]])
				colorsp=cols[ClusterNumber]
				graphics::points(x=xnext,y=ynext,col=colorsp,pch=19,cex=1.25)
				
				labelcluster=ClusterNumber
				
				position=3
				if(!(is.integer(ynext[length(ynext)]))){
					position=1
				}
				graphics::text(xnext,ynext,labelcluster,cex=1.5,pos=position,col="black",font=2)		
				
				
			}
			
			else{
				xprev=xnext
				yprev=ynext
				
				ynext=c()
				colorsp=c()
				
				xnext=m-1
				ynext=length(FoundCl[[m]][[1]])
				colorsp=cols[ClusterNumber]
				
				graphics::points(x=xnext,y=ynext,col=colorsp,pch=19,cex=1.25)
				graphics::segments(x0=xprev,y0=yprev,x1=xnext,y1=ynext,col=colorsp)
				labelcluster=ClusterNumber
				
				position=3
				if(!(is.integer(ynext[length(ynext)]))){
					position=1
				}
				graphics::text(xnext,ynext,labelcluster,cex=1.5,pos=position,col="black",font=2)
				
				#Dissapeared objects
				if(length(FoundCl[[m]]$Dissapeared)!=0){
					
					xdiss=rep(xnext,length(FoundCl[[m]]$Dissapeared))
					ydiss=c()
					colorsd=c()
					labs=c()
					for(d in 1:length(FoundCl[[m]]$Dissapeared)){
						ydiss=c(ydiss,length(FoundCl[[m]]$Dissapeared[[d]][[1]]))
						colorsd=c(colorsd,cols[FoundCl[[m]]$Dissapeared[[d]][[2]]])
						labs=c(labs,FoundCl[[m]]$Dissapeared[[d]][[2]])
					}
					
					if(length(ydiss)>1){
						for(t in 1:(length(ydiss)-1)){
							if(ydiss[t]==ydiss[length(ydiss)]){
								ydiss[length(ydiss)]=ydiss[length(ydiss)]-0.3
							}
						}	
					}
					
					graphics::points(x=xdiss,y=ydiss,col=colorsd,pch=19,cex=1.25)
					labelcluster=labs
					position=3
					if(!(is.integer(ydiss[length(ydiss)]))){
						position=1
					}
					graphics::text(xdiss,ydiss,labelcluster,cex=1.5,pos=position,col="black",font=2)	
					
					for(p in 1:length(ydiss)){
						graphics::segments(x0=xprev,y0=yprev,x1=xdiss[p],y1=ydiss[p],col=colorsd[p])
					}	
					
				}	
				
				
				if(length(FoundCl[[m]]$Joined)!=0){
					
					xjoin=rep(xprev,length(FoundCl[[m]]$Joined))
					yjoin=c()
					colorsj=c()
					labs=c()
					for(d in 1:length(FoundCl[[m]]$Joined)){
						yjoin=c(yjoin,length(FoundCl[[m]]$Joined[[d]][[1]]))
						colorsj=c(colorsj,cols[FoundCl[[m]]$Joined[[d]][[2]]])
						labs=c(labs,FoundCl[[m]]$Joined[[d]][[2]])
					}
					
					if(length(yjoin)>1){
						for(t in 1:(length(yjoin-1))){
							if(yjoin[t]==yjoin[length(yjoin)]){
								yjoin[length(yjoin)]=yjoin[length(yjoin)]-0.3
							}
						}	
					}
					
					graphics::points(x=xjoin,y=yjoin,col=colorsj,pch=19,cex=1.25)
					labelcluster=labs
					position=3
					if(!(is.integer(yjoin[length(yjoin)]))){
						position=1
					}
					graphics::text(xjoin,yjoin,labelcluster,cex=1.5,pos=position,col="black",font=2)	
					
					for(p in 1:length(yjoin)){
						graphics::segments(x0=xjoin[p],y0=yjoin[p],x1=xnext,y1=ynext,col=colorsj[p])
					}	
					
				}	
				
			}
			
		}
		graphics::axis(1,at=c(0.1,seq(1,length(FoundCl)-1)),labels=c(names),las=2,cex=1.5)
		graphics::axis(2,at=seq(0,max(nrcomps)),labels=seq(0,max(nrcomps)),cex=1.5)
		graphics::legend(legendposx,max(nrcluster,nrcomps)+legendposy,pch=c(0),col=c("black"),legend=c("Cluster number"),bty="n",cex=1.2)
		
		plottypeout(plottype)
		
	}	
	
	if(table==TRUE & selectionPlot==TRUE){
		SharedComps=Selection
		Extra=list()
		temp=c()
		for(a in 1:length(Found)){
			SharedComps=intersect(SharedComps,Found[[a]]$Complabels)
		}
		for(a in 1:length(Found)){
			Extra[[a]]=Found[[a]]$Complete.new.cluster[which(!(Found[[a]]$Complete.new.cluster%in%SharedComps))]
			names(Extra)[a]=names[a]
			if(followMaxComps==TRUE){
				names(Extra)[a]=paste(names(Extra)[a],"_",labelcluster[a],sep="")
			}
			temp=c(temp,length(Extra[[a]]))
		}
		
		ExtraOr=Selection[which(!(Selection%in%SharedComps))]
		
		collength=max(length(SharedComps),length(ExtraOr),temp)
		
		if(length(SharedComps)<collength){
			spaces=collength-length(SharedComps)
			SharedComps=c(SharedComps,rep(" ",spaces))
		}
		
		if(length(ExtraOr)<collength){
			spaces=collength-length(ExtraOr)
			ExtraOr=c(ExtraOr,rep(" ",spaces))
		}
		
		for(b in 1:length(Extra)){
			if(length(Extra[[b]])<collength){
				spaces=collength-length(Extra[[b]])
				Extra[[b]]=c(Extra[[b]],rep(" ",spaces))
			}
			
		}		
		table1=data.frame(ExtraOr=ExtraOr,SharedComps=SharedComps)
		for(b in 2:length(Extra)){
			table1[1+b]=Extra[[b]]
			colnames(table1)[1+b]=names(Extra)[b]
		}
		
	}
	else{table1=list()}
	
	names(Found)=names
	Found[[length(Found)+1]]=table1
	names(Found)[length(Found)]="Table"
	return(Found)	
}






