# Project: Integrated Clustering
# 
# Author: lucp8409
###############################################################################

#' Target prediction data
#' 
#' A binary data matrix that contains 477 target predictions for a set of 56 compounds.
#' 
#' @name targetMat
#'
#' @usage data(targetMat)
#'
#' @format An object of class \code{"matrix"}.
#'
#' @keywords datasets
#' @examples
#' data(targetMat)
NULL


#' Fingerprint data
#' 
#' A binary data matrix that contains 250 fingerprints for a set of 56 compounds.
#' 
#' @name fingerprintMat
#'
#' @usage data(fingerprintMat)
#'
#' @format An object of class \code{"matrix"}.
#'
#' @keywords datasets
#' @examples
#' data(fingerprintMat)
NULL


#' Gene expression data
#' 
#' Gene expression data for 2434 genes for a set of 56 compounds.
#' 
#' @name geneMat
#'
#' @usage data(geneMat)
#'
#' @format An object of class \code{"matrix"}.
#'
#' @keywords datasets
#' @examples
#' data(geneMat)
NULL

#' Colour examples
#' 
#' A vector of HEX codes for the colours used in the examples
#' 
#' @name Colors1
#' 
#' @format An object if class \code{"character"}.
#' @keywords datasets
#' @examples
#' data(Colors1)
NULL

#' Information of the genes
#'  
#' Gene info in a data frame
#' 
#' @name GeneInfo
#' 
#' @format A data frame with 3 variables: ENTREZID, SYMBOL and GENENAME
#' @keywords datasets
#' @examples
#' data(GeneInfo)
NULL


#' List of GO Annotations 
#' 
#' A list that contains the GO annotations produced by getGeneSets of the \code{MLP} package for the genes in the geneMat data.
#' 
#' @name GS
#' 
#' @format A data frame with 3 variables: ENTREZID, SYMBOL and GENENAME
#' @keywords datasets
#' @examples
#' data(GS)
NULL


