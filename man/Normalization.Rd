\name{Normalization}
\alias{Normalization}
\title{
A normalization function
}
\description{
If data of different scales are being employed by the user, it is recommended to perform
a normalization to make the data structures comparable. This is performed by the
\code{Normalization} function.
}
\usage{
Normalization(Data,method=c("Quantile","Fisher-Yates","Standardize",
"Range","Q","q","F","f","S","s","R","r"))
}
\arguments{
  \item{Data}{
A data matrix. It is assumed the  rows are corresponding with the objects.
}
  \item{method}{
A method of normalization. Should be one of "Quantile","Fisher-Yates",
"standardize","Range" or any of the first letters of these names.
}
}
\details{
The method "Quantile" refers to the Quantile-Normalization widely used in 
omics data. The "Fisher-Yates" normalization has a similar approach as the Quantile-
Normalization but does not rely on the data, just on the number of rows present in
the data matrix. The "Standardize" method refers to the \code{stdize} function of the
\pkg{pls} package and centers and scales the data matrix. The method "Range" computes
the maximum and minimum value of the matrix and determines the range. Every value is then
reduced by the minimum and divided by the range of the data matrix. The latter normalization will
result in values between 0 and 1.
}
\value{
The returned value is a distance matrix.
}
\author{
Marijke Van Moerbeke
}
\examples{
x=matrix(rnorm(100),ncol=10,nrow=10)
Norm_x=Normalization(x,method="R")
}
\keyword{Normalize}