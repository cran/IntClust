\name{ColorPalette}
\alias{ColorPalette}
\title{
Create a color palette to be used in the plots
}
\description{
In order to facilitate the visualization of the influence of the different
methods on the clustering of the compounds, colors can be used. The function
\code{ColorPalette} is able to pick out as many colors as there are clusters.
This is done with the help of the \code{ColorRampPalette} function of the 
grDevices package
}
\usage{
ColorPalette(colors = c("red", "green"), ncols = 5)
}
\arguments{
  \item{colors}{
A vector containing the colors of choice
}
  \item{ncols}{
The number of colors to be specified. If higher than the number of colors,
it specifies colors in the region between the given colors.
}
}
\value{
A vector containing the hex codes of the chosen colors.
}
\author{
Marijke Van Moerbeke
}
\examples{
Colors1<-ColorPalette(c("cadetblue2","chocolate","firebrick2",
"darkgoldenrod2", "darkgreen","blue2","darkorchid3","deeppink2"), ncols=8)
}
\keyword{ Colors }
