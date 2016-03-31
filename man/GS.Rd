\name{GS}
\alias{GS}
\docType{data}
\title{
List of GO Annotations 
}
\description{
A list that contains the GO annotations produced by getGeneSets of the \code{MLP} package for the
genes in the geneMat data.
}
\usage{data(GS)}
\format{
  The format is:
List of 8734
 $ GO:0000002: chr [1:20] "291" "1763" "1890" "3980" ...
 $ GO:0000003: chr [1:925] "18" "49" "51" "90" ...
 $ GO:0000012: chr [1:8] "3981" "7141" "7515" "23411" ...
 $ GO:0000018: chr [1:51] "604" "641" "940" "958" ...
 $ GO:0000019: chr [1:4] "641" "4292" "4361" "10111"
 $ GO:0000022: chr [1:2] "9055" "9493"
 ...
 $ GO:0000724: chr [1:70] "472" "641" "672" "675" ...
  [list output truncated]
 - attr(*, "species")= chr "Human"
 - attr(*, "geneSetSource")= chr "GOBP"
 - attr(*, "descriptions")= Named chr [1:13226] "mitochondrial genome maintenance" "reproduction" "single strand break repair" "regulation of DNA recombination" ...
  ..- attr(*, "names")= chr [1:13226] "GO:0000002" "GO:0000003" "GO:0000012" "GO:0000018" ...
 - attr(*, "class")= chr [1:2] "geneSetMLP" "list"
}
\keyword{datasets}
