% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SegFunctions.R
\name{Gini2}
\alias{Gini2}
\title{A function to compute Spatial Gini's between group index}
\usage{
Gini2(x)
}
\arguments{
\item{x}{- an object of class matrix (or which can be coerced to that class), 
where each column represents the distribution of a group within 
spatial units. The number of columns should be greater than 1 (at least 2 
groups are required). You should not include a column with total 
population, because this will be interpreted as a group.}
}
\value{
a matrix with between group Gini index
}
\description{
The between group version of Gini index is obtained 
by computing the index for a subpopulation formed by each pair of groups
}
\examples{
x <- segdata@data[ ,1:2]
Gini2(x)
}
\references{
Duncan O. D. and Duncan B. (1955) \emph{A Methodological 
Analysis of Segregation Indexes}. American Sociological Review 41, 
pp. 210-217
}
\seealso{
Other one-group  evenness indices: 
\code{\link{ISDuncan}}, \code{\link{Gini}}, 
\code{\link{Gorard}}, \code{\link{Atkinson}}, 
\code{\link{HTheil}}, '\code{\link{ISWong}}, \code{\link{ISMorrill}},
\code{\link{ISMorrillK}}

Between groups dissimilarity indices: 
\code{\link{DIDuncan}}, \code{\link{DIMorrill}}, 
\code{\link{DIMorrillK}}, \code{\link{DIWong}}
}
