% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SegFunctions.R
\name{DIDuncan}
\alias{DIDuncan}
\title{A function to compute Duncan dissimilarity segregation index}
\usage{
DIDuncan(x)
}
\arguments{
\item{x}{- an object of class matrix (or which can be coerced to that class), 
where each column represents the distribution of a group within 
spatial units. The number of columns should be greater than 1 (at least 2 
groups are required). You should not include a column with total 
population, because this will be interpreted as a group.}
}
\value{
a matrix containing dissimilarity index values
}
\description{
Duncan's dissimilarity index is the segregation index 
most commonly used in the literature. It is derived from Lorenz 
curves as the maximum difference between the segregation curve 
and the diagonal. The index measures the unevenness of a group's 
spatial distribution compared to another group. It can be 
interpreted as the share of the group that would have to move to 
achieve an even distribution compared to another group.
}
\examples{
x <- segdata@data[ ,1:2]
DIDuncan(x) 
}
\references{
Duncan O. D. and Duncan B. (1955) \emph{A Methodological 
Analysis of Segregation Indexes}. American Sociological Review 41, 
pp. 210-217
}
\seealso{
Other one-group  evenness indices: 
\code{\link{ISDuncan}}, \code{\link{Gini}}, \code{\link{Gorard}}, 
\code{\link{Atkinson}}, \code{\link{HTheil}}, 
\code{\link{ISWong}}, \code{\link{ISMorrill}}, \code{\link{ISMorrillK}}

Between groups dissimilarity indices: 
\code{\link{DIMorrill}}, \code{\link{DIMorrillK}}, \code{\link{DIWong}}
}
