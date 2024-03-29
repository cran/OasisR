% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SegFunctions.R
\name{ISDuncan}
\alias{ISDuncan}
\title{A function to compute Duncan & Duncan segregation index}
\usage{
ISDuncan (x)
}
\arguments{
\item{x}{- an object of class matrix (or which can be coerced to that class), 
where each column represents the distribution of a group within 
spatial units. The number of columns should be greater than 1 (at least 2 
groups are required). You should not include a column with total 
population, because this will be interpreted as a group.}
}
\value{
a numeric vector with values of the Duncan's segregation index 
for each group
}
\description{
Duncan's segregation index is one-group form of 
dissimilarity index \code{\link{DIDuncan}} and  
measures the unevenness of a group distribution  
compared to the rest of the population. It can be interpreted
as the share of the group that would have to move to achieve 
an even distribution compared to the rest of the population.
}
\examples{
x <- segdata@data[ ,1:2]
ISDuncan(x) 
}
\references{
Duncan O. D. and Duncan B. (1955) \emph{ 
Residential Distribution and Occupational Stratification}. 
American Journal of Sociology 60 (5), pp. 493-503
}
\seealso{
One-group evenness indices: 
\code{\link{Gini}}, \code{\link{Atkinson}}, \code{\link{Gorard}}, 
\code{\link{HTheil}}, '\code{\link{ISWong}}, \code{\link{ISMorrill}},
\code{\link{ISMorrillK}}

Between groups dissimilarity indices: 
\code{\link{DIDuncan}}, \code{\link{Gini2}}, 
\code{\link{DIMorrill}}, \code{\link{DIMorrillK}}, \code{\link{DIWong}}
}
