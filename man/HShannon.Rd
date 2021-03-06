% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SegFunctions.R
\name{HShannon}
\alias{HShannon}
\title{A function to compute Shannon-Wiener diversity (entropy) index}
\usage{
HShannon(x)
}
\arguments{
\item{x}{- an object of class matrix (or which can be coerced to that class), 
where each column represents the distribution of a group within 
spatial units. The number of columns should be greater than 1 (at least 2 
groups are required). You should not include a column with total 
population, because this will be interpreted as a group.}
}
\value{
Shannon-Wiener diversity index
}
\description{
The Shannon-Wiener diversity index is based on 
the notion of entropy and measures population heterogeneity.
}
\examples{
x <- segdata@data[ ,1:2]
HShannon(x) 
}
\references{
Shannon C. E. (1948) \emph{A mathematical theory 
of communication}. Bell System Technical Journal (27)
}
\seealso{
Social diversity indices: 
\code{\link{NShannon}}, \code{\link{ISimpson}},

Multi-group indices: 
\code{\link{PMulti}}, \code{\link{GiniMulti}}, \code{\link{DMulti}},  
\code{\link{HMulti}}, \code{\link{CMulti}}, \code{\link{RelDivers}}
}
