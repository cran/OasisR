% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/SegFunctions.R
\name{GiniMulti}
\alias{GiniMulti}
\title{A function to compute multigroup Gini index}
\usage{
GiniMulti(x)
}
\arguments{
\item{x}{- an object of class matrix (or that can be coerced
to that class), where each column represents the distribution
of a population group, within spatial units. The number of
columns should be greater than one (at least two population
groups are required).}
}
\description{
Multigroup version of \code{\link{Gini}} index
}
\examples{
x <- slot(GreHSize, 'data')[ ,3:5]
GiniMulti(x)
}
\references{
Reardon S. F. (1998) \emph{Measures of racial
diversity and segregation in multigroup and hierarchical
structured Populations}. Annual meeting of the Eastern
Sociological Society, Philadelphia
}
\seealso{
Other multigroup eveness indices:
\code{\link{HShannon}}, \code{\link{JPielou}},
\code{\link{ISimpson}}, \code{\link{DMulti}}

Other multigroup indices: \code{\link{PIsol}},
\code{\link{RelDivers}}
}
