% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SegFunctions.R
\name{spatmultiseg}
\alias{spatmultiseg}
\title{A function from seg package to compute spatial multi-group segregation indices}
\usage{
spatmultiseg(x, spatobj = NULL, folder = NULL, shape = NULL, ...)
}
\arguments{
\item{x}{- an object of class matrix (or which can be coerced to that class), 
where each column represents the distribution of a group within 
spatial units. The number of columns should be greater than 1 (at least 2 
groups are required). You should not include a column with total 
population, because this will be interpreted as a group.}

\item{spatobj}{- a spatial object (SpatialPolygonsDataFrame) with 
geographic information}

\item{folder}{- a character vector with the folder (directory) 
name indicating where the shapefile is located on the drive}

\item{shape}{- a character vector with the name of the shapefile 
(without the .shp extension) .}

\item{...}{- other parameters of spseg function from seg package.}
}
\value{
A vector with Reardon's spatial multi-group segregation indices:
D* - spatial multi-group dissimilarity index
R* - spatial multi-group relative diversity index
H* - spatial multi-group information theory index
}
\description{
A function adapted from seg package (Hong et al. 2014) 
to compute spatial multi-group segregation indices developed by 
Reardon and O'Sullivan (2004)
}
\examples{
x <- segdata@data[ ,1:2]
foldername <- system.file('extdata', package = 'OasisR')
shapename <- 'segdata'

spatmultiseg(x, spatobj = segdata)

spatmultiseg(x, folder = foldername, shape = shapename) 

}
\references{
Reardon, S. F. and O'Sullivan, D. (2004) 
\emph{Measures of spatial segregation}.
Sociological Methodology, 34, 121-162.

Hong S.Y., O'Sullivan D., Sadahiro Y. (2014) 
\emph{Implementing Spatial Segregation Measures in R'}.
PLoS ONE, 9(11)
}
\seealso{
Multi-group indices: 
\code{\link{PMulti}}, \code{\link{GiniMulti}}, \code{\link{DMulti}},  
\code{\link{HMulti}}, \code{\link{RelDivers}}

Social diversity indices: 
\code{\link{HShannon}}, \code{\link{NShannon}}, 
\code{\link{ISimpson}},
}
