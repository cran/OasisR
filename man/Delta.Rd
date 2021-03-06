% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SegFunctions.R
\name{Delta}
\alias{Delta}
\title{A function to compute Delta index}
\usage{
Delta(x, a = NULL, spatobj = NULL, folder = NULL, shape = NULL)
}
\arguments{
\item{x}{- an object of class matrix (or which can be coerced to that class), 
where each column represents the distribution of a group within 
spatial units. The number of columns should be greater than 1 (at least 2 
groups are required). You should not include a column with total 
population, because this will be interpreted as a group.}

\item{a}{- a numeric vector containing spatial unit areas}

\item{spatobj}{- a spatial object (SpatialPolygonsDataFrame) with 
geographic information}

\item{folder}{- a character vector with the folder (directory) 
name indicating where the shapefile is located on the drive}

\item{shape}{- a character vector with the name of the shapefile 
(without the .shp extension).}
}
\value{
a numeric vector containing Delta index value for 
each group
}
\description{
The Delta index is a specific application of dissimilarity 
index \code{\link{DIDuncan}} which simply measures the dissimilarity
between the spatial distribution of a group and the spatial 
distribution of available land. It can be interpreted as the share of group 
that would have to move to achieve uniform density over all spatial units. 
The function can be used in two ways: to provide an area vector or 
a external geographic information source (spatial object or shape file).
}
\examples{
x <- segdata@data[ ,1:2]
ar <- area(segdata)
foldername <- system.file('extdata', package = 'OasisR')
shapename <- 'segdata'

Delta(x, a = ar) 

Delta(x, spatobj = segdata)

Delta(x, folder = foldername, shape = shapename) 
}
\references{
Duncan O. D., Cuzzoert  and Duncan B. (1961) 
\emph{Problems in analyzing areal data}. Statistical geography, 
Glencoe, Illinois: The free press of Glencoe
}
\seealso{
Absolute Concentration Index: \code{\link{ACO}}

Relative Concentration Index: \code{\link{RCO}}
}
