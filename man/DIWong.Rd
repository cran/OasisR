% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SegFunctions.R
\name{DIWong}
\alias{DIWong}
\title{A function to compute Wongs's dissimilarity index}
\usage{
DIWong(x, b = NULL,  a = NULL, p = NULL, ptype = 'int', variant = 's', 
spatobj = NULL, folder = NULL, shape = NULL)
}
\arguments{
\item{x}{- an object of class matrix (or which can be coerced to that class), 
where each column represents the distribution of a group within 
spatial units. The number of columns should be greater than 1 (at least 2 
groups are required). You should not include a column with total 
totals because this will be interpreted as a group}

\item{b}{- a common boundaries matrix where each element \emph{Bij} 
equals the shared boundary of \emph{i}-th and \emph{j}-th spatial units.}

\item{a}{- a numeric vector containing spatial unit areas}

\item{p}{- a numeric vector containing spatial units perimeters.}

\item{ptype}{- a string variable giving two options for perimeter calculation
when a spatial object or shapefile is provided: 'int' to use only interior
borders of spatial units, and 'all' to use entire borders, including to
the exterior of the area}

\item{variant}{- a character variable that allows to choose the index version: 
variant = 's' for the dissimilarity index adjusted for contiguous spatial units
boundary lengths and perimeter/area ratio (by default) and variant = 'w' 
for the version without perimeter/area ratio}

\item{spatobj}{- a spatial object (SpatialPolygonsDataFrame) with 
geographic information}

\item{folder}{- a character vector with the folder (directory) 
name indicating where the shapefile is located on the drive}

\item{shape}{- a character vector with the name of the shapefile 
(without the .shp extension).}
}
\value{
a matrix containing Wong's dissimilarity index values
}
\description{
Wong's dissimilarity index is a development of 
\code{\link{DIDuncan}}'s which takes into account the interactions 
between spatial units(common boundaries and perimeter/area ratios). 
The function can be used in two ways: to provide spatial data (
boundaries matrix, a perimeter vector and an area vector) 
or a external geographic information source (spatial object or shape file).
}
\examples{
x <- segdata@data[ ,1:2]
bound <- boundaries(segdata)
per <- perimeter(segdata)
ar <- area(segdata)
foldername <- system.file('extdata', package = 'OasisR')
shapename <- 'segdata'

DIWong(x, b = bound, p = per, a = ar) 

DIWong(x, spatobj = segdata, variant = 'w') 

DIWong(x, folder = foldername, shape = shapename, ptype ='all') 
}
\references{
Wong D. W. S. (1993) \emph{Spatial Indices of Segregation}. 
Urban Studies, 30 (3), pp. 559-572.
}
\seealso{
Other one-group  evenness indices: 
\code{\link{ISDuncan}}, \code{\link{Gini}}, \code{\link{Gorard}}, 
\code{\link{Atkinson}}, \code{\link{HTheil}}, 
'\code{\link{ISWong}}, \code{\link{ISMorrill}}, \code{\link{ISMorrillK}}

Between groups dissimilarity indices: 
\code{\link{DIDuncan}}, \code{\link{DIMorrill}}, \code{\link{DIMorrillK}}
}
