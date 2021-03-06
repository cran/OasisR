% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SegFunctions.R
\name{ACEDuncan}
\alias{ACEDuncan}
\title{A function to compute Duncan's Absolute Centralisation Index (ACEDuncan)}
\usage{
ACEDuncan(x, dc = NULL, center = 1, 
spatobj = NULL, folder = NULL, shape = NULL)
}
\arguments{
\item{x}{- an object of class matrix (or which can be coerced to that class), 
where each column represents the distribution of a group within 
spatial units. The number of columns should be greater than 1 (at least 2 
groups are required). You should not include a column with total 
population, because this will be interpreted as a group.}

\item{dc}{- a numeric vector containing the distances between spatial units
centroids and the central spatial unit}

\item{center}{- a numeric value giving the number of the spatial unit that 
represents the center in the table}

\item{spatobj}{- a spatial object (SpatialPolygonsDataFrame) with 
geographic information}

\item{folder}{- a character vector with the folder (directory) 
name indicating where the shapefile is located on the drive}

\item{shape}{- a character vector with the name of the shapefile 
(without the .shp extension).}
}
\value{
a numeric vector containing Duncan's asolute centralisation index 
value for each group
}
\description{
Duncan's absolute centralization index measures the 
proportion of a group that should change its localization to 
achieve the same level of centralization as the rest of the population.
The function can be used in two ways: to provide a vector containing 
the distances between spatial/organizational unit centroids or a external geographic 
information source (spatial object or shape file).
}
\examples{
x <- segdata@data[ ,1:2]
distc<- distcenter(segdata, center = 28)
foldername <- system.file('extdata', package = 'OasisR')
shapename <- 'segdata'

ACEDuncan(x, dc=distc) 

ACEDuncan(x, spatobj = segdata, center = 28) 

ACEDuncan(x, folder = foldername, shape = shapename, center = 28) 

}
\references{
Duncan O. D. and Duncan B. (1955) \emph{A 
Methodological Analysis of Segregation Indexes}. 
American Sociological Review 41, pp. 210-217
}
\seealso{
\code{\link{ACEDuncanPoly}}, \code{\link{ACEDuncanPolyK}},

\code{\link{RCE}}, \code{\link{RCEPoly}}, \code{\link{RCEPolyK}}

\code{\link{ACE}}, \code{\link{ACEPoly}}
}
