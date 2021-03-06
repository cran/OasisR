% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SegFunctions.R
\name{RCL}
\alias{RCL}
\title{A function to compute the relative clustering index (RCL)}
\usage{
RCL(x, d = NULL, fdist = 'e', distin = 'm',  distout = 'm', diagval = '0', 
beta = 1, spatobj = NULL, folder = NULL, shape = NULL)
}
\arguments{
\item{x}{- an object of class matrix (or which can be coerced to that class), 
where each column represents the distribution of a group within 
spatial units. The number of columns should be greater than 1 (at least 2 
groups are required). You should not include a column with total 
population, because this will be interpreted as a group.}

\item{d}{- a matrix of the distances between spatial unit centroids}

\item{fdist}{- the method used for distance interaction matrix: 
e' for inverse exponential function (by default) and 'l' for linear.}

\item{distin}{- input metric conversion, based on  \pkg{bink} package and 
includes conversions from 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.}

\item{distout}{- output metric conversion, based on  \pkg{bink} package and 
includes conversions to 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.}

\item{diagval}{- when providing a spatial object or a shape file, 
the user has the choice of the spatial matrix diagonal definition: 
diagval = '0' (by default) for an null diagonal and diagval = 'a' 
to compute the diagonal as 0.6 * square root (spatial/organizational unitsarea) (White, 1983)}

\item{beta}{- distance decay parameter}

\item{spatobj}{- a spatial object (SpatialPolygonsDataFrame) with 
geographic information}

\item{folder}{- a character vector with the folder (directory) 
name indicating where the shapefile is located on the drive}

\item{shape}{- a character vector with the name of the shapefile 
(without the .shp extension).}
}
\value{
a matrix containing relative clustering index values for each pair of groups
}
\description{
The relative clustering index, RCL, compares the mean 
proximity of a group to the mean proximity of another group. 
The function can be used in two ways: to provide a distance matrix 
or a external geographic information source (spatial object 
or shape file).
}
\examples{
x <- segdata@data[ ,1:2]
ar<-area(segdata)
dist <- distance(segdata)
diag(dist)<-sqrt(ar) * 0.6
foldername <- system.file('extdata', package = 'OasisR')
shapename <- 'segdata'

RCL(x, spatobj = segdata)

RCL(x, folder = foldername, shape = shapename, fdist = 'l') 

RCL(x, spatobj = segdata, diagval ='a')

RCL(x, d = dist, fdist = 'e')

}
\references{
Massey D. S. and Denton N. A. (1988) \emph{The dimensions 
of residential segregation}. Social Forces 67(2),  pp. 281-315.
}
\seealso{
Proximity measures: \code{\link{Pxx}}, 
\code{\link{Pxy}},  \code{\link{Poo}},  \code{\link{SP}}

Clustering Indices: \code{\link{ACL}}
}
