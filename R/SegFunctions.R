

################################################## 

# SEGREGATION INDEXES

################################################## 

# EVENESS INDEXES

################### 


#' A function to compute Duncan & Duncan segregation index
#'
#' @usage ISDuncan (x) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @return a numeric vector containing the Duncan's segregation index value for 
#' each population group 
#' @references Duncan O. D. and Duncan B. (1955) \emph{ 
#' Residential Distribution and Occupational Stratification}. 
#' American Journal of Sociology 60 (5), pp. 493-503
#' @description Duncan's segregation index is one-group form of 
#' dissimilarity index \code{\link{DIDuncan}} and  
#' measures the unevenness of a group spatial distribution  
#' compared to the rest of the population. It can be interpreted
#' as the share of the group that would have to move to achieve 
#' an even distribution, compared to the rest of the population.
#' @examples x <- segdata@data[ ,1:2]
#' ISDuncan(x) 
#' @seealso One-group evenness indices: 
#' \code{\link{Gini}}, \code{\link{Atkinson}}, \code{\link{Gorard}}, 
#' \code{\link{HTheil}}, '\code{\link{ISWong}}, \code{\link{ISMorrill}}
#' @seealso Between groups dissimilarity indices: 
#' \code{\link{DIDuncan}}, \code{\link{Gini2}}, 
#' \code{\link{DIMorrill}}, \code{\link{DIWong}}
#' @export


ISDuncan <- function(x) {
    x <- as.matrix(x)
    result <- vector(length = ncol(x))
    Total <- sum(x)
    pTotal <- colSums(x)/Total
    tx <- rowSums(x)
    px <- x/tx
    for (i in 1:ncol(x)) {
        px[, i] <- tx * abs(px[, i] - pTotal[i])
        result[i] <- sum(px[, i])/(2 * Total * pTotal[i] * (1 - pTotal[i]))
    }
    return(round(result, 4))
}

#' A function to compute Atkinson segregation index
#'
#' @usage Atkinson (x, delta = 0.5) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param delta - a inequality aversion parameter, by default equal to 0.5, 
#' varying from 0 to 1.
#' @return a numeric vector containing the Atkinson's segregation index value for 
#' each population group 
#' @references James, D. and K. E. Taeuber (1985)  \emph{ Measures 
#' of Segregation}. Sociological Methodology 15, pp. 1-32
#' @description The spatial version of Atkinson inequality index is derived
#' from Lorenz curves. It allows to decide wich part of the curve contribute more 
#' to the index, by choosing the value of the shape parameter, delta. 
#' @examples x <- segdata@data[ ,7:8]
#' Atkinson(x) 
#' Atkinson(x, 0.1)
#' Atkinson(x, delta = 0.9)
#' @seealso One-group evenness indices: 
#' \code{\link{ISDuncan}}, \code{\link{Gini}}, \code{\link{Gorard}}, 
#' \code{\link{HTheil}}, '\code{\link{ISWong}}, \code{\link{ISMorrill}}
#' @seealso Between groups dissimilarity indices: 
#' \code{\link{DIDuncan}}, \code{\link{Gini2}}, 
#' \code{\link{DIMorrill}}, \code{\link{DIWong}}
#' @export

Atkinson <- function(x, delta = 0.5) {
    x <- as.matrix(x)
    result <- vector(length = ncol(x))
    Total <- sum(x)
    pTotal <- colSums(x)/Total
    tx <- rowSums(x)
    px <- x/tx
    provi <- px
    for (k in 1:ncol(x)) {
        provi[, k] <- (((1 - px[, k])^(1 - delta) * px[, k]^delta * tx)/(pTotal[k] * Total))
        result[k] <- 1 - (pTotal[k]/(1 - pTotal[k])) * sum(provi[, k])^(1/(1 - delta))
    }
    return(round(result, 4))
}

#' A function to compute Theil's entropy segregation index
#'
#' @usage HTheil (x) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @return a numeric vector containing Theils's entropy index value for 
#' each population group 
#' @references Theil H. (1972)  \emph{Statistical decomposition analysis: with 
#' applications in the social and administrative.} Amsterdam, North-Holland, 337 p.
#' @description The entropy index (also called information index) measures
#' departure from evenness by assessing each spatial unit deviation from the 
#' entropy in the area. 
#' @examples x <- segdata@data[ ,1:2]
#' HTheil(x) 
#' @seealso One-group evenness indices: 
#' \code{\link{ISDuncan}}, \code{\link{Gini}}, \code{\link{Gorard}}, 
#' \code{\link{Atkinson}}, '\code{\link{ISWong}}, \code{\link{ISMorrill}}
#' @seealso Between groups dissimilarity indices: 
#' \code{\link{DIDuncan}}, \code{\link{Gini2}}, \code{\link{DIMorrill}}, \code{\link{DIWong}}
#' @export


HTheil <- function(x) {
    x <- as.matrix(x)
    result <- vector(length = ncol(x))
    E <- matrix(data = 0, nrow = nrow(x), ncol = ncol(x))
    Total <- sum(x)
    pTotal <- colSums(x)/Total
    tx <- rowSums(x)
    px <- x/tx
    Etot <- rep(0, ncol(x))
    for (k in 1:ncol(x)) Etot[k] <- pTotal[k] * log(1/pTotal[k]) + (1 - pTotal[k]) * log(1/(1 - pTotal[k]))
    for (k in 1:ncol(x)) {
        E[, k] <- px[, k] * log(1/px[, k]) + (1 - px[, k]) * log(1/(1 - px[, k]))
        if (any(is.nan(E[, k]))) 
            E[, k][is.nan(E[, k])] <- 0
        result[k] <- sum((tx * (Etot[k] - E[, k]))/(Etot[k] * Total))
    }
    return(round(result, 4))
}



#' A function to compute Morrill's segregation index 
#'
#' @usage ISMorrill(x, c = NULL, queen = TRUE, 
#' spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param c - a standard binary contiguity (adjacency) symmetric matrix where 
#' each element \emph{Cij} equals 1 if \emph{i}-th and \emph{j}-th spatial 
#' units are adjacent, and 0 otherwise.
#' @param queen - logical parameter difining criteria used for contiguity 
#' matrix computation, TRUE (by default) for queen , FALSE for rook 
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)  which contains the geographic information
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return a numeric vector containing the Morrill's segregation index value for 
#' each population group
#' @examples x <- segdata@data[ ,1:2]
#' contiguity <- contig(segdata)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'segdata'
#' 
#' ISMorrill(x, c = contiguity) 
#' 
#' ISMorrill(x, spatobj = segdata)
#' 
#' ISMorrill(x, folder = foldername, shape = shapename) 
#' 
#' @references Morrill B. (1991) \emph{On the measure of geographic 
#' segregation}. Geography research forum, 11, pp. 25-36.
#' @description Morrill's segregation index is a development of 
#' \code{\link{ISDuncan}}'s index which takes into account the 
#' interactions between spatial units (contiguity). 
#' The function can be used in two ways: to provide a contiguity 
#' matrix or a external geographic information source (spatial object 
#' or shape file).
#' @seealso One-group evenness indices: 
#' \code{\link{ISDuncan}}, \code{\link{Gini}}, \code{\link{Gorard}}, 
#' \code{\link{HTheil}}, \code{\link{Atkinson}}, '\code{\link{ISWong}}
#' @seealso Between groups dissimilarity indices: 
#' \code{\link{DIDuncan}}, \code{\link{Gini2}}, 
#' \code{\link{DIMorrill}}, \code{\link{DIWong}}
#' @export


ISMorrill <- function(x, c = NULL, queen = TRUE, spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(c)) 
        c <- contig(spatobj = spatobj, folder = folder, shape = shape, queen = queen)
    IS <- ISDuncan(x)
    result <- vector(length = ncol(x))
    pij <- array(0, dim = c(ncol(x), nrow(x), nrow(x)))
    t <- rowSums(x)
    p <- x/t
    for (k in 1:ncol(x)) {
        for (i in 1:nrow(x)) pij[k, , i] <- p[i, k]
        for (i in 1:nrow(x)) pij[k, i, ] <- abs(pij[k, i, ] - p[i, k])
        matprovi <- c %*% pij[k, , ]
        matprovi <- matprovi * diag(nrow(x))
        result[k] <- IS[k] - sum(matprovi)/sum(c)
    }
    return(round(result, 4))
}


#' A function to compute Wong's segregation index 
#'
#' @usage ISWong(x, b = NULL,  a = NULL, p = NULL, ptype = 'int', variant = 's', 
#' spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' totals because this will be interpreted as a population group
#' @param b - a common boundaries matrix where each element \emph{Bij} 
#' equals the shared boundary of \emph{i}-th and \emph{j}-th spatial units.
#' @param p - a numeric vector containing the perimeters of spatial units
#' @param ptype - a string variable giving two options for perimeter calculation
#' when a spatial object or shapefile is provided: 'int' to use only interior
#' boundaries of spatial units, and 'all' to use entire boundaries, 
#' including the boundaries to the exterior of the area
#' @param a - a numeric vector containing the areas of spatial units
#' @param variant - a character variable that allows to choose the index version: 
#' variant = 's' for the index adjusted for contiguous spatial units 
#' boundary lengths and perimeter/area ratio (by default) and variant = 'w' 
#' for the version based only on shared boundaries length
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension) which contains the geographic information
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return a numeric vector containing the Wong's segregation index value for 
#' each population group
#' @examples x <- segdata@data[ ,1:2]
#' bound <- boundaries(segdata)
#' per <- perimeter(segdata)
#' ar <- area(segdata)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'segdata'
#' 
#' ISWong(x, b = bound, p = per, a = ar) 
#' 
#' ISWong(x, spatobj = segdata, variant = 's', ptype = 'int')
#' 
#' ISWong(x, folder = foldername, shape = shapename, variant = 'w') 
#' 
#' @references Wong D. W. S. (1998) \emph{Measuring multiethnic spatial 
#' segregation}. Urban Geography, 19 (1), pp. 77-87.
#' @description Wong's segregation index is a development of 
#' \code{\link{ISDuncan}}'s which takes into account the interactions 
#' between spatial units (common boundaries and perimeter/area ratio). 
#' The function can be used in two ways: to provide spatial data ( 
#' boundaries matrix, a perimeter vector and an area vector) 
#' or a external geographic information source (spatial object or shape file).
#' @seealso One-group evenness indices: 
#' \code{\link{ISDuncan}}, \code{\link{Gini}}, \code{\link{Gorard}}, 
#' \code{\link{HTheil}}, '\code{\link{Atkinson}}, \code{\link{ISMorrill}}
#' @seealso Between groups dissimilarity indices: 
#' \code{\link{DIDuncan}}, \code{\link{Gini2}}, 
#' \code{\link{DIMorrill}}, \code{\link{DIWong}}
#' @export



ISWong <- function(x, b = NULL, a = NULL, p = NULL, ptype = "int", variant = "s", spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(b)) 
        b <- boundaries(spatobj = spatobj, folder = folder, shape = shape)
    IS <- ISDuncan(x)
    result <- vector(length = ncol(x))
    pij <- array(0, dim = c(ncol(x), nrow(x), nrow(x)))
    pp <- x/rowSums(x)
    w <- b/sum(b)
    # w <- b/ rowSums (b) w <- w/sum(w)
    if (variant == "w") 
        for (k in 1:ncol(x)) {
            for (i in 1:nrow(x)) pij[k, , i] <- pp[i, k]
            for (i in 1:nrow(x)) pij[k, i, ] <- abs(pij[k, i, ] - pp[i, k])
            matprovi <- w %*% (pij[k, , ])
            matprovi <- matprovi * diag(nrow(x))
            result[k] <- IS[k] - sum(matprovi)
        }
    if (variant == "s") {
        if (is.null(a)) 
            a <- area(spatobj = spatobj, folder = folder, shape = shape)
        if (is.null(p)) {
            if (ptype == "all") 
                p <- perimeter(spatobj = spatobj, folder = folder, shape = shape)
            if (ptype == "int") 
                p <- rowSums(b)
        }
        PerAij <- matrix(data = 0, nrow = nrow(x), ncol = nrow(x))
        for (i in 1:nrow(x)) PerAij[, i] <- p[i]/a[i]
        for (i in 1:nrow(x)) PerAij[i, ] <- PerAij[i, ] + p[i]/a[i]
        maxPA <- max(p/a)
        for (k in 1:ncol(x)) {
            for (i in 1:nrow(x)) pij[k, , i] <- pp[i, k]
            for (i in 1:nrow(x)) pij[k, i, ] <- abs(pij[k, i, ] - pp[i, k])
            for (i in 1:nrow(x)) PerAij[, i] <- p[i]/a[i]
            for (i in 1:nrow(x)) PerAij[i, ] <- PerAij[i, ] + p[i]/a[i]
            matprovi <- w %*% (pij[k, , ] * PerAij)
            matprovi <- matprovi * diag(nrow(x))
            result[k] <- IS[k] - sum(matprovi)/(2 * maxPA)
        }
    }
    return(round(result, 4))
}


#' A function to compute Gorard's segregation index 
#'
#' @usage Gorard(x)
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' totals because this will be interpreted as a population group
#' @return a numeric vector containing the Gorard's segregation index value for 
#' each population group
#' @examples x <- segdata@data[ ,1:2]
#' Gorard(x)
#' @references Gorard S. (2000) \emph{Education and Social Justice}. 
#' Cardiff, University of Wales Press
#' @description Gorard's index is an alternative to \code{\link{ISDuncan}}'s 
#' index, which measures the dissimilarity between the spatial  
#' distribution of a group and the total population. 
#' @seealso One-group evenness indices: 
#' \code{\link{ISDuncan}}, \code{\link{Gini}}, \code{\link{Atkinson}}, 
#' \code{\link{HTheil}}, '\code{\link{ISWong}}, \code{\link{ISMorrill}}
#' @seealso Between groups dissimilarity indices: 
#' \code{\link{DIDuncan}}, \code{\link{Gini2}}, 
#' \code{\link{DIMorrill}}, \code{\link{DIWong}}
#' @export


Gorard <- function(x) {
    x <- as.matrix(x)
    result <- vector(length = ncol(x))
    tx <- rowSums(x)
    varTotal <- colSums(x)
    Total <- sum(x)
    for (k in 1:ncol(x)) result[k] <- 0.5 * sum(abs(x[, k]/varTotal[k] - tx/Total))
    return(round(result, 4))
}

#' A function to compute Spatial Gini's segregation index 
#'
#' @usage Gini(x)
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @return a numeric vector containing the Gini's segregation index value for 
#' each population group
#' @examples x <- segdata@data[ ,1:2]
#' Gini(x)
#' @references Duncan O. D. and Duncan B. (1955) \emph{A Methodological 
#' Analysis of Segregation Indexes}. American Sociological Review 41, 
#' pp. 210-217
#' @description The spatial version of the Gini index can be derived from 
#' the Lorenz curve as the area between the segregation curve and the 
#' diagonal. 
#' @seealso Other one-group  evenness indices: 
#' \code{\link{ISDuncan}}, \code{\link{Atkinson}}, \code{\link{Gorard}}, 
#' \code{\link{HTheil}}, '\code{\link{ISWong}}, \code{\link{ISMorrill}}
#' @seealso Between groups dissimilarity indices: 
#' \code{\link{DIDuncan}}, \code{\link{Gini2}}, 
#' \code{\link{DIMorrill}}, \code{\link{DIWong}}
#' @export


Gini <- function(x) {
    x <- as.matrix(x)
    result <- vector(length = ncol(x))
    pij <- array(0, dim = c(ncol(x), nrow(x), nrow(x)))
    t <- rowSums(x)
    p <- x/t
    Total <- sum(x)
    pTotal <- colSums(x)/Total
    for (k in 1:ncol(x)) {
        for (i in 1:nrow(x)) pij[k, , i] <- p[i, k]
        for (i in 1:nrow(x)) pij[k, i, ] <- abs(pij[k, i, ] - p[i, k])
        matprovi <- (t %*% t(t)) * pij[k, , ]
        result[k] <- sum(matprovi)/(2 * Total * Total * pTotal[k] * (1 - pTotal[k]))
    }
    return(round(result, 4))
}


#' A function to compute Spatial Gini's between group index index 
#'
#' @usage Gini2(x)
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @return a matrix with between group Gini index 
#' @examples x <- segdata@data[ ,1:2]
#' Gini(x)
#' @references Duncan O. D. and Duncan B. (1955) \emph{A Methodological 
#' Analysis of Segregation Indexes}. American Sociological Review 41, 
#' pp. 210-217
#' @description The between group version of Gini index is obtained 
#' by computing the index for a subpopulation formed by two groups 
#' @seealso Other one-group  evenness indices: 
#' \code{\link{ISDuncan}}, \code{\link{Gini}}, 
#' \code{\link{Gorard}}, \code{\link{Atkinson}}, 
#' \code{\link{HTheil}}, '\code{\link{ISWong}}, \code{\link{ISMorrill}}
#' @seealso Between groups dissimilarity indices: 
#' \code{\link{DIDuncan}}, \code{\link{DIMorrill}}, \code{\link{DIWong}}
#' @export

Gini2 <- function(x) {
    x <- as.matrix(x)
    result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    for (k1 in 1:(ncol(x) - 1)) for (k2 in (k1 + 1):ncol(x)) {
        xprovi <- x[, c(k1, k2)]
        xprovi <- xprovi[rowSums(xprovi)>0,]
        result[k1, k2] <- Gini(xprovi)[1]
        result[k2, k1] <- result[k1, k2]
    }
    return(round(result, 4))
}


#' A function to compute Duncan dissimilarity segregation index
#'
#' @usage DIDuncan(x) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @return a matrix containing dissimilarity index values
#' @references Duncan O. D. and Duncan B. (1955) \emph{A Methodological 
#' Analysis of Segregation Indexes}. American Sociological Review 41, 
#' pp. 210-217
#' @description Duncan's dissimilarity index is the segregation index 
#' most commonly used in the literature. It is derived from Lorenz 
#' curves as the maximum difference between the segregation curve 
#' and the diagonal. The index measures the unevenness of a group's 
#' spatial distribution compared to another group. It can be 
#' interpreted as the share of the group that would have to move to 
#' achieve an even distribution compared to another group.
#' @examples x <- segdata@data[ ,1:2]
#' DIDuncan(x) 
#' @seealso Other one-group  evenness indices: 
#' \code{\link{ISDuncan}}, \code{\link{Gini}}, \code{\link{Gorard}}, 
#' \code{\link{Atkinson}}, \code{\link{HTheil}}, 
#' \code{\link{ISWong}}, \code{\link{ISMorrill}}
#' @seealso Between groups dissimilarity indices: 
#' \code{\link{DIMorrill}}, \code{\link{DIWong}}
#' @export

DIDuncan <- function(x) {
    x <- as.matrix(x)
    result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    pij <- array(0, dim = c(ncol(x), nrow(x), nrow(x)))
    varTotal <- colSums(x)
    for (k1 in 1:(ncol(x) - 1)) for (k2 in (k1 + 1):ncol(x)) {
        result[k1, k2] <- 0.5 * sum(abs(x[, k1]/varTotal[k1] - x[, k2]/varTotal[k2]))
        result[k2, k1] <- result[k1, k2]
    }
    return(round(result, 4))
}


#' A function to compute Morrill's dissimilarity index
#'
#' @usage DIMorrill(x, c = NULL, queen = TRUE, spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param c - a standard binary contiguity (adjacency) symmetric matrix where 
#' each element \emph{Cij} equals 1 if \emph{i}-th and \emph{j}-th spatial 
#' units are adjacent, and 0 otherwise.
#' @param queen - logical parameter difining criteria used for contiguity 
#' matrix computation, TRUE (by default) for queen , FALSE for rook 
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension)  which contains the geographic information
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return a matrix with Morrill's dissimilarity index values 
#' @references Morrill B. (1991) \emph{On the measure of geographic 
#' segregation}. Geography research forum, 11, pp. 25-36.
#' @description Morrill's dissimilarity index is a development of 
#' \code{\link{DIDuncan}}'s index which takes into account the 
#' interactions between spatial units (contiguity). 
#' The function can be used in two ways: to provide a contiguity 
#' matrix or a external geographic information source (spatial object 
#' or shape file).
#' @examples x <- segdata@data[ ,1:2]
#' contiguity <- contig(segdata)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'segdata'
#' 
#' DIMorrill(x, c = contiguity) 
#' 
#' DIMorrill(x, spatobj = segdata, queen = FALSE)
#' 
#' DIMorrill(x, folder = foldername, shape = shapename) 
#' @seealso Other one-group  evenness indices: 
#' \code{\link{ISDuncan}}, \code{\link{Gini}}, \code{\link{Gorard}}, 
#' \code{\link{Atkinson}}, \code{\link{HTheil}}, 
#' \code{\link{ISWong}}, \code{\link{ISMorrill}}
#' @seealso Between groups dissimilarity indices: 
#' \code{\link{DIDuncan}}, \code{\link{DIWong}}
#' @export


DIMorrill <- function(x, c = NULL, queen = TRUE, spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    if (is.null(c)) 
        c <- contig(spatobj = spatobj, folder = folder, shape = shape, queen = queen)
    DI <- DIDuncan(x)
    pij <- array(0, dim = c(ncol(x), nrow(x), nrow(x)))
    for (k1 in 1:(ncol(x) - 1)) for (k2 in (k1 + 1):ncol(x)) {
        for (i in 1:nrow(x)) pij[k1, , i] <- x[i, k1]/(x[i, k1] + x[i, k2])
        for (i in 1:nrow(x)) pij[k1, i, ] <- abs(pij[k1, i, ] - x[i, k1]/(x[i, k1] + x[i, k2]))
        pij[k1, , ][is.nan(pij[k1, , ])] <- 0
        matprovi <- c %*% pij[k1, , ]
        matprovi <- matprovi * diag(nrow(x))
        result[k1, k2] <- DI[k1, k2] - sum(matprovi)/sum(c)
        result[k2, k1] <- result[k1, k2]
    }
    return(round(result, 4))
}




#' A function to compute Wongs's dissimilarity index
#'
#' @usage DIWong(x, b = NULL,  a = NULL, p = NULL, ptype = 'int', variant = 's', 
#' spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' totals because this will be interpreted as a population group
#' @param b - a common boundaries matrix where each element \emph{Bij} 
#' equals the shared boundary of \emph{i}-th and \emph{j}-th spatial units.
#' @param p - a numeric vector containing the perimeters of spatial units
#' @param ptype - a string variable giving two options for perimeter calculation
#' when a spatial object or shapefile is provided: 'int' to use only interior
#' borders of spatial units, and 'all' to use entire borders, including to
#' the exterior of the area
#' @param a - a numeric vector containing the areas of spatial units
#' @param variant - a character variable that allows to choose the index version: 
#' variant = 's' for the dissimilarity index adjusted for contiguous spatial units 
#' boundary lengths and perimeter/area ratio (by default) and variant = 'w' 
#' for the version without perimeter/area ratio
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension) which contains the geographic information
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return a matrix containing Wong's dissimilarity index values 
#' @references Wong D. W. S. (1993) \emph{Spatial Indices of Segregation}. 
#' Urban Studies, 30 (3), pp. 559-572.
#' @description Wong's dissimilarity index is a development of 
#' \code{\link{DIDuncan}}'s which takes into account the interactions 
#' between spatial units (common boundaries and perimeter/area ratios). 
#' The function can be used in two ways: to provide spatial data ( 
#' boundaries matrix, a perimeter vector and an area vector) 
#' or a external geographic information source (spatial object or shape file).
#' @examples x <- segdata@data[ ,1:2]
#' bound <- boundaries(segdata)
#' per <- perimeter(segdata)
#' ar <- area(segdata)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'segdata'
#' 
#' DIWong(x, b = bound, p = per, a = ar) 
#' 
#' DIWong(x, spatobj = segdata, variant = 'w') 
#' 
#' DIWong(x, folder = foldername, shape = shapename, ptype ='all') 
#' @seealso Other one-group  evenness indices: 
#' \code{\link{ISDuncan}}, \code{\link{Gini}}, \code{\link{Gorard}}, 
#' \code{\link{Atkinson}}, \code{\link{HTheil}}, 
#' '\code{\link{ISWong}}, \code{\link{ISMorrill}}
#' @seealso Between groups dissimilarity indices: 
#' \code{\link{DIDuncan}}, \code{\link{DIMorrill}}
#' @export


DIWong <- function(x, b = NULL, a = NULL, p = NULL, ptype = "int", variant = "s", spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(b)) 
        b <- boundaries(spatobj = spatobj, folder = folder, shape = shape)
    result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    DI <- DIDuncan(x)
    pij <- array(0, dim = c(ncol(x), nrow(x), nrow(x)))
    w <- b/sum(b)
    # w <- b/ rowSums (b) w <- w/sum(w)
    if (variant == "w") 
        for (k1 in 1:(ncol(x) - 1)) for (k2 in (k1 + 1):ncol(x)) {
            for (i in 1:nrow(x)) pij[k1, , i] <- x[i, k1]/(x[i, k1] + x[i, k2])
            for (i in 1:nrow(x)) pij[k1, i, ] <- abs(pij[k1, i, ] - x[i, k1]/(x[i, k1] + x[i, k2]))
            pij[k1, , ][is.nan(pij[k1, , ])] <- 0
            matprovi <- w %*% (pij[k1, , ])
            matprovi <- matprovi * diag(nrow(x))
            result[k1, k2] <- DI[k1, k2] - sum(matprovi)
            result[k2, k1] <- result[k1, k2]
        }
    if (variant == "s") {
        if (is.null(a)) 
            a <- area(spatobj = spatobj, folder = folder, shape = shape)
        if (is.null(p)) {
            if (ptype == "all") 
                p <- perimeter(spatobj = spatobj, folder = folder, shape = shape)
            if (ptype == "int") 
                p <- rowSums(b)
        }
        PerAij <- matrix(data = 0, nrow = nrow(x), ncol = nrow(x))
        for (i in 1:nrow(x)) PerAij[, i] <- p[i]/a[i]
        for (i in 1:nrow(x)) PerAij[i, ] <- PerAij[i, ] + p[i]/a[i]
        maxPA <- max(p/a)
        for (k1 in 1:(ncol(x) - 1)) for (k2 in (k1 + 1):ncol(x)) {
            for (i in 1:nrow(x)) pij[k1, , i] <- x[i, k1]/(x[i, k1] + x[i, k2])
            for (i in 1:nrow(x)) pij[k1, i, ] <- abs(pij[k1, i, ] - x[i, k1]/(x[i, k1] + x[i, k2]))
            pij[k1, , ][is.nan(pij[k1, , ])] <- 0
            matprovi <- w %*% (pij[k1, , ] * PerAij)
            matprovi <- matprovi * diag(nrow(x))
            result[k1, k2] <- DI[k1, k2] - sum(matprovi)/(2 * maxPA)
            result[k2, k1] <- result[k1, k2]
        }
    }
    return(round(result, 4))
}


#################### 

# EXPOSITION INDEXES

#################### 

#' A function to compute Bell's isolation index (xPx)
#'
#' @usage xPx(x, exact = FALSE) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param exact - a logical variable to specifiy the index version: 
#' exact = FALSE (by default) for the approximate version of the index, 
#' and exact = TRUE for the exact version
#' @return a numeric vector containing the isolation index value for 
#' each population group
#' @references Bell W. (1954) \emph{A probability model for the 
#' measurement of ecological segregation}. Social Forces 32(4), 
#' pp. 357-364
#' @description The isolation index, xPx, is an exposure index 
#' that measures the probability that two members of a group share 
#' the same spatial unit. This index can be calculated using the 
#' approximate or the exact method (see Bell, 1954).
#' @examples x <- segdata@data[ ,7:8]
#' xPx(x) 
#' xPx(x, exact = TRUE)
#' @seealso Isolation indices: 
#' \code{\link{Eta2}},  \code{\link{DPxx}}
#' @seealso Interaction indices: 
#' \code{\link{xPy}}, \code{\link{DPxy}}
#' @export

xPx <- function(x, exact = FALSE) {
    x <- as.matrix(x)
    result <- vector(length = ncol(x))
    varTotal <- colSums(x)
    t <- rowSums(x)
    if (exact == FALSE) 
        for (k in 1:ncol(x)) result[k] <- sum(x[, k]/varTotal[k] * x[, k]/t)
    if (exact == TRUE) 
        for (k in 1:ncol(x)) result[k] <- sum(x[, k]/varTotal[k] * (x[, k] - 1)/(t - 1))
    return(round(result, 4))
}


#' A function to compute the distance-decay isolation index (DPxx)
#'
#' @usage DPxx(x, d = NULL, distin = 'm',  distout = 'm', diagval = '0', 
#' spatobj = NULL, folder = NULL, shape = NULL) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param d - a matrix of the distances between spatial unit centroids
#' @param distin - input metric conversion, based on  \pkg{bink} package and 
#' includes conversions from 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.
#' @param distout - output metric conversion, based on  \pkg{bink} package and 
#' includes conversions to 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.
#' @param diagval - when providing a spatial object or a shape file, 
#' the user has the choice of the spatial matrix diagonal definition: 
#' diagval = '0' (by default) for an null diagonal and diagval = 'a' 
#' to compute the diagonal as 0.6 * square root (spatial units area) (White, 1983) 
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension) which contains the geographic information
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return a numeric vector containing the distance-decay isolation index 
#' value for each population group
#' @references Morgan, B. S. (1983) \emph{A Distance-Decay Based Interaction 
#' Index to Measure Residential Segregation}. Area 15(3),  pp. 211-217.
#' @description The distance decay isolation index, DPxx, is a spatial
#' adaptation of isolation index \code{\link{xPx}}. 
#' The function can be used in two ways: to provide a distance matrix 
#' or a external geographic information source (spatial object or shape file).
#' @examples x <- segdata@data[ ,1:2]
#' ar <- area(segdata)
#' dist <- distance(segdata)
#' diag(dist)<-sqrt(ar) * 0.6
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'segdata'
#' 
#' DPxx(x, d = dist)
#' 
#' DPxx(x, spatobj = segdata, diagval = 'a')
#' 
#' DPxx(x, folder = foldername, shape = shapename, diagval = '0') 
#' @seealso Isolation indices: 
#' \code{\link{xPx}},  \code{\link{Eta2}}
#' @seealso Interaction indices: 
#' \code{\link{xPy}}, \code{\link{DPxy}}
#' @export


DPxx <- function(x, d = NULL, distin = "m", distout = "m", diagval = "0", spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    result <- vector(length = ncol(x))
    varTotal <- colSums(x)
    t <- rowSums(x)
    if (is.null(d)) 
        d <- distance(spatobj = spatobj, folder = folder, shape = shape, distin = distin, distout = distout, diagval = diagval)
    dd <- exp(-d)
    K1 <- dd %*% t
    K2 <- dd * t
    K <- K2/as.vector(K1)
    for (k in 1:ncol(x)) {
        X1 <- K %*% (x[, k]/t)
        X2 <- (x[, k]/varTotal[k]) * X1
        result[k] <- sum(X2)
    }
    return(round(result, 4))
}



#' A function to compute adjusted isolation index (Eta2)
#'
#' @usage Eta2(x) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @return a numeric vector containing the adjusted isolation index value for 
#' each population group
#' @references Bell W. (1954) \emph{A probability model for the 
#' measurement of ecological segregation}. Social Forces 32(4), 
#' pp. 357-364
#' @references Duncan O. D. and Duncan B. (1955) \emph{ 
#' Residential Distribution and Occupational Stratification.}. 
#' American Journal of Sociology 60 (5), pp. 493-503
#' @description The adjusted isolation index is the standardized 
#' version of the isolation index, \code{\link{xPx}}, which 
#' controls for the effect of total population structure. Using 
#' the approximate version of xPx, the adjusted index is equal 
#' to Eta2 (the square of the correlation ratio) which, in the 
#' case of the binomial variable, is identical to the square of 
#' the mean square contingency coefficient phi. It can be used 
#' as a segregation score and varies from 0 (minimum segregation) 
#' to 1 (maximum segregation).
#' @examples x <- segdata@data[ ,1:2]
#' Eta2(x) 
#' @seealso Isolation indices: 
#' \code{\link{xPx}},  \code{\link{DPxx}}
#' @seealso Interaction indices: 
#' \code{\link{xPy}}, \code{\link{DPxy}}
#' @export



Eta2 <- function(x) {
    x <- as.matrix(x)
    result <- vector(length = ncol(x))
    varTotal <- colSums(x)
    Total <- sum(x)
    pTotal <- colSums(x)/Total
    t <- rowSums(x)
    
    for (k in 1:ncol(x)) {
        result[k] <- sum(x[, k]/varTotal[k] * x[, k]/t)
        result[k] <- (result[k] - pTotal[k])/(1 - pTotal[k])
    }
    return(round(result, 4))
}

#' A function to compute interaction index (xPy)
#'
#' @usage xPy(x, exact = FALSE) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param exact - a logical variable to specifiy the index version: 
#' exact = FALSE (by default) for the approximate version of the index, 
#' and exact = TRUE for the exact version
#' @return a matrix with interaction index values
#' @references Bell W. (1954) \emph{A probability model for the 
#' measurement of ecological segregation}. Social Forces 32(4),
#'  pp. 357-364
#' @description The interaction index, xPy, is an exposure 
#' between groups index which measures the probability that a member 
#' of a group shares the same spatial unit with a member of another 
#' group. The index can be calculated with the approximate or exact 
#' method (see Bell, 1954).
#' @examples x <- segdata@data[ ,1:2]
#' xPy(x) 
#' @seealso Isolation indices: 
#' \code{\link{xPx}}, \code{\link{Eta2}},  \code{\link{DPxx}}
#' @seealso Distance decay interaction index: \code{\link{DPxy}}
#' @export



xPy <- function(x, exact = FALSE) {
    x <- as.matrix(x)
    result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    varTotal <- colSums(x)
    t <- rowSums(x)
    if (exact == FALSE) 
        for (k1 in 1:ncol(x)) for (k2 in 1:ncol(x)) result[k1, k2] <- sum(x[, k1]/varTotal[k1] * x[, k2]/t)
    if (exact == T) 
        for (k1 in 1:ncol(x)) for (k2 in 1:ncol(x)) if (k1 != k2) {
            result[k1, k2] <- sum(x[, k1]/varTotal[k1] * x[, k2]/(t - 1))
        } else {
            result[k1, k2] <- sum(x[, k1]/varTotal[k1] * (x[, k2] - 1)/(t - 1))
        }
    
    return(round(result, 4))
}





#' A function to compute the distance-decay interaction index (DPxy)
#'
#' @usage DPxy(x, d = NULL, distin = 'm',  distout = 'm', diagval = '0', 
#' spatobj = NULL, folder = NULL, shape = NULL) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param d - a matrix of the distances between spatial unit centroids
#' @param distin - input metric conversion, based on  \pkg{bink} package and 
#' includes conversions from 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.
#' @param distout - output metric conversion, based on  \pkg{bink} package and 
#' includes conversions to 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.
#' @param diagval - when providing a spatial object or a shape file, 
#' the user has the choice of the spatial matrix diagonal definition: 
#' diagval = '0' (by default) for an null diagonal and diagval = 'a' 
#' to compute the diagonal as 0.6 * square root (spatial units area) (White, 1983) 
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension) which contains the geographic information
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return a numeric vector containing the distance-decay isolation index 
#' value for each population group
#' @references Morgan, B. S. (1983) \emph{An Alternate Approach to the 
#' Development of a Distance-Based Measure of Racial Segregation}. 
#' American Journal of Sociology 88,  pp. 1237-1249.
#' @description The distance decay interaction index, DPxy, is a 
#' spatial adaptation of interaction index \code{\link{xPy}}. 
#' The function can be used in two ways: to provide a distance matrix 
#' or a external geographic information source (spatial object or shape file).
#' @examples x <- segdata@data[ ,1:2]
#' ar <- area(segdata)
#' dist <- distance(segdata)
#' diag(dist)<-sqrt(ar) * 0.6
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'segdata'
#' 
#' DPxy(x, d = dist)
#' 
#' DPxy(x, spatobj = segdata, diagval = 'a')
#' 
#' DPxy(x, folder = foldername, shape = shapename, diagval = '0') 
#' @seealso Isolation indices: 
#' \code{\link{xPx}}, \code{\link{Eta2}},  \code{\link{DPxx}}
#' @seealso Interaction index: \code{\link{xPy}}
#' @export


DPxy <- function(x, d = NULL, distin = "m", distout = "m", diagval = "0", spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    t <- rowSums(x)
    varTotal <- colSums(x)
    if (is.null(d)) 
        d <- distance(spatobj = spatobj, folder = folder, shape = shape, distin = distin, distout = distout, diagval = diagval)
    result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    varTotal <- colSums(x)
    dd <- exp(-d)
    K1 <- dd %*% t
    K2 <- dd * t
    K <- K2/as.vector(K1)
    
    for (k1 in 1:ncol(x)) for (k2 in 1:ncol(x)) {
        X1 <- K %*% (x[, k2]/t)
        X2 <- (x[, k1]/varTotal[k1]) * X1
        result[k1, k2] <- sum(X2)
    }
    return(round(result, 4))
}



####################### 

# CONCENTRATION INDEXES

####################### 


#' A function to compute Delta index
#'
#' @usage Delta(x, a = NULL, spatobj = NULL, folder = NULL, shape = NULL) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param a - a numeric vector containing the areas of spatial units
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension) which contains the geographic information
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return a numeric vector containing Delta index value for 
#' each population group
#' @references Duncan O. D., Cuzzoert  and Duncan B. (1961) 
#' \emph{Problems in analyzing areal data}. Statistical geography, 
#' Glencoe, Illinois: The free press of Glencoe
#' @description The Delta index is a specific application of dissimilarity
#'  index \code{\link{DIDuncan}} which simply measures the dissimilarity 
#' between the spatial distribution of a population group and the spatial 
#' distribution of available area. It can be interpreted as the share of group 
#' that would have to move to achieve uniform density over all spatial units. 
#' The function can be used in two ways: to provide an area vector or 
#' a external geographic information source (spatial object or shape file).
#' @examples x <- segdata@data[ ,1:2]
#' ar <- area(segdata)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'segdata'
#' 
#' Delta(x, a = ar) 
#' 
#' Delta(x, spatobj = segdata)
#' 
#' Delta(x, folder = foldername, shape = shapename) 
#' @seealso Absolute Concentration Index: \code{\link{ACO}}
#' @seealso Relative Concentration Index: \code{\link{RCO}}
#' @export


Delta <- function(x, a = NULL, spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(a)) 
        a <- area(spatobj = spatobj, folder = folder, shape = shape)
    result <- vector(length = ncol(x))
    varTotal <- colSums(x)
    areaTotal <- sum(a)
    for (k in 1:ncol(x)) result[k] <- 0.5 * sum(abs(x[, k]/varTotal[k] - a/areaTotal))
    return(round(result, 4))
}



#' A function to compute Absolute Concentration index (ACO)
#'
#' @usage ACO(x, a = NULL, spatobj = NULL, folder = NULL, shape = NULL) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param a - a numeric vector containing the areas of spatial units
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension) which contains the geographic information
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return a numeric vector containing Absolute Concentration index value for 
#' each population group
#' @references Massey D. S. and Denton N. A. (1988) \emph{
#' The dimensions of residential segregation}. 
#' Social Forces 67(2),  pp. 281-315.
#' @description The absolute concentration index, ACO, computes 
#' the total area inhabited by a group, and compares the result 
#' to the minimum and maximum possible areas that could be 
#' inhabited by that group in the study area. 
#' #' The function can be used in two ways: to provide an area vector or a 
#' external geographic information source (spatial object or shape file). 
#' @examples x <- GreHSize@data[ ,3:5]
#' ar <- area(GreHSize)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' 
#' ACO(x, a = ar) 
#' 
#' ACO(x, spatobj = GreHSize)
#' 
#' ACO(x, folder = foldername, shape = shapename) 
#' 
#' @seealso Delta Index: \code{\link{Delta}}
#' @seealso Relative Concentration Index: \code{\link{RCO}}
#' @export


ACO <- function(x, a = NULL, spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(a)) 
        a <- area(spatobj = spatobj, folder = folder, shape = shape)
    varTotal <- colSums(x)
    xprovi <- as.data.frame(cbind(x, a))
    xprovi <- xprovi[order(xprovi$a), ]
    xprovi$Total <- rowSums(xprovi) - xprovi$a
    areaTotal <- sum(a)
    result <- vector(length = ncol(x))
    n1 <- vector(length = ncol(x))
    n2 <- vector(length = ncol(x))
    T1 <- vector(length = ncol(x))
    T2 <- vector(length = ncol(x))
    # t <- rowSums(xprovi)
    
    for (k in 1:ncol(x)) {
        T1[k] <- 0
        i <- 0
        while (T1[k] < varTotal[k]) {
            i <- i + 1
            T1[k] <- T1[k] + xprovi$Total[i]
        }
        n1[k] <- i
        T2[k] <- 0
        i <- nrow(xprovi) + 1
        while (T2[k] < varTotal[k]) {
            i <- i - 1
            T2[k] <- T2[k] + xprovi$Total[i]
        }
        n2[k] <- i
    }
    
    for (k in 1:ncol(x)) {
        vartemp1 <- sum(xprovi[, k] * xprovi$a/varTotal[k])
        vartemp2 <- sum(xprovi$Total[1:n1[k]] * xprovi$a[1:n1[k]]/T1[k])
        vartemp3 <- sum(xprovi$Total[n2[k]:nrow(xprovi)] * xprovi$a[n2[k]:nrow(xprovi)]/T2[k])
        result[k] <- 1 - (vartemp1 - vartemp2)/(vartemp3 - vartemp2)
    }
    
    return(round(result, 4))
}


#' A function to compute Relative Concentration index (RCO)
#'
#' @usage RCO(x, a = NULL, spatobj = NULL, folder = NULL, shape = NULL) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param a - a numeric vector containing the areas of spatial units
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension) which contains the geographic information
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return a matrix containing relative concentration index values
#' @references Massey D. S. and Denton N. A. (1988) \emph{
#' The dimensions of residential segregation}. 
#' Social Forces 67(2),  pp. 281-315.
#' @description The relative concentration index, RCO, measures 
#' the share of space occupied by a group compared to another group.
#' The function can be used in two ways: to provide an area vector or a 
#' external geographic information source (spatial object or shape file).
#' @examples x <- GreHSize@data[ ,3:5]
#' ar <- area(GreHSize)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'GreHSize'
#' 
#' RCO(x, a = ar) 
#' 
#' RCO(x, spatobj = GreHSize)
#' 
#' RCO(x, folder = foldername, shape = shapename) 
#' 
#' @seealso one-group concentration indices: 
#' \code{\link{Delta}},  \code{\link{ACO}}
#' @export


RCO <- function(x, a = NULL, spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(a)) 
        a <- area(spatobj = spatobj, folder = folder, shape = shape)
    result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    varTotal <- colSums(x)
    xprovi <- as.data.frame(cbind(x, a))
    xprovi <- xprovi[order(xprovi$a), ]
    xprovi$Total <- rowSums(xprovi) - xprovi$a
    areaTotal <- sum(a)
    n1 <- vector(length = ncol(x))
    n2 <- vector(length = ncol(x))
    T1 <- vector(length = ncol(x))
    T2 <- vector(length = ncol(x))
    for (k in 1:ncol(x)) {
        T1[k] <- 0
        i <- 0
        while (T1[k] < varTotal[k]) {
            i <- i + 1
            T1[k] <- T1[k] + xprovi$Total[i]
        }
        n1[k] <- i
        T2[k] <- 0
        i <- nrow(xprovi) + 1
        while (T2[k] < varTotal[k]) {
            i <- i - 1
            T2[k] <- T2[k] + xprovi$Total[i]
        }
        n2[k] <- i
    }
    for (k1 in 1:ncol(x)) for (k2 in 1:ncol(x)) {
        vartemp1 <- sum(xprovi[, k1] * xprovi$a/varTotal[k1])
        vartemp1b <- sum(xprovi[, k2] * xprovi$a/varTotal[k2])
        vartemp2 <- sum(xprovi$Total[1:n1[k1]] * xprovi$a[1:n1[k1]]/T1[k1])
        vartemp3 <- sum(xprovi$Total[n2[k2]:nrow(xprovi)] * xprovi$a[n2[k2]:nrow(xprovi)]/T2[k2])
        result[k1, k2] <- ((vartemp1/vartemp1b) - 1)/((vartemp2/vartemp3) - 1)
    }
    return(round(result, 4))
}

####################### 

# CLUSTERING INDEXES

####################### 


#' A function to compute Absolute Clustering Index (ACL)
#'
#' @usage ACL(x, spatmat = 'c', c = NULL, queen = TRUE, distin = 'm',  
#' distout = 'm', diagval = '0',  spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param spatmat - the method used for spatial calculations: 'c' for the 
#' contiguity matrix (by default) or any other user spatial interaction matrix 
#' and 'd' for the inverse exponential function of the distance. 
#' @param c - a modified binary contiguity (adjacency) symmetric matrix where 
#' each element \emph{Cij} equals 1 if \emph{i}-th and \emph{j}-th spatial 
#' units are adjacent or identical, and 0 otherwise.
#' @param queen - logical parameter difining criteria used for contiguity 
#' matrix computation, TRUE (by default) for queen , FALSE for rook 
#' @param distin - input metric conversion, based on  \pkg{bink} package and 
#' includes conversions from 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.
#' @param distout - output metric conversion, based on  \pkg{bink} package and 
#' includes conversions to 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.
#' @param diagval - when providing a spatial object or a shape file, 
#' the user has the choice of the spatial matrix diagonal definition: 
#' diagval = '0' (by default) for an null diagonal and diagval = 'a' 
#' to compute the diagonal as 0.6 * square root (spatial units area) (White, 1983) 
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension) which contains the geographic information
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return a numeric vector containing Absolute Clustering index value for 
#' each population group
#' @references Massey D. S. and Denton N. A. (1988) \emph{
#' The dimensions of residential segregation}. 
#' Social Forces 67(2),  pp. 281-315.
#' @description The absolute clustering index, ACL, expresses the 
#' average number of a population group's members in nearby spatial 
#' units, as a proportion of the total population in those spatial units. 
#' The spatial interactions can be expressed as a contiguity matrix 
#' (with diagonal equal to 1), as an inverse exponential function of the 
#' distance between spatial units centres (with diagonal equal to 0.6 of the 
#' square root of each spatial units area) or other user specified interaction matrix. 
#' The function can be used in two ways: to provide a spatial interactions matrix 
#' (a contiguity matrix or a distance matrix) or a external 
#' geographic information source (spatial object or shape file).
#' @examples x <- segdata@data[ ,1:2]
#' contiguity <- contig(segdata)
#' diag(contiguity) <- 1
#' ar<-area(segdata)
#' dist <- distance(segdata)
#' dist <- exp(-dist)
#' diag(dist)<-sqrt(ar) * 0.6
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'segdata'
#' 
#' ACL(x, c = contiguity) 
#' 
#' ACL(x, spatobj = segdata)
#' 
#' ACL(x, spatmat = 'd', folder = foldername, shape = shapename) 
#'  
#' ACL(x, spatobj = segdata, spatmat = 'd', diagval = 'a')
#' 
#' ACL(x, c = dist, spatmat = 'd')
#'
#' @seealso Proximity measures: \code{\link{Pxx}}, 
#' \code{\link{Pxy}}, \code{\link{Poo}}, \code{\link{SP}}
#' @seealso Relative Clustering Index: \code{\link{RCL}}
#' @export


ACL <- function(x, spatmat = "c", c = NULL, queen = TRUE, distin = "m", distout = "m", diagval = "0", spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (spatmat == "c" & is.null(c)) {
        c <- contig(spatobj = spatobj, folder = folder, shape = shape, queen = queen)
        diag(c) <- 1
    }
    if (spatmat == "d" & is.null(c)) 
        c <- distance(spatobj = spatobj, folder = folder, shape = shape, distin = distin, distout = distout, diagval = diagval)
    if (spatmat == "d" )  c <- exp(-c)
    result <- vector(length = ncol(x))
    varTotal <- colSums(x)
    t <- as.vector(rowSums(x))
    for (k in 1:ncol(x)) {
        vartemp1 <- sum((c %*% x[, k]) * x[, k]/varTotal[k])
        vartemp2 <- (sum(c) * varTotal[k])/(nrow(x)^2)
        vartemp3 <- sum((c %*% t) * x[, k]/varTotal[k])
        result[k] <- (vartemp1 - vartemp2)/(vartemp3 - vartemp2)
    }
    return(round(result, 4))
}


#' A function to compute the mean proximity between members of a group (Pxx)
#'
#' @usage Pxx(x, d = NULL, fdist = 'e', distin = 'm',  distout = 'm', diagval = '0', 
#' spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param d - a matrix of the distances between spatial unit centroids
#' @param fdist - the method used for distance interaction matrix: 
#' e' for inverse exponential function (by default) and 'l' for linear.
#' @param distin - input metric conversion, based on  \pkg{bink} package and 
#' includes conversions from 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.
#' @param distout - output metric conversion, based on  \pkg{bink} package and 
#' includes conversions to 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.
#' @param diagval - when providing a spatial object or a shape file, 
#' the user has the choice of the spatial matrix diagonal definition: 
#' diagval = '0' (by default) for an null diagonal and diagval = 'a' 
#' to compute the diagonal as 0.6 * square root (spatial units area) (White, 1983) 
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension) which contains the geographic information
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return a numeric vector containing Pxx index value for 
#' each population group
#' @references White M. J. (1983) \emph{The Measurement of Spatial 
#' Segregation}. American Journal of Sociology, 88, p. 1008-1019
#' @description  Mean proximity, Pxx, computes the mean distance 
#' between the members of a group. The distance matrix can be expressed as  
#' a linear or as an inverse exponential function of the distance between 
#' spatial units centroids.The function can be used in two ways: to providea 
#' a distance matrix  or a external geographic information source (spatial 
#' object or shape file).
#' @examples x <- segdata@data[ ,1:2]
#' ar<-area(segdata)
#' dist <- distance(segdata)
#' diag(dist)<-sqrt(ar) * 0.6
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'segdata'
#' 
#' Pxx(x, spatobj = segdata)
#' 
#' Pxx(x, folder = foldername, shape = shapename, fdist = 'l') 
#' 
#' Pxx(x, spatobj = segdata, diagval ='a')
#' 
#' Pxx(x, d = dist, fdist = 'e')
#' 
#' @seealso Proximity measures: 
#' \code{\link{Pxy}}, \code{\link{Poo}}, \code{\link{SP}}
#' @seealso Clustering Indices: 
#' \code{\link{ACL}}, \code{\link{RCL}}
#' @export



Pxx <- function(x, d = NULL, fdist = "e", distin = "m", distout = "m", diagval = "0", spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(d)) 
        d <- distance(spatobj = spatobj, folder = folder, shape = shape, distin = distin, distout = distout, diagval = diagval)
    if (fdist == "e") 
        d <- exp(-d)
    varTotal <- colSums(x)
    result <- vector(length = ncol(x))
    for (k in 1:ncol(x)) result[k] <- sum((d %*% x[, k]) * x[, k]/(varTotal[k])^2)
    return(round(result, 4))
}




#' A function to compute the mean proximity between 
#' persons without regard to group (Poo)
#'
#' @usage Poo(x, d = NULL, fdist = 'e', distin = 'm',  distout = 'm', diagval = '0', 
#' itype = 'multi', spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param d - a matrix of the distances between spatial unit centroids
#' @param fdist - the method used for distance interaction matrix: 
#' e' for inverse exponential function (by default) and 'l' for linear.
#' @param distin - input metric conversion, based on  \pkg{bink} package and 
#' includes conversions from 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.
#' @param distout - output metric conversion, based on  \pkg{bink} package and 
#' includes conversions to 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.
#' @param diagval - when providing a spatial object or a shape file, 
#' the user has the choice of the spatial matrix diagonal definition: 
#' diagval = '0' (by default) for an null diagonal and diagval = 'a' 
#' to compute the diagonal as 0.6 * square root (spatial units area) (White, 1983) 
#' @param itype - a character string defining the index type:
#' itype = 'multi' (by default) for the multigroup index (White, 1986)
#' or itype = 'between' for the between groups version (White, 1983)
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension) which contains the geographic information
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return a numeric vector containing Poo index value for 
#' each population group
#' @references White M. J. (1983) \emph{The Measurement of Spatial 
#' Segregation}. American Journal of Sociology, 88, p. 1008-1019
#' @references  White, M. J. (1986) \emph{Segregation and Diversity Measures 
#' in Population Distribution}E. Population Index 52(2): 198-221.
#' @description Mean proximity, Poo, computes the mean distance 
#' between the individuals in the area with no regard for group.
#' The function can be used in two ways: to provide a distance matrix 
#' or a external geographic information source (spatial object 
#' or shape file) 
#' @examples x <- segdata@data[ ,1:2]
#' ar<-area(segdata)
#' dist <- distance(segdata)
#' diag(dist)<-sqrt(ar) * 0.6
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'segdata'
#' 
#' Poo(x, spatobj = segdata)
#' 
#' Poo(x, folder = foldername, shape = shapename, fdist = 'l') 
#' 
#' Poo(x, spatobj = segdata, diagval ='a')
#' 
#' Poo(x, d = dist, fdist = 'e') 
#'
#' @seealso Proximity measures: \code{\link{Pxx}}, 
#' \code{\link{Pxy}},  \code{\link{SP}}
#' @seealso Clustering Indices: 
#' \code{\link{ACL}}, \code{\link{RCL}}
#' @export



Poo <- function(x, d = NULL, fdist = "e", distin = "m", distout = "m", diagval = "0", itype = "multi", spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    tx <- rowSums(x)
    varTotal <- colSums(x)
    if (is.null(d)) 
        d <- distance(spatobj = spatobj, folder = folder, shape = shape, distin = distin, distout = distout, diagval = diagval)
    if (fdist == "e") 
        d <- exp(-d)
    if (itype == "between") {
        result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
        for (k1 in 1:ncol(x)) for (k2 in 1:ncol(x)) result[k1, k2] <- sum((d %*% (x[, k1] + x[, k2])) * (x[, k1] + x[, k2])/(varTotal[k1] + varTotal[k2])^2)
    }
    if (itype == "multi") 
        result <- sum((d %*% tx) * tx/(sum(tx))^2)
    return(round(result, 4))
}

#' A function to compute the mean proximity between 
#' persons of different groups (Pxy)
#'
#' @usage Pxy(x, d = NULL, fdist = 'e', distin = 'm',  distout = 'm', diagval = '0', 
#' spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param d - a matrix of the distances between spatial unit centroids
#' @param fdist - the method used for distance interaction matrix: 
#' e' for inverse exponential function (by default) and 'l' for linear.
#' @param distin - input metric conversion, based on  \pkg{bink} package and 
#' includes conversions from 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.
#' @param distout - output metric conversion, based on  \pkg{bink} package and 
#' includes conversions to 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.
#' @param diagval - when providing a spatial object or a shape file, 
#' the user has the choice of the spatial matrix diagonal definition: 
#' diagval = '0' (by default) for an null diagonal and diagval = 'a' 
#' to compute the diagonal as 0.6 * square root (spatial units area) (White, 1983) 
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension) which contains the geographic information
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return a matrix containing Pxy index values for each pair of groups
#' @references White M. J. (1983) \emph{The Measurement of 
#' Spatial Segregation}. American Journal of Sociology, 88, p. 1008-1019
#' @description Mean proximity, Pxy, computes the mean distance 
#' between the members of different groups. 
#' The function can be used in two ways: to provide a distance matrix 
#' or a external geographic information source (spatial object 
#' or shape file).
#' @examples x <- segdata@data[ ,1:2]
#' ar<-area(segdata)
#' dist <- distance(segdata)
#' diag(dist)<-sqrt(ar) * 0.6
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'segdata'
#' 
#' Pxy(x, spatobj = segdata)
#' 
#' Pxy(x, folder = foldername, shape = shapename, fdist = 'l') 
#' 
#' Pxy(x, spatobj = segdata, diagval ='a')
#' 
#' Pxy(x, d = dist, fdist = 'e')
#'
#' @seealso Proximity measures: \code{\link{Pxx}}, 
#' \code{\link{Poo}},  \code{\link{SP}}
#' @seealso Clustering Indices: 
#' \code{\link{ACL}}, \code{\link{RCL}}
#' @export


Pxy <- function(x, d = NULL, fdist = "e", distin = "m", distout = "m", diagval = "0", spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(d)) 
        d <- distance(spatobj = spatobj, folder = folder, shape = shape, distin = distin, distout = distout, diagval = diagval)
    if (fdist == "e") 
        d <- exp(-d)
    varTotal <- colSums(x)
    result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    for (k1 in 1:ncol(x)) for (k2 in 1:ncol(x)) result[k1, k2] <- sum((d %*% (x[, k2])) * (x[, k1]))/varTotal[k1]/varTotal[k2]
    return(round(result, 4))
}


#' A function to compute the spatial proximity index (SP)
#'
#' @usage SP(x, d = NULL, fdist = 'e', distin = 'm',  distout = 'm', diagval = '0', 
#' itype = 'multi', spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param d - a matrix of the distances between spatial unit centroids
#' @param fdist - the method used for distance interaction matrix: 
#' e' for inverse exponential function (by default) and 'l' for linear.
#' @param distin - input metric conversion, based on  \pkg{bink} package and 
#' includes conversions from 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.
#' @param distout - output metric conversion, based on  \pkg{bink} package and 
#' includes conversions to 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.
#' @param diagval - when providing a spatial object or a shape file, 
#' the user has the choice of the spatial matrix diagonal definition: 
#' diagval = '0' (by default) for an null diagonal and diagval = 'a' 
#' to compute the diagonal as 0.6 * square root (spatial units area) (White, 1983) 
#' @param itype - a character string defining the index type:
#' itype = 'multi' (by default) for the multigroup index (White, 1986)
#' itype = 'between' for the between groups version (White, 1983) or
#' itype = 'one' for the one-group version (Apparicio et al, 2008)
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension) which contains the geographic information
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return a matrix containing spatial proximity index values for each pair of groups
#' @references White M. J. (1983) \emph{The Measurement of Spatial 
#' Segregation}. American Journal of Sociology, 88, p. 1008-1019.
#' @references  White, M. J. (1986) \emph{Segregation and Diversity Measures 
#' in Population Distribution}E. Population Index 52(2): 198-221.
#' @references  Apparicio, P., V. Petkevitch and M. Charron (2008): \emph{Segregation 
#' Analyzer: A C#.Net application for calculating residential segregation indices}, 
#' Cybergeo: European Journal of Geography, 414, 1-27.
#' @description The spatial proximity index, SP, compares the clustering 
#' level (mean proximity) of a group compared to another group. 
#' The function can be used in two ways: to provide a distance matrix 
#' or a external geographic information source (spatial object 
#' or shape file).
#' @examples x <- segdata@data[ ,1:2]
#' ar<-area(segdata)
#' dist <- distance(segdata)
#' diag(dist)<-sqrt(ar) * 0.6
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'segdata'
#' 
#' SP(x, spatobj = segdata)
#' 
#' SP(x, folder = foldername, shape = shapename, fdist = 'l', itype = 'between') 
#' 
#' SP(x, spatobj = segdata, diagval ='a', itype = 'one')
#' 
#' SP(x, d = dist, fdist = 'e')
#'
#' @seealso Proximity measures: \code{\link{Pxx}}, 
#' \code{\link{Pxy}},  \code{\link{Poo}}
#' @seealso Clustering Indices: 
#' \code{\link{ACL}}, \code{\link{RCL}}
#' @export



SP <- function(x, d = NULL, fdist = "e", distin = "m", distout = "m", diagval = "0", itype = "multi", spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    varTotal <- colSums(x)
    tx <- rowSums(x)
    if (is.null(d)) 
        d <- distance(spatobj = spatobj, folder = folder, shape = shape, distin = distin, distout = distout, diagval = diagval)
    if (fdist == "e") 
        d <- exp(-d)
    if (itype == "between") {
        Poo1 <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
        for (k1 in 1:ncol(x)) for (k2 in 1:ncol(x)) Poo1[k1, k2] <- sum((d %*% (x[, k1] + x[, k2])) * (x[, k1] + x[, k2])/(varTotal[k1] + varTotal[k2])^2)
    }
    if (itype == "multi" || itype == "one") Poo1 <- sum((d %*% tx) * tx/(sum(tx))^2)
    Pxx1 <- vector(length = ncol(x))
    for (k in 1:ncol(x)) Pxx1[k] <- sum((d %*% x[, k]) * x[, k]/(varTotal[k])^2)
    if (itype == "multi") 
        result <- sum(Pxx1 * varTotal)/(sum(x) * Poo1)
    if (itype == "one") 
        result <- Pxx1/Poo1
    if (itype == "between") {
        result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
        for (k1 in 1:ncol(x)) for (k2 in 1:ncol(x)) result[k1, k2] <- (varTotal[k1] * Pxx1[k1] + varTotal[k2] * Pxx1[k2])/((varTotal[k1] + varTotal[k2]) * Poo1[k1, k2])
    }
    return(round(result, 4))
}




#' A function to compute the relative clustering index (RCL)
#'
#' @usage RCL(x, d = NULL, fdist = 'e', distin = 'm',  distout = 'm', diagval = '0', 
#' spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param d - a matrix of the distances between spatial unit centroids
#' @param fdist - the method used for distance interaction matrix: 
#' e' for inverse exponential function (by default) and 'l' for linear.
#' @param distin - input metric conversion, based on  \pkg{bink} package and 
#' includes conversions from 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.
#' @param distout - output metric conversion, based on  \pkg{bink} package and 
#' includes conversions to 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.
#' @param diagval - when providing a spatial object or a shape file, 
#' the user has the choice of the spatial matrix diagonal definition: 
#' diagval = '0' (by default) for an null diagonal and diagval = 'a' 
#' to compute the diagonal as 0.6 * square root (spatial units area) (White, 1983) 
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension) which contains the geographic information
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return a matrix containing relative clustering index values for each pair of groups
#' @references Massey D. S. and Denton N. A. (1988) \emph{The dimensions 
#' of residential segregation}. Social Forces 67(2),  pp. 281-315.
#' @description The relative clustering index, RCL, compares the mean 
#' proximity of a group to the mean proximity of another group. 
#' The function can be used in two ways: to provide a distance matrix 
#' or a external geographic information source (spatial object 
#' or shape file).
#' @examples x <- segdata@data[ ,1:2]
#' ar<-area(segdata)
#' dist <- distance(segdata)
#' diag(dist)<-sqrt(ar) * 0.6
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'segdata'
#' 
#' RCL(x, spatobj = segdata)
#' 
#' RCL(x, folder = foldername, shape = shapename, fdist = 'l') 
#' 
#' RCL(x, spatobj = segdata, diagval ='a')
#' 
#' RCL(x, d = dist, fdist = 'e')
#' 
#' @seealso Proximity measures: \code{\link{Pxx}}, 
#' \code{\link{Pxy}},  \code{\link{Poo}},  \code{\link{SP}}
#' @seealso Clustering Indices: \code{\link{ACL}}
#' @export


RCL <- function(x, d = NULL, fdist = "e", distin = "m", distout = "m", diagval = "0", spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(d)) 
        d <- distance(spatobj = spatobj, folder = folder, shape = shape, distin = distin, distout = distout, diagval = diagval)
    if (fdist == "e") 
        d <- exp(-d)
    varTotal <- colSums(x)
    Pxx1 <- vector(length = ncol(x))
    for (k in 1:ncol(x)) Pxx1[k] <- sum((d %*% x[, k]) * x[, k]/(varTotal[k])^2)
    result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    for (k1 in 1:ncol(x)) for (k2 in 1:ncol(x)) {
        if (fdist == "e") 
            result[k1, k2] <- Pxx1[k1]/Pxx1[k2] - 1
        if (fdist == "l") 
            result[k1, k2] <- 1 - Pxx1[k1]/Pxx1[k2]
    }
    return(round(result, 4))
}

######################## 

# CENTRALISATION INDEXES

######################## 

#' A function to compute Duncan's Absolute Centralisation Index (ACEDuncan)
#'
#' @usage ACEDuncan(x, dc = NULL, center = 1, 
#' spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param dc - a numeric vector containing the distances between spatial units 
#' centroids and the central spatial unit
#' @param center - a numeric value giving the number of the spatial unit that 
#' represents the centre in the table
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension) which contains the geographic information
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return a numeric vector containing Duncan's asolute centralisation index 
#' value for each population group
#' @references Duncan O. D. and Duncan B. (1955) \emph{A 
#' Methodological Analysis of Segregation Indexes}. 
#' American Sociological Review 41, pp. 210-217
#' @description Duncan's absolute centralization index measures the 
#' proportion of a group that should change its localization to 
#' achieve the same level of centralization as the rest of the population.
#' The function can be used in two ways: to provide a vector containing 
#' the distances between spatial units centroids or a external geographic 
#' information source (spatial object or shape file).
#' @examples x <- segdata@data[ ,1:2]
#' ar<-area(segdata)
#' distc<- distcenter(segdata, center = 45)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'segdata'
#' 
#' ACEDuncan(x, dc=distc) 
#' 
#' ACEDuncan(x, spatobj = segdata, center = 45) 
#' 
#' ACEDuncan(x, folder = foldername, shape = shapename, center = 45) 
#'
#' @seealso Relative Centralisation Index: \code{\link{RCE}}
#' @seealso Absolute Centralisation Index: \code{\link{ACE}}
#' @export


ACEDuncan <- function(x, dc = NULL, center = 1, spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(dc)) 
        dc <- distcenter(spatobj = spatobj, folder = folder, shape = shape, center)
    result <- vector(length = ncol(x))
    varTotal <- colSums(x)
    t <- rowSums(x)
    prop <- varTotal/sum(varTotal)
    xprovi <- cbind(x, t, dc)
    xprovi <- xprovi[order(xprovi[, ncol(xprovi)]), ]
    xprovi <- as.data.frame(xprovi)
    xxprovi <- xprovi[, c(1, ncol(xprovi) - 1, ncol(xprovi))]
    for (k in 1:ncol(x)) {
        xxprovi[, 1] <- xprovi[k]
        varTotalb <- colSums(xxprovi)
        XI1 <- cumsum(xxprovi[, 1])[1:(nrow(xxprovi) - 1)]/varTotalb[1]
        XI <- cumsum(xxprovi[, 1])[2:nrow(xxprovi)]/varTotalb[1]
        YI1 <- cumsum(xxprovi[, 2])[1:(nrow(xxprovi) - 1)]/varTotalb[2]
        YI <- cumsum(xxprovi[, 2])[2:nrow(xxprovi)]/varTotalb[2]
        result[k] <- (XI1 %*% YI - XI %*% YI1)/(1 - prop[k])
    }
    return(round(result, 4))
}


#' A function to compute Relatice Centralisation Index (RCE)
#'
#' @usage RCE(x, dc = NULL, center = 1, spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param dc - a numeric vector containing the distances between spatial units 
#' centroids and the central spatial unit
#' @param center - a numeric value giving the number of the spatial unit that 
#' represents the centre in the table
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension) which contains the geographic information
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return a matrix containing relative centralisation index values 
#' @references Duncan O. D. and Duncan B. (1955) \emph{A 
#' Methodological Analysis of Segregation Indexes}. 
#' American Sociological Review 41, pp. 210-217
#' @description The relative centralisation index measures the 
#' proportion of a group that should change its localization to 
#' achieve the same level of centralization as another group.
#' The function can be used in two ways: to provide a vector containing 
#' the distances between spatial units centroids or a external geographic 
#' information source (spatial object or shape file).
#' @examples x <- segdata@data[ ,1:2]
#' distc<- distcenter(segdata, center = 45)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'segdata'
#' 
#' RCE(x, dc=distc) 
#' 
#' RCE(x, spatobj = segdata, center = 45) 
#' 
#' RCE(x, folder = foldername, shape = shapename, center = 45) 
#'
#' @seealso Duncan's Absolute Centralisation Index: \code{\link{ACEDuncan}}
#' @seealso Absolute Centralisation Index: \code{\link{ACE}}
#' @export


RCE <- function(x, dc = NULL, center = 1, spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(dc)) 
        dc <- distcenter(spatobj = spatobj, folder = folder, shape = shape, center)
    result <- matrix(data = 0, nrow = ncol(x), ncol = ncol(x))
    varTotal <- colSums(x)
    xprovi <- cbind(x, dc)
    xprovi <- xprovi[order(xprovi[, ncol(xprovi)]), ]
    xprovi <- as.data.frame(xprovi)
    for (k1 in 1:ncol(x)) for (k2 in 1:ncol(x)) {
        XI1 <- cumsum(xprovi[, k1])[1:(nrow(xprovi) - 1)]/varTotal[k1]
        XI <- cumsum(xprovi[, k1])[2:nrow(xprovi)]/varTotal[k1]
        YI1 <- cumsum(xprovi[, k2])[1:(nrow(xprovi) - 1)]/varTotal[k2]
        YI <- cumsum(xprovi[, k2])[2:nrow(xprovi)]/varTotal[k2]
        result[k1, k2] <- XI1 %*% YI - XI %*% YI1
    }
    return(round(result, 4))
}



#' A function to compute Massey Absolute Centralisation Index (ACE)
#'
#' @usage ACE(x, a = NULL, dc = NULL, center = 1, 
#' spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param a - a numeric vector containing the areas of spatial units
#' @param dc - a numeric vector containing the distances between spatial units 
#' centroids and the central spatial unit
#' @param center - a numeric value giving the number of the spatial unit that 
#' represents the centre in the table
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the name of the shapefile 
#' (without the .shp extension) which contains the geographic information
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return a numeric vector containing asolute centralisation index value for 
#' each population group
#' @references Massey D. S. and Denton N. A. (1988) \emph{
#' The dimensions of residential segregation}. 
#' Social Forces 67(2),  pp. 281-315.
#' @description The absolute centralization index measures a group
#' spatial distribution compared to the distribution of land area 
#' around the city centre. The function can be used in two ways: to provide 
#' an area vector and a vector containing the distances between spatial units 
#' centroids and  the central spatial unit or a external geographic information 
#' source (spatial object or shape file).
#' 
#' @examples x <- segdata@data[ ,1:2]
#' ar<-area(segdata)
#' distc<- distcenter(segdata, center = 45)
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'segdata'
#' 
#' ACE(x, a = ar, dc=distc) 
#' 
#' ACE(x, spatobj = segdata, center = 45) 
#' 
#' ACE(x, folder = foldername, shape = shapename, center = 45) 
#' 
#' @seealso Duncan's Absolute Centralisation Index: \code{\link{ACEDuncan}}
#' @seealso Relative Centralisation Index: \code{\link{RCE}}
#' @export


ACE <- function(x, a = NULL, dc = NULL, center = 1, spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    if (is.null(a)) 
        a <- area(spatobj = spatobj, folder = folder, shape = shape)
    if (is.null(dc)) 
        dc <- distcenter(spatobj = spatobj, folder = folder, shape = shape, center)
    result <- vector(length = ncol(x))
    varTotal <- colSums(x)
    t <- rowSums(x)
    prop <- varTotal/sum(varTotal)
    xprovi <- cbind(x, a, dc)
    xprovi <- xprovi[order(xprovi[, ncol(xprovi)]), ]
    xprovi <- as.data.frame(xprovi)
    for (k in 1:ncol(x)) {
        XI1 <- cumsum(xprovi[, k])[1:(nrow(xprovi) - 1)]/varTotal[k]
        AI <- cumsum(xprovi$a)[2:nrow(xprovi)]/sum(xprovi$a)
        XI <- cumsum(xprovi[, k])[2:nrow(xprovi)]/varTotal[k]
        AI1 <- cumsum(xprovi$a)[1:(nrow(xprovi) - 1)]/sum(xprovi$a)
        result[k] <- XI1 %*% AI - XI %*% AI1
    }
    return(round(result, 4))
}




######################## 

# MULTIGROUP INDEXES

######################## 


#' A function to compute Shannon-Wiener diversity (entropy) index
#'
#' @usage HShannon(x) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @return Shannon-Wiener diversity index 
#' @references Shannon C. E. (1948) \emph{A mathematical theory 
#' of communication}. Bell System Technical Journal (27) 
#' @description The Shannon-Wiener diversity index is based on 
#' the notion of entropy and measures population heterogeneity.
#' @examples x <- segdata@data[ ,1:2]
#' HShannon(x) 
#' @seealso  Social diversity indices: 
#' \code{\link{NShannon}}, \code{\link{ISimpson}}, 
#' @seealso Multigroup indices: 
#' \code{\link{PMulti}}, \code{\link{GiniMulti}}, \code{\link{DMulti}},  
#' \code{\link{HMulti}}, \code{\link{CMulti}}, \code{\link{RelDivers}}
#' @export


HShannon <- function(x) {
    x <- as.matrix(x)
    pTotal <- colSums(x)/sum(x)
    result <- -sum(pTotal * log(pTotal))
    return(round(result, 4))
}


#' A function to compute Shannon-Wiener diversity normalized index
#'
#' @usage NShannon(x) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @return Shannon-Wiener normalized diversity index 
#' @references Shannon C. E. (1948) \emph{A mathematical theory 
#' of communication}. Bell System Technical Journal (27) 
#' @description The Shannon-Wiener diversity index is based on 
#' the notion of entropy and measures population heterogeneity.
#' @examples x <- segdata@data[ ,1:2]
#' NShannon(x) 
#' @seealso  Other multigroup eveness indices: 
#' \code{\link{HShannon}}, \code{\link{ISimpson}}, 
#' \code{\link{GiniMulti}}, \code{\link{DMulti}}, \code{\link{HMulti}}, 
#' \code{\link{CMulti}}
#' @seealso Other multigroup indices: \code{\link{PMulti}}, 
#' \code{\link{RelDivers}}
#' @export


NShannon <- function(x) {
    x <- as.matrix(x)
    pTotal <- colSums(x)/sum(x)
    result <- -sum(pTotal * log(pTotal))
    result <- result/log(ncol(x))
    return(round(result, 4))
}

#' A function to compute Simpson's interaction index
#'
#' @usage ISimpson(x) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @return Simpson's interaction index 
#' @references Simpson E. H. (1949) \emph{Measurement of diversity}. 
#' Nature 163:688 
#' @description Simpson's interaction index measures the probability 
#' that randomly selected individuals are not in the same group. 
#' @examples x <- segdata@data[ ,1:2]
#' ISimpson(x) 
#' @seealso  Social diversity indices: 
#' \code{\link{HShannon}}, \code{\link{NShannon}}, 
#' @seealso Multigroup indices: 
#' \code{\link{PMulti}}, \code{\link{GiniMulti}}, \code{\link{DMulti}},  
#' \code{\link{HMulti}}, \code{\link{CMulti}}, \code{\link{RelDivers}}
#' @export



ISimpson <- function(x) {
    x <- as.matrix(x)
    pTotal <- colSums(x)/sum(x)
    result <- sum(pTotal * (1 - pTotal))
    return(round(result, 4))
}



#' A function to compute multigroup Gini index
#'
#' @usage GiniMulti(x) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @references Reardon S. F. (1998) \emph{Measures of racial 
#' diversity and segregation in multigroup and hierarchical 
#' structured Populations}. Annual meeting of the Eastern 
#' Sociological Society, Philadelphia 
#' @description GiniMulti is a multigroup version of 
#' the \code{\link{Gini}} index 
#' @examples x <- segdata@data[ ,1:2]
#' GiniMulti(x) 
#' @seealso Multigroup indices: 
#' \code{\link{PMulti}}, \code{\link{GiniMulti}},   
#' \code{\link{HMulti}}, \code{\link{CMulti}}, \code{\link{RelDivers}}
#' @seealso  Social diversity indices: 
#' \code{\link{HShannon}}, \code{\link{NShannon}}, 
#' \code{\link{ISimpson}}, 
#' @export



GiniMulti <- function(x) {
    x <- as.matrix(x)
    pTotal <- colSums(x)/sum(x)
    II <- sum(pTotal * (1 - pTotal))
    vartemp <- vector(length = ncol(x))
    t <- rowSums(x)
    Total <- sum(x)
    p <- x/t
    pij <- array(0, dim = c(ncol(x), nrow(x), nrow(x)))
    for (k in 1:ncol(x)) {
        for (i in 1:nrow(x)) pij[k, , i] <- p[i, k]
        for (i in 1:nrow(x)) pij[k, i, ] <- abs(pij[k, i, ] - p[i, k])
        vartemp[k] <- sum(t * (t %*% pij[k, , ]))
    }
    result <- sum(vartemp)/(2 * Total * Total * II)
    return(round(result, 4))
}





#' A function to compute multigroup dissimilarity index
#'
#' @usage DMulti(x) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @return Multigroup dissimilarity index 
#' @references Sakoda J. N. (1981) \emph{A generalized Index of 
#' dissimilarity}. Demography,18, 245-250 
#' @description Multigroup dissimilarity, DMulti, is a multigroup 
#' version of  Duncan's dissimilarity index (\code{\link{DIDuncan}})
#' @examples x <- segdata@data[ ,1:2]
#' DMulti(x) 
#' @seealso Multigroup indices: 
#' \code{\link{PMulti}}, \code{\link{GiniMulti}},   
#' \code{\link{HMulti}}, \code{\link{CMulti}}, \code{\link{RelDivers}}
#' @seealso  Social diversity indices: 
#' \code{\link{HShannon}}, \code{\link{NShannon}}, 
#' \code{\link{ISimpson}}, 
#' @export 



DMulti <- function(x) {
    x <- as.matrix(x)
    result <- 0
    pTotal <- colSums(x)/sum(x)
    II <- sum(pTotal * (1 - pTotal))
    t <- rowSums(x)
    Total <- sum(x)
    p <- x/t
    for (k in 1:ncol(x)) result <- result + t %*% abs(p[, k] - pTotal[k])
    result <- result/(2 * Total * II)
    return(round(result, 4))
}




#' A function to compute multigroup normalised exposure (PMulti)
#'
#' @usage PMulti(x) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @return Multigroup normalised isolation index 
#' @references  James, F. J. (1986) \emph{A New Generalized 'Exposure-Based' 
#' Segregation Index}. Sociological Methods and Research, 14, pp. 301-316
#' @references  Reardon S. F. and G. Firebaugh (2002) \emph{Measures of 
#' Multigroup Segregation}. Sociological Methodology, 32(1), pp 33-67
#' @description The multigroup normalised isolation index, PMulti, 
#' is a multigroup version of the isolation index (\code{\link{xPx}})
#' @examples x <- segdata@data[ ,1:2]
#' PMulti(x) 
#' @seealso Multigroup indices: 
#' \code{\link{GiniMulti}}, \code{\link{DMulti}},  
#' \code{\link{HMulti}}, \code{\link{CMulti}}, \code{\link{RelDivers}}
#' @seealso  Social diversity indices: 
#' \code{\link{HShannon}}, \code{\link{NShannon}}, 
#' \code{\link{ISimpson}}, 
#' @export 



PMulti <- function(x) {
    x <- as.matrix(x)
    result <- 0
    pTotal <- colSums(x)/sum(x)
    t <- rowSums(x)
    Total <- sum(x)
    p <- x/t
    for (k in 1:ncol(x)) result <- result + t %*% ((p[, k] - pTotal[k]) * (p[, k] - pTotal[k])/(1 - pTotal[k]))
    result <- result/Total
    return(round(result, 4))
}


#' A function to compute multigroup relative diversity index
#'
#' @usage RelDivers(x) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @return Multigroup relative diversity index 
#' @references Carlson S. M. (1992) \emph{Trends in race/sex 
#' occupational inequality:  conceptual and measurement issues}. 
#' Social Problems, 39, p. 269-290
#' @description The relative diversity index is a multigroup 
#' index based on Simpson's interaction index \code{\link{ISimpson}} 
#' @examples x <- segdata@data[ ,1:2]
#' RelDivers(x) 
#' @seealso Multigroup indices: 
#' \code{\link{PMulti}}, \code{\link{GiniMulti}}, \code{\link{DMulti}},  
#' \code{\link{HMulti}}, \code{\link{CMulti}}
#' @seealso  Social diversity indices: 
#' \code{\link{HShannon}}, \code{\link{NShannon}}, 
#' \code{\link{ISimpson}}, 
#' @export 



RelDivers <- function(x) {
    x <- as.matrix(x)
    result <- 0
    pTotal <- colSums(x)/sum(x)
    II <- sum(pTotal * (1 - pTotal))
    t <- rowSums(x)
    Total <- sum(x)
    p <- x/t
    for (k in 1:ncol(x)) result <- result + t %*% ((p[, k] - pTotal[k])^2)
    result <- result/(Total * II)
    return(round(result, 4))
}

#' A function to compute multigroup entropy segregation index
#'
#' @usage HMulti(x) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @return Multigroup entropy segregation index
#' @references Theil H. (1972)  \emph{Statistical decomposition analysis: with 
#' applications in the social and administrative.} Amsterdam, North-Holland, 337 p.
#' @description The multigroup version of Theil's entropy index \code{\link{HTheil}} 
#' @examples x <- segdata@data[ ,1:2]
#' HMulti(x) 
#' @seealso Multigroup indices: 
#' \code{\link{PMulti}}, \code{\link{GiniMulti}}, \code{\link{DMulti}},  
#' \code{\link{CMulti}}, \code{\link{RelDivers}}
#' @seealso  Social diversity indices: 
#' \code{\link{HShannon}}, \code{\link{NShannon}}, 
#' \code{\link{ISimpson}}, 
#' @export 



HMulti <- function(x) {
    x <- as.matrix(x)
    E <- vector("list", ncol(x))
    Total <- sum(x)
    pTotal <- colSums(x)/Total
    tx <- rowSums(x)
    px <- x/tx
    Etot <- sum(pTotal * log(1/pTotal))
    result <- 0
    for (k in 1:ncol(x)) {
        E[[k]] <- tx * px[, k] * log(px[, k]/pTotal[k])
        E[[k]][is.na(E[[k]])] <- 0
        result <- result + sum(E[[k]])
    }
    result <- result/(Etot * Total)
    return(round(result, 4))
}


#' A function to compute multigroup squared coefficient of variation index
#'
#' @usage CMulti(x) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @return Multigroup entropy segregation index
#' @references Reardon S. F. and Firebaugh G. (2002) \emph{Measures of multigroup 
#' segregation}. Sociological Methodology, 32, pp. 33-67.
#' @description The index can be interpreted as a measure of the variance of the 
#' spatial representation of the groups accros spatial unite, or as a normalized 
#' chi-squared measure of association between groups and units. 
#' @examples x <- segdata@data[ ,1:2]
#' CMulti(x) 
#' @seealso Multigroup indices: 
#' \code{\link{PMulti}}, \code{\link{GiniMulti}}, \code{\link{DMulti}},  
#' \code{\link{HMulti}}, \code{\link{RelDivers}}
#' @seealso  Social diversity indices: 
#' \code{\link{HShannon}}, \code{\link{NShannon}}, 
#' \code{\link{ISimpson}}, 
#' @export 



CMulti <- function(x) {
    x <- as.matrix(x)
    Total <- sum(x)
    pTotal <- colSums(x)/Total
    tx <- rowSums(x)
    px <- x/tx
    pTotal2 <- matrix(rep(pTotal, nrow(x)), nrow = nrow(x), ncol = ncol(x), byrow = T)
    result <- sum((tx/Total) * ((px - pTotal2)^2)/((ncol(x) - 1) * pTotal2))
    return(round(result, 4))
}



######################## 

# LOCAL INDEXES

######################## 

#' A function to compute location quotients (LQs)
#'
#' @usage LQ(x) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @return a matrix of location quotiens 
#' @references Isard W. (1960) \emph{Methods of regional analysis: 
#' an introduction to regional science}. The MIT Press, Cambridge
#' @description Location quotients compare the relative part of a 
#' group in a particular spatial unit, to the relative part of that 
#' same group in the area.
#' @examples x <- segdata@data[ ,1:2]
#' LQ(x) 
#' @seealso Other local indices \code{\link{LShannon}}
#' \code{\link{HLoc}}, \code{\link{LSimpson}}   
#' @export

LQ <- function(x) {
    x <- as.matrix(x)
    result <- x
    Total <- sum(x)
    t <- rowSums(x)
    for (i in 1:ncol(x)) result[, i] <- (x[, i]/t)/(sum(x[, i])/sum(t))
    return(round(result, 4))
}




#' A function to compute local diversity index
#'
#' @usage HLoc(x) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @return a numeric vector containing  diversity index value for 
#' each population group
#' @references Theil H. (1972) \emph{Statistical Decomposition Analysis}. 
#' North-Holland, Amsterdam
#' @description Local diversity index, HLoc, is a local 
#' adaptation of Pielou's normalized diversity index \code{\link{NShannon}}.
#' @examples x <- segdata@data[ ,1:2]
#' HLoc(x) 
#' @seealso Other local indices \code{\link{LQ}}
#' \code{\link{LShannon}}, \code{\link{LSimpson}}   
#' @export


HLoc <- function(x) {
    x <- as.matrix(x)
    Total <- sum(x)
    t <- rowSums(x)
    result <- x/t
    result <- cbind(result, 0)
    for (i in 1:nrow(result)) {
        n <- 0
        for (j in 1:(ncol(result) - 1)) if (result[i, j] > 0) {
            result[i, ncol(result)] <- result[i, ncol(result)] + result[i, j] * log(result[i, j])
            n <- n + 1
        }
        result[i, ncol(result)] <- -result[i, ncol(result)]/log(n)
    }
    result <- result[, ncol(result)]
    result[is.na(result)] <- 0
    return(round(result, 4))
}


#' A function to compute Shannon-Wiener local diversity (entropy) index
#'
#' @usage LShannon(x) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @return Local Shannon-Wiener diversity index 
#' @references Shannon C. E. (1948) \emph{A mathematical theory 
#' of communication}. Bell System Technical Journal (27) 
#' @description The Shannon-Wiener diversity index is based on 
#' the notion of entropy and measures population heterogeneity.
#' @examples x <- segdata@data[ ,1:2]
#' LShannon(x) 
#' @seealso Other local indices: \code{\link{LQ}}, 
#' \code{\link{HLoc}}, \code{\link{LSimpson}}   
#' @export


LShannon <- function(x) {
  x <- as.matrix(x)
  Total <- sum(x)
  t <- rowSums(x)
  result <- x/t
  result <- cbind(result, 0)
  for (i in 1:nrow(result)) {
    n <- 0
    for (j in 1:(ncol(result) - 1)) if (result[i, j] > 0) {
      result[i, ncol(result)] <- result[i, ncol(result)] + result[i, j] * log(result[i, j])
      n <- n + 1
    }
    result[i, ncol(result)] <- -result[i, ncol(result)]
  }
  result <- result[, ncol(result)]
  result[is.na(result)] <- 0
  return(round(result, 4))
}


#' A function to compute local Simpson's index
#'
#' @usage LSimpson (x) 
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @return Local Simpson's interaction index 
#' @references Simpson E. H. (1949) \emph{Measurement of diversity}. 
#' Nature 163:688 
#' @description Local Simpson's interaction index measures the probability 
#' that randomly selected individuals are not in the same group in 
#' each spatial unit. 
#' @examples x <- segdata@data[ ,1:2]
#' LSimpson (x) 
#' @seealso Other local indices: \code{\link{LQ}}, 
#' \code{\link{HLoc}}, \code{\link{LShannon}}   
#' @export


LSimpson <- function(x) {
  x <- as.matrix(x)
  p <- x/rowSums(x)
  result <- rowSums(p*(1-p))
  result[is.na(result)] <- 0
  return(round(result, 4))
}




################################################## 

# MONTE CARLO SIMULATIONS

################################################## 



#' A function to test segregation indices using Monte Carlo simulations
#'
#' @usage MCTest(x, fun, simtype = "rand", proba = NULL, delta = 0.5, ptype = "int", 
#' variant = "s", distin = "m", distout = "m", itype = "multi", exact = F, 
#' diagval = "0", spatmat = "c", c = NULL, queen = TRUE, b = NULL, p = NULL, 
#' a = NULL, d = NULL, dc = NULL, fdist = "e", center = 1, nsim = 99, 
#' spatobj = NULL, folder = NULL, shape = NULL)
#' @param x - an object of class matrix (or which can be coerced to that class), 
#' where each column represents the distribution of a population group, within 
#' spatial units. The number of columns should be greater than 1 (at least 2 
#' population groups are required). You should not include a column with total 
#' population in each unit, because this will be interpreted as a group.
#' @param fun - a character vector with the segregation function 
#' to be tested (only one-group or multigroup indexes)
#' @param simtype - a character vector with type of simulation. 
#' If simtype='perm', the function will use the permutation test. 
#' If simtype='rand', the population is realocated in the spatial 
#' units randomly. If simtype='area', the location probabilities 
#' are proportional to the spatial units' area. For the default 
#' value type='user', the function expects you to introduce a 
#' probability location vector.
#' @param proba - a vector with location probabilities. By 
#' default proba = NULL being calculated depending on the 
#' simulation type. The parameter is necessary when a user method 
#' is specified.
#' @param nsim - the number of simulations
#' @param delta - a specific Atikinson inequality aversion parameter,  
#' by default equal to 0.5,varying from 0 to 1.
#' @param variant - a character variable that allows to choose the Wong's index 
#' version: variant = 's' for the index adjusted for contiguous spatial units 
#' boundary lengths and perimeter/area ratio (by default) and variant = 'w' 
#' for the version based only on shared boundaries length
#' @param itype - a character string defining the type of some proximity indices:
#' itype = 'multi' (by default) for the multigroup index (White, 1986) or
#' itype = 'one' for the one-group version (Apparicio et al, 2008)
#' @param exact - a logical variable to specifiy the version of interaction indices 
#' exact = FALSE (by default) for the approximate version of the index, 
#' and exact = TRUE for the exact version
#' @param spatmat - the method used for some spatial calculations: 'c' for the 
#' contiguity matrix (by default) or any other user spatial interaction matrix 
#' and 'd' for the inverse exponential function of the distance. 
#' @param c - a standard binary contiguity (adjacency) symmetric matrix where 
#' each element \emph{Cij} equals 1 if \emph{i}-th and \emph{j}-th spatial 
#' units are adjacent, and 0 otherwise.
#' @param queen - logical parameter difining criteria used for contiguity 
#' matrix computation, TRUE (by default) for queen , FALSE for rook 
#' @param d - a matrix of the distances between spatial unit centroids
#' @param dc - a numeric vector containing the distances between spatial units 
#' centroids and the central spatial unit
#' @param distin - input metric conversion, based on  \pkg{bink} package and 
#' includes conversions from 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.
#' @param distout - output metric conversion, based on  \pkg{bink} package and 
#' includes conversions to 'm', 'km', 'inch', 'ft', 'yd', 'mi', 'naut_mi', etc.
#' @param fdist - the method used for distance interaction matrix: 
#' e' for inverse exponential function (by default) and 'l' for linear.
#' @param diagval - when providing a spatial object or a shape file, 
#' the user has the choice of the spatial matrix diagonal definition: 
#' diagval = '0' (by default) for an null diagonal and diagval = 'a' 
#' to compute the diagonal as 0.6 * square root (spatial units area) (White, 1983) 
#' @param b - a common boundaries matrix where each element \emph{Bij} 
#' equals the shared boundary of \emph{i}-th and \emph{j}-th spatial units.
#' @param p - a numeric vector containing the perimeters of spatial units
#' @param ptype - a string variable giving two options for perimeter 
#' calculation, when needed: 'int' to use only interior
#' boundaries of spatial units, and 'all' to use entire boundaries, 
#' including the boundaries to the exterior of the area
#' @param a - a numeric vector containing the areas of spatial units
#' @param center - a numeric value giving the number of the spatial unit that 
#' represents the centre in the table
#' @param folder - a character vector with the folder (directory) 
#' name indicating where the shapefile with the geographic information 
#' is located.
#' @param shape - a character vector with the shape file name
#' @param spatobj - a spatial object (SpatialPolygonsDataFrame) containing 
#' geographic information
#' @return A list with thre elements: index's name, a Summary of the simulation
#' and index simulated distribution
#' @references Tivadar M., Schaeffer Y, Torre A. and Bray F. (2014) 
#' \emph{OASIS - un Outil d'Analyse de la Segregation et des Inegalites 
#' Spatiales}.  Cybergeo : European Journal of Geography, GeOpenMod, 
#' document 699
#' @description The Monte Carlo tests for one-group or multigroup 
#' segregation indexes.
#' @examples x <- segdata@data[ ,1:2]
#' foldername <- system.file('extdata', package = 'OasisR')
#' shapename <- 'segdata'
#' ar <- area(segdata)
#' probavector<-ar/sum(ar)
#' 
#' MCTest (x, fun='ISMorrill', simtype = 'perm', spatobj = segdata)
#' 
#' MCTest (x, fun='ISMorrill', simtype = 'rand', spatobj = segdata)
#' 
#' MCTest (x, fun='ISMorrill', simtype = 'area', spatobj = segdata)
#' 
#' MCTest (x, fun='ISMorrill', simtype='user', proba = probavector, 
#' spatobj = segdata)
#' 
#' @seealso \code{\link{MCPlot}}
#' @importFrom stats density
#' @importFrom graphics plot segments mtext
#' @export


MCTest <- function(x, fun, simtype = "rand", proba = NULL, delta = 0.5, ptype = "int", variant = "s", distin = "m", distout = "m", itype = "multi", exact = F, diagval = "0", 
    spatmat = "c", c = NULL, queen = TRUE, b = NULL, p = NULL, a = NULL, d = NULL, dc = NULL, fdist = "e", center = 1, nsim = 99, spatobj = NULL, folder = NULL, shape = NULL) {
    x <- as.matrix(x)
    assign("func", eval(parse(text = fun)))
    xdistr <- vector("list", nsim)
    if (simtype == "perm") 
        for (k in 1:nsim) {
            xdistr[[k]] <- matrix(nrow = nrow(x), ncol = ncol(x))
            neworder <- sample(c(1:nrow(x)), size = nrow(x))
            for (i in 1:nrow(x)) for (j in 1:ncol(x)) xdistr[[k]][i, j] <- x[neworder[i], j]
        }
    if (simtype == "rand") 
        proba <- rep(1/nrow(x), nrow(x))
    if (simtype == "area") {
        if (is.null(a)) 
            a <- area(spatobj = spatobj, folder = folder, shape = shape)
        proba <- a/sum(a)
    }
    if (simtype != "perm") 
        for (i in 1:nsim) {
            xdistr[[i]] <- matrix(0, nrow = nrow(x), ncol = ncol(x))
            for (k in 1:ncol(x)) {
                xprovi <- table(sample(nrow(x), size = sum(x[, k]), replace = T, prob = proba))
                dimv <- as.numeric(dimnames(xprovi)[[1]])
                xdistr[[i]][dimv, k] <- xprovi
            }
            xdistr[[i]] <- xdistr[[i]][rowSums(xdistr[[i]]) > 0, ]
        }
    nvar <- ncol(x)
    if (fun == "HMulti" || fun == "PMulti" || fun == "GiniMulti" || fun == "DMulti" || fun == "RelDivers" || fun == "CMulti") 
        nvar <- 1
    if (fun == "SP" && itype == "multi") 
        nvar <- 1
    if (fun == "Poo" && itype == "multi") 
        nvar <- 1
    IndTest <- matrix(nrow = nvar, ncol = 5)
    IndTest <- as.data.frame(IndTest)
    names(IndTest) <- c("Var", fun, "Simulated", "Rank", "P.Value")
    IndTest$Var <- 1:nvar
    resim <- matrix(nrow = nvar, ncol = nsim)
    if (fun == "ISDuncan" || fun == "Gorard" || fun == "Gini" || fun == "HTheil" || fun == "Eta2" || fun == "GiniMulti" || fun == "DMulti" || fun == "PMulti" || fun == "RelDivers" || 
        fun == "HMulti" || fun == "CMulti") {
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]])
        IndTest[, 2] <- func(x)
    }
    if (fun == "Atkinson") {
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]], delta)
        IndTest[, 2] <- func(x, delta)
    }
    if (fun == "xPx") {
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]], exact = exact)
        IndTest[, 2] <- func(x, exact = exact)
    }
    if (fun == "ISMorrill") {
        if (is.null(c)) 
            c <- contig(spatobj = spatobj, folder = folder, shape = shape, queen = queen)
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]], c, queen = queen)
        IndTest[, 2] <- func(x, c)
    }
    if (fun == "ACL") {
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]], c = c, spatmat = spatmat, distin = distin, distout = distout, diagval = diagval, spatobj = spatobj, folder = folder, 
            shape = shape)
        IndTest[, 2] <- func(x, c = c, spatmat = spatmat, distin = distin, distout = distout, diagval = diagval, spatobj = spatobj, folder = folder, shape = shape)
    }
    if (fun == "ISWong") {
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]], b = b, a = a, p = p, ptype = ptype, variant = variant, spatobj = spatobj, folder = folder, shape = shape)
        IndTest[, 2] <- func(x, b = b, a = a, p = p, ptype = ptype, variant = variant, spatobj = spatobj, folder = folder, shape = shape)
    }
    if (fun == "Delta" || fun == "ACO") {
        if (is.null(a)) 
            a <- area(spatobj = spatobj, folder = folder, shape = shape)
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]], a)
        IndTest[, 2] <- func(x, a)
    }
    
    if (fun == "DPxx") {
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]], d = d, distin = distin, distout = distout, diagval = diagval, spatobj = spatobj, folder = folder, shape = shape)
        IndTest[, 2] <- func(x, d = d, distin = distin, distout = distout, diagval = diagval, spatobj = spatobj, folder = folder, shape = shape)
    }
    
    if (fun == "Pxx") {
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]], d = d, fdist = fdist, distin = distin, distout = distout, diagval = diagval, spatobj = spatobj, folder = folder, 
            shape = shape)
        IndTest[, 2] <- func(x, d = d, fdist = fdist, distin = distin, distout = distout, diagval = diagval, spatobj = spatobj, folder = folder, shape = shape)
    }
    if (fun == "Poo") {
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]], d = d, fdist = fdist, itype = "multi", distin = distin, distout = distout, diagval = diagval, spatobj = spatobj, 
            folder = folder, shape = shape)
        IndTest[, 2] <- func(x, d = d, fdist = fdist, itype = "multi", distin = distin, distout = distout, diagval = diagval, spatobj = spatobj, folder = folder, shape = shape)
    }
    if (fun == "SP") {
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]], d = d, fdist = fdist, itype = itype, distin = distin, distout = distout, diagval = diagval, spatobj = spatobj, folder = folder, 
            shape = shape)
        IndTest[, 2] <- func(x, d = d, fdist = fdist, itype = itype, distin = distin, distout = distout, diagval = diagval, spatobj = spatobj, folder = folder, shape = shape)
    }
    if (fun == "ACE") {
        if (is.null(a)) 
            a <- area(spatobj = spatobj, folder = folder, shape = shape)
        if (is.null(dc)) 
            dc <- distcenter(spatobj = spatobj, folder = folder, shape = shape, center)
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]], a, dc)
        IndTest[, 2] <- func(x, a, dc)
    }
    if (fun == "ACEDuncan") {
        if (is.null(dc)) 
            dc <- distcenter(spatobj = spatobj, folder = folder, shape = shape, center)
        for (k in 1:nsim) resim[, k] <- func(xdistr[[k]], dc)
        IndTest[, 2] <- func(x, dc)
    }
    for (i in 1:nvar) IndTest[i, 3] <- mean(resim[i, ])
    for (i in 1:nvar) {
        ncontor <- 0
        for (k in 1:nsim) if (IndTest[i, 2] > resim[i, k]) 
            ncontor <- ncontor + 1
        IndTest[i, 4] <- ncontor + 1
        IndTest[i, 5] <- max((nsim - ncontor)/(nsim + 1), 1/(nsim + 1))
    }
    result <- list(fun, IndTest, resim)
    names(result) <- c("Index", "Summary", "Distribution")
    return(result)
}




#' A function to plot the results of Monte Carlo simulations
#'
#' @usage MCPlot(MCTest = NULL, var = 1, dens = NULL, ind = NULL, pval = NULL, 
#' indexname = 'Index', coldist = 'red', colind = 'blue', legend = T, 
#' legendpos = 'top', cex.legend = 1, bty = 'o')
#' @param MCTest - a MCTest object prodused with \code{\link{MCTest}} function
#' to be tested (only one-group or multigroup indeces)
#' @param var - the number of the group to be plot (if several groups are simulated)
#' @param dens - a vector with the simulatated distribution of the index
#' @param ind - index value
#' @param pval - pseudo p-value
#' @param indexname - a string with the name of the index
#' @param coldist - color used to plot the simulated distribution 
#' @param colind - color used to plot the index
#' @param legendpos - a character string giving the legend's position: 
#' "bottomright", "bottom", "bottomleft", "left", "topleft", "top", 
#' "topright", "right" and "center".
#' @param cex.legend - a numerical value giving the amount by which 
#' plotting text and symbols in legend should be magnified relative to the default. 
#' @param bty - a character string which determines the type of box 
#' of the legend. If bty is one of "o" (the default), "l", "7", "c", 
#' "u", or "]" the resulting box resembles the corresponding upper 
#' case letter. A value of "n" suppresses the box.
#' @param legend - logical parameter, to control the legend's plots
#' @return A plot with Monte Carlo results
#' @references Tivadar M., Schaeffer Y, Torre A. and Bray F. (2014) 
#' \emph{OASIS - un Outil d'Analyse de la Segregation et des Inegalites 
#' Spatiales}.  Cybergeo : European Journal of Geography, GeOpenMod, 
#' document 699
#' @description Plot of Monte Carlo simulations results. The function can
#' be used in two ways: buy providing a MCTest object, using \code{\link{MCTest}} 
#' or a simulated distribution vector, a value and a name of the index
#' @examples x <- segdata@data[ ,1:2]
#' test <- MCTest (x, fun='ISMorrill', simtype = 'perm', spatobj = segdata)
#' 
#' MCPlot(test, cex.legend = 0.8)
#' 
#' MCPlot(dens = test$Distribution[1,], ind = test$Summary$ISMorrill[1], 
#' pval = test$Summary$P.Value[1], indexname = test$Index, cex.legend = 0.8)
#' 
#' @seealso \code{\link{MCTest}} 
#' @export


MCPlot <- function(MCTest  = NULL, var = 1, dens = NULL, ind = NULL, pval = NULL, indexname = "Index", 
                   coldist = "red", colind = "blue", legend = T, legendpos = "top", 
                   cex.legend = 1, bty = "o") {
    k <- var
    if (!is.null(MCTest)) {
        indexname <- MCTest$Index
        dens <- MCTest$Distribution[k, ]
        ind <- MCTest$Summary[k, 2]
        pval <- MCTest$Summary[k, 5]
    }
    nsim <- length(dens)
    if (is.null(pval)) {
        ncontor <- 0
        for (kk in 1:nsim) if (ind > dens[kk]) 
            ncontor <- ncontor + 1
        pval <- max((nsim - ncontor)/(nsim + 1), 1/(nsim + 1))
    }
    dens2 <- density(dens)
    liminf <- min(dens2$x, ind)
    limsup <- max(dens2$x, ind)
    ylimit = c(0, ((max(dens2$y)) + 2))
    esper <- mean(dens)
    plot(dens2, xlim = c(liminf, limsup), ylim = ylimit, lwd = 2, col = coldist, main = "", xlab = "Index values")
    segments(esper, 0, esper, max(dens2$y), lwd = 1, col = coldist)
    segments(ind, 0, ind, max(dens2$y), lwd = 3, col = colind)
    mtext(paste0("Monte Carlo Test: ", indexname), side = 3, font = 2, line = 2)
    mtext(paste(indexname, "=", ind, " P-Value =", pval, " NSim =", nsim), side = 3, font = 1, line = 1)
    if (legend) 
        legend(legendpos, c("Simulated distribution", "Simulated expectation", indexname), 
               col = c(coldist, coldist, colind), lty = 1, lwd = c(2, 1, 2), 
               bty = bty, cex = cex.legend)
}
