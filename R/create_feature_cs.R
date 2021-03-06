#' Create a Landscape Feature cost surface
#'
#' Creates a Landscape Feature Cost Surface representing the attraction/repulsion of a feature in the landscape. See Llobera (2000) for theoretical discussion in its application
#'
#' @param raster \code{RasterLayer} (raster package). The Resolution, Extent, and Spatial Reference System of the provided RasterLayer is used when creating the resultant Barrier Cost Surface
#'
#' @param locations \code{SpatialPoints*} (sp package). Location of Features within the landscape
#'
#' @param x \code{numeric vector}. Values denoting the attraction/repulsion of the landscape features within the landscape. Each value in the vector is assigned to each ring of cells moving outwards from supplied locations
#'
#' @param neighbours \code{numeric} value. Number of directions used in the Least Cost Path calculation. See Huber and Church (1985) for methodological considerations when choosing number of neighbours. Expected numeric values are 4, 8, 16, 32, 48 or a matrix object. Default is numeric value 16
#'
#' @return \code{TransitionLayer} (gdistance package) numerically expressing the attraction/repulsion of a feature in the landscape. The resultant \code{TransitionLayer} can be incorporated with other \code{TransitionLayer} through Raster calculations.
#'
#' @references Llobera, M. (2000). Understanding movement: a pilot model towards the sociology of movement. In: Lock G (ed) Beyond the map. Archaeology and spatial technologies. (pp 66-84). Amsterdam: IOS Press/Ohmsha.
#'
#' @author Joseph Lewis
#'
#' @import rgdal
#' @import rgeos
#' @import sp
#' @import raster
#' @import gdistance
#'
#' @export
#'
#' @examples
#' r <- raster::raster(system.file('external/maungawhau.grd', package = 'gdistance'))
#' loc1 = cbind(2667670, 6479000)
#' loc1 = sp::SpatialPoints(loc1)
#'
#' num <- seq(200, 1, length.out = 20)
#'
#' feature <- create_feature_cs(raster = r, locations = loc1, x = num)

create_feature_cs <- function(raster, locations, x, neighbours = 16) {
    message("note: create_feature_cs expects planar coordinates")
    
    if (!inherits(raster, "RasterLayer")) {
        stop("raster argument is invalid. Expecting a RasterLayer object")
    }
    
    if (!inherits(locations, c("SpatialPoints", "SpatialPointsDataFrame"))) {
        stop("locations argument is invalid. Expecting a SpatialPoints* object")
    }
    
    if (!inherits(x, "numeric")) {
        stop("x argument is invalid. Expecting a numeric vector object")
    }
    
    if (any(!neighbours %in% c(4, 8, 16, 32, 48)) & (!inherits(neighbours, "matrix"))) {
        stop("neighbours argument is invalid. Expecting 4, 8, 16, 32, 48, or matrix object")
    }
    
    if (inherits(neighbours, "numeric")) {
        if (neighbours == 32) {
            neighbours <- neighbours_32
            
        } else if (neighbours == 48) {
            neighbours <- neighbours_48
        }
        
    }
    
    locations <- as(locations, "SpatialPoints")
    
    r <- raster::rasterize(locations, raster)
    
    r <- raster::distance(r)
    
    mround <- function(x, base) {
        base * round(x/base)
    }
    
    r <- mround(r, res(r)[1]/100)
    
    r_mod <- matrix(cbind(unique(r)[1:length(x)], unique(r)[1:length(x)] + res(r)[1]/100, x), ncol = 3)
    
    r_mod <- r_mod[stats::complete.cases(r_mod), ]
    
    r[r > unique(r)[length(x)]] <- NA
    
    rc <- raster::reclassify(r, r_mod, include.lowest = TRUE, right = FALSE)
    
    rc[is.na(rc)] <- 1
    
    rc <- gdistance::transition(rc, min, neighbours, symm = TRUE)
    
    return(rc)
    
}
