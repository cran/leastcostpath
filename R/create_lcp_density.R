#' create_lcp_density
#'
#' Creates a cumulative least cost path raster
#'
#' Cumulatively combines Least Cost Paths in order to identify routes of preferential movement within the landscape.
#'
#' @param lcps \code{SpatialLines} or \code{SpatialLinesDataFrame}. Least Cost Paths
#'
#' @param raster \code{RasterLayer} (raster package). This is used to derive the resolution, extent, and spatial reference system to be used when calculating the cumulative least cost path raster
#'
#' @param rescale if TRUE raster values scaled to between 0 and 1. Default is FALSE
#'
#' @return \code{RasterLayer} (raster package). The resultant object is the cumulatively combined Least Cost Paths. This identifies routes of preferential movement within the landscape.
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
#'
#'r <- raster::raster(nrow=50, ncol=50,  xmn=0, xmx=50, ymn=0, ymx=50, crs='+proj=utm')
#'
#'r[] <- stats::runif(1:length(r))
#'
#'slope_cs <- create_slope_cs(r, cost_function = 'tobler')
#'
#'x1 <- c(seq(1,10), seq(11,25), seq(26,30))
#'y1 <- c(seq(1,10), seq(11,25), seq(26,30))
#'line1 <- sp::SpatialLines(list(sp::Lines(sp::Line(cbind(x1,y1)), ID='a')))
#'
#'x2 <- c(seq(1,10), seq(11,25), seq(26, 30))
#'y2 <- c(seq(1,10), seq(11,25), rep(25, 5))
#'line2 <- sp::SpatialLines(list(sp::Lines(sp::Line(cbind(x2,y2)), ID='b')))
#'
#'lcp_network <- rbind(line1, line2)
#'
#'cumulative_lcps <- create_lcp_density(lcps = lcp_network, raster = r, rescale = FALSE)

create_lcp_density <- function(lcps, raster, rescale = FALSE) {
    
    if (!inherits(lcps, c("SpatialLines", "SpatialLinesDataFrame"))) {
        stop("lcps expects a SpatialLines or SpatialLinesDataFrame object")
    }
    
    if (!inherits(raster, "RasterLayer")) {
        stop("raster expects a RasterLayer object")
    }
    
    lcp_pts <- methods::as(lcps, "SpatialPoints")
    
    cumulative_pts <- raster::rasterize(x = lcp_pts, y = raster, fun = "count")
    
    if (rescale) {
        
        rasterRescale <- function(r) {
            ((r - raster::cellStats(r, "min"))/(raster::cellStats(r, "max") - raster::cellStats(r, "min")))
        }
        
        cumulative_pts <- rasterRescale(cumulative_pts)
    }
    
    cumulative_pts[is.na(cumulative_pts)] <- 0
    
    return(cumulative_pts)
    
}