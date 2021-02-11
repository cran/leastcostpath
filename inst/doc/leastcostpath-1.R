## ----set-options, echo=FALSE, cache=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
options(width = 500)

## ----libraries, echo = TRUE, message= FALSE, warning= FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(rgdal)
library(rgeos)
library(sp)
library(raster)
library(gdistance)
library(leastcostpath)

## ----raster, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dem <- raster::raster(system.file('external/maungawhau.grd', package = 'gdistance'))

dem_extent <- as(raster::extent(dem), 'SpatialPolygons')

raster::crs(dem_extent) <- raster::crs(dem)

raster::plot(dem)
raster::plot(dem_extent, add = T, border = "red")

## ----cs, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# cost functions currently implemented within leastcostpath

cfs <- c("tobler", "tobler offpath", "irmischer-clarke male", "irmischer-clarke offpath male", "irmischer-clarke female", "irmischer-clarke offpath female",
        "modified tobler", "wheeled transport", "herzog", "llobera-sluckin", "campbell 2019")

# neighbours can be 4, 8, 16, 32, or 48. 8 used for illustration purposes, but greater the number the better the cost surface / LCP approximates reality. 
neigh <- 4

slope_cs <- leastcostpath::create_slope_cs(dem = dem, cost_function = "tobler", neighbours = neigh)

plot(raster(slope_cs), col = grey.colors(100))

## ----barrier, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# create barrier of altitude values above 180

altitude <- dem > 180
# convert value 0 to NA - this ensures that the areas below 180 are not viewed as barriers. That is, if the raster value is NA then create_barrier_cs() will assume that the cell is background and will therefore assign the background argument value 
altitude[altitude == 0] <- NA

# the values NOT NA will be assigned the field argument value. If 0 (default) then movement within the area will be completely inhibited.
altitude_cs <- leastcostpath::create_barrier_cs(raster = altitude, barrier = altitude, neighbours = neigh, field = 0, background = 1)

# multiplying the two cost surfaces ensures that barriers continue to completely inhibit movement (i.e. slope_cs values * 0 = 0)
slope_altitude_cs <- slope_cs * altitude_cs

plot(raster(slope_altitude_cs), col = grey.colors(100))

## ----lcp = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
A <- sp::SpatialPoints(cbind(2667930, 6479466))
B <- sp::SpatialPoints(cbind(2667668, 6478818))

# if cost function is anisotropic (e.g. Tobler's Hiking Function) then LCP from A-to-B and B-to-A may be different. To create LCP in both directions set the directional argument to FALSE (default)
lcp <- leastcostpath::create_lcp(cost_surface = slope_altitude_cs, origin = A, destination = B, directional = FALSE)

plot(raster(slope_altitude_cs), col = grey.colors(100))
plot(A, add = T, col = "black")
plot(B, add = T, col = "black")
# LCP from A-to-B
plot(lcp[1,], add = T, col = "red")
# LCP from B-to-A
plot(lcp[2,], add = T, col = "blue")

## ----cc, TRUE,  fig.height = 6, fig.width = 6, warning = FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# calculates accumulated cost surface from A and from B. These are then averaged. rescale argument TRUE (not default) rescales the accumulated cost surface to 1
cc <- leastcostpath::create_cost_corridor(cost_surface = slope_altitude_cs, origin = A, destination = B, rescale = TRUE)

plot(cc, col = heat.colors(100))
plot(A, add = T, col = "white", pch = 16)
plot(B, add = T, col = "white", pch = 16)

## ----fete, TRUE,  fig.height = 6, fig.width = 6, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
random_locs <- sp::spsample(x = dem_extent, n = 25, type = 'random')

fete_lcp <- leastcostpath::create_FETE_lcps(cost_surface = slope_altitude_cs, locations = random_locs, cost_distance = FALSE, parallel = FALSE)

plot(raster(slope_altitude_cs), col = grey.colors(100))
plot(fete_lcp, add = T, col = "red")
plot(random_locs, add = T, col = "black", pch = 16)

## ----lcp_d, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# rasterises the LCPs and calculates the number of times an LCP crosses a raster cell. If rescale argument TRUE (not default) then LCP density scaled to 1. If rasterize_as_points argument TRUE (default) then LCP vertices points are rasterized. This is quicker than rasterizing lines (rasterize_as_lines argument FALSE) but will contain 'gaps' where no LCP vertex point is present 
fete_lcp_density <- leastcostpath::create_lcp_density(lcps = fete_lcp, raster = dem, rescale = TRUE, rasterize_as_points = TRUE)

fete_lcp_density[fete_lcp_density == 0] <- NA

plot(fete_lcp_density, col = heat.colors(100))

## ----p_lcp, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# least cost paths are calculated using input data that contains errors. The most notable is the Digital Elevation Model (DEM) which is used to model movement difficulty within the landscape. By incorporating DEM error (in this case vertical error), the error is propagated throughout the LCP modelling process. This results in the input data better approximating the real world, with the LCP being more reliable as a result  

# for illustration purposes, we will assume that the vertical error (root-mean-square-error) is +/- 0.5m.
RMSE <- 0.5
n <- 100
lcps <- list()

for (i in 1:n) {
  
  print(i)

lcps[[i]] <- leastcostpath::create_lcp(cost_surface = leastcostpath::create_slope_cs(dem = leastcostpath::add_dem_error(dem = dem, rmse = RMSE, type = "autocorrelated"), cost_function = "tobler", neighbours = neigh) * altitude_cs, origin = A, destination = B, directional = FALSE)

}

lcps <- do.call(rbind, lcps)

lcps_A_B <- leastcostpath::create_lcp_density(lcps = lcps[lcps$direction == "A to B",], raster = dem, rescale = TRUE, rasterize_as_points = TRUE)

lcps_A_B[lcps_A_B == 0] <- NA

lcps_A_B <- (lcps_A_B/n)

lcps_B_A <- leastcostpath::create_lcp_density(lcps = lcps[lcps$direction == "B to A",], raster = dem, rescale = FALSE, rasterize_as_points = TRUE)

lcps_B_A[lcps_B_A == 0] <- NA

lcps_B_A <- (lcps_B_A/n)

plot(raster(slope_altitude_cs), legend = FALSE, col = grey.colors(100))
plot(lcps_A_B, add = T, col = heat.colors(100))
plot(A, add = T, col = "white", pch = 16)
plot(B, add = T, col = "white", pch = 16)

plot(raster(slope_altitude_cs), legend = FALSE, col = grey.colors(100))
plot(lcps_B_A, add = T, col = heat.colors(100))
plot(A, add = T, col = "white", pch = 16)
plot(B, add = T, col = "white", pch = 16)

