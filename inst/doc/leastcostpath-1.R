## ----libraries, echo = TRUE, message= FALSE, warning= FALSE-------------------
library(rgdal)
library(rgeos)
library(sp)
library(raster)
library(gdistance)
library(leastcostpath)

## ----raster, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE-----
r <- raster::raster(system.file('external/maungawhau.grd', package = 'gdistance'))

r_extent <- as(raster::extent(r), 'SpatialPolygons')

crs(r_extent) <- crs(r)

plot(r)
plot(r_extent, add = T, border = "red")

## ----slope_cs, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE----
cs <- create_slope_cs(dem = r, cost_function = 'tobler', neighbours = 16)

## ----lcp, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE--------

loc1 = cbind(2667876, 6479424)
loc1 = sp::SpatialPoints(loc1)

loc2 = cbind(2667677, 6478737)
loc2 = sp::SpatialPoints(loc2)

lcp <- create_lcp(cost_surface = cs, origin = loc1, destination = loc2, directional = FALSE)

plot(raster(cs))
plot(lcp[1,], add = T, col = "red")
plot(lcp[2,], add = T, col = "blue")

## ----slope_traverse_cs, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE----
cs <- create_slope_cs(dem = r, cost_function = 'tobler', neighbours = 16) %>%
  "*" (create_traversal_cs(dem = r, neighbours = 16))

## ----lcp_2, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE------

loc1 = cbind(2667876, 6479424)
loc1 = sp::SpatialPoints(loc1)

loc2 = cbind(2667677, 6478737)
loc2 = sp::SpatialPoints(loc2)

lcp <- create_lcp(cost_surface = cs, origin = loc1, destination = loc2, directional = FALSE)

plot(raster(cs))
plot(lcp[1,], add = T, col = "red")
plot(lcp[2,], add = T, col = "blue")

## ----slope_traverse_feature_cs, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE----
feature_loc = cbind(2667652, 6478997)
feature_loc = sp::SpatialPoints(feature_loc)

x <- seq(100, 1, length.out = 20)

cs <- create_slope_cs(dem = r, cost_function = 'tobler', neighbours = 16) %>%
  "*" (create_traversal_cs(dem = r, neighbours = 16)) %>%
  "*" (create_feature_cs(raster = r, locations = feature_loc, x))

## ----lcp_3, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE------

loc1 = cbind(2667876, 6479424)
loc1 = sp::SpatialPoints(loc1)

loc2 = cbind(2667677, 6478737)
loc2 = sp::SpatialPoints(loc2)

lcp <- create_lcp(cost_surface = cs, origin = loc1, destination = loc2, directional = FALSE)

plot(raster(cs))
plot(feature_loc, add = T, col = "black")
plot(lcp[1,], add = T, col = "red")
plot(lcp[2,], add = T, col = "blue")

## ----slope_traverse_feature_cc, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE----
feature_loc = cbind(2667652, 6478997)
feature_loc = sp::SpatialPoints(feature_loc)

x <- seq(100, 1, length.out = 20)

cs <- create_slope_cs(dem = r, cost_function = 'tobler', neighbours = 16) %>%
  "*" (create_traversal_cs(dem = r, neighbours = 16)) %>%
  "*" (create_feature_cs(raster = r, locations = feature_loc, x))

## ----cc, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE---------

loc1 = cbind(2667876, 6479424)
loc1 = sp::SpatialPoints(loc1)

loc2 = cbind(2667677, 6478737)
loc2 = sp::SpatialPoints(loc2)

cc <- create_cost_corridor(cs, loc1, loc2)

plot(cc)

## ----fete, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE-------
locs <- sp::spsample(r_extent,n=25,'regular')

lcp_network <- create_slope_cs(dem = r, cost_function = 'tobler', neighbours = 16) %>%
  "*" (create_traversal_cs(dem = r, neighbours = 16)) %>%
  create_FETE_lcps(cost_surface = ., locations = locs, cost_distance = FALSE, parallel = FALSE)

plot(r)
plot(locs, add = T)
plot(lcp_network, add = T, col = "red")

## ----lcp_d, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE------
lcp_network_density <- create_lcp_density(lcps = lcp_network, raster = r, rescale = TRUE)

plot(lcp_network_density)

## ----lcp_kd, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE-----
fwModel <- raster::focalWeight(lcp_network_density, max(res(r)) * 2, type='circle')

fwModel[fwModel>0] <- 1

lcp_network_kd <- raster::focal(lcp_network_density, w=fwModel ,fun=sum , na.rm=TRUE)

plot(lcp_network_kd)

## ----ccp, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE--------
locs = cbind(2667652, 6478997)
locs = sp::SpatialPoints(locs)

lcp_network <- create_slope_cs(dem = r, cost_function = 'tobler', neighbours = 16) %>%
  "*" (create_traversal_cs(dem = r, neighbours = 16)) %>%
  create_CCP_lcps(cost_surface = ., location = locs, distance = 100, radial_points = 50, cost_distance = FALSE, parallel = FALSE)

plot(r)
plot(locs, add = T)
plot(lcp_network, add = T, col = "red")

## ----banded_lcp, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE----
locs = cbind(2667652, 6478997)
locs = sp::SpatialPoints(locs)

lcp_network <- create_slope_cs(dem = r, cost_function = 'tobler', neighbours = 16) %>%
  "*" (create_traversal_cs(dem = r, neighbours = 16)) %>%
  create_banded_lcps(cost_surface = ., location = locs, min_distance = 5, max_distance = 250, radial_points = 50, cost_distance = FALSE, parallel = FALSE)

plot(r)
plot(locs, add = T)
plot(lcp_network, add = T, col = "red")

## ----lcp_network, echo = TRUE,  fig.height = 6, fig.width = 6, warning = FALSE----
locs <- sp::spsample(r_extent,n=10,'random')

lcp_network <- create_slope_cs(dem = r, cost_function = 'tobler', neighbours = 16) %>%
  "*" (create_traversal_cs(dem = r, neighbours = 16)) %>%
  create_lcp_network(., locations = locs, nb_matrix = cbind(c(1, 4, 2, 1), c(2, 2, 4, 3)), cost_distance = FALSE, parallel = FALSE)

plot(r)
plot(locs, add = T)
plot(lcp_network, add = T, col = "red")

