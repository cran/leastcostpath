% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rasterise.R
\name{rasterise}
\alias{rasterise}
\title{converts conductanceMatrix to SpatRaster}
\usage{
rasterise(x)
}
\arguments{
\item{x}{\code{conductanceMatrix}}
}
\value{
\code{spatRaster}
}
\description{
converts conductanceMatrix to SpatRaster
}
\examples{

r <- terra::rast(system.file("extdata/SICILY_1000m.tif", package="leastcostpath"))

slope_cs <- create_slope_cs(x = r, cost_function = "tobler", neighbours = 4)

cs_rast <- rasterise(slope_cs)
}
\author{
Joseph Lewis
}
