
# Set options ...
# Ignore some gdal warnings:
options("rgdal_show_exportToProj4_warnings"="none")

# Load required libraries
library(raster)
library(rgdal)
library(rgeos)
library(sf)

# ---

# Country code for which admin boundary will be retrieved:
country_code <- 'PER'

# Cell size used for sampling, in meters:
cell_size <- 8000
# Shape file directory:
shp_dir <- '.'
# Shape file containing biome data to load:
shp_file <- 'PERU_WWF'

# How many times to loop / sample the data:
loop_count <- 1000
# Sample size for each loop:
sample_size <- 50

# Output CSV file:
csv_file = 'biome_areas.csv'

# ---

# Function to get EPSG code for UTM zone from a coordinate, where longitude is
# first coordinate and latitude is second coordinate.
# More details : https://stackoverflow.com/a/9188972
ll2utm <- function(coord) {
  utm_zone = (floor((coord[1] + 180) / 6) %% 60) + 1
  if (coord[2] > 0) {
    epsg_code = utm_zone + 32600
  } else {
    epsg_code = utm_zone + 32700
  }
  return(epsg_code)
}

# ---

# Get admin boundary using getData:
aoi_shp <- getData(country = country_code, level = 0)
# Get the centroid coordinates of the boundary, using gCentroid.
# This is used to work out the UTM zone:
aoi_cntr <- coordinates(gCentroid(aoi_shp))
# Get the CRS for the UTM zone:
aoi_crs <- st_crs(
  ll2utm(aoi_cntr))

# CRS proj4 sring:
aoi_proj <- aoi_crs$proj4string
# Transform boundary to UTM CRS:
aoi_shp <- spTransform(aoi_shp, CRSobj=aoi_proj)
# Convert to sf object:
aoi_sf <- st_as_sf(aoi_shp)

# Read in the shapefile data:
data_shp <- readOGR(shp_dir, shp_file)

#remove intersections 
data_shp <- gBuffer(data_shp, byid=TRUE, width=0)

# Transform to match boundary crs:
data_shp <- spTransform(data_shp, CRSobj=aoi_proj)
# Convert to sf object to enable use of intersection function:
data_sf <- st_as_sf(data_shp)

# Create a union of all polygons in data, so we can get all grid polygons
# which are contained entirely within any of the data types:
data_union <- gUnaryUnion(data_shp)
# Convert to sf for intersecting:
data_union_sf <- st_as_sf(data_union)

# Create a cell of requested size:
grid_cell <- c(cell_size, cell_size)

# Create a grid of cells of this size, which fills the bounds of the shapefile
# data:
bbox_grid <- makegrid(bbox(data_shp), cellsize=grid_cell)
# Convert the grid to SpatialPoints:
grid_spo <- SpatialPoints(bbox_grid, proj4string=CRS(aoi_proj))

# Delete bbox_grid and free the memory:
remove(bbox_grid)
gc()

# Convert grid points to SpatialPixels ...
grid_spi <- SpatialPixels(grid_spo)

# Delete grid_spo and free the memory:
remove(grid_spo)
gc()

# Get the BIOME ids which are in use:
biome_ids <- sort(unique(data_shp$BIOME))
# Create matrix for storing results:
biome_areas <- matrix(ncol=length(biome_ids), nrow=loop_count)
# Set the column names to BIOME_ids:
colnames(biome_areas) <- biome_ids
# Set inital values to 0:
biome_areas[] <- 0

# If you want to get the same results every time the code is run, set a
# seed value:
# set.seed(123)

# For specified number of loops:
for (i in seq_len(loop_count)) {
  # While the sample_count is less than requested size:
  sample_count <- 0
  while (sample_count < sample_size) {
    # Sample a single spatial pixel:
    sample_spi <- sample(grid_spi, 1)
    # Then convert to SpatialPolygons:
    sample_polys <- as(sample_spi, "SpatialPolygons") 
    # Convert to sf object for intersection function:
    sample_sf <- st_as_sf(sample_polys)
    # Check if this polygon intersects with data:
    sample_intersects <- st_intersects(data_union_sf, sample_sf) 
    # If not, move on:
    if (length(sample_intersects[[1]]) == 0) {
      next
    }
    # Intersect the sample with the shapefile data:
    intersect <- st_intersection(sample_sf, data_sf)
    # Add the area to the intersect data:
    intersect$area <- st_area(intersect$geometry)
    # For each BIOME id, get sum of areas for this sample:
    for (biome_id in biome_ids) {
      # Sum areas for this biome id:
      biome_area <- sum(
        intersect$area[intersect$BIOME == biome_id]
      )
      # Add the value:
      biome_areas[i, paste(biome_id)] <- biome_areas[i, paste(biome_id)] +
        as.double(biome_area)
    }
    # Increment sample_count:
    sample_count <- sample_count + 1
  }
}

# Write area data to csv:
write.csv(biome_areas, csv_file)