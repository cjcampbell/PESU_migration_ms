
# Load packages.
if(packageVersion("isocat") < "0.2.4.9000") warning("Install development version of isocat.")
# devtools::install_github("cjcampbell/isocat")

# Load libraries required for analyses.
library(isocat)
library(dplyr)
library(sf)

# Load libraries relied on for a few analyses.
library(assignR)
library(chron)
library(geosphere)
library(ggmap)
library(ggplot2)
library(lubridate)
library(lwgeom)
library(measurements)
library(purrr)
library(readxl)
library(rgdal)
library(smatr)
library(stringr)
library(tidyr)

# Make an object to help navigate the subdirectories.
my_dir_path <- "/Users/cjcampbell/PESU_migration_project_directory/PESU_migration_ms"
wd <- list()
wd$R       <- file.path( my_dir_path, "R" )
wd$bin     <- file.path( my_dir_path, "bin" )
wd$data    <- file.path( my_dir_path, "data" )
wd$iucn    <- file.path( my_dir_path, "data", "iucn" )
wd$figs    <- file.path( my_dir_path, "figs" )
wd$out     <- file.path( my_dir_path, "out" )

# Check for presence of subdirectories. Create if needed.
invisible({
  lapply(wd, function(i) if( dir.exists(i) != 1 ) dir.create(i) )
})


# Define extent of spatial analysis.
# Unit == meters
my_extent_aea <- raster::extent(
  -16e5, 34e5,
  -30e5, 16e5
)

# CRS for aea projection:
myCRS <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs +towgs84=0,0,0"
