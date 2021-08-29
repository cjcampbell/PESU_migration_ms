# Load objects -------------------------------------------------------------------

mydata_transformed <- readRDS( file.path(wd$bin, "mydata_transformed.rds") )
maps_df    <- readRDS( file.path(wd$bin, "maps_df.rds") )

# Assemble lat/lon coordinates --------------------------------------------

# Because of how distance is most accurately measured (e.g., by distGeo), I'm
# going to convert my surface coordinates (currently in equal area projections)
# to lat/lon.

# Get unique coordinates to convert
coords_aea <- maps_df %>%
  dplyr::select(x,y) %>%
  distinct()
coords_dd <- coords_aea %>%
  SpatialPointsDataFrame(coords = .[, 1:2], proj4string = CRS(myCRS)) %>%
  st_as_sf() %>%
  st_transform(crs = 4326) %>%
  st_coordinates() %>%
  as.data.frame %>%
  dplyr::rename(x_dd = X, y_dd = Y)

mydata_xyz <- cbind(coords_aea, coords_dd) %>%
  right_join(maps_df, by = c("x", "y")) %>%
  full_join(mydata_transformed, by = "ID")
saveRDS(mydata_xyz, file = file.path(wd$bin, "mydata_xyz.rds"))


# Functions ----------------------------------------------------------------
# Function that extracts distance and bearing relating each potential cell of
# geographic origin with individual's sample site.

getDistanceDirection <- function(
  rowNumber, dataframe, ID, fromLat, toLat, fromLon, toLon, getDistance = TRUE,
  getDirection = TRUE, roundTo = 2
){
  
  p1 <- c( dataframe[ rowNumber, fromLon ], dataframe[ rowNumber, fromLat ] )
  p2 <- c( dataframe[ rowNumber, toLon ],   dataframe[ rowNumber, toLat ] )
  myResults <- list()
  if( getDistance == TRUE ){
    dist_km   <- round( geosphere::distGeo(p1, p2) / 1000 , roundTo) #Convert to km, round.
    myResults <- cbind(myResults, dist_km)
  }
  if( getDirection == TRUE ){
    theta_from_site   <- round( geosphere::bearing(p2, p1), roundTo)
    theta_from_origin <- round( geosphere::bearing(p1, p2), roundTo)
    myResults         <- cbind(myResults, theta_from_site, theta_from_origin)
  }
  return(myResults)
  
}

# Apply -------------------------------------------------------------------

# Measure distance and direction from every cell to each sampling location.

# Set to 'FALSE' is you don't want to run in parallel.
nclusters <- parallel::detectCores() 

# First make df with unique combinations of from/to variables.
mydata_FromTo <- mydata_xyz %>%
  dplyr::select(decimalLatitude, decimalLongitude, x_dd, y_dd) %>%
  dplyr::distinct()

distDirection_list <- pbmcapply::pbmclapply( # mclapply, alternatively
  FUN = getDistanceDirection, mc.cores = nclusters,
  1:nrow(mydata_FromTo),
  dataframe =  mydata_FromTo,
  fromLat = "decimalLatitude", toLat = "y_dd", fromLon = "decimalLongitude", toLon = "x_dd",
  getDistance = TRUE, getDirection = TRUE
)

# Here I do some hackey stuff to minimize file size a bit...
wd$tmp_dist <- file.path(wd$bin,"tmp_dist")
if(!dir.exists(wd$tmp_dist) ) dir.create(wd$tmp_dist)
if(length(list.files(wd$tmp_dist)) > 1 ) stop("Files already exist in this directory!")

chunks <- 50
start <- floor( seq(1, nrow(mydata_FromTo), length.out = chunks)[1:(chunks-1)] )
stop <- c( start[2:(chunks-1)]-1, nrow(mydata_FromTo))

getDistDirection_chunked <- function(start, stop){
  
  mdf <- pbmcapply::pbmclapply( # mclapply, alternatively
    FUN = getDistanceDirection, mc.cores = nclusters,
    start:stop,
    dataframe =  mydata_FromTo,
    fromLat = "decimalLatitude", toLat = "y_dd", fromLon = "decimalLongitude", toLon = "x_dd",
    getDistance = TRUE, getDirection = TRUE
  ) %>%
    lapply(as.data.frame) %>%
    plyr::ldply() %>%
    data.frame(mydata_FromTo[start:stop, ], .)
  saveRDS(mdf, file = file.path(wd$tmp_dist, paste0("distDir_", start,"_", stop,".rds")))
  print( paste0( signif(stop/nrow(mydata_FromTo)*100 , 2), "% of the way done"))
}

mapply(getDistDirection_chunked, start = start, stop = stop)

mydata_distDir_loaded <- list.files(wd$tmp_dist, full.names = TRUE) %>%
  lapply(readRDS) %>%
  plyr::ldply()
saveRDS(mydata_distDir_loaded, file = file.path(wd$bin, "mydata_distDir_loaded.rds"))

# Stitch together chunks of mydata_distDir_loaded into mydata_xyz.
if(!exists("mydata_xyz")) mydata_xyz <- readRDS( file.path(wd$bin, "mydata_xyz.rds") )
if(!exists("mydata_distDir_loaded")) mydata_distDir_loaded <- readRDS( file.path(wd$bin, "mydata_distDir_loaded.rds") )

system.time({
  mydata_distDir <- mydata_distDir_loaded %>%
    dplyr::mutate(
      theta_from_site = unlist(theta_from_site),
      theta_from_origin = unlist(theta_from_origin)
    ) %>% 
    data.table::merge.data.table(., mydata_xyz, by = c("decimalLatitude", "decimalLongitude", "x_dd", "y_dd") ) 
})

saveRDS(mydata_distDir, file = file.path(wd$bin, "mydata_distDir.rds"))


# Find prob-of-origin at sample site --------------------------------------

# To threshold conclusions about min dist traveled, find the probability-of-
# origin for each individual.

# Find PoO at sample site for all individuals.
mydata_transformed <- readRDS(file.path(wd$bin, "mydata_transformed.rds"))

assignmentRasts <- file.path(wd$bin) %>%
  list.files(
    path = ., pattern = "Combined.*grd", full.names = TRUE, recursive = TRUE
  ) %>%
  lapply(raster::stack) %>%
  raster::stack()

probs_at_site_list <- pbmcapply::pbmclapply(
  1:nrow(mydata_transformed),
  FUN = function(i){
    thisRast <- assignmentRasts[[ mydata_transformed[i,"ID"] ]]
    pt <- SpatialPoints(coords = cbind( mydata_transformed[i,"metersLongitude"], mydata_transformed[i,"metersLatitude"] ) )
    
    probVal_raw <- raster::extract(thisRast, pt)
    probVal_OR <- isocat::oddsAtSamplingLocation(thisRast, Lat = pt@coords[[2]], Lon = pt@coords[[1]])
    probVal_quant <- isocat::quantileAtSamplingLocation(thisRast, Lat = pt@coords[[2]], Lon = pt@coords[[1]])
    
    data.frame(ID = mydata_transformed[i,"ID"], probVal_raw, probVal_OR, probVal_quant)
  },
  mc.cores = 4
)

probs_at_site_df <- probs_at_site_list %>% plyr::ldply()
saveRDS(probs_at_site_df, file = file.path(wd$bin, "probs_at_site_df.rds"))

mydata_PoOatSampleSite <- probs_at_site_df %>%
  full_join(mydata_transformed) %>%
  tidyr::pivot_longer(cols = starts_with("probVal_"), names_to = "method", values_to = "valAtSampleSite")


# Make odds-simulation surfaces ---------------------------------------

# Parametric bootstrapping
oddsvals <- mydata_PoOatSampleSite %>%
  filter(Season == "Summer", method == "probVal_OR") %>%
  #dplyr::filter(dist2centroid <= 100) %>% 
  dplyr::select(valAtSampleSite) %>%
  unlist %>% 
  na.omit

set.seed(42)
sampledoddsvals <- sample(oddsvals, size = 1e5, replace = T)


# Conduct odds-simulation to make surfaces using odds of
# known_origin individuals.

if(!exists("mydata_distDir")) {
  mydata_distDir <- readRDS( file.path(wd$bin, "mydata_distDir.rds") )
}

distDir_OR <- mydata_distDir %>% 
  dplyr::filter(method == "raw") %>% 
  dplyr::select(-method) %>% 
  dplyr::group_by(ID) %>% 
  dplyr::mutate(
    maxVal = max(value),
    dist_km = unlist(dist_km)
    ) %>% 
  dplyr::mutate(OR = (value/(1-value))/(maxVal/(1-maxVal)) )

distDir_ORSim <- distDir_OR 
distDir_ORSim$OR_sim <- pbmcapply::pbmclapply(
  1:nrow(distDir_ORSim), 
  mc.cores = 4,
  function(i) {
  sum( sampledoddsvals <= distDir_OR$OR[i] ) / length(sampledoddsvals)
  }) %>% 
  unlist

saveRDS(distDir_ORSim, file = file.path(wd$bin, "distDir_ORSim.rds"))

minDistDeets0 <- pbmcapply::pbmclapply(
  seq(0,1,by=0.01), 
  mc.cores = parallel::detectCores(),
  function(j) {
    distDir_ORSim %>% 
      dplyr::mutate(threshold = j) %>% 
      dplyr::filter(OR_sim >= j) %>% 
      dplyr::group_by(ID, threshold) %>% 
      dplyr::summarise(minDist = min(dist_km, na.rm = T))
  }) %>% 
  bind_rows()
minDistDeets <- full_join(minDistDeets0, mydata)
saveRDS(minDistDeets, file = file.path(wd$bin, "minDistDeets.rds"))

