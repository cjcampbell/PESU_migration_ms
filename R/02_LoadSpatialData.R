
# Setup -------------------------------------------------------------------
download_GADM <- FALSE

# This script assumes that there are candidate isoscapes downloaded from isoMAP
# somewhere in the wd$data directory.
reload_isoscapes <- FALSE

# This script assumes that there are IUCN rangemaps somewhere in the wd$iucn
# directory.
reload_IUCN_rangemaps <- FALSE

# Load GADM data ----------------------------------------------------------
if(download_GADM == TRUE){

  library(rmapshaper)
  message("Loading GADM data...")

  # First, make NoAm_sf (an "sf" object) fit extent object (my_extent).
  # my_extent_aea_st <- st_bbox(st_transform(st_as_sfc(st_bbox(my_extent, crs = 4326)), myCRS))
  # saveRDS(my_extent_aea_st, file = file.path(wd$bin, "my_extent_aea_st.rds"))

  # Get GADM data to state level.
  USA <- raster::getData('GADM', path = wd$bin, country='USA', level=0)
  MEX <- raster::getData('GADM', path = wd$bin, country='MEX', level=0)
  CAN <- raster::getData('GADM', path = wd$bin, country='CAN', level=0)
  GTM <- raster::getData('GADM', path = wd$bin, country='GTM', level=0)

  # Prepare to remove areas outside of desired extent.
  # ## Remove Hawaii
  USA_1 <- raster::getData('GADM', path = wd$bin, country='USA', level=1)
  hawaii <- USA_1[USA_1@data$NAME_1 == "Hawaii",] # Select Hawaii
  hawaii_simpl <- hawaii %>%
    sf::st_as_sf() %>%
    st_transform(crs = myCRS) %>%
    st_simplify(preserveTopology = FALSE, dTolerance = 5000) %>%
    st_buffer(dist = 1e6)
  # ## Remove water bodies.
  USA_2 <- raster::getData('GADM', path = wd$bin, country='USA', level=2)
  MEX_2 <- raster::getData('GADM', path = wd$bin, country='MEX', level=2)
  CAN_2 <- raster::getData('GADM', path = wd$bin, country='CAN', level=2)
  GTM_2 <- raster::getData('GADM', path = wd$bin, country='GTM', level=2)
  waterbodies <- lapply(list(USA_2,MEX_2,CAN_2,GTM_2), function(x){
    x[(x$ENGTYPE_2) == "Water body",]
  }) %>%
    do.call(rbind, .) %>%
    sf::st_as_sf() %>%
    st_transform(crs = myCRS) %>%
    st_simplify(preserveTopology = FALSE, dTolerance = 5000)

  # Combine into one polygon, convert to sf object.
  NoAm <- raster::bind(
    MEX, USA, CAN, GTM#, BLZ, SLV, HND, NIC
  ) %>% 
    sf::st_as_sf(.) %>%
    st_transform(crs = myCRS) %>%
    st_simplify(preserveTopology = FALSE, dTolerance = 5e3) %>%
    st_difference(., hawaii_simpl) # Remove Hawaii
  
  saveRDS(NoAm, file = file.path(wd$bin, "NoAm.rds"))
  
  NoAm_boundary_aea <- NoAm %>% 
    st_buffer(dist = 5e4) %>%
    rmapshaper::ms_erase(., waterbodies)  # Remove water bodies

  saveRDS(NoAm_boundary_aea, file = file.path(wd$bin, "NoAm_boundary_aea.rds"))

} else message("Not redownloading GADM Data...")


# Load isoscapes ----------------------------------------------------------

if(reload_isoscapes == TRUE){
  message("reloading isoscapes...")

  NoAm_boundary_aea <- readRDS( file.path(wd$bin, "NoAm_boundary_aea.rds") )

  # Function to extend/crop/mask by above.
  ECM <- function(rasterLayer){
    rasterLayer %>%
      raster::projectRaster(., crs = myCRS) %>%
      raster::extend( ., my_extent_aea ) %>%
      raster::crop(   ., my_extent_aea ) %>%
      raster::resample(., refIsoscape) %>%
      raster::mask(   ., NoAm_boundary_aea  )
  }

  # Load then write assignR's isoscape. Could call directly, but this works okay.
  assignRpathPattern <- "assignR_d2h_world"
  GGS_H <- assignR::d2h_world
  if( !dir.exists(file.path(wd$data, assignRpathPattern)) ) {
    dir.create(file.path(wd$data, assignRpathPattern))
  }
  writeRaster(GGS_H[[1]], overwrite = TRUE,
              filename = file.path(wd$data, assignRpathPattern, "isoscape.tif"))
  writeRaster(GGS_H[[2]], overwrite = TRUE,
              filename = file.path(wd$data, assignRpathPattern, "sd.tif"))

  # Load then crop and mask environmental data rasters.
  loadAdjustIsoscapeTiff <- function( directory, path_pattern, isoscape_pattern,
                                      sd_pattern, refIsoscape){
    l <- list()
    l$directory          <- directory
    l$path_pattern       <- path_pattern
    l$isoscape_pattern   <- isoscape_pattern
    l$sd_pattern         <- sd_pattern
    l$reference_isoscape_name <- names(refIsoscape)

    l$isoscape <- list.files(
      directory, recursive = TRUE, pattern = isoscape_pattern, full.names = TRUE
    ) %>%
      grep(path_pattern, ., value = TRUE) %>%
      raster::raster(.) %>%
      ECM(.)

    l$sd <- list.files(
      directory, recursive = TRUE, pattern = sd_pattern, full.names = TRUE
    ) %>%
      grep(path_pattern, ., value = TRUE) %>%
      raster::raster(.) %>%
      ECM(.)

    return(l)
  }

  suppressWarnings({
    
    refIsoscape <- file.path(wd$data, assignRpathPattern, "isoscape.tif") %>%
      raster::raster(.) %>%
      raster::projectRaster(., crs = myCRS) %>%
      raster::extend( ., my_extent_aea ) %>%
      raster::crop(   ., my_extent_aea )

    my_isoscapes <- mapply(
      FUN = loadAdjustIsoscapeTiff,
      path_pattern = c("70052", "70055", "70064","73715", "66098", "66100", assignRpathPattern),
      isoscape_pattern = c(rep("predkrig.tiff$", 6), "isoscape.tif$"),
      sd_pattern = c(rep("stdkrig.tiff$", 6), "sd.tif$"),
      MoreArgs = list(
        directory = wd$data, 
        refIsoscape = refIsoscape
        ),
      SIMPLIFY = FALSE
    )
    
  })

  # Check if everything turned out okay.
  lapply(my_isoscapes, function(i) c( i$isoscape, i$sd) ) %>%
    unlist %>%
    lapply(., compareRaster, x = .[[1]] ) %>%
    unlist %>%
    all %>%
    {if(.!=TRUE) stop("Something is wrong! CompareRaster yeilds FALSE results.")}

  # Save.
  save(my_isoscapes, file = file.path(wd$bin, "my_isoscapes.RData"))
} else message("Not reloading isoscapes, loading saved version...")


# Load and buffer IUCN Rangemaps -----------------------------------------------------
if(reload_IUCN_rangemaps == TRUE){
  load(file.path(wd$bin, "my_isoscapes.RData"), verbose = TRUE)
  NoAm_boundary_aea <- readRDS( file.path(wd$bin, "NoAm_boundary_aea.rds") )

  # Crop Polygon above the tropics.
  NoAm_aboveTropics <- st_multipoint(rbind(c(-110, 23.4366), c(-50,50)), crs) %>% 
    st_sfc() %>% 
    st_set_crs(., "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>% 
    st_transform(myCRS) %>% 
    st_bbox() %>% 
    st_as_sfc() %>%
    st_intersection(., NoAm_boundary_aea)
  
  # Load candidate rangemaps.
  # Convert to simple features object and reproject to myCRS.
  range_PESU <- list.dirs(wd$iucn) %>% 
    grep(pattern = "species_data", value = T) %>% 
    rgdal::readOGR(dsn = ., layer = "data_0") %>% 
    st_as_sf(crs = 4326) %>%
    st_transform(crs = myCRS) %>%
    st_simplify(preserveTopology = TRUE, dTolerance = 5000) %>%
    st_make_valid() %>%
    st_buffer(dist = 10e3) %>% 
    st_intersection(., NoAm_aboveTropics) %>%
    st_combine()
  
  # Check that all points fall in the polygon.
  pts <- mydata %>%
    dplyr::select(metersLongitude, metersLatitude) %>%
    sp::SpatialPoints(., proj4string = crs(myCRS)) %>%
    st_as_sf()
  myDists <- sf::st_distance(range_PESU, pts)
  if(sum(as.numeric(myDists)) > 0) stop("Some points don't fall on the polygon! Uh-oh!")

  # Convert buffered rangemaps to rasters with appropriate
  ex_rast <- my_isoscapes[[7]]$isoscape
  ex_rast[] <- 1
  range_raster <- raster::mask(
    ex_rast, 
    mask = as_Spatial(range_PESU), 
    updatevalue = NA
    ) %>%
    raster::mask(., NoAm_boundary_aea) %>%
    raster::crop(., my_extent_aea)
  
  save( range_raster, file = file.path(wd$bin, "range_raster.Rdata" ) )
} else message("Not reloading IUCN rangemaps...")
