source('~/PESU_migration_project_directory/PESU_migration/R/00_Setup.R')
source('~/PESU_migration_project_directory/PESU_migration/R/01_loadIsotopeData.R')

# Load direction stats.
if(!exists("distDir_ORSim") ) distDir_ORSim <- readRDS( file.path(wd$bin, "distDir_ORSim.rds") )
dirSummary <- distDir_ORSim %>% 
  dplyr::select(-Sex) %>% 
  dplyr::mutate(theta_rounded = round(theta_from_origin, 0)) %>% 
  group_by(ID, theta_rounded, .drop = F) %>% 
  dplyr::summarise(
    mean_rawVal = mean(value, na.rm = T),
    mean_OR = mean(OR,  na.rm = T )
    ) %>% 
  ungroup()

dir_deets <- dirSummary %>% 
  dplyr::select(-mean_OR) %>% 
  pivot_wider(
    names_from = theta_rounded,
    names_prefix = "probOfOriginAtAngle_" ,
    values_from = mean_rawVal
    )

dir_deets_SHORTER <- distDir_ORSim %>% 
  dplyr::mutate(
    direction = case_when(
      theta_from_origin >= -45 & theta_from_origin < 45 ~ "North",
      theta_from_origin >= 45 & theta_from_origin < 135 ~ "East",
      theta_from_origin >= -135 & theta_from_origin < -45 ~ "West",
      TRUE ~ "South"
    )
    ) %>% 
  dplyr::group_by(ID, direction) %>% 
  dplyr::summarise(
    ave_prob = mean(value, na.rm = T)
    ) %>% 
  pivot_wider(
    values_from = ave_prob,
    names_from = direction
  )

# Load distance stats.
if(!exists("minDistDeets") ) minDistDeets <- readRDS( file.path(wd$bin, "minDistDeets.rds") )
dist_deets <- minDistDeets %>% 
  dplyr::select(-Sex) %>% 
  dplyr::select(ID, threshold, minDist) %>% 
  dplyr::filter(threshold %in% c( 0.25, 0.32, 0.5, 0.68, 0.75) ) %>% 
  arrange(ID, threshold) %>% 
  pivot_wider(
    names_from = threshold,
    names_prefix = "probDistTraveledAtThreshold_" ,
    values_from = minDist
  )


# Combine -----------------------------------------------------------------
mydata_transformed <- readRDS( file.path(wd$bin, "mydata_transformed.rds") )

PESUResults <- mydata_transformed %>% 
  dplyr::select(
    ID, verbatimID, contains("decimal"), contains("meters"),
    everything()
    ) %>% 
  full_join(dist_deets) %>% 
  full_join(dir_deets_SHORTER) %>% 
  full_join(dir_deets) 

# Which direction is most likely origin?

directionCols <- which( names(PESUResults) %in% c("East", "North", "South", "West") )
PESUResults$mostLikelyDirection <- names(PESUResults)[directionCols][max.col(PESUResults[,directionCols])]


# Which theta/degree has most likely origin?
angleCols <- grep(pattern = "probOfOriginAtAngle_", names(PESUResults) )

PESUResults <- mutate_at(PESUResults, angleCols, ~replace(., is.na(.), 0)) # replace NA with 0
PESUResults$mostLikelyAngle <- names(PESUResults)[angleCols][max.col(PESUResults[,angleCols])] %>% 
  gsub("probOfOriginAtAngle_","",.) %>% 
  as.numeric()

# Finally combine with clustering info and pull out some summaries.
myResults <- PESUResults %>% 
  dplyr::mutate(
    mostLikelyDiagonal = case_when(
      mostLikelyAngle >= 0 &  mostLikelyAngle < 90 ~ "NE",
      mostLikelyAngle >= -90 &  mostLikelyAngle < 0 ~ "NW",
      mostLikelyAngle >= -180 &  mostLikelyAngle < -90 ~ "SW",
      mostLikelyAngle >= 90 &  mostLikelyAngle < 180 ~ "SE"
    ),
    mostLikelyDiagonal = factor(mostLikelyDiagonal, levels = c("NE", "SE", "SW", "NW")),
    NorthOrSouth = case_when(
      mostLikelyDiagonal %in% c("NW", "NE") ~ "North",
      mostLikelyDiagonal %in% c("SW", "SE") ~ "South"),
    distanceTraveled = case_when(
      probDistTraveledAtThreshold_0.25 < 100 ~ "Resident",
      probDistTraveledAtThreshold_0.25 >= 100 & 
        probDistTraveledAtThreshold_0.25 <1000 ~ "Regional",
      probDistTraveledAtThreshold_0.25 >= 1000 ~"Long-distance"
    ),
    distanceTraveled = factor(distanceTraveled, levels = c("Resident", "Regional", "Long-distance")),
    distDir = case_when(
      distanceTraveled == "Resident" ~ "Resident",
      distanceTraveled == "Regional" & NorthOrSouth == "North" ~ "Regional-N",
      distanceTraveled == "Regional" & NorthOrSouth == "South" ~ "Regional-S",
      distanceTraveled == "Long-distance" & NorthOrSouth == "North" ~ "Long-distance-N",
      distanceTraveled == "Long-distance" & NorthOrSouth == "South" ~ "Long-distance-S"
    ),
    Region_long = case_when(
      Region == "NW" ~ "Northwest",
      Region == "NC" ~ "North-central"
    ),
    Region_long = factor(Region_long, levels = c("Northwest", "North-central")),
    # Specify season_lab:
    seasonLab = paste(Season, analysisLab, sep = "_"),
    seasonLab = factor(seasonLab, levels=c("Summer_UWO", "Summer_CASIF", "Winter_CASIF"))
  ) %>% 
  # Remove angle data...
  dplyr::select_at(vars(!starts_with("probOfOriginAtAngle"))) %>% 
  dplyr::select_at(vars(!starts_with("propSimulations")))

write.csv(myResults, file = file.path(wd$bin, "myResults.csv"), row.names = F)
