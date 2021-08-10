
# Find location of csv file we want to load.
which_file <- list.files( path = wd$data, pattern = "dD_analyzed_obscuredCoords", full.names = TRUE)

# Load the file.
mydata0 <- read.csv( which_file, header = T,  na.strings=c("","NA") )

mydata1 <- mydata0%>% 
# Convert to projection.
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  st_transform(crs = myCRS) %>%
  st_coordinates() %>%
  as.data.frame %>%
  dplyr::rename(
    metersLatitude  = Y,
    metersLongitude = X
  ) %>%
  cbind(mydata0) %>%
  dplyr::select(ID, ends_with("tude"), everything() )

# Optionally, check it out and make sure it looks okay.
# View(mydata)

# Transform data based on information on standards ------------------------

adjust_OldSA.1_H_1 <- function(d2H) {
  # Function that serves as a wrapper for assignR::refTrans for a specific 
  # standard scale.
  # 
  # WSU samples were analyzed with respect to the following standards: "cow hoof, chicken feathers and whale baleen"
  # 
  # Arbitrary sd (currently not used in this wrapper). 
  s <- data.frame("d2H" = d2H, "d2H.sd" = rep(2), "d2H_cal" = rep("OldEC.1_H_1"))

  d1 <- assignR::refTrans(s,  marker = "d2H")
  return(d1$data$d2H)
  
}

mydata2 <- mydata1 %>% 
  dplyr::mutate(
    fur_adjusted = case_when(
      analysisLab == "CASIF" ~ fur,
      analysisLab == "UWO" ~ adjust_OldSA.1_H_1(fur)
    )
  ) 

mydata <- mydata2 %>% 
  dplyr::mutate(
    cave = case_when(
      Region == "Summer" ~ as.character(NA) ,
      !is.na(Region) ~ gsub("^.*\\.","", ID)
    )
  )
