# Load  ----------------------------------------------------------

if(!exists("my_isoscapes")) load(file.path(wd$bin, "my_isoscapes.RData"), verbose = TRUE)
if(!exists("mydata")) source(file.path(wd$R, "01_loadIsotopedata.R"))


# Extract vals from isoscapes. -----------------------------------

mydata_isoVals <- lapply(my_isoscapes, function(which_isoscape){
  isoVals <- lapply(1:nrow(mydata), function(n){
    coords = mydata[ n, c("metersLongitude", "metersLatitude") ]
    raster::extract(which_isoscape$isoscape, coords)
  })
  test <- data.frame(mydata, rename_me = unlist(isoVals))
  names(test)[ncol(test)] <- paste( "precip_val", which_isoscape$path_pattern, sep = "_" )
  return(test)
}) %>%
  purrr::reduce(left_join, by = names(mydata)) %>%
  dplyr::mutate_if(is.factor, as.character)

saveRDS(mydata_isoVals, file = file.path(wd$bin, "mydata_isoVals.Rds"))


# Fitting Transfer functions ---------------------------------------------------

mydata_isoVals_knownOrigin <- mydata_isoVals %>%
  dplyr::filter(Season == "Summer")

mydata_isoVals_knownOrigin %>% group_by(Location) %>% dplyr::summarise(n=n())


# Write regression functions ---------------------------------------------------

# Regression function.
doSMA_Regression <- function(data, y.data, x.data, ...){
  
  if(class(y.data) == "character"){ y.data <- data[ , y.data] }
  if(class(x.data) == "character"){ x.data <- data[ , x.data] }
  if(length(unique(y.data)) <= 2 | length(unique(x.data)) <= 2){
    stop("Not enough data to run regression") }
  
  mymodel <- smatr::sma(y.data ~ x.data)
  
  m <- list(...)
  m$intercept <- coef(mymodel)[1]
  m$slope     <- coef(mymodel)[2]
  m$sdRes     <- sd(residuals(mymodel))
  m$n         <- mymodel$groupsummary$n
  m$p         <- mymodel$groupsummary$pval
  m$r2        <- mymodel$groupsummary$r2
  m$model      <- mymodel
  return(m)
}

regWithResampling <- function( data, xCol, colsToResampleBy, savedirectory, min_n){
  
  if(!dir.exists(savedirectory)) dir.create(savedirectory)
  
  sampledData <- data.frame()
  while( nrow(sampledData) < min_n ){
    sampledData  <- data %>%
      group_by_at(vars(all_of(colsToResampleBy))) %>%
      dplyr::sample_n(1)
  }
  
  ydat <- unlist( sampledData[ , "fur_adjusted"] )
  xdat <- unlist( sampledData[ , xCol] )
  
  mymodel <- smatr::sma(ydat~xdat)
  
  modelData <- cbind(mymodel$data, resid(mymodel))
  modelData$predicted_y_from_x <- ( mymodel$coef[[1]][2,1] * modelData[ , "xdat"] ) + mymodel$coef[[1]][1,1]
  modelData$predicted_x_from_y <- ( modelData[ , "ydat"] - mymodel$coef[[1]][1,1] ) / mymodel$coef[[1]][2,1]
  
  saveRDS(modelData, file = file.path(paste0( tempfile(tmpdir = savedirectory), ".rds")) )
  
  return(
    cbind(
      coef(mymodel)[1],               # intercept
      coef(mymodel)[2],               # slope
      sd(residuals(mymodel)),         # sdRes
      mymodel$groupsummary$n,         # n
      mymodel$groupsummary$pval,      # pval
      mymodel$groupsummary$r2         # r2
    )
  )
}

runBootstrappedSMA <- function(mydataframe, nClusters, R, saveDir = file.path(wd$bin, "tmp_boot"), from, to, min_n = 10){
  
  if(!dir.exists(saveDir)) dir.create(saveDir)

  lapply( grep(pattern = "precip_val", names(mydataframe), value = T ),
    function(isoColName){
      
      myBoot <- boot::boot(
        statistic = regWithResampling, R = R, sim = "parametric",
        data = mydataframe,
        colsToResampleBy = c("decimalLatitude", "decimalLongitude"),
        xCol = isoColName,
        savedirectory = file.path(saveDir, paste(isoColName, sep = "_")),
        min_n = min_n,
        parallel = "multicore", ncpus = nClusters
      )
      
      # Assemble and write variance outputs.
      file.path(saveDir, paste(isoColName, sep = "_")) %>%
        list.files(pattern = "rds", full.names = T) %>%
        lapply(readRDS) %>%
        do.call(rbind, .) %>%
        data.frame(., from = from, to = to, isoscape = isoColName) %>%
        saveRDS(., file = file.path(saveDir, paste(isoColName, from, to, "boot_var_output.rds", sep = "_")) )
      
      unlink( file.path(saveDir, paste(isoColName, sep = "_") ), recursive = T)
      
      # Get and tidy parameter estimate outputs.
      myBoot$t %>%
        as.data.frame() %>%
        {names(.) <- c("intercept", "slope", "sdRes", "n", "pval", "r2"); .} %>%
        lapply(function(x) data.frame(
          mean = mean(x),
          sd = sd(x),
          min = min(x),
          max = max(x)
        )) %>%
        plyr::ldply() %>%
        dplyr::rename(param = '.id') %>%
        data.frame(isoscape = gsub("precip_val_", "", isoColName), ., from = from, to = to)
    }) %>% plyr::ldply()
}


# Fit regression models. -------------------------------------------------------

myR <- 5000
nClusters <- 4

# directory to add coefficient data:
saveDir <- file.path(wd$bin, "tmp_boot")

system.time({
  sma_results_litDates_raw <- runBootstrappedSMA(
    mydataframe = mydata_isoVals_knownOrigin, R = myR,
    nClusters = nClusters, from = "lit", to = "lit",
    saveDir = saveDir
  )
})

sma_results <- sma_results_litDates_raw %>%
  tidyr::pivot_wider(
    names_from = param,
    values_from = c("mean", "sd", "min", "max")
  ) %>%
  dplyr::mutate_if(is.factor, as.character)

saveRDS(sma_results, file = file.path(wd$bin, "sma_results.rds"))


# Select transfer functions to apply. ----------------------------

sma_selected <- sma_results %>%
  dplyr::mutate_if(is.factor, as.character) %>%
  dplyr::arrange(desc(mean_r2)) %>%
  slice(1)


# Calculate residuals -----------------------------------------------------

# Make temporary model.
mod0 <- smatr::sma(precip_val_66098~fur_adjusted,data = mydata_isoVals_knownOrigin)
# Replace coefficients with those estimated.
mod0$coef[[1]][ ,1] <- c(sma_selected$mean_intercept, sma_selected$mean_slope)
myresiduals <- residuals(mod0)
df_resids <- data.frame(ID = mydata_isoVals_knownOrigin$ID, known_origin_resid = myresiduals)

# Apply transfer functions, define molt status  ---------------------------------

mydata_transformed <- mydata %>%
  dplyr::mutate_if(is.factor, as.character) %>%
  dplyr::mutate(
    dDprecip =  (fur_adjusted - sma_selected$mean_intercept) / sma_selected$mean_slope,
    sdResid = sma_selected$mean_sdRes
  ) %>% 
  full_join(df_resids, by = "ID")

saveRDS(mydata_transformed, file = file.path(wd$bin, "mydata_transformed.rds"))


# Pull out selected isoscape ----------------------------------------------

bestFitIso <- my_isoscapes[[sma_selected$isoscape]]
save(bestFitIso, file = file.path(wd$bin, "bestFitIso.rdata"))

