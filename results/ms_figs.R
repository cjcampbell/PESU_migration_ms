# Setup -------------------------------------------------------------------

plotEM <- TRUE

library(tidyverse)
library(raster)
library(ggsn)
library(gridExtra)
library(kableExtra)
library(lubridate)
library(ggstatsplot)
library(ggpubr)
library(patchwork)
library(viridisLite)

# Make an object to help navigate the subdirectories.
wd <- list()
wd$bin  <- file.path( "/Users/cjcampbell/PESU_migration_project_directory/PESU_migration_ms", "bin" )
wd$data <- file.path( "/Users/cjcampbell/PESU_migration_project_directory/PESU_migration_ms", "data" )
wd$R    <- file.path( "/Users/cjcampbell/PESU_migration_project_directory/PESU_migration_ms", "R" )
wd$figs <- file.path( "/Users/cjcampbell/PESU_migration_project_directory/PESU_migration_ms", "figs" )

# Set plot themes.
theme_set(theme_light())

# Define extent of spatial analysis.
# This helps keep multiple spatial layers playing nicely with one another.
my_extent <- raster::extent(-178, -30, 10, 80)


if(!exists("maps_df")) maps_df <- readRDS(file.path(wd$bin, "maps_df.rds"))
source(file.path(wd$R, "00_Setup.R"))
source(file.path(wd$R, "01_loadIsotopeData.R"))

## Load results ##

# Load dD analysis results
myResults <- read.csv( file.path(wd$bin, "myResults.csv") )

# Load sma results
mydata_isoVals <- readRDS( file.path(wd$bin, "mydata_isoVals.Rds") )
load(file.path(wd$bin, "bestFitIso.rdata"))
sma_results <- readRDS( file.path(wd$bin, "sma_results.rds") )

# Load odds ratio sim surfaces
if(!exists("distDir_ORSim")) {
  distDir_ORSim <- readRDS(file.path(wd$bin, "distDir_ORSim.rds")) %>% ungroup }

# Transfer function plots ---------------------------------------------------

plotTransferFunctions <- function(myColors, mean.color = "darkred") {
  
  ## Dists of fur vals#
  p_dfur0 <- myResults %>% 
    ggbetweenstats(
      data = ,
      y = fur_adjusted,
      x = seasonLab,
      type = "r",
      centrality.point.args = list(size = 5, color = mean.color),
      ggplot.component	= list(
        ggplot2::scale_color_manual(values = myColors),
        ggplot2::scale_y_continuous(
          name = expression(paste(delta^2 ~ H[fur], " (", "\u2030", ", VSMOW)") ),
          limits = c(-65,35)
        ),
        ggplot2::xlab("Season Sampled / Analysis Laboratory")
      ),
      pairwise.display = "all", results.subtitle = FALSE
    )
  
  
  p_dfur00 <- ggplot_build(p_dfur0)
  # Remove caption.
  p_dfur00$plot$labels$caption <- NULL
  # Gotta get hackey to change the shapes...
  p_dfur00$data[[1]][ , "shape"] <- if_else(p_dfur00$data[[1]][ , "group"] == 1, 15, 16)
  p_dfur <- ggplot_gtable(p_dfur00) 
  
  
  ## regression Plots ##
  sma_selected <- sma_results %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::arrange(desc(mean_r2)) %>%
    slice(1)
  
  mySubtitle <- bquote( paste(
    atop(
      widehat(delta^2~H[fur])~'='~beta[0]+beta[1]~delta^2~H[precip]~
        ", "~
        
        bar(beta[0])~'*'~'='~.(round(sma_selected$mean_intercept,2))~
        "\u00B1"~.(round(sma_selected$sd_intercept,2))~
        ", "~
        bar(beta[1])~'*'~'='~.(round(sma_selected$mean_slope,2))~
        "\u00B1"~.(round(sma_selected$sd_slope,2))~
        ", "
      ,
      bar(R^2)~'*'~'='~.(round(sma_selected$mean_r2,2))~
        "\u00B1"~.(round(sma_selected$sd_r2,2))~
        ", "~
        bar(n)~'*='~.(sma_selected$mean_n)~
        "\u00B1"~.(round(sma_selected$sd_n,2))~
        ", B=5000"
    ))
  )
  
  p <- mydata_isoVals %>% 
    dplyr::filter(Season == "Summer") %>% 
    ggplot() +
    geom_point(
      aes(
        x = get(paste0("precip_val_", bestFitIso$path_pattern)), 
        y = fur_adjusted,
        color = analysisLab, shape = analysisLab
      ),
      size = 3, stroke = 0, alpha = 0.5
    )  +
    geom_abline(
      intercept = sma_selected$mean_intercept, 
      slope = sma_selected$mean_slope
    ) +
    scale_x_continuous(
      name = expression(paste(delta^2~H[precip]("\u2030", VSMOW))),
      limits = c(-60,-10)
    ) +
    scale_y_continuous(
      name = expression(paste(delta^2~H[fur]("\u2030", VSMOW))),
      limits = c(-65, 10)
    ) +
    scale_color_manual( name = "Analysis\nLaboratory", values = myColors , breaks = c("CASIF", "UWO")) +
    scale_shape_manual( name = "Analysis\nLaboratory" , values =c(15, 16) , breaks = c("CASIF", "UWO")) +
    theme(
      legend.position = c(0.16,0.76),
      plot.subtitle = element_text(size = 10)
    )
  
  p_d <- ggExtra::ggMarginal(p, type = "density", groupFill = TRUE)
  
  p_reg <- arrangeGrob(
    grobs = list(
      p_d, 
      textGrob(
        label = mySubtitle, hjust=0.5, gp = gpar(fontsize = 8) )
    ), 
    heights = c(5,1)
  )
  
  # Source data plot
  p1 <- arrangeGrob(grobs = list(
    grid.arrange(p_dfur , top = grid::textGrob(label = 'A', x = 0, hjust = 0, gp = gpar(fontsize = 16))) , 
    grid.arrange(p_reg  , top = grid::textGrob(label = 'B', x = 0, hjust = 0, gp = gpar(fontsize = 16)))
    ), ncol = 2)
  return(p1)
}

# Color version:
p1_col <- plotTransferFunctions(myColors = c("#440154FF", "#21908CFF", "#FDE725FF"))
ggsave(p1_col, filename = file.path(wd$figs, "sample_dDVals_plots_COLOR.png"),
       units = "in", height = 4, width = 8 )


# Black and white:
p1_bw <- plotTransferFunctions(myColors = grey.colors(5)[1:4], mean.color = "grey20")
ggsave(p1_bw, filename = file.path(wd$figs, "sample_dDVals_plots_bw.png"),
       units = "in", height = 4, width = 8 )

# Boxplots ----------------------------------------------------------------

load(file.path(wd$bin,"modelResults.Rdata"))
myThresh <- 0.25
myinterval <- -qnorm((1-0.95)/2)  # 95%

allModeldf_tidy <- allModeldf %>% 
  dplyr::mutate(
    Variable = case_when(
      Variable == "Hib. Size (Small)" ~ "Colony Size (Small)",
      Variable == "Hib. Region (N-C)" ~ "Karst Region (N-C)"
    )
    )

makeBoxPlots <- function(myColors, mean.color = "darkred", coeffColors = c("red", "grey50")){
  suppressMessages({suppressWarnings({
    set.seed(42)
    # Min dist by sex
    box_sex <- myResults %>% 
      filter(Season == "Winter") %>% 
      ggbetweenstats(
        y = !!(paste0("probDistTraveledAtThreshold_", myThresh)),
        x = Sex,
        type = "r",
        centrality.point.args = list(size = 5, color = mean.color),
        ggplot.component	= list(
          ggplot2::scale_color_manual(values = myColors),
          scale_y_continuous( name = "Minimum Distance Moved (km)") ,
          scale_x_discrete(position = "top") ,
          xlab("Sex"),
          ggpubr::theme_pubclean(),
          theme(legend.position = "none")
        ),
        pairwise.display = "all", results.subtitle = TRUE
      )
    subtit <- box_sex$labels$subtitle
    box_sex$labels$subtitle <- NULL
    box_sex <- grid.arrange(box_sex, bottom = textGrob(subtit, x = 0.55, y = 1, gp = gpar(fontsize = 11)))
    
    
    # Min dist by hibernaculum size
    box_size <- myResults %>% 
      filter(Season == "Winter") %>% 
      ggbetweenstats(
        y = !!(paste0("probDistTraveledAtThreshold_", myThresh)),
        x = Size,
        type = "r",
        centrality.point.args = list(size = 5, color = mean.color),
        ggplot.component	= list(
          ggplot2::scale_color_manual(values = myColors),
          scale_y_continuous( name = "Minimum Distance Moved (km)"),
          scale_x_discrete(position = "top") ,
          xlab("Colony Size"),
          ggpubr::theme_pubclean(),
          theme(legend.position = "none")
          ),
        pairwise.display = "all", results.subtitle = TRUE
      )
    subtit <- box_size$labels$subtitle
    box_size$labels$subtitle <- NULL
    box_size <- grid.arrange(box_size, bottom = textGrob(subtit, x = 0.55, y = 1, gp = gpar(fontsize = 11)))
    
    # Min dist by region
    box_region <- myResults %>% 
      filter(Season == "Winter") %>% 
      ggbetweenstats(
        y = !!(paste0("probDistTraveledAtThreshold_", myThresh)),
        x = Region_long,
        type = "r",
        centrality.point.args = list(size = 5, color = mean.color),
        ggplot.component	= list(
          ggplot2::scale_color_manual(values = myColors),
          xlab("Karst Region"),
          scale_y_continuous( name = "Minimum Distance Moved (km)"),
          scale_x_discrete(position = "top") ,
          ggpubr::theme_pubclean(),
          theme(legend.position = "none")
        ),
        pairwise.display = "all", results.subtitle = TRUE
      )
    subtit <- box_region$labels$subtitle
    box_region$labels$subtitle <- NULL
    box_region <- grid.arrange(box_region, bottom = textGrob(subtit, x = 0.55, y = 1, gp = gpar(fontsize = 11)))
    
    # Model Outputs
    coeff_plot <- ggplot(allModeldf_tidy, aes(colour = modelName, linetype = modelName)) + 
      geom_vline(xintercept = 0, colour = "grey20", lty = 5) +
      geom_pointrange(aes(y = Variable, x = Coefficient, xmin = Coefficient - SE*myinterval,
                          xmax = Coefficient + SE*myinterval),
                      lwd = 1/2, position = position_dodge(width = 1/2),
                      shape = 21, fill = "WHITE") +
      ggpubr::theme_pubr() +
      scale_color_manual(values = coeffColors) +
      scale_linetype_manual(values = c(1,1)) +
      theme(legend.position = c(0.15,0.1),  legend.box = "vertical", legend.title = element_blank())
    
    boxplots <- arrangeGrob(grobs = list(
      grid.arrange(box_region , top = grid::textGrob(label = 'A', x = 0, hjust = 0, gp = gpar(fontsize = 16))) , 
      grid.arrange(box_size   , top = grid::textGrob(label = 'B', x = 0, hjust = 0, gp = gpar(fontsize = 16))) , 
      grid.arrange(box_sex    , top = grid::textGrob(label = 'C', x = 0, hjust = 0, gp = gpar(fontsize = 16))) , 
      grid.arrange(coeff_plot , top = grid::textGrob(label = 'D', x = 0, hjust = 0, gp = gpar(fontsize = 16))) 
      ), ncol = 2)
    
  })})
  return(boxplots)
}

box_colors <- makeBoxPlots(myColors = viridis(3)[1:2], coeffColors = c("grey60", "maroon2"))
Sys.sleep(10)
ggsave(box_colors, filename = file.path(wd$figs, "minDistBoxplots_color.png"),
       units = "in", width = 12, height = 9)

box_bw <- makeBoxPlots(myColors = grey.colors(10)[c(2,5)], mean.color = "grey20", coeffColors = c("grey60", "black"))
Sys.sleep(10)
ggsave(box_bw, filename = file.path(wd$figs, "minDistBoxplots_bw.png"),
       units = "in", width = 12, height = 9)


# Compass rose and boxplots -----------------------------------------------

makeCompassRose <- function(myColors, add_n_labels = FALSE) {

  set.seed(420)
  
  # Shading
  shader <- data.frame(
    xmin = rep(-180, 3), 
    xmax = rep(180, 3), 
    ymin = c(log10(0), log10(100),log10(1000)),
    ymax = c(log10(100), log10(1000), Inf),
    dist = factor(c("Resident", "Regional", "Long-distance"), levels = c("Resident", "Regional", "Long-distance"))
    )
  
  # Plot.
  p_circlePoints <- myResults %>% filter(Season == "Winter") %>% 
    ggplot() +
    # Shading.
    geom_rect( data = shader,
               aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill = dist),
               alpha = 0.33
    ) +
    # Lines between shaded areas.
    geom_hline(yintercept = shader$ymin, colour = "grey50", size = 0.2) +
    # Lines along where panel.grid.major.x should go if it didn't suck so much and overlap with axis labels.
    geom_vline(xintercept = c(-180,0), color = "grey50", size = 0.2) +
    geom_vline(xintercept = c(-90,90), color = "grey50", size = 0.4) +
    # Tick marks.
    geom_segment( color = "grey50", x = -89, xend = -91, y = 1, yend = 1) +
    geom_text(size = 5, color = "grey50", x= -91, y = 1, label = "10 km", angle = 45, vjust = 0, hjust = 0) +
    
    geom_segment(color = "grey50", x = -89, xend = -91, y = 2, yend = 2) +
    geom_text(size = 5, color = "grey50", x= -91, y = 2, label = "100 km", angle = 45, vjust = 0, hjust = 0) +
    
    geom_segment(color = "grey50", x = -89, xend = -91, y = 3, yend = 3) +
    geom_text(size = 5, color = "grey50", x= -91, y = 3, label = "1000 km", angle = 45, vjust = 0, hjust = 0) +
    # Plot points
    geom_jitter( aes(
      x = mostLikelyAngle, 
      y = log10(probDistTraveledAtThreshold_0.25) 
      ), width = 8, height = 0, size = 3, stroke = 0, alpha = 0.5
    ) +
    # Make it look nice.
    scale_x_continuous(
      name = NULL,
      breaks = seq(-180,179,90),
      limits = c(-180,180), 
      expand = c(0,0),
      labels = c("S", "W", "N", "E")
    ) +
    scale_y_continuous(
      breaks = log10(c(0,100,1000)), 
      limits = c(0,log10(5000)),
      expand = c(0,0)
    ) +
    coord_polar(theta = "x", start = pi) +
    scale_fill_manual(
      values = myColors,
      name = "Min. Distance\nTraveled") +
    facet_wrap(~Region_long) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 12.5),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.ticks.x.top = element_line(),
      strip.background = element_blank(),
      strip.text = element_text(size = 15),
      legend.position = c(0.51, -.03), legend.justification = c(0.5, 0),
      legend.background = element_rect(fill = NA),
      plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")
    )
  
  if(add_n_labels == TRUE) {
    m2 <- myResults %>% filter(Season == "Winter") %>% 
      dplyr::group_by(Region, mostLikelyDiagonal, distanceTraveled, .drop = T) %>% 
      dplyr::summarise(n = n()) %>% 
      mutate(
        x = case_when(
          mostLikelyDiagonal == "NE" ~ 90,
          mostLikelyDiagonal == "SE" ~ 90,
          mostLikelyDiagonal == "NW" ~ -90,
          mostLikelyDiagonal == "SW" ~ -90,
        ),
        y = case_when(
          distanceTraveled == "Resident" ~ 0.8,
          distanceTraveled == "Regional"   ~ 2.1,
          distanceTraveled == "Long-distance"  ~ 3.1
        ),
        vjust = case_when(
          mostLikelyDiagonal == "NE" ~ -0.3,
          mostLikelyDiagonal == "SE" ~ 1.3,
          mostLikelyDiagonal == "NW" ~ -0.3,
          mostLikelyDiagonal == "SW" ~ 1.3
        )
      )
    
    p_circlePoints_n <- p_circlePoints + 
      geom_text(
        data = m2, 
        aes(label = paste0("n = ", n), x= x, y = y, vjust=vjust),
        angle = 0, hjust = 0,  color = "gray30"
      )
    return(p_circlePoints_n)
  } else {
    return(p_circlePoints)
  }
}

## Stacked bar plots

makeCompassBarPlots <- function(myColors = c("#f9de59", "#e8a628", "#f98365"), add_n_labels = FALSE) {
  p_circlePoints <- makeCompassRose(myColors, add_n_labels = add_n_labels)
  
  tt <- myResults %>% 
    filter(Season == "Winter") %>% 
    dplyr::mutate(
      NS = if_else(distanceTraveled == "Resident", "Indeterminate", as.character(NorthOrSouth)),
      NS = factor(NS, levels = c("Indeterminate", "North", "South"))) %>% 
    dplyr::group_by(Region_long, distDir, distanceTraveled, NS, .drop = F) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr::mutate(mylab = case_when(
      distanceTraveled == "Resident" ~ paste0("All Origins\nn = ", n),
      distanceTraveled != "Resident" ~ paste0(NS,"erly Origins\nn = ", n)
    ))
  
  tt_Northies <- filter(tt, n != 0, NS=="North") %>% 
    dplyr::mutate(y = case_when(
      distDir == "Regional-N" & Region_long == "Northwest" ~ 34,
      distDir == "Long-distance-N" & Region_long == "Northwest" ~ 1,
    ))
  
  p_bar <- tt %>% 
    ggplot() +
    aes(fill = distanceTraveled, y = n, x = distanceTraveled, label = paste0("n=", n), group = NS) +
    geom_bar(fill = "white", position="stack", stat = "identity", width = 1) +
    geom_bar(color = "black", position="stack", stat = "identity", width = 1, alpha = 0.33) +
    facet_wrap(~Region_long, ncol = 2) +
    geom_text(data = tt_Northies,
              aes(label = mylab, y=y), size = 3.5, hjust = 0.5, vjust = -0.5,
              lineheight = .8) +
    geom_segment(
      data = tt_Northies, 
      aes(y=y-(n/2), yend = y+1.8, xend=distanceTraveled),
      color = "grey50") +
    geom_text(
      data = filter(tt, n != 0, NS!="North") ,
      aes(label = mylab), stat = "identity", position = position_stack(0.5),
      size = 3.5, hjust = 0.5, vjust = 0,lineheight = .8) +
    scale_fill_manual(values = myColors) +
    scale_y_continuous(name = "Sample Size", limits = c(0,40)) +
    scale_x_discrete(name = "Minimum Distance Traveled") +
    theme_pubclean() +
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")
    ) +
    guides(fill = "none")
  
  p_circleBar <- ggpubr::ggarrange(p_circlePoints, p_bar, ncol = 1, common.legend = F, heights = c(0.7,0.5))
  return(p_circleBar)
}

p_circleBar_color <- makeCompassBarPlots(myColors = c("#f9de59", "#e8a628", "#f98365"), add_n_labels = FALSE)
ggsave(p_circleBar_color, file = file.path(wd$figs, "p_circleBar_color.png"),
       units = "in", width = 9, height = 7)

p_circleBar_bw <- makeCompassBarPlots(myColors = grey.colors(3, start = 0.9, end = 0.4), add_n_labels = FALSE)
ggsave(p_circleBar_bw, file = file.path(wd$figs, "p_circleBar_bw.png"),
       units = "in", width = 9, height = 8)

# Summary Origin Maps -----------------------------------------------------

plotSummaryOriginMaps <- function(lowCol = "grey90", highCol = "black", pointcol = "black", textAnnotationSize = 3.5) {
  if(!exists("na")){
    na <- rnaturalearth::ne_countries(country = c("United States of America", "Canada", "Mexico"), returnclass = "sf", scale = "large") %>% 
      st_transform(crs = myCRS)
  }
  
  colorbarName <-  expression(paste("Cumulative Origins of ", italic("n"), " Individuals") )
  
  ## North-Central ##
  ORSim_Northcentral <- distDir_ORSim %>% 
    dplyr::select(ID, x, y, ID, OR_sim, dist_km) %>% 
    right_join(myResults, by = "ID") %>% 
    mutate(over0.25 = if_else(OR_sim >= 0.25, 1, 0)) %>% 
    filter(Season == "Winter", Region == "NC")
  
  r <- ORSim_Northcentral %>%
    dplyr::filter(distDir %in% c("Regional-N", "Long-distance-N", "Regional-S")) %>% 
    dplyr::group_by(x,y) %>%
    dplyr::summarise(summed = sum(over0.25)) %>%
    ungroup %>% 
    as.data.frame() 
  r$summed[r$summed==0] <- NA
  
  plot_params <-
    list(
      scale_color_gradient(
        name = colorbarName, low = lowCol, high = highCol, na.value = NA,
        limits = c(0,40),
        guide = guide_colorbar(
          nbin = 300, barwidth=10, frame.colour=NA, frame.linewidth=1, 
          ticks.colour="black",  direction="horizontal", title.position="top",
          title.hjust = 0.5
        )
      ),
      scale_fill_gradient( 
        name = colorbarName , low = lowCol, high = highCol, na.value = NA,
        limits = c(0,40),
        guide = guide_colorbar(
          nbin = 300, barwidth=10, frame.colour=NA, frame.linewidth=1, 
          ticks.colour="black",  direction="horizontal", title.position="top",
          title.hjust = 0.5
        )
      ) ,
      theme(
        legend.position = "bottom",
        axis.title = element_blank(),
        plot.subtitle = element_text(hjust = 0.5)
      ),
      coord_sf(
        xlim = c(min(r$x), max(r$x)*1.1),
        ylim = c(min(r$y), max(r$y)),
        crs = myCRS
      ) 
    )
  
  u <- ggplot() +
    geom_sf(na, mapping = aes(), fill= "white") +
    geom_tile( data = r, mapping = aes(x=x,y=y, fill = summed), alpha = 0.8) +
    geom_sf(na, mapping = aes(), fill= NA) +
    plot_params +
    scale_color_gradient(
      name = "Number of\nIndividuals" , low = "grey90", high = "black", na.value = NA,
      limits = c(0,40),
      guide = guide_colorbar(
        nbin = 300, barwidth=10, frame.colour=NA, frame.linewidth=1, 
        ticks.colour="black",  direction="horizontal", title.position="top")
    ) 
  
  # Find closest positions.
  minCellLocations <- ORSim_Northcentral %>% 
    filter(over0.25 == 1) %>% 
    group_by(ID) %>% 
    arrange(dist_km) %>% slice(1)
  
  northies <- minCellLocations %>% 
    filter(NorthOrSouth == "North", distanceTraveled != "Resident")
  
  southies <- minCellLocations %>% 
    filter(NorthOrSouth == "South", distanceTraveled != "Resident") %>% 
    group_by(distanceTraveled) %>% 
    arrange(dist_km) %>% slice(1)
  
  samplingLocations_NC <- ORSim_Northcentral %>% 
    dplyr::select(metersLongitude, metersLatitude) %>% 
    distinct
  
  # Top annotate location
  northiesLabelLocation <- c(2e6-2e5, 0)
  southiesLabelLocation <- c(2e6-2e5, -1e6)
  
  map_NC <- u +
    geom_point(samplingLocations_NC, mapping = aes(x=metersLongitude, y = metersLatitude), color = pointcol) +
    labs(subtitle = "North-central Florida Hibernacula") +
    # Southies
    annotate(
      geom = "text", 
      x=southiesLabelLocation[1]+5e4, y = southiesLabelLocation[2], color = "black",
      hjust = 0, size = textAnnotationSize,
      label = "Sixteen bats moved\n>100km from southerly\nsummering grounds."
    ) +
    annotate(
      geom = "curve", x=southiesLabelLocation[1], y = southiesLabelLocation[2], 
      xend = southies$x[1], yend = southies$y[1],
      arrow = arrow(length = unit(0.03, "npc") ) )

  
  
  ### NorthWest!!! ###
  ORSim_Northwest <- distDir_ORSim %>% 
    dplyr::select(ID, x, y, ID, OR_sim, dist_km) %>% 
    right_join(myResults, by = "ID") %>% 
    mutate(over0.25 = if_else(OR_sim >= 0.25, 1, 0)) %>% 
    filter(Season == "Winter", Region == "NW")
  r <- ORSim_Northwest %>%
    dplyr::filter(distDir %in% c("Regional-N", "Long-distance-N", "Regional-S")) %>% 
    dplyr::group_by(x,y) %>%
    dplyr::summarise(summed = sum(over0.25)) %>%
    ungroup %>% 
    as.data.frame() 
  r$summed[r$summed==0] <- NA
  
  u <- ggplot() +
    geom_sf(na, mapping = aes(), fill= "white") +
    geom_tile( data = r, mapping = aes(x=x,y=y, fill = summed), alpha = 0.8) +
    geom_sf(na, mapping = aes(), fill= NA) +
    plot_params +
    scale_color_gradient( # I have no idea why I need to do this twice...
      name = "Number of\nIndividuals" , low = "grey90", high = "black", na.value = NA,
      limits = c(0,40),
      guide = guide_colorbar(
        nbin = 300, barwidth=10, frame.colour=NA, frame.linewidth=1, 
        ticks.colour="black",  direction="horizontal", title.position="top")
    ) 
  
  # Find closest positions.
  minCellLocations <- ORSim_Northwest %>% 
    filter(over0.25 == 1) %>% 
    group_by(ID) %>% 
    arrange(dist_km) %>% slice(1)
  
  northies <- minCellLocations %>% 
    filter(NorthOrSouth == "North", distanceTraveled != "Resident")
  
  southies <- minCellLocations %>% 
    filter(NorthOrSouth == "South", distanceTraveled != "Resident") %>% 
    group_by(distanceTraveled) %>% 
    arrange(dist_km) %>% slice(1)
  
  samplingLocations_NW <- ORSim_Northwest %>% 
    dplyr::select(metersLongitude, metersLatitude) %>% 
    distinct

  map_NW <- u  +
    geom_point(samplingLocations_NW, mapping = aes(x=metersLongitude, y = metersLatitude), color = pointcol) +
    labs(subtitle = "Northwest Florida hibernacula") +
    annotate(
      geom = "text", 
      x=northiesLabelLocation[1]+5e4, y = northiesLabelLocation[2], color = "black",
      hjust = 0, size = textAnnotationSize,
      label = "Three individuals moved\n>100km from northerly\nsummering grounds."
    ) +
    annotate(
      geom = "curve", x=northiesLabelLocation[1], y = northiesLabelLocation[2], 
      xend = northies$x[1], yend = northies$y[1],
      arrow = arrow(length = unit(0.03, "npc") ) ) +
    annotate(
      geom = "curve", x=northiesLabelLocation[1], y = northiesLabelLocation[2], 
      xend = northies$x[2], yend = northies$y[2],
      arrow = arrow(length = unit(0.03, "npc") ) ) +
    annotate(
      geom = "curve", x=northiesLabelLocation[1], y = northiesLabelLocation[2], 
      xend = northies$x[3], yend = northies$y[3],
      arrow = arrow(length = unit(0.03, "npc") ) ) +
    # Southies
    annotate(
      geom = "text", 
      x=southiesLabelLocation[1]+5e4, y = southiesLabelLocation[2], color = "black",
      hjust = 0, size = textAnnotationSize,
      label = "Thirty-two bats moved\n>100km from southerly\nsummering grounds."
    ) +
    annotate(
      geom = "curve", x=southiesLabelLocation[1], y = southiesLabelLocation[2], 
      xend = southies$x[1], yend = southies$y[1],
      arrow = arrow(length = unit(0.03, "npc") ) )
  #map_NW
  p_originMap <- ggpubr::ggarrange(map_NW, map_NC, ncol = 2, common.legend = T, legend = "bottom") %>%
    annotate_figure(top = "Summer Origins of Regional- and Long-distance Migrators")

}

p_originMap_col <- plotSummaryOriginMaps(lowCol = "#BFE1B0", highCol = "#0A2F51", pointcol = "#A30808")
ggsave(p_originMap_col, filename = file.path(wd$figs, "p_originMap_col.png"),
       units = "in", width = 12, height = 6)


p_originMap_bw <- plotSummaryOriginMaps(lowCol = "grey90", highCol = "black")
ggsave(p_originMap_bw, filename = file.path(wd$figs, "p_originMap_bw.png"),
       units = "in", width = 12, height = 6)


# Pull out stats for results section -----------------------------------

myResults %>% 
  filter(Season == "Summer") %>% 
  ggplot() +
  geom_sf(na, mapping = aes(), fill= "white") +
  geom_point(aes(x=metersLongitude, y = metersLatitude, color = analysisLab))

myResults %>% 
  filter(Season == "Summer") %>% 
  group_by(analysisLab, Location) %>% 
  dplyr::summarise(n=n(), min_dD = min(fur_adjusted))

myResults %>% filter(Season == "Summer") %>%  
  group_by(Location) %>% dplyr::summarise(n=n())
                   
myResults %>% 
  filter(Season == "Winter") %>% 
  group_by(cave, Location) %>% 
  dplyr::summarise(n=n(), min_dD = min(fur_adjusted))

myResults %>% filter(Season == "Winter") %>% 
  group_by(cave, Location) %>% dplyr::summarise(n=n()) %>% ungroup %>% 
  {summary(.$n)}

myResults %>% filter(Season == "Winter") %>% 
  {summary(.$fur_adjusted)}

myResults %>% 
  filter(Season == "Winter") %>% 
  group_by(Sex) %>% 
  dplyr::summarise(n=n(), min_dD = min(fur_adjusted))


# Distance traveled
myResults %>% group_by(Season) %>% dplyr::summarise(mean=mean(probDistTraveledAtThreshold_0.25), sd = sd(probDistTraveledAtThreshold_0.25), median = median(probDistTraveledAtThreshold_0.25))

myResults %>% filter(Season == "Winter") %>% group_by(distanceTraveled) %>% dplyr::summarise(n=n())

myResults %>% filter(Season == "Winter", NorthOrSouth == "North", distanceTraveled != "Resident") %>% dplyr::summarise(probDistTraveledAtThreshold_0.25 = probDistTraveledAtThreshold_0.25)
