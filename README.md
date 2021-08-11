# PESU_migration

Code and data used to generate results and figures of "Subtropical tri-colored bats (Perimyotis subflavus) migrate poleward in autumn" manuscript. If you clone this repository, results and figures are largely replicable-- note that the data we have included contains coordinates to caves that, as sensitive information, have been obscured post-analysis.

## This repository contains:
### /R/
Scripts used to read, tidy, and standardize data and run analyses.
| Script                                 | Description |
| -----------                            | ----------- |
| 00_Setup.R                             | Set up workspace, define extents and CRS |
| 01_loadIsotopeData.R                   | Import data, tidy, transform to reflect uniform standards |
| 02_LoadSpatialData.R                   | Load and standardize administrative area, isoscape, and range map datasets. |
| 03_FitApplyTransferFunctions.R         | Compare candidate isoscapes by finding isoscape predicted precipitation values at sampling sites for reference individuals, fit regression functions, estimate model parameters, select top model |
| 04_MakeProbabilityOfOriginMaps.R       | Create probability-of-origin maps for each individual |
| 05_ComparePoOMaps.R                    | Compare maps to one another for similarity |
| 06_MetricsConnectingPoOToSampleSite.R  | Estimate distance and directions traveled by individuals|
| 08_CombineOutputs.R                    | Combine results into a tidy dataframe |
| 09_modelFitting.R                      | Linear modeling to test predictive power of variables on distance traveled by individuals |

### /results/
#### ms_figs.R
Code used to generate figures used in manuscript.


### /data/
#### dD_analyzed_obscuredCoords.csv 
Contains individual bat collection data and &#948;D results. This includes samples collected as part of this study (analysisLab == "CASIF") and for Fraser et al. 2021. (analysisLab == "UWO")<sup>1</sup> [https://doi.org/10.1371/journal.pone.0031419]


| Column name      | Description |
| ----------- | ----------- |
| ID      | Individual ID Code       |
| fur   | &#948;D analysis results (‰, scale dependant on analysisLab)        |
| Date | Date individual sampled (%d/%m/%y for analysisLab == "CASIF", %yday-%y for analysisLab == "UWO" |
| Location | stateProvince of sampling location |
| Size | Tri-colored bat population of hibernacula sampling site, if applicable. Binned into "large" and "small". |
| decimalLatitude | Latitude of sample site, in decimal degrees (WGS84). Note that coordinates for Florda-sampled overwintering bats were obscured prior to release of these data. |
| decimalLongitude | Longitude of sample site, in decimal degrees (WGS84). Note that coordinates for Florda-sampled overwintering bats were obscured prior to release of these data. |
| Sex | Sex of bat |
| Season | Sampling season, binned to "Summer" and "Winter" (See <sup>2</sup> for information on molt periods.) |
| Region | Region within Florida containing the hibernacula where wintering bats were sampled |
| analysisLab | Stable isotope analysis lab where samples were analyzed for &#948;D values. Important because different standards were used at each lab -- see text of this ms and also <sup>3</sup>. |
| coordinatePrecision | Decimal representation of the precision of coordinates from decimalLatitude and decimalLongitude<sup>4</sup>. |

#### /iucn/
Contains IUCN redlist range map information, downloaded from https://www.iucnredlist.org/. Citation: Solari, S. 2018. IUCN Red List of Threatened Species: Perimyotis subflavus. https://dx.doi.org/10.2305/IUCN.UK.2018-2.RLTS.T17366A22123514.en

#### Isoscape data
- 66098_caitjcampbell_JJA_NoAm_Map_1980_2009/
- 66100_caitjcampbell_Annual_NoAm_Map_H_1980_2010/
- 70047_lisa.smith_annual_ele_lat_lat2/
- 70052_lisa.smith_annual_Lat_elev_map/ 
- 70055_lisa.smith_molt_lat_elev_map/
- 70064_lisa.smith_molt/
- 73715_dnelson99_May_aug_ele_lat_lat2/

Candidate isoscapes generated from [IsoMAP](https://isomap.rcac.purdue.edu/isomap/)

## Session Info
All analyses were conducted with the following configuration:

```
R version 4.0.2 (2020-06-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Mojave 10.14.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] viridisLite_0.4.0  patchwork_1.1.1    ggpubr_0.4.0       ggstatsplot_0.8.0 
 [5] kableExtra_1.3.4   gridExtra_2.3      ggsn_0.5.0         forcats_0.5.1     
 [9] readr_1.4.0        tibble_3.1.2       tidyverse_1.3.1    tidyr_1.1.3       
[13] stringr_1.4.0      smatr_3.4-8        rgdal_1.5-23       readxl_1.3.1      
[17] purrr_0.3.4        measurements_1.4.0 lwgeom_0.2-6       lubridate_1.7.10  
[21] ggmap_3.0.0        ggplot2_3.3.5      geosphere_1.5-10   chron_2.3-56      
[25] assignR_2.1.1      sf_1.0-0           dplyr_1.0.7        isocat_0.2.6      
[29] raster_3.4-13      sp_1.4-5          

loaded via a namespace (and not attached):
  [1] utf8_1.2.1                tidyselect_1.1.1          gmp_0.6-2                
  [4] maptools_1.1-1            munsell_0.5.0             effectsize_0.4.5         
  [7] codetools_0.2-18          units_0.7-2               miniUI_0.1.1.1           
 [10] withr_2.4.2               colorspace_2.0-2          knitr_1.33               
 [13] rstudioapi_0.13           ipmisc_6.0.2              ggsignif_0.6.2           
 [16] pbmcapply_1.5.0           emmeans_1.6.2-1           RgoogleMaps_1.4.5.3      
 [19] coda_0.19-4               vctrs_0.3.8               generics_0.1.0           
 [22] TH.data_1.0-10            xfun_0.24                 BWStest_0.2.2            
 [25] R6_2.5.0                  BayesFactor_0.9.12-4.2    bitops_1.0-7             
 [28] cachem_1.0.5              reshape_0.8.8             assertthat_0.2.1         
 [31] promises_1.2.0.1          scales_1.1.1              multcomp_1.4-17          
 [34] ggExtra_0.9               gtable_0.3.0              multcompView_0.1-8       
 [37] sandwich_3.0-1            rlang_0.4.11              MatrixModels_0.5-0       
 [40] zeallot_0.1.0             systemfonts_1.0.2         PMCMRplus_1.9.0          
 [43] splines_4.0.2             rstatix_0.7.0             broom_0.7.8              
 [46] abind_1.4-5               modelr_0.1.8              backports_1.2.1          
 [49] httpuv_1.6.1              rmapshaper_0.4.5          tools_4.0.2              
 [52] ellipsis_0.3.2            proxy_0.4-26              WRS2_1.1-2               
 [55] jsonvalidate_1.1.0        Rcpp_1.0.7                plyr_1.8.6               
 [58] rnaturalearth_0.1.0       classInt_0.4-3            viridis_0.6.1            
 [61] pbapply_1.4-3             correlation_0.6.1         zoo_1.8-9                
 [64] haven_2.4.1               ggrepel_0.9.1             fs_1.5.0                 
 [67] crul_1.1.0                magrittr_2.0.1            data.table_1.14.0        
 [70] openxlsx_4.2.4            reprex_2.0.0              mvnfast_0.2.7            
 [73] mvtnorm_1.1-2             hms_1.1.0                 mime_0.11                
 [76] evaluate_0.14             xtable_1.8-4              rio_0.5.27               
 [79] jpeg_0.1-8.1              pairwiseComparisons_3.1.6 compiler_4.0.2           
 [82] V8_3.4.2                  KernSmooth_2.23-20        crayon_1.4.1             
 [85] htmltools_0.5.1.1         mc2d_0.1-21               later_1.2.0              
 [88] DBI_1.1.1                 SuppDists_1.1-9.5         kSamples_1.2-9           
 [91] dbplyr_2.1.1              MASS_7.3-54               boot_1.3-28              
 [94] Matrix_1.3-4              car_3.0-11                cli_3.0.0                
 [97] parallel_4.0.2            insight_0.14.2            pkgconfig_2.0.3          
[100] statsExpressions_1.1.0    foreign_0.8-81            xml2_1.3.2               
[103] paletteer_1.3.0           foreach_1.5.1             svglite_2.0.0            
[106] geojsonlint_0.4.0         webshot_0.5.2             estimability_1.3         
[109] rvest_1.0.0               digest_0.6.27             parameters_0.14.0        
[112] httpcode_0.3.0            rmarkdown_2.9             cellranger_1.1.0         
[115] curl_4.3.2                shiny_1.6.0               gtools_3.9.2             
[118] rjson_0.2.20              lifecycle_1.0.0           jsonlite_1.7.2           
[121] carData_3.0-4             fansi_0.5.0               pillar_1.6.1             
[124] lattice_0.20-44           fastmap_1.1.0             httr_1.4.2               
[127] survival_3.2-11           glue_1.4.2                bayestestR_0.10.0        
[130] zip_2.2.0                 png_0.1-7                 iterators_1.0.13         
[133] class_7.3-19              stringi_1.6.2             performance_0.7.2        
[136] rematch2_2.1.2            memoise_2.0.0             Rmpfr_0.8-4              
[139] e1071_1.7-7   
```



## Cited

<sup>1</sup>. Fraser, E. E., McGuire, L. P., Eger, J. L., Longstaffe, F. J., & Fenton, M. B. (2012). Evidence of latitudinal migration in tri-colored bats, Perimyotis subflavus. PLoS One, 7(2), e31419. https://doi.org/10.1371/journal.pone.0031419

<sup>2</sup>. Fraser, E. E., Longstaffe, F. J., & Fenton, M. B. (2013). Moulting matters: the importance of understanding moulting cycles in bats when using fur for endogenous marker analysis. Canadian Journal of Zoology, 91(8), 533-544. https://doi.org/10.1139/cjz-2013-0072

<sup>3</sup>. Magozzi, Sarah, Clement P. Bataille, Keith A. Hobson, Michael B. Wunder, John D. Howa, Andrea Contina, Hannah B. Vander Zanden, and Gabriel J. Bowen. "Calibration chain transformation improves the comparability of organic hydrogen and oxygen stable isotope data." Methods in Ecology and Evolution 12, no. 4 (2021): 732-747. https://doi.org/10.1111/2041-210X.13556

<sup>4</sup> [DarwinCore#coordinatePrecision](https://dwc.tdwg.org/terms/#dwc:coordinatePrecision)
