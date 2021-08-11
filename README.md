# PESU_migration

Code and data used to generate results and figures of "Subtropical tri-colored bats (Perimyotis subflavus) migrate poleward in autumn" manuscript. If you clone this repository, results and figures are largely replicable-- note that the data we have included contains coordinates to caves that, as sensitive information, have been obscured post-analysis.

# This repository contains:
## /data/
### dD_analyzed_obscuredCoords.csv 
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

### /iucn/
Contains IUCN redlist range map information, downloaded from https://www.iucnredlist.org/. Citation: Solari, S. 2018. IUCN Red List of Threatened Species: Perimyotis subflavus. https://dx.doi.org/10.2305/IUCN.UK.2018-2.RLTS.T17366A22123514.en

### Isoscape data
- 66098_caitjcampbell_JJA_NoAm_Map_1980_2009/
- 66100_caitjcampbell_Annual_NoAm_Map_H_1980_2010/
- 70047_lisa.smith_annual_ele_lat_lat2/
- 70052_lisa.smith_annual_Lat_elev_map/ 
- 70055_lisa.smith_molt_lat_elev_map/
- 70064_lisa.smith_molt/
- 73715_dnelson99_May_aug_ele_lat_lat2/

Candidate isoscapes generated from [IsoMAP](https://isomap.rcac.purdue.edu/isomap/)




# Cited

<sup>1</sup>. Fraser, E. E., McGuire, L. P., Eger, J. L., Longstaffe, F. J., & Fenton, M. B. (2012). Evidence of latitudinal migration in tri-colored bats, Perimyotis subflavus. PLoS One, 7(2), e31419. https://doi.org/10.1371/journal.pone.0031419

<sup>2</sup>. Fraser, E. E., Longstaffe, F. J., & Fenton, M. B. (2013). Moulting matters: the importance of understanding moulting cycles in bats when using fur for endogenous marker analysis. Canadian Journal of Zoology, 91(8), 533-544. https://doi.org/10.1139/cjz-2013-0072

<sup>3</sup>. Magozzi, Sarah, Clement P. Bataille, Keith A. Hobson, Michael B. Wunder, John D. Howa, Andrea Contina, Hannah B. Vander Zanden, and Gabriel J. Bowen. "Calibration chain transformation improves the comparability of organic hydrogen and oxygen stable isotope data." Methods in Ecology and Evolution 12, no. 4 (2021): 732-747. https://doi.org/10.1111/2041-210X.13556

<sup>4</sup> [DarwinCore#coordinatePrecition](https://dwc.tdwg.org/terms/#dwc:coordinatePrecision)
