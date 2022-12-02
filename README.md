# Global effects of land-use intensity on local pollinator biodiversity

[![DOI](https://zenodo.org/badge/346324215.svg)](https://zenodo.org/badge/latestdoi/346324215) [![Generic badge](https://img.shields.io/badge/License-MIT-green.svg)](https://shields.io/)

This repository contains all the scripts used for the analysis carried out in the below study:

> **Millard _et al_. (2021), Global effects of land-use intensity on local pollinator biodiveristy. _Nature Communications_ DOI: https://doi.org/10.1038/s41467-021-23228-3

The scripts for the main figures are the following: 
* Figure 1, '01a_PREDICTS_compilation_figures'
* Figure 2, '02_pollinator_analysis.R'
* Figure 3, '02c_cropland_geographic_region_response.R'
* Figure 4, '02a_taxonomic-breakdown_analysis_anthropogenic.R'
* Figure 5, '03a_rate_earth-stat_intensity-models.R'

To reproduce the analysis and figures in these scripts download the pollinator PREDICTS subset from either Figshare here (https://figshare.com/articles/dataset/Global_effects_of_land-use_intensity_on_local_pollinator_biodiversity/12815669/2), or from the 'outputs' folder of this GitHub repo. 

For Figure 5 download the fertiliser application rate data from Earthstat here (http://www.earthstat.org/nutrient-application-major-crops/). For Figure 5 you'll also need to use the script '01b_fertiliser_aggregation.R' to aggregate total fertiliser application rate for all crops. The fertiliser raster output from this script ('fertiliser_application_rate_1.tif') is then read in at the top of '03a_rate_earth-stat_intensity-models.R'.
