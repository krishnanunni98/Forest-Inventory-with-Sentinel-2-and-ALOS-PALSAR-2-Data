# Forest Inventory with Sentinel-2 and ALOS PALSAR-2 Data

## Project Overview

This project combines Sentinel-2 optical imagery and ALOS PALSAR-2 radar data for forest inventory in the Juupajoki/Orivesi area. The workflow includes visual interpretation of satellite images, Random Forest-based timber volume estimation, and species-specific volume prediction using KNN in R. The lab is explicitly designed to connect multisensor remote sensing with practical forest inventory modelling.

## Objectives

- Compare Sentinel-2 and ALOS PALSAR-2 imagery with a base map.
- Evaluate HH and HV radar bands and the PALSAR-2 land cover product.
- Estimate timber volume using Random Forest models.
- Compare models built with Sentinel-2 only, PALSAR-2 only, and both datasets combined.
- Generate wall-to-wall volume prediction maps.
- Estimate species-specific volumes using the KNN method.

## Data

- S2_Juupajoki.tif — Sentinel-2 optical image resampled to 10 m resolution.
- PALSAR-2_HH_Juupajoki.tif — HH polarization radar band.
- PALSAR-2_HV_Juupajoki.tif — HV polarization radar band.
- PALSAR-2_LC_Juupajoki.tif — radar-derived land cover classification.
- Juupajoki_field_data.csv — field plots from the Juupajoki sub-area.
- Orivesi_field+RS_data.csv — additional field plots and precomputed remote-sensing variables from the rest of the inventory area.
- Juupajoki_NFI_Volume.tif — official NFI-based volume map used as reference and forest mask.
  
## Methodology

### 1. Visual Interpretation
- Compared Sentinel-2 band composites against the base map.
- Inspected HH and HV radar bands for forest structure and land-cover differences.
- Evaluated the PALSAR-2 land cover classification for forest, non-forest, and water classes.

### 2. Random Forest Timber Volume Modelling
- Extracted Sentinel-2 and PALSAR-2 predictors for field plots.
- Built three Random Forest models:
- s2p2 using both Sentinel-2 and PALSAR-2 predictors,
- s2 using Sentinel-2 only,
- p2 using PALSAR-2 only.
- Evaluated the models using variable importance plots, scatter plots, RMSE, and relative RMSE.
- Used the fitted s2p2 model to produce a wall-to-wall volume map for the study area.

### 3. Species-Specific Volume Prediction
- Applied KNN-based modelling with Sentinel-2 predictors.
- Estimated pine, spruce, and deciduous volume components.
- Used yaImpute with the MSN distance metric and k = 5 as described in the lab instructions.

## Main Outputs

- Visual comparison maps for Sentinel-2, PALSAR-2 HH/HV, and land cover.
- Variable-importance plots for the three Random Forest models.
- Predicted-vs-observed scatter plots.
- RMSE and relative RMSE values for each model.
- Wall-to-wall RF volume map.
- Mean volume comparison against the NFI reference map.
- Species-specific volume estimation outputs from the KNN workflow.

## Key Results

- The combined s2p2 model performed best, with RMSE 85.68 and relative RMSE 0.4596. The Sentinel-2-only model performed very similarly, while the PALSAR-2-only model was weaker, with RMSE 130.45 and relative RMSE 0.6998. - The report also notes that the mean volume predicted by the s2p2 model was 152.1918 m³/ha, compared with 145.6553 m³/ha from the NFI reference map.

## Skills Demonstrated

- Optical and radar image interpretation
- Multisensor forest inventory
- Random Forest modelling in R
- KNN-based volume prediction
- Predictor extraction from geospatial raster data
- Model validation using RMSE and scatter plots
- Wall-to-wall raster prediction
- Comparison against reference inventory products
