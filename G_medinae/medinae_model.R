################################################################################
# Ecological Niche Modeling (ENM) Workflow for Genieridium medinae
#
# Developed by Diego E. Martínez-Revelo
# date: 2026-03-14
#
# Description:
# This R script implements a complete workflow for species distribution modeling
# and ecological niche estimation for Genieridium medinae using Maxent (via 
# maxnet/ENMeval). The workflow includes:
#   1. Loading and preprocessing WorldClim bioclimatic variables for Colombia
#   2. Defining the calibration area (M) based on occurrence records
#   3. Extracting environmental predictors and assessing multicollinearity
#   4. Calibrating, tuning, and evaluating Maxent models
#   5. Selecting the optimal model based on AICc
#   6. Calculating variable importance (coefficients & permutation)
#   7. Projecting climatic suitability within calibration and expanded areas
#   8. Estimating the potential climatic niche (fundamental niche)
#   9. Integrating land cover (Hansen forest cover) to estimate realized habitat
#  10. Calculating habitat area, potential loss, and plotting results
#
# Notes:
# - The script uses terra and sf packages for raster/vector operations.
# - Occurrence records are read from Excel; background points are sampled 
#   within the accessible area (M).
# - Thresholding for potential niche is based on the 10th percentile of 
#   suitability at occurrences.
# - Forested habitat is defined as cells with >=70% cover ("adequate forest cover").
#
# Intended use: reproducible SDM and niche analyses, suitable for publication,
# conservation assessment, or further spatial analysis.
#
################################################################################
setwd("F:/Doctorado/encuentro_ciencias")

library(readxl)
library(dplyr)
library(terra)
library(sf)
library(rnaturalearth)
library(geodata)
library(ggplot2)
library(viridis)
# ------------------------------------------------------------------------------
# WorldClim 2.1 bioclimatic variables (historical), 30 arc-second resolution
#-------------------------------------------------------------------------------
#
# Load administrative boundaries of Colombia and DEM
  col <- ne_states(country = "colombia", returnclass = "sf")
  dem_col <- elevation_30s(country = "COL", path = tempdir())
#
# Path to the directory containing WorldClim BIO variables (.tif files)
  path_bio <- "F:/Capas/worldclim/wc2.1_30s_bio"
#
# List all GeoTIFF files in the directory
  files <- list.files(path_bio, pattern = "\\.tif$", full.names = TRUE)
#
# Extract BIO variable numbers from filenames (e.g., bio_1, bio_2, ..., bio_19)
  step1 <- gsub(".*_bio_", "", files)
  step2 <- as.numeric(gsub("\\.tif", "", step1))
#
# Reorder files numerically to ensure correct stacking order (BIO1 → BIO19)
  files <- files[order(step2)]
#
# Load all bioclimatic layers as a SpatRaster stack
  bioclim <- rast(files)
#
# Print Coordinate Reference System (CRS) information for verification
  cat("Bioclimatic CRS:", crs(bioclim, proj = TRUE), "\n")
  cat("Colombia CRS:", st_crs(col)$input, "\n")
#
# Crop raster stack to Colombia's bounding box (fast operation)
  bioclim_col_crop <- crop(bioclim, col)
#
# Mask raster using Colombia's polygon
  bioclim_col <- mask(bioclim_col_crop, col)
#
# Optionally save the processed raster stack for later use
# saveRDS(bioclim_col, "bioclim_col.rds")
#
# Reload the preprocessed raster stack (avoids repeating heavy operations)
  bioclim_col <- readRDS("./poster_2/bioclim_col.rds")
#
# Plot BIO1 (Annual Mean Temperature) for Colombia
  plot(bioclim_col[[1]], main = "BIO1 - Annual Mean Temperature (Colombia)")
  plot(st_geometry(col), add = TRUE, border = "black", lwd = 1.2)
#
# ------------------------------------------------------------------------------
# Calibration area (M) based on occurrence records
  
# Read occurrence data from Excel file
  records <- read_excel("./poster_2/G_medinae.xlsx", sheet = "General")
  records <- records %>% select("Species", "Lat", "Long")
  records_sf <- st_as_sf(records, coords = c("Long", "Lat"), crs = 4326)  

# Albers Equal Area projection (meters) 
  AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  
# Reproject to Albers Equal Area
  records_aea <- st_transform(records_sf, crs = AEAstring)
  
# Create lines connecting points to represent potential dispersal corridors
  occurrence_line <- records_aea %>%
        st_union() %>%
        st_cast("LINESTRING")
  
# Create a 50 km buffer around the line, defines accessible area (M)
  m_buffer <- st_buffer(occurrence_line, dist = 50000)
  m_calibration <- st_transform(m_buffer, crs = 4326)
  
# Plot M and records
  plot(dem_col, axes = FALSE, box = TRUE, main = "Calibration area (M)")
  plot(m_calibration, col = adjustcolor("grey70", alpha.f = 0.5), border = "grey60", lwd = 2, add = TRUE, f.alpha = 0.5)
  plot(st_geometry(records_sf), pch = 21, bg = "red", add = TRUE)
#-------------------------------------------------------------------------------
# Bioclimatic variables within the calibration area (M)
  
# Crop bioclimatic layers to the bounding box of M
  bioclim_M_crop <- crop(bioclim_col, m_calibration)
  
# Convert sf polygon to terra vector format for masking
  m_calibration_vect <- vect(m_calibration)
  
# Mask layers using the exact shape of M (continuous area)
# Cells outside the accessible area are set to NA
  bioclim_M <- mask(bioclim_M_crop, m_calibration_vect)
  
# Optionally save final calibration stack
# saveRDS(bioclim_M, "./poster_2/bioclim_M.rds")
  
# Reload masked layers
  bioclim_M <- readRDS("./poster_2/bioclim_M.rds")
  
# Clean layer names to simple format: bio_1 ... bio_19
# WorldClim files often include long or inconsistent names
  names(bioclim_M) <- paste0("bio_", 1:19)
  
# Quick visual check of one variable (e.g., BIO12 - Annual Precipitation)
  plot(bioclim_M[[12]], main = "BIO12 within calibration area (M)")
  plot(st_geometry(col), add = TRUE, border = "black", lwd = 1.2)
  plot(st_geometry(records_sf), add = TRUE, pch = 21, bg = "red", co = "black")
#  
# ------------------------------------------------------------------------------
# Multicollinearity assessment and predictor selection
  
  library(usdm)
  library(corrplot)
  
 # Extract raster pixel values within calibration area (M)
  env_values_M <- as.data.frame(values(bioclim_M, na.rm = TRUE))
  
# Spearman correlation matrix
  cor_mat <- cor(env_values_M, method = "spearman")
  corrplot(cor_mat,
           method = "circle",
           type = "upper",
           tl.col = "black",
           tl.srt = 45,
           title = "Spearman Correlation - Bioclimatic Variables (M)",
           mar = c(0, 0, 1, 0))
  
# VIF analysis
  vif_res <- vif(env_values_M)
  print(vif_res)
  
  vif_select <- vifstep(env_values_M, th = 10)
  print(vif_select)
  
# Selected predictor names
  selected_vars <- vif_select@results$Variables
  print(selected_vars)
  
# Raster stack with selected predictors
  bioclim_selected <- bioclim_M[[selected_vars]]
  
# Extract environmental values at occurrence points
  env_at_occ <- terra::extract(bioclim_selected, records_sf)
  
# Dataset for modeling
  model_data <- cbind(st_coordinates(records_sf), env_at_occ[, -1])
#
# -----------------------------------------------------------------------------
#  Ecological niche model calibration using Maxent (maxnet) via ENMeval
  
  library(ENMeval)
  
# Set seed for reproducibility of background sampling
  set.seed(456)
  
# Generate random background points within calibration area (M)
# Sampling is constrained to cells with environmental data (non-NA)
  bg_points <- spatSample(bioclim_selected,
                          size = 10000,
                          method = "random",
                          na.rm = TRUE,
                          as.points = TRUE)
  
# Extract coordinates of occurrence records
  occ_coords <- st_coordinates(records_sf)
  
# Extract coordinates of background points
  bg_coords <- crds(bg_points)
  
# Rename columns to required format for ENMeval (longitude, latitude)
  colnames(occ_coords) <- c("long", "lat")
  colnames(bg_coords)  <- c("long", "lat")
  
  # Model tuning and evaluation
    enm_results <- ENMevaluate(
        # Occurrence coordinates (presence-only data)
        occs = occ_coords,
        # Environmental predictors (selected after VIF filtering)
        envs = bioclim_selected,
        # Background points representing available environment in M
        bg = bg_coords,
        # Maxent implementation in R (maxnet package)
        algorithm = "maxnet",
        # Spatial partitioning using geographic blocks
        # Recommended for spatially structured occurrence data
        partitions = "block",
        # Model tuning parameters
        tune.args = list(
        # Feature classes controlling response complexity
        fc = c("L", "LQ", "LH", "LQH"),
        # Regularization multiplier controlling model smoothness
        rm = seq(0.5, 2, by = 0.5)
        )
       )
# saveRDS(resultados_enmeval, "./poster_2/resultados_enmeval.rds")
#
# ------------------------------------------------------------------------------
# Model selection, variable importance, and habitat suitability prediction
      
# Load ENMeval results object
  enm_results <- readRDS("./poster_2/enm_results.rds")
      
# Model selection based on AICc
      
# Extract evaluation results table
  results_table <- eval.results(enm_results)
      
# Order models by delta AICc (lower = better)
  results_table <- results_table[order(results_table$delta.AICc), ]
      
# Select optimal model (index corresponds to best model in results)
  optimal_model <- enm_results@models[[11]]
      
  # ---------------------------------------------
  # Variable importance based on model coefficients (Maxnet)
  # ---------------------------------------------
  # This calculates the importance of each variable by summing the absolute values
  # of its coefficients in the Maxnet model (including hinge and linear features).
  # Interpretation:
  #   - Represents the "mathematical weight" that the model assigns to each predictor.
  #   - Useful to understand the internal structure of the model.
  #   - Note: Large coefficients do not always mean the variable strongly influences
  #     model predictions if the variable has low variation or interacts with others.
      
      # Extract beta coefficients
      beta_values <- optimal_model$betas
      
      # Clean coefficient names to recover original predictor names
      # (e.g., hinge(bio_2):... → bio_2)
      clean_names <- gsub(".*\\((.*)\\):.*", "\\1", names(beta_values))
      clean_names <- gsub("linear\\((.*)\\)", "\\1", clean_names)
      
      coef_df <- data.frame(
        Variable = clean_names,
        Coefficient = as.numeric(beta_values)
      )
      
      # Aggregate absolute coefficients by original variable
      # (a predictor may appear multiple times due to hinge features)
      var_importance <- aggregate(abs(Coefficient) ~ Variable,
                                  data = coef_df,
                                  sum)
      
      colnames(var_importance)[2] <- "Relative_Importance"
      
      # Sort variables from most to least important
      var_importance <- var_importance[order(-var_importance$Relative_Importance), ]
      
      print(var_importance)
      
      # --- Plot variable importance
      
      barplot(sort(var_importance$Relative_Importance, decreasing = TRUE),
              names.arg = var_importance$Variable,
              las = 2,
              col = "steelblue",
              main = "Variable Importance: Genieridium medinae",
              ylab = "Relative Importance")
      
  # ---------------------------------------------
  # Manual permutation importance (Maxnet)
  # ---------------------------------------------
  # This calculates the importance of each variable by permuting its values across
  # the dataset and measuring the drop in correlation between the original and
  # permuted predictions.
  # Interpretation:
  #   - Represents the "predictive importance" of each variable for the model outputs.
  #   - Indicates how much the model predictions rely on that variable in practice.
  #   - Unlike coefficient-based importance, it accounts for non-linearities and
  #     interactions, and is normalized as a percentage of total importance.
      
      # Presence (1) + background (0)
      pa <- c(rep(1, nrow(occ_coords)),
              rep(0, nrow(bg_coords)))
      # Environmental values for presences + background
      env_occ <- terra::extract(bioclim_selected, occ_coords)[, -1]
      env_bg  <- terra::extract(bioclim_selected, bg_coords)[, -1]
      env_values <- rbind(env_occ, env_bg)
      env_values <- na.omit(env_values)
      # Fit full model
      full_model <- maxnet(pa, env_values)
      # Predict baseline suitability
      pred_full <- predict(full_model, env_values, type = "cloglog")
      # Container
      perm_importance <- numeric(ncol(env_values))
      # Permute each variable
      for (i in seq_len(ncol(env_values))) {
        permuted <- env_values
        permuted[, i] <- sample(permuted[, i])
        pred_perm <- predict(full_model, permuted, type = "cloglog")
        # Importance = drop in correlation
        perm_importance[i] <- 1 - cor(pred_full, pred_perm)}
      
      names(perm_importance) <- colnames(env_values)
      perm_importance <- 100 * perm_importance / sum(perm_importance)
      
      print(perm_importance)
      
      barplot(sort(perm_importance, decreasing = TRUE),
              las = 2,
              col = "darkgreen",
              main = "Permutation Importance",
              ylab = "Importance (%)")      
#      
# ------------------------------------------------------------------------------      
# Climatic suitability prediction
      
      # Extract predicted suitability raster from ENMeval output
      suitability_map <- enm_results@predictions[[11]]
      
      # Plot suitability map with occurrence points
      plot(suitability_map,
           main = "Climatic Suitability: Genieridium medinae (LH, RM 1.5)")
      
      points(occ_coords,
             pch = 21,
             bg = "red",
             col = "white",
             cex = 0.8)
# ------------------------------------------------------------------------------
# Expanded projection area (M) for model extrapolation
      
    # Expand occurrence points by 150 km buffer
      expanded_buffer <- st_buffer(occurrence_line, dist = 150000)
      # Reproject to WGS84 (EPSG:4326)
      m_projection <- st_transform(expanded_buffer, crs = 4326)
    # Crop raster to extent of expanded area
      bioclim_m_projection_crop <- crop(bioclim_col, m_projection)
    # Convert sf polygon to terra vector format for masking
      m_projection_vect <- vect(m_projection)
    # Mask raster using the polygon
      bioclim_M_projection <- mask(bioclim_m_projection_crop, m_projection_vect)
    # Clean layer names
      names(bioclim_M_projection) <- paste0("bio_", 1:19)
  
# Predict climatic suitability in the expanded area using the trained Maxnet model
  expanded_suitability <- predict(bioclim_M_projection, optimal_model, 
                                  type = "cloglog", na.rm = TRUE)

# Plot expanded suitability map with occurrence points
  plot(expanded_suitability, main = "Climatic suitability: Expanded area (150 km buffer)")
  points(occ_coords, pch = 21, bg = "red", col = "white", cex = 0.8)
#
# ------------------------------------------------------------------------------
# Final Habitat Estimation: Potential climatic niche vs Realized Distribution

  # Extract raster values at occurrence points
  occ_values <- terra::extract(expanded_suitability, occ_coords)
  
  # Take the last column (the raster value) and remove NAs
  value_col <- ncol(occ_values)
  occ_values_present <- na.omit(occ_values[, value_col])
  
  # Calculate 10th percentile threshold
  threshold_10 <- quantile(occ_values_present, probs = 0.1)
  
  # Create binary suitability map (1 = suitable, 0 = unsuitable)
  suitability_binary <- expanded_suitability >= threshold_10
  
  print(paste("Threshold (10th percentile):", round(threshold_10, 4)))
  
  # Plot binary map to check filtered areas
  plot(suitability_binary, main = "Potential Climatic Niche (10% threshold)", col = c("grey90", "forestgreen"))
  points(occ_coords, pch = 21, bg = "red", cex = 0.5)
  
  # Load Hansen forest cover
  hansen_raster <- rast("F:/Capas/World/Hansen/coverGeoTIF.tif")
  hansen_crs <- crs(hansen_raster)
  
  # Reproject expanded projection area to Hansen CRS
  expanded_m_reproj <- project(m_projection_vect, hansen_crs)
  
  # Crop and mask Hansen using the expanded polygon
  hansen_crop <- crop(hansen_raster, expanded_m_reproj)
  hansen_mask <- mask(hansen_crop, expanded_m_reproj)
  
  # Reproject raster back to WGS84 to match SDM
  hansen_4326 <- project(hansen_mask, crs(suitability_binary))
  
  # Resample to SDM raster resolution (1km)
  hansen_resampled <- resample(hansen_4326, suitability_binary, method = "bilinear")
  
  # Define forested habitat (cover >= 70% "adequate forest cover")
  forest_mask <- hansen_resampled >= 70
  
  # Realized habitat = suitable climate + forest cover
  realized_habitat <- suitability_binary[[1]] * forest_mask
  
  # Calculate pixel area in km2
  pixel_area <- cellSize(suitability_binary[[1]], unit = "km")
  
  # Potential climatic area (sum of suitable pixels)
  km2_climatic <- sum(values(suitability_binary[[1]] * pixel_area), na.rm = TRUE)
  
  # Realized habitat area
  km2_realized <- sum(values(realized_habitat[[1]] * pixel_area), na.rm = TRUE)
  
  # Percentage of habitat loss
  habitat_loss_pct <- 100 - (km2_realized * 100 / km2_climatic)
  
  # Print report
  cat("----------------------------------------------\n",
      "FINAL REPORT: Genieridium medinae\n",
      "----------------------------------------------\n",
      "Potential Climatic Area: ", round(km2_climatic, 2), " km2\n",
      "Forested Habitat Area: ", round(km2_realized, 2), " km2\n",
      "Estimated Habitat Loss: ", round(habitat_loss_pct, 2), "%\n",
      "----------------------------------------------\n")
  
  # Load Colombia administrative map
  colombia_map <- vect("F:/Capas/Colombia/Colombia/COL_adm1.shp")
  
  # Make 0s NA so they are fully transparent
  niche_plot <- suitability_binary[[1]]
  niche_plot[niche_plot == 0] <- NA

  habitat_plot <- realized_habitat[[1]]
  habitat_plot[habitat_plot == 0] <- NA
 
  # Set layout for side-by-side plotting
  par(mfrow = c(1, 2))
  
  # Plot potential climatic niche over DEM
  plot(dem_col, axes = FALSE, main = "Potential climatic niche")  # DEM base
  plot(niche_plot[[1]], col = "orange", alpha = 0.8, legend = FALSE, add = TRUE)
  plot(colombia_map, add = TRUE, border = "gray40")
  points(occ_coords, pch = 21, bg = "red", col = "white", cex = 1)
  
  # Plot realized habitat over DEM
  plot(dem_col, axes = FALSE, main = "Realized habitat (Forest cover)")  # DEM base
  plot(habitat_plot[[1]], col = "orange", alpha = 0.8, legend = FALSE, add = TRUE)
  plot(colombia_map, add = TRUE, border = "gray40")
  points(occ_coords, pch = 21, bg = "red", col = "white", cex = 1)