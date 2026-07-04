
#========================================================================
# TASK 2. ESTIMATION OF TIMBER VOLUME WITH THE RANDOM FOREST (RF) METHOD
#========================================================================

rm(list = ls()) #Clean R's memory

# install.packages("sf")
# install.packages("terra")
# install.packages("exactextractr") 
# install.packages("randomForest")
# install.packages("yaImpute")
# install.packages("gower")

# Load the installed libraries
require(sf)
require(terra)
require(randomForest)
require(yaImpute)
require(gower)
require(exactextractr)

# Set your own working directory
setwd("C:\\Users\\Documents\\ope\\kkjk\\2024\\harkat24\\S2P2")
# setwd("W:\\ope\\kkjk\\2024\\harkat24\\S2P2")

# Load custom functions
source("S2_radar_functions.r")

#=============================================================
# Extract Sentinel-2 for Juupajoki plots 
#=============================================================

# Read in field plots
juu_plots <- read.csv("Juupajoki_field_data.csv", sep=",", dec=".")
tail(juu_plots)

# Set coordinate system
juu_plots <- st_as_sf(juu_plots, coords = c("xcoordinate","ycoordinate"))
st_crs(juu_plots) <- 3067

# Read Sentinel-2 and set coordinate system
sentinel2 <- rast("S2_Juupajoki.tif")
crs(sentinel2) <- "EPSG:3067"

# Plot satellite image and add the field plot locations
plotRGB(sentinel2, r=7, g=3, b=2, stretch="lin") 
plot(juu_plots, col="yellow", pch=16, cex=.75, add=T)

# Set plot radius based on sample plot type
juu_plots$rad <- 9
juu_plots$rad[juu_plots$sampleplottype==2] <- 5.64
juu_plots$rad[juu_plots$sampleplottype==4] <- 12.62

# Convert plot centers to polygons based on the set radius
buffered <- st_buffer(juu_plots, dist = juu_plots$rad)

# Extract as data frame, get a normalized weight for each pixel within a polygon
image_df <- exact_extract(sentinel2, buffered, fun="mean", weights="area",
                          max_cells_in_memory = 4e+07) 

# Rename columns and view last rows
names(image_df) <- c("blue","green","red","re1","re2","re3","nir","nirnarrow", 
                     "swir1","swir2")
tail(image_df)

# Merge with plot data
juu_s2 <- cbind(juu_plots, image_df)

# View last rows
tail(juu_s2)

#=============================================================
# Extract PALSAR-2 for Juupajoki plots 
#=============================================================

# Read PALSAR-2 bands
hh <- rast("PALSAR-2_HH_Juupajoki.tif")
hv <- rast("PALSAR-2_HV_Juupajoki.tif")

# Merge bands into a single raster brick
palsar2 <- c(hh,hv)
crs(palsar2) <- "EPSG:3067"

# Plot HH band
plotRGB(palsar2, r=1, g=1, b=1, stretch="lin") 
plot(juu_plots, col="yellow", pch=16, cex=.75, add=T)

# Extract radar data from a fixed 20 m radius around each plot to decrease noise
juu_plots$rad <- 20
buffered <- st_buffer(juu_plots, dist = juu_plots$rad)

# Extract as data frame
image_df <- exact_extract(palsar2, buffered, fun="mean", weights="area") 

# Rename columns 
names(image_df) <- c("hh", "hv")

# Merge with previous data
juu <- cbind(juu_s2, image_df)

# Remove geometry column
juu <- st_set_geometry(juu, NULL)

# view last rows
tail(juu)

#=============================================================
# Merge with Orivesi data sets
#=============================================================

# Read plots from the remaining Orivesi area
ori <- read.csv("Orivesi_field+RS_data.csv", sep=",", dec=".")

# Check that the column names match
names(ori) == names(juu)

# Merge Juupajoki and Orivesi plots into new data frame d
d <- rbind(juu, ori)

#=============================================================
# Construct random forest models
#=============================================================

# Close previous plots and reset graphics
dev.off() 

# View column names
names(d)

# Sentinel-2 model
s2 <- randomForest(volume ~ blue + green + red + re1 + re2 + re3 + nir + nirnarrow + 
                     swir1 + swir2, data = d)

plot(predict(s2), d$volume, xlim=c(0,800), ylim=c(0,800), xlab="Predicted volume, m3/ha", 
     ylab="Observed volume, m3/ha", main="Random forest with Sentinel-2 only", abline(0,1))

rmse(predict(s2), d$volume)
relrmse(predict(s2), d$volume)

varImpPlot(s2,type=2)

# PALSAR-2 model
p2 <- randomForest(volume ~ hh + hv, data = d)

plot(predict(p2), d$volume, xlim=c(0,800), ylim=c(0,800), xlab="Predicted volume, m3/ha", 
     ylab="Observed volume, m3/ha", main="Random forest with PALSAR-2 only", abline(0,1))

rmse(predict(p2), d$volume)
relrmse(predict(p2), d$volume)

varImpPlot(p2,type=2)

# Sentinel-2 + PALSAR-2 model
s2p2 <- randomForest(volume ~ blue + green + red + re1 + re2 + re3 + nir + nirnarrow + 
                       swir1 + swir2 + hh + hv, data = d)

plot(predict(s2p2), d$volume, xlim=c(0,800), ylim=c(0,800), xlab="Predicted volume, m3/ha", 
     ylab="Observed volume, m3/ha", main="Random forest with Sentinel-2 and PALSAR-2", 
     abline(0,1))

rmse(predict(s2p2), d$volume)
relrmse(predict(s2p2), d$volume)

varImpPlot(s2p2, type=2)

#=============================================================
# Process map data
#=============================================================

# Resample PALSAR-2 to same resolution as Sentinel-2
palsar2 <- resample(palsar2, sentinel2, method="bilinear") 

# Stack S2 and  layers into a single raster
sentinel2 <- c(sentinel2, palsar2) 

# Visualize
plotRGB(sentinel2, r=11, g=4, b=3, stretch="lin")

# Convert raster stack into data frame
imagedataframe <- as.data.frame(sentinel2) 

# Add column names for bands
names(imagedataframe) <- c("blue", "green", "red", "re1", "re2", "re3", "nir", 
                           "nirnarrow", "swir1", "swir2", "hh", "hv")

# Read NFI volume raster and set coordinate reference system
nfivol <- rast("Juupajoki_NFI_Volume.tif")
crs(nfivol) <- "EPSG:3067"

# Set up a forest mask: areas that are not nodata are forest
formask <- !is.na(nfivol)

# Plot the mask
plot(formask)

# Convert also forest mask into data frame
maskdataframe <- as.data.frame(formask)


#=============================================================
# Predict with the earlier model object s2p2
#=============================================================

# Apply the previous rf model with also radar bands
pred <- predict(s2p2, imagedataframe)

# Replace the predicted values for non-forest areas with nodata
pred[! maskdataframe$Layer_1] <- NA

# Construct a new map
# Make a copy of the first image band
volmap <- sentinel2[[1]] 

# Replace original values by model predictions
values(volmap) <- pred 

# Show the  map
plot(volmap) 

# Write the map as a tif file
writeRaster(volmap,"RF_Volume.tif",overwrite=T)


#=============================================================
# Comparison
#=============================================================

# Initialize plotting
dev.off()
par(mfrow=c(1,2))

### RF map ###

v_s2p2 <- pred[! is.na(pred)] #Volume vector without NA's
mean(v_s2p2)
hist(v_s2p2, main="RF volume")

### NFI map ###

v_nfi <- nfivol[! is.na(nfivol)] #Volume vector without NA's
mean(v_nfi)
hist(v_nfi, main="NFI volume")

### Extra: compare map pixels with each other and plot with hexbin library ###

rmse(v_nfi, v_s2p2)
bias(v_nfi, v_s2p2)

# library(hexbin)
# bin <- hexbin(v_nfi, v_s2p2)
# plot(bin)
# dev.off()

#===============================================================================
# TASK 3: SPECIES-SPECIFIC VOLUME ESTIMATION WITH THE K NEAREST NEIGHBOR METHOD 
#===============================================================================

# The code continues from where the previous task finished

# Define a vector with response variable names
yvar <- c("pinevol", "sprucevol", "decidvol") 

# Define a weight for each response variable for variable selection
weights <- c(1,1,1) 

# Define the column range in d that contains the predictor variables 
colmin <- 8  # 8 = predictors start from 8th column (blue)
colmax <- 19 # 19 = predictors end at 19 column (swir2)

# KNN parameters
prednum <- 5        # How many predictors to search for
KNNK <- 5           # How many nearest neighbors for KNN
met <- "msn"        # Distance metric for knn. Type ?yai for more info.
wm <- "dstWeighted" # Weighting of neighbors = inverse distance .

# Simulated annealing parameters
t_ini <- 0.2        # Initial temperature at simulated annealing
n_iter <- 10000     # How many iterations in optimization; larger is better but takes longer

# Extract response variables into their own data frame "ytrain" 
ytrain<-as.data.frame(d[,yvar]) 

# Rename columns
names(ytrain)<-yvar 


#=============================================================
# Variable selection by simulated annealing
#=============================================================

t <- t_ini                                     # Initial temperature
s <- sample(colmin:colmax, prednum, replace=F) # Initial variables

e <- yaisel(s) # Run KNN  for initial variables

ebest <- e     # Save initial mean rmse
sbest <- s     # Save initial variable combination

k <- 0         # Initialize iteration counter

while(k < n_iter){ 
  
  # sdot  = new experimental solution
  # s     = current solution to be improved,
  # sbest = best solution ever found
  
  sdot <- PickNewSolution(s, k, n_iter) # New candidate variables
  edot <- yaisel(sdot)                  # KNN result for the new candidate variables
  
  # Implement the simulated annealing algorithm
  if(exp((-(edot-e))/t) > runif(1)){
    e <- edot
    s <- sdot
  }
  if(edot < ebest){
    ebest <- edot
    sbest <- sdot
  }
  t <- max(0, -.2/.8*k/n_iter+.2)      # Cool temperature
  k <- k+1
  
}

names(d)[sbest] # Print selected variable names 
ebest           # Print the mean rmse with these variables

#=============================================================
# Fitting a KNN model with the selected variables
#=============================================================

# Type a vector with the selected variable names here
xvar <- c("green",  "nir", "nirnarrow", "swir1", "swir2")

xtrain <- d[, xvar] # Extract x-variables into their own raster

tail(xtrain) # View last rows

# Rename rows to avoid problems later on
row.names(ytrain) <- 1:nrow(ytrain)
row.names(xtrain) <- 1:nrow(xtrain)

# Train the model
knn <- yai(y=ytrain, x=xtrain, method = met, k=KNNK); 

# Get and view predicted values
pred <- impute(knn, k=KNNK, method=wm) 

tail(pred)  # Variables with .o = observed, without .o = predicted

# Calculate relative RMSEs for all species!
# Draw scatter plots with estimated volume on x and observed volume on y axis
# for all species! Remember to add an 1:1 line.

dev.off()

# Pine
plot(pred$pinevol, pred$pinevol.o, xlim = c(0, 600), ylim = c(0, 600), 
     main = "Pine volume, m3 / ha"); abline(0,1)
rmse(pred$pinevol.o, pred$pinevol)
relrmse(pred$pinevol.o, pred$pinevol)

# Spruce
plot(pred$sprucevol, pred$sprucevol.o, xlim = c(0, 600), ylim = c(0, 600), 
     main = "Spruce volume, m3 / ha"); abline(0,1)
rmse(pred$sprucevol.o, pred$sprucevol)
relrmse(pred$sprucevol.o, pred$sprucevol)

# Deciduous
plot(pred$decidvol, pred$decidvol.o, xlim = c(0, 600), ylim = c(0, 600), 
     main = "Deciduous volume, m3 / ha"); abline(0,1)
rmse(pred$decidvol.o, pred$decidvol)
relrmse(pred$decidvol.o, pred$decidvol)


#=============================================================
# Species mapping
#=============================================================

# Use data extracted previously from the images
df <- imagedataframe[, xvar]

# Change the row names to start where the training data row names ended
row.names(df) <- as.numeric(row.names(df)) + nrow(d)

# Finding the nn references for the test data based on the previous model
knn2 <- newtargets(knn, newdata=df)

# Prediction
pred  <- impute(knn2, vars=yvars(knn2), method=wm, k=KNNK)

tail(pred)

# Replace the predicted values for non-forest areas with nodata
pred$pinevol[! maskdataframe$Layer_1] <- NA
pred$sprucevol[! maskdataframe$Layer_1] <- NA
pred$decidvol[! maskdataframe$Layer_1] <- NA

# Construct a new map

# Make copies of an image band
pinevmap <- sentinel2[[1]] 
sprucevmap <- sentinel2[[1]]
decidvmap <- sentinel2[[1]] 

# Replace values
values(pinevmap) <- pred$pinevol 
values(sprucevmap) <- pred$sprucevol 
values(decidvmap) <- pred$decidvol 

# Stack species layers
s <- c(pinevmap, sprucevmap, decidvmap)

# Write stacked raster as .tif
writeRaster(s, "Species_stack.tif", overwrite=T) 
