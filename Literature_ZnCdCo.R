library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf) 
library(marmap)
library(terra)
library(orsifronts)
library(patchwork)
library(tidyverse)
library(viridis)
library(ggnewscale)

setwd("/Users/eggboy/Dropbox/Science/Data/Growth Rates") #setwd

# import the growth rates
ZnCdCo_oceans <- read.csv("/Users/eggboy/Dropbox/Science/Data/Growth Rates/Lit_dZN_concentrations.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
  mutate(Location = as.factor(Location), Reference = as.factor(Reference),
         Year = as.numeric(Year), Month = as.factor(Month), Season = as.factor(Season), Depth = as.numeric(Depth),
         dZn = as.numeric(dZn), dCo = as.numeric(dCo), dCd = as.numeric(dCd),
         Lat = as.numeric(Lat), Long = as.numeric(Long))

# Filter, group, and summarize the data
Surface_SO <- ZnCdCo_oceans %>%
  filter(Depth < 35) %>%
  group_by(Lat, Long, Reference, Location, Year, Season, Month) %>%
  summarize(
    surface_dZn = mean(dZn, na.rm = TRUE),
    surface_dCo = mean(dCo, na.rm = TRUE),
    surface_dCd = mean(dCd, na.rm = TRUE)
  ) %>%
  mutate(
    dCo_dZn_ratio = surface_dCo / (surface_dZn * 1000),
    dCd_dZn_ratio = surface_dCd / (surface_dZn * 1000)
  )


# Transform data points
Surface_sf <- st_as_sf(Surface_SO, coords = c("Long", "Lat"), crs = 4326)  # Set initial CRS
Surface_robinson <- st_transform(Surface_sf, crs = "+proj=robin +lon_0")

# Transform to Antarctic Polar Stereographic (EPSG:3031)
Surface_sf_polar <- st_transform(Surface_sf, crs = 3031)

#loading the world map
world <- ne_countries(scale = "large", returnclass = "sf")
world_robinson <- st_transform(world, crs = "+proj=robin +lon_0")
# Transform to Antarctic Polar Stereographic (EPSG:3031)
world_polar <- st_transform(world, crs = 3031) 

# 1. Retrieve bathymetry data
bathymetry <- getNOAA.bathy(lon1 = -180, lon2 = 180, lat1 = -90, lat2 = -20, 
                            resolution = 4)
# 2. Convert bathy object to a data frame with longitude, latitude, and depth
bathymetry_df <- as.xyz(bathymetry)
# 3. Create a terra raster from the data frame
bathymetry_raster <- rast(bathymetry_df, 
                          type = "xyz",     # Specifies data format
                          crs = "EPSG:4326") # Assign WGS84 CRS
# 4. Plot to verify
plot(bathymetry_raster, main = "Bathymetry Raster (terra)")
# Reproject the raster to EPSG:3031 (Antarctic Polar Stereographic)
bathymetry_polar <- project(bathymetry_raster, "EPSG:3031")
bathymetry_robinson <- project(bathymetry_raster, "+proj=robin +lon_0")

# Create a high-resolution global raster with desired extent
global_raster <- rast(
  nrows = 1800, ncols = 3600,  # Adjust resolution
  xmin = -180, xmax = 180,
  ymin = -90, ymax = 90,
  crs = "EPSG:4326"  # WGS 84
)

# Resample the bathymetry raster to match the global raster
bathymetry_resampled <- resample(bathymetry_raster, global_raster, method = "bilinear")

# Reproject to Robinson projection
bathymetry_robinson <- project(
  bathymetry_resampled,
  "+proj=robin +lon_0=0"
)
# Plot reprojected raster
plot(bathymetry_polar)
plot(bathymetry_robinson)

# 1. Convert the polar raster to a data frame
bathymetry_polar_df <- as.data.frame(bathymetry_polar, xy = TRUE)  # Extract coordinates and values
colnames(bathymetry_polar_df) <- c("x", "y", "Depth")       
bathymetry_polar_df$Depth <- ifelse(bathymetry_polar_df$Depth > 0, NA, bathymetry_polar_df$Depth) # Set positive Depth values (land) to NA

bathymetry_robinson_df <- as.data.frame(bathymetry_robinson, xy = TRUE)  # Extract coordinates and values
colnames(bathymetry_robinson_df) <- c("x", "y", "Depth")       
bathymetry_robinson_df$Depth <- ifelse(bathymetry_robinson_df$Depth > 0, NA, bathymetry_robinson_df$Depth) # Set positive Depth values (land) to NA


# creating the graticule.
graticule <- st_graticule(crs = 3031, lon = seq(-180, 180, by = 30), lat = seq(-80, -30, by = 10))
graticule_rob <- st_graticule(crs = "+proj=robin +lon_0", lon = seq(-180, 180, by = 30), lat = seq(-90, 90, by = 30))

#getting the park fronts from orsifronts package
parkfronts_sf <- st_as_sf(parkfronts, coords = c("longitude", "latitude"), crs = 3031)

# Inset map for Zn (Ross Sea)
ross_inset <- ggplot() +
  geom_sf(data = graticule, color = "#909090", linetype = "solid", size = 0.2) +
  geom_sf(data = Surface_sf_polar %>% 
            filter(!is.na(surface_dZn)), 
          aes(color = surface_dZn), size = 2, alpha = 1) +
  geom_sf(data = world_polar, fill = "#eeeeee", color = "#909090") +
  scale_colour_viridis_c(option = "magma", 
                         limits = c(0, 10),
                         name = "dCd/dZn", direction = -1) +
  coord_sf(xlim = c(-750000, 600000), ylim = c(-2000000, -900000), crs = st_crs(3031)) +
  theme_bw() +
  theme(axis.title = element_blank(), axis.text = element_blank(), 
        axis.ticks = element_blank())



# ---- interpolation for Zn first------------------ ----------
library(MBA)

# Extract the coordinates and values from the Surface_sf_polar object, for Summer data
surface_summer_coords <- st_coordinates(Surface_sf_polar %>% filter(Month %in% c("November", "December", "January", "February")))

# Filter and prepare the data for summer while removing NAs
surface_summer_filtered <- Surface_sf_polar %>%
  filter(Month %in% c("November", "December", "January", "February"), !is.na(surface_dZn))

# Extract coordinates after filtering
surface_summer_coords <- st_coordinates(surface_summer_filtered)  # Ensure only filtered data is used

# Prepare the final data
surface_summer_data <- surface_summer_filtered %>%
  ungroup() %>%                              # Remove any grouping
  st_drop_geometry() %>%                     # Remove geometry column
  bind_cols(as.data.frame(surface_summer_coords)) %>%  # Add filtered coordinates
  select(x = X, y = Y, surface_dZn)          # Select relevant columns

# Run the interpolation using mba.surf
interpolation_result <- mba.surf(
  xyz = surface_summer_data,           # The input data (coordinates and values)
  no.X = 1000,                   # Number of grid points in X direction
  no.Y = 1000,                   # Number of grid points in Y direction
  extend = F                 # Option to extend the grid beyond data points
)

# Check the interpolation result
interpolation_result$xyz.est$z

# Extract the components
x_coords <- interpolation_result$xyz.est$x
y_coords <- interpolation_result$xyz.est$y
z_values <- interpolation_result$xyz.est$z

# Convert to a tidy data frame
interpolated_df <- expand.grid(x = x_coords, y = y_coords) %>%
  mutate(surface_dZn = as.vector(z_values))  # Flatten the matrix

# CREATING A MASK OF SORTS 
# Convert grid points to an sf object
interpolated_sf <- st_as_sf(interpolated_df, coords = c("x", "y"), crs = st_crs(Surface_sf_polar))

# Subset the summer data points
summer_points <- Surface_sf_polar %>% filter(Month %in% c("November", "December", "January", "February"))

# Calculate pairwise distances and take the minimum for each grid point
distances <- st_distance(interpolated_sf, summer_points)  # Pairwise distances
interpolated_sf$distance <- apply(distances, 1, min)      # Minimum distance for each grid point

# Filter grid points within a specified threshold (e.g., 300 km)
threshold <- 500000  # Distance in meters
interpolated_filtered_sf <- interpolated_sf %>% filter(distance <= threshold)

# Convert back to a data frame for plotting and handle surface_dZn conditions
interpolated_filtered_df <- interpolated_filtered_sf %>%
  st_coordinates() %>%
  as.data.frame() %>%
  bind_cols(interpolated_filtered_sf %>% st_drop_geometry()) %>%
  mutate(surface_dZn = ifelse(is.na(surface_dZn) | surface_dZn < 0, 0, surface_dZn)) %>%  # Replace NA or <0 with 0
  rename(x = X, y = Y)

# da plot
SO_Zn_megaplot <- ggplot() +
  geom_raster(data = bathymetry_polar_df, aes(x = x, y = y, fill = Depth), alpha = 0.6) +
  scale_fill_gradient(low = "black", high = "grey", 
                      limits = c(min(bathymetry_polar_df$Depth, na.rm = TRUE), 0), 
                      name = "Depth (m)") +
  new_scale_fill() +
  geom_sf(data = graticule, color = "#909090", linetype = "solid", size = 0.2) +
  geom_tile(data = interpolated_filtered_df, aes(x = x, y = y, fill = surface_dZn)) +
  scale_fill_viridis(option = "turbo", name = "dZn (nM)", 
                     #limits = c(0, 5), 
                     direction = -1) +
  geom_sf(data = world_polar, fill = "#eeeeee", color = "#909090") +
  geom_sf(data = parkfronts_sf, aes(geometry = geometry), 
          color = c("white", "black", "white", "white", "white")[parkfronts_sf$front], 
          linetype = c("dashed", "solid", "solid", "solid", "dashed")[parkfronts_sf$front], 
          size = 0.8) +
  geom_sf(data = Surface_sf_polar %>% filter(Month %in% c("November", "December", "January", "February")), 
          color = "black", size = 0.2, alpha = 1) +
  coord_sf(xlim = c(-6550000, 6550000), ylim = c(-6550000, 6550000), crs = st_crs(3031)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave("dZn_megaplot_SO_summer.svg", SO_Zn_megaplot, width = 5.4, height = 5)


# ----- now for Cadmium --------------------------- --------------------
# Extract the coordinates and values from the Surface_sf_polar object, for Summer data
surface_summer_coords <- st_coordinates(Surface_sf_polar %>% filter(Month %in% c("November", "December", "January", "February"), !is.na(surface_dCd)))

# Filter and prepare the data for summer while removing NAs
surface_summer_filtered <- Surface_sf_polar %>%
  filter(Month %in% c("November", "December", "January", "February"), !is.na(surface_dCd))

# Extract coordinates after filtering
surface_summer_coords <- st_coordinates(surface_summer_filtered)  # Ensure only filtered data is used

# Prepare the final data
surface_summer_data <- surface_summer_filtered %>%
  ungroup() %>%                              # Remove any grouping
  st_drop_geometry() %>%                     # Remove geometry column
  bind_cols(as.data.frame(surface_summer_coords)) %>%  # Add filtered coordinates
  select(x = X, y = Y, surface_dCd)          # Select relevant columns

# Run the interpolation using mba.surf
interpolation_result <- mba.surf(
  xyz = surface_summer_data,           # The input data (coordinates and values)
  no.X = 1000,                   # Number of grid points in X direction
  no.Y = 1000,                   # Number of grid points in Y direction
  extend = F                 # Option to extend the grid beyond data points
)

# Check the interpolation result
interpolation_result$xyz.est$z

# Extract the components
x_coords <- interpolation_result$xyz.est$x
y_coords <- interpolation_result$xyz.est$y
z_values <- interpolation_result$xyz.est$z

# Convert to a tidy data frame
interpolated_df <- expand.grid(x = x_coords, y = y_coords) %>%
  mutate(surface_dCd = as.vector(z_values))  # Flatten the matrix

# CREATING A MASK OF SORTS 
# Convert grid points to an sf object
interpolated_sf <- st_as_sf(interpolated_df, coords = c("x", "y"), crs = st_crs(Surface_sf_polar))

# Subset the summer data points
summer_points <- Surface_sf_polar %>% filter(Month %in% c("November", "December", "January", "February"))

# Calculate pairwise distances and take the minimum for each grid point
distances <- st_distance(interpolated_sf, summer_points)  # Pairwise distances
interpolated_sf$distance <- apply(distances, 1, min)      # Minimum distance for each grid point

# Filter grid points within a specified threshold (e.g., 300 km)
threshold <- 500000  # Distance in meters
interpolated_filtered_sf <- interpolated_sf %>% filter(distance <= threshold)

# Convert back to a data frame for plotting and handle surface_dCd conditions
interpolated_filtered_df <- interpolated_filtered_sf %>%
  st_coordinates() %>%
  as.data.frame() %>%
  bind_cols(interpolated_filtered_sf %>% st_drop_geometry()) %>%
  mutate(surface_dCd = ifelse(is.na(surface_dCd) | surface_dCd < 0, 0, surface_dCd)) %>%  # Replace NA or <0 with 0
  rename(x = X, y = Y)


# da plot
SO_Cd_megaplot <- ggplot() +
  geom_raster(data = bathymetry_polar_df, aes(x = x, y = y, fill = Depth), alpha = 0.6) +
  scale_fill_gradient(low = "black", high = "grey", 
                      limits = c(min(bathymetry_polar_df$Depth, na.rm = TRUE), 0), 
                      name = "Depth (m)") +
  new_scale_fill() +
  geom_sf(data = graticule, color = "#909090", linetype = "solid", size = 0.2) +
  geom_tile(data = interpolated_filtered_df, aes(x = x, y = y, fill = surface_dCd)) +
  scale_fill_viridis(option = "turbo", name = "dCd (pM)", 
                     #limits = c(0, 100), 
                     direction = -1) +
  geom_sf(data = world_polar, fill = "#eeeeee", color = "#909090") +
  geom_sf(data = parkfronts_sf, aes(geometry = geometry), 
          color = c("white", "black", "white", "white", "white")[parkfronts_sf$front], 
          linetype = c("dashed", "solid", "solid", "solid", "dashed")[parkfronts_sf$front], 
          size = 0.8) +
  geom_sf(data = Surface_sf_polar %>% filter(Month %in% c("November", "December", "January", "February")), 
          color = "black", size = 0.2, alpha = 1) +
  coord_sf(xlim = c(-6550000, 6550000), ylim = c(-6550000, 6550000), crs = st_crs(3031)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave("dCd_megaplot_SO_summer.svg", SO_Cd_megaplot, width = 5.4, height = 5)


# ----- now for Cobalt? ---------------------------
# Extract the coordinates and values from the Surface_sf_polar object, for Summer data
surface_summer_coords <- st_coordinates(Surface_sf_polar %>% filter(Month %in% c("November", "December", "January", "February"), !is.na(surface_dCo)))

# Filter and prepare the data for summer while removing NAs
surface_summer_filtered <- Surface_sf_polar %>%
  filter(Month %in% c("November", "December", "January", "February"), !is.na(surface_dCo))

# Extract coordinates after filtering
surface_summer_coords <- st_coordinates(surface_summer_filtered)  # Ensure only filtered data is used

# Prepare the final data
surface_summer_data <- surface_summer_filtered %>%
  ungroup() %>%                              # Remove any grouping
  st_drop_geometry() %>%                     # Remove geometry column
  bind_cols(as.data.frame(surface_summer_coords)) %>%  # Add filtered coordinates
  select(x = X, y = Y, surface_dCo)          # Select relevant columns

# Run the interpolation using mba.surf
interpolation_result <- mba.surf(
  xyz = surface_summer_data,           # The input data (coordinates and values)
  no.X = 1000,                   # Number of grid points in X direction
  no.Y = 1000,                   # Number of grid points in Y direction
  extend = F                 # Option to extend the grid beyond data points
)

# Check the interpolation result
interpolation_result$xyz.est$z

# Extract the components
x_coords <- interpolation_result$xyz.est$x
y_coords <- interpolation_result$xyz.est$y
z_values <- interpolation_result$xyz.est$z

# Convert to a tidy data frame
interpolated_df <- expand.grid(x = x_coords, y = y_coords) %>%
  mutate(surface_dCo = as.vector(z_values))  # Flatten the matrix

# CREATING A MASK OF SORTS 
# Convert grid points to an sf object
interpolated_sf <- st_as_sf(interpolated_df, coords = c("x", "y"), crs = st_crs(Surface_sf_polar))

# Subset the summer data points
summer_points <- Surface_sf_polar %>% filter(Month %in% c("November", "December", "January", "February"), !is.na(surface_dCo))

# Calculate pairwise distances and take the minimum for each grid point
distances <- st_distance(interpolated_sf, summer_points)  # Pairwise distances
interpolated_sf$distance <- apply(distances, 1, min)      # Minimum distance for each grid point

# Filter grid points within a specified threshold (e.g., 300 km)
threshold <- 500000  # Distance in meters
interpolated_filtered_sf <- interpolated_sf %>% filter(distance <= threshold)

# Convert back to a data frame for plotting and handle surface_dCo conditions
interpolated_filtered_df <- interpolated_filtered_sf %>%
  st_coordinates() %>%
  as.data.frame() %>%
  bind_cols(interpolated_filtered_sf %>% st_drop_geometry()) %>%
  mutate(surface_dCo = ifelse(is.na(surface_dCo) | surface_dCo < 0, 0, surface_dCo)) %>%  # Replace NA or <0 with 0
  rename(x = X, y = Y)


# da plot
SO_Co_megaplot <- ggplot() +
  geom_raster(data = bathymetry_polar_df, aes(x = x, y = y, fill = Depth), alpha = 0.6) +
  scale_fill_gradient(low = "black", high = "grey", 
                      limits = c(min(bathymetry_polar_df$Depth, na.rm = TRUE), 0), 
                      name = "Depth (m)") +
  new_scale_fill() +
  geom_sf(data = graticule, color = "#909090", linetype = "solid", size = 0.2) +
  geom_tile(data = interpolated_filtered_df, aes(x = x, y = y, fill = surface_dCo)) +
  scale_fill_viridis(option = "turbo", name = "dCo (pM)", 
                     limits = c(0, 80), 
                     direction = -1) +
  geom_sf(data = world_polar, fill = "#eeeeee", color = "#909090") +
  geom_sf(data = parkfronts_sf, aes(geometry = geometry), 
          color = c("white", "black", "white", "white", "white")[parkfronts_sf$front], 
          linetype = c("dashed", "solid", "solid", "solid", "dashed")[parkfronts_sf$front], 
          size = 0.8) +
  geom_sf(data = Surface_sf_polar %>% filter(Month %in% c("November", "December", "January", "February"), !is.na(surface_dCo)), 
          color = "black", size = 0.2, alpha = 1) +
  coord_sf(xlim = c(-6550000, 6550000), ylim = c(-6550000, 6550000), crs = st_crs(3031)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

ggsave("dCo_megaplot_SO_summer.svg", SO_Co_megaplot, width = 5.4, height = 5)


# ---- now for the Robinson whole world perspective for the Co and Cd ratios ----------------
# Filter and prepare the data for dCo/dZn ratio
surface_dCo_dZn_filtered <- Surface_robinson %>%
  filter(!is.na(dCo_dZn_ratio))  # Filter out NA values

# Extract coordinates after filtering
surface_dCo_dZn_coords <- st_coordinates(surface_dCo_dZn_filtered)  # Ensure alignment with filtered data

surface_dCo_dZn_data <- surface_dCo_dZn_filtered %>%
  ungroup() %>%                              # Remove any grouping
  st_drop_geometry() %>%                     # Remove geometry column
  bind_cols(as.data.frame(surface_dCo_dZn_coords)) %>%  # Add filtered coordinates
  dplyr::select(X, Y, dCo_dZn_ratio) %>%     # Select only coordinates and dCo_dZn_ratio
  rename(x = X, y = Y) 

# Run the interpolation using mba.surf
interpolation_result <- mba.surf(
  xyz = surface_dCo_dZn_data,  # The input data (coordinates and values)
  no.X = 1000,                # Number of grid points in X direction
  no.Y = 1000,                # Number of grid points in Y direction
  extend = TRUE              # Does the grid extend beyond the data points?
)

# Extract the components of the interpolation
x_coords <- interpolation_result$xyz.est$x
y_coords <- interpolation_result$xyz.est$y
z_values <- interpolation_result$xyz.est$z  # Interpolated dCo/dZn_ratio values

# Convert the interpolated result to a tidy data frame
interpolated_df <- expand.grid(x = x_coords, y = y_coords) %>%
  mutate(dCo_dZn_ratio = as.vector(z_values))  # Flatten the matrix of interpolated values

# Convert the interpolated grid to an sf object for spatial processing
interpolated_sf <- st_as_sf(interpolated_df, coords = c("x", "y"), crs = st_crs(Surface_robinson))

# Calculate pairwise distances between grid points and original data points
distances <- st_distance(interpolated_sf, surface_dCo_dZn_filtered)  # Pairwise distances
interpolated_sf$distance <- apply(distances, 1, min)                 # Minimum distance for each grid point

# Apply a threshold (e.g., 500,000 meters)
threshold <- 500000  # Distance in meters
interpolated_filtered_sf <- interpolated_sf %>%
  filter(distance <= threshold)  # Retain points within the threshold distance

# Convert back to a data frame for plotting and handle dCo_dZn_ratio conditions
interpolated_filtered_df <- interpolated_filtered_sf %>%
  st_coordinates() %>%  # Extract coordinates into a matrix
  as.data.frame() %>%   # Convert to a data frame
  bind_cols(interpolated_filtered_sf %>% st_drop_geometry()) %>%  # Combine with original data minus geometry
  mutate(dCo_dZn_ratio = ifelse(is.na(dCo_dZn_ratio) | dCo_dZn_ratio < 0, 0, dCo_dZn_ratio)) %>%  # Replace NA or <0 with 0
  rename(x = X, y = Y)  # Rename coordinate columns for clarity

CoZn_worldplot <- ggplot() +
  geom_raster(data = bathymetry_robinson_df, aes(x = x, y = y, fill = Depth), alpha = 0.6) +
  scale_fill_gradient(low = "black", high = "grey", 
                      limits = c(min(bathymetry_robinson_df$Depth, na.rm = TRUE), 0), 
                      name = "Depth (m)") +
  new_scale_fill() +
  geom_sf(data = graticule_rob, color = "#909090", linetype = "solid", size = 0.2) +  # Graticule for geographic lines
  geom_tile(data = interpolated_filtered_df %>% st_drop_geometry(), aes(x = x, y = y, fill = dCo_dZn_ratio)) +  # Interpolated grid
  scale_fill_viridis_c(option = "viridis", name = "dCo/dZn Ratio", 
                       limits = c(0, max(interpolated_filtered_sf$dCo_dZn_ratio, na.rm = TRUE)),  # Adjust limits dynamically
                       direction = -1) +
  geom_sf(data = world_robinson, fill = "#eeeeee", color = "#909090") +  # World map layer
  geom_sf(data = surface_dCo_dZn_filtered, color = "black", size = 0.2, alpha = 1) +  # Original data points for reference
  coord_sf(crs = "+proj=robin +lon_0=0", expand = TRUE) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 16))

ggsave("dCo_dZn_ratioplot.svg", CoZn_worldplot, width = 6.5, height = 5)

# Filter and prepare the data for dCd/dZn ratio
surface_dCd_dZn_filtered <- Surface_robinson %>%
  filter(!is.na(dCd_dZn_ratio))  # Filter out NA values

# Extract coordinates after filtering
surface_dCd_dZn_coords <- st_coordinates(surface_dCd_dZn_filtered)  # Ensure alignment with filtered data

surface_dCd_dZn_data <- surface_dCd_dZn_filtered %>%
  ungroup() %>%                              # Remove any grouping
  st_drop_geometry() %>%                     # Remove geometry column
  bind_cols(as.data.frame(surface_dCd_dZn_coords)) %>%  # Add filtered coordinates
  dplyr::select(X, Y, dCd_dZn_ratio) %>%     # Select only coordinates and dCo_dZn_ratio
  rename(x = X, y = Y) 

# Run the interpolation using mba.surf
interpolation_result <- mba.surf(
  xyz = surface_dCd_dZn_data,  # The input data (coordinates and values)
  no.X = 1000,                # Number of grid points in X direction
  no.Y = 1000,                # Number of grid points in Y direction
  extend = TRUE              # Prevent the grid from extending beyond the data points
)

# Extract the components of the interpolation
x_coords <- interpolation_result$xyz.est$x
y_coords <- interpolation_result$xyz.est$y
z_values <- interpolation_result$xyz.est$z  # Interpolated dCo/dZn_ratio values

# Convert the interpolated result to a tidy data frame
interpolated_df <- expand.grid(x = x_coords, y = y_coords) %>%
  mutate(dCd_dZn_ratio = as.vector(z_values))  # Flatten the matrix of interpolated values

# Convert the interpolated grid to an sf object for spatial processing
interpolated_sf <- st_as_sf(interpolated_df, coords = c("x", "y"), crs = st_crs(Surface_robinson))

# Calculate pairwise distances between grid points and original data points
distances <- st_distance(interpolated_sf, surface_dCd_dZn_filtered)  # Pairwise distances
interpolated_sf$distance <- apply(distances, 1, min)                 # Minimum distance for each grid point

# Apply a threshold (e.g., 500,000 meters)
threshold <- 500000  # Distance in meters
interpolated_filtered_sf <- interpolated_sf %>%
  filter(distance <= threshold)  # Retain points within the threshold distance

# Convert back to a data frame for plotting and handle dCo_dZn_ratio conditions
interpolated_filtered_df <- interpolated_filtered_sf %>%
  st_coordinates() %>%  # Extract coordinates into a matrix
  as.data.frame() %>%   # Convert to a data frame
  bind_cols(interpolated_filtered_sf %>% st_drop_geometry()) %>%  # Combine with original data minus geometry
  mutate(dCd_dZn_ratio = ifelse(is.na(dCd_dZn_ratio) | dCd_dZn_ratio < 0, 0, dCd_dZn_ratio)) %>%  # Replace NA or <0 with 0
  rename(x = X, y = Y)  # Rename coordinate columns for clarity

  CdZn_worldplot <- ggplot() +
    geom_raster(data = bathymetry_robinson_df, aes(x = x, y = y, fill = Depth), alpha = 0.6) +
    scale_fill_gradient(low = "black", high = "grey", 
                        limits = c(min(bathymetry_robinson_df$Depth, na.rm = TRUE), 0), 
                        name = "Depth (m)") +
    new_scale_fill() +
  geom_sf(data = graticule_rob, color = "#909090", linetype = "solid", size = 0.2) +
  geom_tile(data = interpolated_filtered_df %>% 
              st_drop_geometry() %>% 
              filter(dCd_dZn_ratio <= 3),  # Filter data for the first range
            aes(x = x, y = y, fill = dCd_dZn_ratio)) +
  scale_fill_viridis_c(option = "viridis", name = "dCd/dZn (0-3)", limits = c(0, 3), direction = -1) +
  ggnewscale::new_scale_fill() +  # Add a new fill scale
  geom_tile(data = interpolated_filtered_df %>% 
              st_drop_geometry() %>% 
              filter(dCd_dZn_ratio > 3 & dCd_dZn_ratio <= 8),  # Filter data for the second range
            aes(x = x, y = y, fill = dCd_dZn_ratio)) +
  scale_fill_viridis_c(option = "rocket", name = "dCd/dZn (3-8)", limits = c(3, 8), direction = 1) +
  geom_sf(data = world_robinson, fill = "#eeeeee", color = "#909090") +  # World map layer
  geom_sf(data = surface_dCd_dZn_filtered, color = "black", size = 0.2, alpha = 1) +  # Original data points for reference
  coord_sf(crs = "+proj=robin +lon_0", expand = TRUE) +  # Robinson projection
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 16))
  
  ggsave("dCd_dZn_ratioplot.svg", CdZn_worldplot, width = 6.5, height = 5)
  
# -------------- attempting data import from ODV weighted gridding average --------------
# Load the NetCDF file
raster_dCddZn <- read.csv("/Users/eggboy/Dropbox/Science/Data/Growth Rates/dCddZn_ODVraster.csv", fileEncoding="UTF-8-BOM", header = TRUE)
  
# Convert data to sf object
dCddZn_sf <- st_as_sf(raster_dCddZn, coords = c("Long", "Lat"), crs = 4326) 
dCddZn_robinson <- st_transform(dCddZn_sf, crs = "+proj=robin +lon_0=0")

# Extract coordinates and flatten the data
dCddZn_flat <- dCddZn_robinson %>%
  st_drop_geometry() %>%
  cbind(st_coordinates(dCddZn_robinson))  # Add x and y columns

# Filter out rows with NA values in Estimated.dCd_dZn_ratio
dCddZn_flat <- dCddZn_flat %>%
  filter(!is.na(Estimated.dCd_dZn_ratio)) %>%  
  as.data.frame()


CdZn_worldplot <- ggplot() +
  geom_sf(data = graticule, color = "#909090", linetype = "solid", size = 0.2) +  # Graticule
  geom_tile(data = dCddZn_flat, aes(x = X, y = Y, fill = Estimated.dCd_dZn_ratio)) +  # Interpolated grid
  scale_fill_viridis(option = "viridis",
                       name = "dCd/dZn Ratio",
                       direction = -1) +
  geom_sf(data = world_robinson, fill = "#eeeeee", color = "#909090") +  # World map
  coord_sf(crs = "+proj=robin +lon_0=0", expand = TRUE) +  # Robinson projection
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 16))


# ---- Locations of strains --------------
strain_locations <- read.csv("/Users/eggboy/Dropbox/Science/Data/Growth Rates/CloneLocations.csv", fileEncoding="UTF-8-BOM", header = TRUE) %>%
  mutate(Species = as.factor(Species), Ref = as.factor(Ref),
         Clone = as.factor(Clone), Lat = as.numeric(Lat), Long = as.numeric(Long))

# Transform strain locations into an sf object with geographic coordinates (WGS84)
strain_sf <- strain_locations %>%
  filter(!is.na(Lat) & !is.na(Long)) %>%  # Exclude rows with missing coordinates
  st_as_sf(coords = c("Long", "Lat"), crs = 4326) %>%  # Convert to sf object
  st_transform(crs = 3031)  # Transform to polar projection (EPSG:3031)

# Create labels with Species and Clone in brackets
strain_sf <- strain_sf %>%
  mutate(label = paste(Species, "(", Clone, ")", sep = ""))

# Plot the map
clone_locations <- ggplot() +
 geom_raster(data = bathymetry_polar_df, aes(x = x, y = y, fill = Depth), alpha = 0.6) +
 scale_fill_gradient(low = "black", high = "grey", 
                     limits = c(min(bathymetry_polar_df$Depth, na.rm = TRUE), 0), 
                     name = "Depth (m)") +
  geom_sf(data = world_polar, fill = "#eeeeee", color = "#909090") +  # Landmass
  geom_sf(data = graticule, color = "#909090", linetype = "solid", size = 0.2) +  # Graticule
  geom_sf(data = strain_sf, size = 2, shape = 23, colour = "black", fill = "white") +  # Species points
  geom_text(data = st_drop_geometry(strain_sf), aes(x = st_coordinates(strain_sf)[, 1],
                                                    y = st_coordinates(strain_sf)[, 2],
                                                    label = label), size = 3, nudge_y = 100000) +  # Labels
  coord_sf(xlim = c(-6550000, 6550000), ylim = c(-6550000, 6550000), crs = st_crs(3031)) +
  theme_bw() +
  theme(axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank())

ggsave("clone_locations.svg", clone_locations, width = 5.4, height = 5)
