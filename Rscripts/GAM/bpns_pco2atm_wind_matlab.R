if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("rhdf5")
library(rhdf5)


# Define the path to your .mat file
mat_file_path <- "atm_pco2_19822020sansmaskcotier.mat"  # Replace with the actual path to your .mat file
mat_file_path2 <- "ERA5windfirst025_19822020.mat"  # Replace with the actual path to your .mat file

# List all groups and datasets in the file
h5ls(mat_file_path2)

# Read a specific dataset (replace 'dataset_name' with the actual name)
data <- h5read(mat_file_path, "atm")
data2 <- h5read(mat_file_path2, "windfirst")
dim(data2)

# find row and column for stations
lon_min = -179.875
lon_max = 179.875  
lat_min = -89.875
lat_max = 89.875 
resolution = 0.25

lon = seq(lon_min, lon_max + resolution, resolution) 
lat = seq(lat_min, lat_max + resolution, resolution) 

lat_bpns <- which.min(abs(lat - 51.57989))
lon_bpns <- which.min(abs(lon - 2.993217))
lat_nas <- which.min(abs(lat - 45.698))
lon_nas <- which.min(abs(lon - 13.708))

# extract time series
# bpns
BPNS_pco2atm <- data[,lat_bpns,lon_bpns]
BPNS_wind <- data2[,lat_bpns,lon_bpns]
year <- c(rep(1982:2020,12))
year <- year[order(year)]
month <- c(rep(1:12,39))
BPNS <- data.frame("year" = year, "month" = month, "pco2atm" = BPNS_pco2atm, "wind" = BPNS_wind)

write.csv(BPNS, "bpns_pco2atm_windspeed.csv",row.names = F)


# northern adriatic
nas_pco2atm <- data[,(lat_nas-3),lon_nas] # lat_nas - 3 was closest point with values and is 0.75° to the west
nas_wind <- data2[,lat_nas,lon_nas]
year <- c(rep(1982:2020,12))
year <- year[order(year)]
month <- c(rep(1:12,39))
nas <- data.frame("year" = year, "month" = month, "pco2atm" = nas_pco2atm, "wind" = nas_wind)

write.csv(nas, "nas_pco2atm_windspeed.csv",row.names = F)
