# Load packages for creating GAM models
packages <- list("mgcv",
                 "lubridate",
                 "ggplot2",
                 "stats",
                 "xts",
                 "lwdataexplorer",
                 "plyr",
                 "dplyr")
lapply(packages, library, character.only=TRUE)

# Load functions for creating GAM models
source(paste0(wd,"Rscripts/GAM/functions_GAM_v1_SP.R"))
source(paste0(wd,"Rscripts/GAM/LW_DataAccess_functions_for_GAM.R"))

  select_vars <- function(df_c){
  
    vars <- c("Station","Time","Air pressure (hPa)","Air temperature (°C)","Average wind speed @ 10m (m/s)",
    "Salinity (PSU)","Water temperature (°C)","xCO2 air (ppm)", "xCO2 sea (ppm)",
    "CO2-atm concentration (ppm)","CO2-seawater concentration (ppm)")
    vars_available <- vars[which(vars %in% colnames(df_c))]
    df_c1 <- df_c %>% 
                select(all_of(vars_available))
    df_c1 <- subset(df_c1, Station == "Buoy at C-Power")
  return(df_c1)
  }

get_carbon_Thornton_buoy_data <- function(start, stop, station_selection, station_region){
  
 #### Code used to get the data from the Life Watch Portal ####
  
  ################## ABIOTIC ##################
  
  #Select station data, dates and stations of interest. 
  df_carbon <- NULL
  for (year in c(2018,2019,2020,2021,2022,2023)){
  df_all <- lwdataexplorer::getBuoyData(paste0(year,"-01-01"), paste0(year,"-12-31"), "All")
  
  df_carbon1 <- subset(df_all, Station == "Buoy at C-Power")
  df_carbon1 <- df_carbon1 %>% 
    select(Station,Time,`Air pressure (hPa)`,`Air temperature (°C)`,`Average wind speed @ 10m (m/s)`,
           `Salinity (PSU)`,`Water temperature (°C)`,
           `CO2-atm concentration (ppm)`,`CO2-seawater concentration (ppm)`)
  
  
  df_carbon <- rbind(df_carbon, df_carbon1)
  
  }
 
df_2018 <- lwdataexplorer::getBuoyData("2018-01-01", "2018-12-31", "All")
df_2019 <- lwdataexplorer::getBuoyData("2019-01-01", "2019-12-31", "All")
df_2020 <- lwdataexplorer::getBuoyData("2020-01-01", "2020-12-31", "All")
df_2021 <- lwdataexplorer::getBuoyData("2021-01-01", "2021-12-31", "All")
df_2022 <- lwdataexplorer::getBuoyData("2022-01-01", "2022-12-31", "All")
df_2023 <- lwdataexplorer::getBuoyData("2023-01-01", "2023-12-31", "All")
df_2024 <- lwdataexplorer::getBuoyData("2024-01-01", "2024-12-31", "All")


  df_2018_2 <- select_vars(df_2018)
  df_2019_2 <- select_vars(df_2019)
  df_2020_2 <- select_vars(df_2020)
  df_2021_2 <- select_vars(df_2021)
  df_2022_2 <- select_vars(df_2022)
  df_2023_2 <- select_vars(df_2023)
  df_2024_2 <- select_vars(df_2024)
  
  df_pCO2A <- df_carbon[which(!is.na(df_carbon$`CO2-atm concentration (ppm)`)),]
  df_pCO2W <- df_carbon[which(!is.na(df_carbon$`CO2-seawater concentration (ppm)`)),]
  }
 
 # Calculate the carbon flux
 FC_buoy_validation <- function(u10,Temp,Sal,pCO2W,pCO2A){
   kn = ((0.222*u10^2+0.333*u10)/3600)*100 # Nightingale
   
   solub <- exp(-58.0931+(90.5069*(100/(273.15+Temp)) + (22.2940*log((273.15+Temp)/100)) +
                            Sal*(0.027766-(0.025888*((273.15+Temp)/100)) + (0.0050578*((273.15+Temp)/100)^2))))
   
   dpCO2_val <- pCO2W*10^(-6) - pCO2A*10^(-6)
   
   FC_val <- -1 * kn*solub*dpCO2_val*24*10000 #mmol/(m2*d)     -1 to switch uptake and release CO²
   
   return(FC_val)
 }

req_vars_co2 <- function(start, stop, station_selection, station_region ) {
#function to retrieve required variables (ph,sal and temp) to estimate PCO2w

# Load functions for creating GAM models
source(paste0(wd,"Rscripts/GAM/functions_GAM_v1_SP.R"))
source(paste0(wd,"Rscripts/GAM/LW_DataAccess_functions_for_GAM.R", encoding = "UTF-8"))

#devtools::install_github("lifewatch/lwdataexplorer")
  packages <- list("mgcv",
          "lubridate",
          "stats",
          "xts",
          "lwdataexplorer",
          "dplyr")
  lapply(packages, library, character.only=TRUE)

if (area == "bpns"){
#Select station data, dates and stations of interest. 
df_all <- getStationData("2012-01-01", stop,
                               stations = station_selection,
                               categories="all")
  
  # Column names to lower case
  colnames(df_all) <- tolower(colnames(df_all))

# Create a new column 'areas' and map stations to the 3 areas:
  # Select only columns of interest
  df_all <- df_all %>%
        select( `temperature(degc)`, 'salinity(psu)','ph', time )
  colnames(df_all) <- c('Temp','sal','ph','Date')

# run gam function for nutrients ------------------------------------------
  gam_input <- NULL
  for(var in colnames(df_all[,-4])) {
  df_gam <- gam_nutrients_input(df_all, var, startdate, stopdate)
  gam_input[["Date"]] <- df_gam$Date
  gam_input[[var]] <- df_gam[,var]
}
predictions <- as.data.frame(gam_input)
predictions$station <- station_code
predictions$month <- month(predictions$Date)
predictions$year <- year(predictions$Date) 

required_vars_co2_df <- subset(predictions, Date > as.Date(start) - 1 & Date < as.Date(stop) + 1)
} else {
    # run gam function for carbon ------------------------------------------
  gam_C_input <- NULL
  for(var in colnames(carbon_input[,-4])) {
  df_C_gam <- gam_nutrients_input(carbon_input, var, start, stop)
  gam_C_input[["Date"]] <- df_C_gam$Date
  gam_C_input[[var]] <- df_C_gam[,var]
    }
predictions <- as.data.frame(gam_C_input)
predictions$station <- "C1"
predictions$month <- month(predictions$Date)
predictions$year <- year(predictions$Date)
required_vars_co2_df <- subset(predictions, Date > as.Date(start) - 1 & Date < as.Date(stop) + 1)
    
}


pco2atm_wind <- read.csv(paste0(wd,"Input data/",pco2atm_wind_file))
if (pco2atm_wind$year[length(pco2atm_wind$year)] != year(stop) ){
    missing_years <- year(stop) - pco2atm_wind$year[length(pco2atm_wind$year)] 
    for ( i in 1:(missing_years) ){
        adding_years <- pco2atm_wind[which(pco2atm_wind$year == pco2atm_wind$year[length(pco2atm_wind$year)]),]
        adding_years$year <- unique(tail(pco2atm_wind$year)) + 1
        pco2atm_wind <- rbind(pco2atm_wind, adding_years)
        
        }
}

required_vars_co2_df$month <- month(required_vars_co2_df$Date)

# Combine the data frames based on "year" and "month"
required_vars_co2_df <- required_vars_co2_df %>%
  inner_join(pco2atm_wind, by = c("year", "month"))

    
return(required_vars_co2_df)
}
