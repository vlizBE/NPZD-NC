get_needed_information <- function(){
  area <<- readline(prompt = "What is the area you want to explore? Belgian part of the North Sea, 
                  Northern Adriatic Sea or another place? \n possible answers: bpns, nas, other ")
  
  if (area == "bpns"){
    # Load LifeWatch (LW) data with LW package
    station_list <- c("LW01", "LW02", "435", "W07bis", "W08", "W09", "W10", "421", "130", "700",
                        "780", "330", "230", "710", "ZG02", "120", "215")
    station_coord <- data.frame("lat" = c(2.256, 2.556, 2.7903, 3.0125, 2.35, 2.7, 2.4167, 
                                      2.45, 2.9054, 3.221, 3.0573, 2.8091, 2.8504, 3.1383, 2.5007, 2.7025, 2.6108),
                            "long" = c(51.5687 , 51.8, 51.5807, 51.588, 51.4583, 51.75,
                                       51.6833, 51.4805, 51.2706, 51.377, 51.4714, 51.4341, 51.3087, 51.4412, 51.3352, 51.1861, 51.2749))
    station_region <<- readline("What is your region of interest? nearshore, midshore or offshore?  ")
    
    if (station_region != "offshore") {
        station_code <<- readline("What is the code of your station? For example 130 ")
        if (station_code %in% station_list){
               station_lat <<- station_coord$lat[which(station_list %in% station_code)]       # coordinates station: latitude 
               station_lon <<- station_coord$long[which(station_list %in% station_code)]       # coordinates station: longitude
     } else {
        print("Station was not found? Please check.")
        }
        } else {
        station_code <<- c("LW01", "LW02", "435", "W07bis", "W08", "W09", "W10", "421")
        if (station_code[1] %in% station_list){
               station_lat <<- station_coord$lat[which(station_list %in% station_code[1])]       # coordinates station: latitude 
               station_lon <<- station_coord$long[which(station_list %in% station_code[1])]       # coordinates station: longitude
         } else {
        print("Station was not found? Please check.")
        }
        }
    
    # period of your time series
    startdate <<- as.Date(readline("What is first date of your data? format YYYY-MM-DD  "))
    stopdate <<- as.Date(readline("What is last date of your data? format YYYY-MM-DD  "))
  }
    
  if (area == "nas"){
    station_region <<- readline("What is your region of interest? nearshore, midshore or offshore?  ")
    station_code <<- "C1"
    startdate <<- as.Date(readline("What is first date of your data? format YYYY-MM-DD  "))
    stopdate <<- as.Date(readline("What is last date of your data? format YYYY-MM-DD  "))
    station_lat <<- 44.123 # coordinates station: latitude 
    station_lon <<- 12.57
    input_filename <<- readline("Please provide the name of your file containing the raw input (downloaded from EMODnet). For example Eutrophication_Med_timeseries_2024_unrestricted.txt  ")
      
  }
    
  if (area == "other"){
    station_region <<- readline("What is your region of interest? nearshore, midshore or offshore?  ")
    station_code <<- readline("What is the code of your station? For example 130 ")
    startdate <<- as.Date(readline("What is first date of your data? format YYYY-MM-DD  "))
    stopdate <<- as.Date(readline("What is last date of your data? format YYYY-MM-DD  "))
    station_lat <<- as.numeric(readline("Please provide coordinates for your station: Latitude:  ")) # coordinates station: latitude 
    station_lon <<- as.numeric(readline("Please provide coordinates for your station: Longitude:  ")) # coordinates station: longitude
    input_filename <<- readline("Please provide the name of your file containing the raw input (downloaded from EMODnet). For example Eutrophication_Med_timeseries_2024_unrestricted.txt  ")
  }

  # Nutrients sampling depth
  if (area == 'bpns') {
    sample_depth <<- 3 
        } else {
        sample_depth <<-  as.numeric(readline("At what depth were the nutrients sampled? e.g. 5m depth"))
    }

  # pco2atm and wind data for input
    pco2atm_wind_file <<- readline("What is the name of your file containing pCO2atmosphere and wind data? For example 'bpns_pco2atm_windspeed.csv'")
  # carbon data for validation
  carbon_validation_available <<- readline("Do you have pCO2seawater data for validation? yes or no")

  if (carbon_validation_available == "yes") {
      carbon_validation_folder <<- readline("What is the name of your folder containing pCO2seawater data? For example 'pco2w_validation_bpns'")
    } else {
        carbon_validation_folder <<- NA
        pco2w_validation[1,1 ] <<-  data.frame("Date" = NA, "pco2w" = NA)
    }
}

number_of_iterations <- function(){
  numSimulations <<- as.numeric(readline(prompt = "How many iterations do you want to run?"))
}

gam_nutrients_input <- function(df,var,start,stop){
  df <- df[which(!is.na(df[,var])),]
  
  #omit outliers
  lower_bound <- 0
  upper_bound <- median(df[,var]) + 3 * mad(df[,var], constant = 1)
  outlier_ind <- which(df[,var] < lower_bound | df[,var] > upper_bound)
  df <- df[!c(1:nrow(df)) %in% outlier_ind,]
  
  df_day <- yday(df$Date)
  df_year <- year(df$Date)
  
  models <- list()

  if(length(unique(df_year)) < 3){
    print(paste("Time series of",var,"is not long enough for GAM modelling"))
  }
  if(length(unique(df_year)) == 3){
    models[[1]] <- test_model(unlist(df[var]), df_day, df_year, 3, 3)
    mod_df_select <-models[[1]]
    } 
  if(length(unique(df_year)) == 4){
    models[[1]] <- test_model(unlist(df[var]), df_day, df_year, 3, 3)
    models[[2]] <- test_model(unlist(df[var]), df_day, df_year, 4, 4)
    AICs <- c(models[[1]]$aic,
              models[[2]]$aic)
    mod_df_select <- models[[which.min(AICs)]]
  }
  
  if(length(unique(df_year)) >=5){
    models[[1]] <- test_model(unlist(df[var]), df_day, df_year, 3, 3)
    models[[2]] <- test_model(unlist(df[var]), df_day, df_year, 4, 4)
    models[[3]] <- test_model(unlist(df[var]), df_day, df_year, 5, 5)
    AICs <- c(models[[1]]$aic,
              models[[2]]$aic,
              models[[3]]$aic)
    
    mod_df_select <- models[[which.min(AICs)]]
  }
  

  predictions <- data.frame(matrix(nrow = as.numeric(difftime(as.Date(stop) , as.Date(start))) + 1, ncol = 4))
  
  names(predictions)[1] <- "day"
  names(predictions)[2] <- "year"
  names(predictions)[3] <- var
  names(predictions)[4] <- "Date"
  
  # Generate a sequence of dates between the two times
  predictions$Date <- seq(from = as.Date(start), to = as.Date(stop), by = "day")
  predictions$day <- yday(predictions$Date)
  predictions$year <- year(predictions$Date)
  
  # predict nutrient values
  predictions[3]<-predict.gam(mod_df_select, newdata = predictions[,c(1,2)], type = "response") 
  
  predictions[which(predictions[3] < 0),3] <- 0.001
  
  
  return(predictions)
  }

get_inputdata <- function(station_region, station_selection, startdate, stopdate){
  ## GAM models
 if (area == "bpns"){
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


  
  
  days_run = as.numeric(difftime(as.Date(stopdate) , as.Date(startdate))) + 1 # how many days are in the run = diff (start - stop) + 1
 
 
  
        
station_data <- get_abiotic_pigment(start = startdate, stop = stopdate, station_selection = station_selection, station_region = station_region)
station_data <- station_data %>% 
     select(c("NO2", "PO4", "SiO4", "Temperature", "NH4", "NO3", "DIN", "Date"))
# run gam function for nutrients ------------------------------------------
  gam_input <- NULL
  for(var in colnames(station_data[,-8])) {
  df_gam <- gam_nutrients_input(station_data, var, startdate, stopdate)
  gam_input[["Date"]] <- df_gam$Date
  gam_input[[var]] <- df_gam[,var]
}
predictions <- as.data.frame(gam_input)
predictions$station <- station_code
colnames(predictions) <- c("Date","no2", "po4", "sio4", "Temp", "nh4", "no3", "DIN", "station")
# adjust negative DIN values to small concentration.
if (length(which(predictions$DIN < 0 )) != 0) {
neg_values_DIN <- which(predictions$DIN < 0 )

min_values <- station_data$DIN[order(station_data$DIN, decreasing=F)]
min_values <- na.omit(min_values)[1:20]
min_DIN <- min(min_values)
max_DIN <- max(min_values)

predictions$DIN[neg_values_DIN] <- runif(length(neg_values_DIN), min_DIN, max_DIN)
}

# adjust negative values to small concentration.
if (length(which(predictions$po4 < 0 )) != 0) {
  neg_values_po4 <- which(predictions$po4 < 0 )
  
  min_values <- station_data$PO4[order(station_data$PO4, decreasing=F)]
  min_values <- na.omit(min_values)[1:20]
  min_po4 <- min(min_values)
  max_po4 <- max(min_values)
  
  predictions$po4[neg_values_po4] <- runif(length(neg_values_po4), min_po4, max_po4)
}

# adjust negative values to small concentration.
if (length(which(predictions$sio4 < 0 )) != 0) {
  neg_values_sio4 <- which(predictions$sio4 < 0 )
  
  min_values <- station_data$SiO4[order(station_data$SiO4, decreasing=F)]
  min_values <- na.omit(min_values)[1:20]
  min_sio4 <- min(min_values)
  max_sio4 <- max(min_values)
  
  predictions$sio4[neg_values_sio4] <- runif(length(neg_values_sio4), min_sio4, max_sio4)
}

     
} else {
lines <- readLines(paste0(wd,"Input data/",input_filename))
skiprows <- grep("^Cruise", lines)[1] - 1  # Find the first line starting with "ID"
raw_input <- read.table(paste0(wd,"Input data/",input_filename),skip = skiprows, header = TRUE, sep = '\t')

raw_input_vars <- raw_input %>%
  select(c('Station.name',
           'yyyy.mm.ddThh.mm.ss.sss',
           'Longitude..degrees_east.',
           'Latitude..degrees_north.',
           'Depth..m.',
           'ITS.90.water.temperature..degrees.C.',
           'Water.body.salinity..per.mille.',
           'Water.body.pH..pH.units.',
           'Water.body.phosphate..umol.l.',
           'Water.body.silicate..umol.l.',
           'Water.body.dissolved.inorganic.nitrogen..DIN...umol.l.',
           'Water.body.chlorophyll.a..mg.m.3.'))
raw_input_vars$Date <- as.Date(substr(raw_input_vars$yyyy.mm.ddThh.mm.ss.sss,1,10))

input_data<- raw_input_vars %>%
  select(c('Water.body.dissolved.inorganic.nitrogen..DIN...umol.l.', 
           'ITS.90.water.temperature..degrees.C.','Water.body.phosphate..umol.l.',
           'Water.body.silicate..umol.l.', Date
           )) 
  colnames(input_data) <- c("DIN", "Temp", "po4", "sio4", "Date")
  input_data <- subset(input_data, Date > "2010-12-31")

     #prepare already a data set for carbon function
carbon_input <<-raw_input_vars %>%
  select(c('ITS.90.water.temperature..degrees.C.',
           'Water.body.salinity..per.mille.',
           'Water.body.pH..pH.units.', Date
           ))
    colnames(carbon_input) <<- c("temp", "sal", "ph", "Date")
    carbon_input <<- subset(carbon_input, Date > "2010-12-31")

# prepare chla validation data
validation_data <-raw_input_vars %>%
  select(c(Date,
           'Water.body.chlorophyll.a..mg.m.3.'
           ))
    colnames(validation_data) <- c("Date", "chla")
    validation_data <- subset(validation_data, Date > "2010-12-31")
    validation_data <- subset(validation_data, Date > startdate & Date < stopdate)
    validation_data <<- na.omit(validation_data)

     
# run gam function for nutrients ------------------------------------------
  gam_input <- NULL
  for(var in colnames(input_data[,-5])) {
  df_gam <- gam_nutrients_input(input_data, var, startdate, stopdate)
  gam_input[["Date"]] <- df_gam$Date
  gam_input[[var]] <- df_gam[,var]
}
predictions <- as.data.frame(gam_input)
predictions$station <- station_code
     
}
#Save predictions as csv
write.csv(predictions, file = paste0(wd,"Input data/inputData_npzd_station", station_region[1],".csv"),row.names = F)

return(predictions)
}
