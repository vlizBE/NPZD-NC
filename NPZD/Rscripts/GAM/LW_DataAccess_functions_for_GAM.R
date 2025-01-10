#################################################################
# title: "Preparing input data GAM and NPZ model"               #
# author: "Patricia Cabrera, Salvador Fernandez, Steven Pint"   #
# date: "06/12/2021"                                            #
#################################################################

#packages  
library(lwdataexplorer)
library(plyr)
library(dplyr)
library(lubridate)

##Stations
#Select station data, dates: 2014-2017, and stations of interest. 

get_abiotic_pigment <- function(start, stop, station_selection, station_region){
  
#### Code used to get the data from the Life Watch Portal ####

################## ABIOTIC ##################
 
  #Select station data, dates and stations of interest. 
  df_all <- getStationData(start, stop,
                               stations = station_selection,
                               categories="all")
  
  # Column names to lower case
  colnames(df_all) <- tolower(colnames(df_all))
  
  # Create a new column 'areas' and map stations to the 3 areas:
  # Calculate Dissolved Inorganic Nitrogen (DIN) values based on NH4 + NO3 + NO2 and add it to a new column
  # Select only columns of interest
  df_all <- df_all %>%
    mutate(
  areas = recode(station, !!!setNames(rep(station_region,length(station_selection)), station_selection)),
      din = `ammonium_nh4(µmol n_nh4/l)` + `nitrate_no3(µmol n_no3/l)` + `nitrite_no2(µmol n_no2/l)`
    ) %>% 
    select(
      station, time, latitude, longitude, `ammonium_nh4(µmol n_nh4/l)`, `nitrate_nitrite(µmol n_no3-no2/l)`, `nitrate_no3(µmol n_no3/l)`,
      `nitrite_no2(µmol n_no2/l)`, `phosphate_po4(µmol p_po4/l)`, `silicate_sio4(µmol si_sio4/l)`,`temperature(degc)`, areas, din
    )
  
################## PIGMENTS ##################
  #Get Chla data
  df_pigments <- getStationData(start, stop,
                               categories = "Pigments",
                               stations = station_selection)

  # Column names to lower case
  colnames(df_pigments) <- tolower(colnames(df_pigments))
  
  # Select a few columns
  df_pigments <- df_pigments %>%
    select(station, time, `chlorophyll_a(µg/l)`)
  
  #joint the 2 df based on Date and fill with NA ?????
  abiotic_pigment <- left_join(df_all, df_pigments, by = c("station", "time"))

  #select data needed for GAM models
  df_GAM <- abiotic_pigment %>%
    select(station, time, `nitrite_no2(µmol n_no2/l)`, `phosphate_po4(µmol p_po4/l)`, `silicate_sio4(µmol si_sio4/l)`, 
           `temperature(degc)`, `ammonium_nh4(µmol n_nh4/l)`, `nitrate_no3(µmol n_no3/l)`, areas, `chlorophyll_a(µg/l)`, din)

  #rename columns df_GAM
  colnames(df_GAM) <- c("Station", "Date", "NO2", "PO4", "SiO4", "Temperature", "NH4", "NO3", "Region", "Chla", "DIN")

  #select data needed for NPZ model
  abiotic_pigment_2020 <- abiotic_pigment %>%
    select(station, time, latitude,	longitude, `ammonium_nh4(µmol n_nh4/l)`, `chlorophyll_a(µg/l)`, `nitrate_nitrite(µmol n_no3-no2/l)`, 
           `nitrate_no3(µmol n_no3/l)`, `nitrite_no2(µmol n_no2/l)`, `phosphate_po4(µmol p_po4/l)`, `silicate_sio4(µmol si_sio4/l)`, 
           `temperature(degc)`)
  
  #rename columns df_GAM
  colnames(abiotic_pigment_2020) <- c("Station",	"Time",	"Latitude",	"Longitude",	"Ammonium_NH4",	"Chlorophyll_a",
                                      "Nitrate_Nitrite",	"Nitrate_NO3", "Nitrite_NO2",	"Phosphate_PO4",	
                                      "Silicate_SiO4",	"Temperature")
  
  # write .csv file for NPZ model
  write.csv(abiotic_pigment_2020, file = paste0(wd,"Input data/abiotic_pigment_2020.csv"),row.names = F)

return(df_GAM)
}

#########################################################################################################################
get_zooplankton <- function(start, stop, taxa_select, station_selection){

  #Get zooplankton data
  df_zoo <- getZooscanData(start, stop)

  # Columns to lower case
  colnames(df_zoo) <- tolower(colnames(df_zoo))

  #We are only interested in the most abundant taxa in the BPNS: Calanoida, Noctiluca, Harpacticoida and Appendicularia
  df_zoo <- df_zoo %>%
    filter(taxon %in% taxa_select) %>%
    # select our station(s) of interest
    subset(station %in% station_selection) %>%
    # select the columns that are needed
    select(station, longitude, latitude, tripaction, time, taxon, density)

  # change date format
  df_zoo$time <- strftime(df_zoo$time, format = "%Y-%m-%d")
  
  # add columns with month, year, source and two columns in front (NA's)
  df_zoo$month <- month(df_zoo$time)
  df_zoo$year <- year(df_zoo$time)
  df_zoo$source <- "Zooscan"
  df_zoo$empty <- NA
  df_zoo$x <- NA
  
  df_zoo <- df_zoo %>%
    # select the columns that are needed
    select(empty, x, station, longitude, latitude, tripaction, time, taxon, density, month, year, source)
  
    colnames(df_zoo) <- c("x", "X", "Station",	"Longitude",	"Latitude",	"Tripaction",	"Date",
                         "Taxon",	"Density",	"Month",	"Year",	"Source")
  
  
  # write .csv file for NPZ model
  write.csv(df_zoo, file = paste0(wd,"Input data/zooscan.csv"), row.names = F)

return(df_zoo)
}

# Sea surface temperature MVB
get_MVB_SST <- function(start, stop){

  # use for loop to fill datagaps with SST from other station than 'Westhinder - Measuring pile'
  station_list <- c('Westhinder - Measuring pile', "Thorntonbank South","Westhinder - Buoy")
  MVB <- NA
  
  for (i in station_list){
    #### Code used to get the MVB data from the Life Watch Portal ####
    MVB1 <- getMvbData(
      start,
      stop,
      "Sea water temperature",
      stations = i,
      "day",
      "avg",
      usr = readline("give me a user!!: "), 
      pwd = readline("give me a password!!: "),
      params = FALSE
    )
  
    colnames(MVB1) <- tolower(colnames(MVB1))
    # select date and SST
    MVB1 <- MVB1 %>%
      select(time, `sea water temperature`)
    # Rename and date format
    colnames(MVB1) <- c("Date", "SST")
    MVB1$Date <- as.Date(MVB1$Date)
    # if Westhinder just copy df, otherwise rbind
    if (i == 'Westhinder - Measuring pile'){
      MVB <- MVB1
    } else {
      MVB <- rbind(MVB, MVB1[which(MVB1$Date %in% datagaps),])
    }
  
    # check which date are still missing data
    date_range <- seq(min(MVB$Date), max(MVB$Date), by = 1) 
    datagaps <- date_range[!date_range %in% MVB$Date] 
    # set new start and stop date to extract data
    start <- min(datagaps)
    stop <- max(datagaps)
  }

  # get extra data from buoy at Thornton bank to fill data gaps
  buoy <- getBuoyData(start, stop, "Buoy at C-Power", params = FALSE)
  colnames(buoy) <- tolower(colnames(buoy))
  # select date and SST
  buoy <- buoy %>%
    select(time, `water temperature (°c)`)
  # Rename and date format
  colnames(buoy) <- c("Date", "SST")
  buoy$Date <- as.Date(buoy$Date)
  # aggregate by day
  buoy <- aggregate(SST ~ Date, buoy, mean)
  # add new data to MVB df
  MVB <- rbind(MVB, buoy[which(buoy$Date %in% datagaps),])
  # check for data gaps left
  datagaps <- date_range[!date_range %in% MVB$Date]
  # fill last data gaps by copying SST from previous day
  copy <- data.frame("Date" = NA, "SST" = NA)
  for (j in c(1:length(datagaps))){
    copy$Date <- datagaps[j]
    copy$SST <- MVB$SST[which(MVB$Date == (datagaps[j]-1))]
    MVB <- rbind(MVB, copy)
  }
  
  # Last data gap check
  datagaps <- date_range[!date_range %in% MVB$Date]
  if (length(datagaps) == 0) {
    MVB <- arrange(MVB, Date)
    return(MVB)
  } else {
    print("There is a data gap!")
  }
  
}
    
