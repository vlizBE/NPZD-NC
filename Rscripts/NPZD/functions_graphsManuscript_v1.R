library(ggplot2)
library(xts)
library(reshape2)
library(lubridate)
library(stats)
library(dplyr)
library(viridis)
theme_set(theme_bw())

getBestSim <- function(id_station){
  
  #Read the file with the errors calculated for all simulations
  
  error_log <- read.csv(paste0("~/workspace/NPZ/Output/Final results/NPZ/station", id_station,
                               "_iter2/errors_parameters_iter2.csv"))
  error_log <- error_log[which(!is.na(error_log$total_rmse)),]
  
  #order simulations from lowest to highest error
  order_error <- error_log[order(error_log$total_rmse),]
  
  #select best simulations (10% or 5%)
  percentage <- 0.10
  numSim <- trunc(percentage * nrow(order_error))
  selectedSim <- order_error[1:numSim,]
  
  #This function extracts the details of the best 10% best simulations
  simulations_best_rmse <- bestSimulations(selectedSim$ID_sim, id_station)
  simulations_best_rmse <- simulations_best_rmse[,c(2,3,4,5,6)]
  
  #This function calculates the mean value for each month for chlorophyll-a abundances
  bestSim_month_chla <- reduceSimulations(simulations_best_rmse[,c(3,4,5)])
  bestSim_month_chla$Station <- id_station
  
  #This function calculates the mean value for each month for zooplankton abundances
  bestSim_month_zoo <- reduceSimulations(simulations_best_rmse[,c(2,4,5)])
  bestSim_month_zoo$Station <- id_station
  
  #Two elements in the list, the first the best Chlorophyll-a values, and the second the zooplankton abundances
  list_bestSim <- list(bestSim_month_chla, bestSim_month_zoo)
  
  return(list_bestSim)
  
}

reduceSimulations <- function(all_sim){
  
  #Create a new date to refer to the first day of the month as data will be averaged per month
  all_sim$year <- year(all_sim$Date)
  all_sim$roundDate <- floor_date(all_sim$Date, unit = "month")
  
  #Calculate average values per month mantaining the different ID of each simulation
  simulations_by_month <- aggregate( .~ID+roundDate, data = all_sim, mean )
  
  #Remove the original column date that after all transformations does not contain relevant information anymore
  simulations_by_month <- subset(simulations_by_month, select = -c(Date))
  
  return(simulations_by_month)
  
}

bestSimulations <- function(bestID, station_int){
  #Extract the number of the 10% best simualtions
  count_ID <- length(bestID)
  
  #Read the input data for the correct date for each day
  readDate  <- read.csv(paste0("~/workspace/NPZ/Input data/NPZ/inputData_npzd_station", station_int, ".csv"))
  
  bestSim<-NA
  
  #Extract the details of each simulation 
  
  for(i in 1:count_ID){
    idSim <- bestID[i]
    resultsSim <-read.csv(paste0("~/workspace/NPZ/Output/Final results/NPZ/station", station_int,
                                 "_iter2/detailed_simulation_", idSim, ".csv"))

    #Select only the columns of interest (final results, no limitation factors)
    temp_sim <-resultsSim[,c(2,12,13,17)]
    #Attach the id simulation to differentiate the 10% best simulations from each other to calculate the variance
    temp_sim$ID <- idSim
    temp_sim$Date <- as.Date(readDate$Date, format = "%m/%d/%Y")
    #Remove dummy years
    temp_sim<- temp_sim[-c(1:1096),]
    #Bind all results together to be able to plot the variance of the 10% best simulations
    bestSim <- rbind(bestSim,temp_sim)
  }
  
  bestSim <- bestSim[-c(1),]
  bestSim$Date <- as.Date(bestSim$Date)
  bestSim$month <- month(bestSim$Date)
  #Add that these are the simulated data, to compare in the graph afterwards with observed data
  bestSim$type <- "simulated"
  
  return(bestSim)
}
