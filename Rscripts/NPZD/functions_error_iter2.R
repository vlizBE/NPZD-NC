library(ggplot2)
library(xts)
library(reshape2)
library(lubridate)
library(stats)
library(dplyr)

calculateErrors <- function(num_station, chla, zoo_aggregated, num_sim){
  
  cont <- 1
  possible_parameters_1 <- read.csv(paste0("~/workspace/NPZ/Output/Final results/NPZ/station", num_station, "_iter2/possible_parameters_spring_iter2.csv"))
  possible_parameters_1 <- possible_parameters_1[,-c(1)]
  
  possible_parameters_2 <- read.csv(paste0("~/workspace/NPZ/Output/Final results/NPZ/station", num_station, "_iter2/possible_parameters_summer_iter2.csv"))
  possible_parameters_2 <- possible_parameters_2[,-c(1)]
  
  errors_parameters <- createErrorFile(num_sim)
  
  readDate <- read.csv(paste0("~/workspace/NPZ/Input data/NPZ/inputData_npzd_station", num_station, ".csv"))
  
  
  if(num_station == "offshore"){
    zoo_aggregated <- aggregate(zoo_aggregated$n_mol, by=list(zoo_aggregated$Date), FUN="median",  data = zoo_aggregated)  
    chla <- aggregate(chla$Chlorophyll_a, by=list(chla$Time), FUN="median",  data = chla)  
  }
  
  for(i in 1:num_sim){
    
    output <- read.csv(paste0("~/workspace/NPZ/Output/Final results/NPZ/station", num_station, "_iter2/detailed_simulation_", i, ".csv"))
    output$Date <- as.Date(readDate$Date, format = "%m/%d/%Y")
    
    #remove dummy years (one leap year: 2012)
    output <- output[-c(1:1096),] 
    
    if(!anyNA(output)){
      
      ################
      #calculate error
      print(paste("Num sim:", i))
      
      if(num_station == "offshore"){
        temp_error<-calculateRMSE(output, possible_parameters_1[i,12],  possible_parameters_2[i,12], zoo_aggregated, chla)
        
      }else{
        temp_error<-calculateRMSE(output, possible_parameters_1[i,12], possible_parameters_2[i,12], zoo_aggregated[,c(2,3)], chla[,c(2,3)])
        
      }
      
      errors_parameters$ID_sim[cont] <- i
      errors_parameters$rmse_zoo[cont] <- temp_error[[1]]
      errors_parameters$rmse_phyto[cont] <- temp_error[[2]]
      errors_parameters$total_rmse[cont] <- temp_error[[1]] + temp_error[[2]]
      
      #####
      #attach set of parameters
      
      errors_parameters[cont,c(5:18)] <-possible_parameters_1[i,]
      errors_parameters[cont,c(19:32)] <-possible_parameters_2[i,]
      
      cont <- cont + 1
      
    }
  }
  
  write.csv(errors_parameters, paste0("~/workspace/NPZ/Output/Final results/NPZ/station", num_station,"_iter2/errors_parameters_iter2.csv"))
  
  return(errors_parameters)
  
}


createErrorFile <- function(numSim){
  errors_parameters <- data.frame(matrix(nrow = numSim, ncol = 32))
  names(errors_parameters)[1] <- "ID_sim"
  names(errors_parameters)[2] <- "rmse_phyto"
  names(errors_parameters)[3] <- "rmse_zoo"
  names(errors_parameters)[4] <- "total_rmse"
  
  names(errors_parameters)[5] <- "ID_1"
  names(errors_parameters)[6] <- "maxUptake_1"
  names(errors_parameters)[7] <- "ksPAR_1"
  names(errors_parameters)[8] <- "ksDIN_1"
  names(errors_parameters)[9] <- "ksP_1"
  names(errors_parameters)[10] <- "ksSi_1"
  names(errors_parameters)[11] <- "maxGrazing_1"
  names(errors_parameters)[12] <- "ksGrazing_1"
  names(errors_parameters)[13] <- "pFaeces_1"
  names(errors_parameters)[14] <- "excretionRate_1"
  names(errors_parameters)[15] <- "mortalityRate_1"
  names(errors_parameters)[16] <- "ChlNratio_1"
  names(errors_parameters)[17] <- "Tobs_1"
  names(errors_parameters)[18] <- "kd_1"
  
  names(errors_parameters)[19] <- "ID_2"
  names(errors_parameters)[20] <- "maxUptake_2"
  names(errors_parameters)[21] <- "ksPAR_2"
  names(errors_parameters)[22] <- "ksDIN_2"
  names(errors_parameters)[23] <- "ksP_2"
  names(errors_parameters)[24] <- "ksSi_2"
  names(errors_parameters)[25] <- "maxGrazing_2"
  names(errors_parameters)[26] <- "ksGrazing_2"
  names(errors_parameters)[27] <- "pFaeces_2"
  names(errors_parameters)[28] <- "excretionRate_2"
  names(errors_parameters)[29] <- "mortalityRate_2"
  names(errors_parameters)[30] <- "ChlNratio_2"
  names(errors_parameters)[31] <- "Tobs_2"
  names(errors_parameters)[32] <- "kd_2"
  
  return(errors_parameters)
}

calculateRMSE <- function(resultsSimulation1, chla_conversion1, chla_conversion2, zoo_obs, chla_obs){
  
  #get only zoo, phyto and chl_a results
  resultsSimulation1 <- resultsSimulation1[,c(12,13,17,19)]
  
  ######################
  #Calculate RMSE zooplankton
  names(zoo_obs)[1]<-"Date_obs"
  names(zoo_obs)[2]<-"zoo_obs"
  
  zoo_obs$Date_obs <- as.Date(zoo_obs$Date_obs)
  
  temp_zoo <- merge(resultsSimulation1, zoo_obs , by.x = "Date", by.y= "Date_obs" ,  all.x = FALSE)
  
  rmse_sim_zoo <- sqrt(mean((temp_zoo$zoo - temp_zoo$zoo_obs)^2))
  
  ######################
  #Calculate error Chl_a in units of mmolN
  names(chla_obs)[1]<-"Date_obs"
  names(chla_obs)[2]<-"chla_obs"
  
  temp_chla <- merge(resultsSimulation1, chla_obs , by.x = "Date", by.y= "Date_obs" ,  all.x = FALSE)
  
  temp_chla$doy <- yday(temp_chla$Date)
  
  temp_chla$phyto_obs <- ifelse( temp_chla$doy > 172 , temp_chla$chla_obs/chla_conversion2 , temp_chla$chla_obs/chla_conversion1)    
  
  rmse_sim_phyto <- sqrt(mean((temp_chla$phyto - temp_chla$phyto_obs)^2))

  list_results <- list(rmse_sim_zoo, rmse_sim_phyto)
  print(list_results)
  
  return(list_results)
  
}


createBestSim <- function(ID_minError, station_int, chla_obs, zoo_obs){
  
  countID <- length(ID_minError)
  
  readDate  <- read.csv(paste0("~/workspace/NPZ/Input data/NPZ/inputData_npzd_station", station_int, ".csv"))
  
  bestSim<-NA
  
  for(i in 1:countID){
    idSim <- ID_minError[i]
    print(idSim)
    resultsSim <-read.csv(paste0("~/workspace/NPZ/Output/Final results/NPZ/station", station_int, "_iter2/detailed_simulation_", idSim, ".csv"))
    temp_sim <-resultsSim[,c(2,12,13,17)]
    temp_sim$ID <- idSim
    temp_sim$Date <- as.Date(readDate$Date, format = "%m/%d/%Y")
    
    #Remove dummy years
    temp_sim<- temp_sim[-c(1:1096),]
    bestSim <- rbind(bestSim,temp_sim)
  }
  
  bestSim <- bestSim[-c(1),]
  bestSim$Date <- as.Date(bestSim$Date, format = "%Y-%m-%d")
  
  
  #####Chl_a
  bestSim_mean_chla <- aggregate( chl_a ~ Date, data = bestSim, mean )
  
  figureBest_chla <- ggplot() + 
    geom_line(data = bestSim, aes(x=as.Date(Date), y = chl_a),alpha = 0.1) +
    geom_line(data = bestSim_mean_chla, aes(x=as.Date(Date), y = chl_a),colour = "blue", size = 1) +
    geom_point(data = chla_obs, aes(x = as.Date(Time), y= Chlorophyll_a), colour = "red")+
    labs(x = "Time", y = "Chl_a (mg Chl_a/m3)") + 
    scale_x_date(date_breaks = "6 month", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=11), axis.text.x=element_text(angle=90, hjust=1))
  
  ##### Zoo
  bestSim_mean_zoo <- aggregate( zoo ~ Date, data = bestSim, mean )
  
  figureBest_zoo <- ggplot() + 
    geom_line(data = bestSim, aes(x=as.Date(Date), y = zoo),alpha = 0.1) +
    geom_line(data = bestSim_mean_zoo, aes(x=as.Date(Date), y = zoo),colour = "blue", size = 1) +
    geom_point(data = zoo_obs, aes(x = as.Date(Date), y= n_mol), colour = "red")+
    labs(x = "Time", y = "Zooplankton (mmol N/m3)") + 
    scale_x_date(date_breaks = "6 month", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=11), axis.text.x=element_text(angle=90, hjust=1))
  
  
  #####Combined
  bestSim_mean <- aggregate( phyto ~ Date, data = bestSim, mean )
  
  figureBest_combined <- ggplot() + 
    geom_line(data = bestSim, aes(x=as.Date(Date), y = phyto),alpha = 0.1, colour = "blue") +
    geom_line(data = bestSim_mean, aes(x=as.Date(Date), y = phyto),colour = "blue", size = 1) +
    geom_line(data = bestSim, aes(x=as.Date(Date), y = zoo),alpha = 0.1) +
    geom_line(data = bestSim_mean_zoo, aes(x=as.Date(Date), y = zoo),colour = "black", size = 1) +
    labs(x = "Time", y = "Plankton (mmol N/m3)") + 
    scale_x_date(date_breaks = "6 month", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=11), axis.text.x=element_text(angle=90, hjust=1))
  
  ############################
  ###logarithmic scale
  
  #####Chl_a
  figureBest_chla_log <- ggplot() + 
    geom_line(data = bestSim, aes(x=as.Date(Date), y = log(chl_a)),alpha = 0.1) +
    geom_line(data = bestSim_mean_chla, aes(x=as.Date(Date), y = log(chl_a)),colour = "blue", size = 1) +
    geom_point(data = chla_obs, aes(x = as.Date(Time), y= log(Chlorophyll_a)), colour = "red")+
    labs(x = "Time", y = "Chl_a (log mg Chl_a/m3)") + 
    scale_x_date(date_breaks = "6 month", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=11), axis.text.x=element_text(angle=90, hjust=1))
  
  ##### Zoo
  figureBest_zoo_log <- ggplot() + 
    geom_line(data = bestSim, aes(x=as.Date(Date), y = log(zoo)),alpha = 0.1) +
    geom_line(data = bestSim_mean_zoo, aes(x=as.Date(Date), y = log(zoo)),colour = "blue", size = 1) +
    geom_point(data = zoo_obs, aes(x = as.Date(Date), y= log(n_mol)), colour = "red")+
    labs(x = "Time", y = "Zooplankton (log mmol N/m3)") + 
    scale_x_date(date_breaks = "6 month", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=11), axis.text.x=element_text(angle=90, hjust=1))
  
  
  #####Combined
  figureBest_combined_log <- ggplot() + 
    geom_line(data = bestSim, aes(x=as.Date(Date), y = log(phyto)),alpha = 0.1, colour = "blue") +
    geom_line(data = bestSim_mean, aes(x=as.Date(Date), y = log(phyto)),colour = "blue", size = 1) +
    geom_line(data = bestSim, aes(x=as.Date(Date), y = log(zoo)),alpha = 0.1) +
    geom_line(data = bestSim_mean_zoo, aes(x=as.Date(Date), y = log(zoo)),colour = "black", size = 1) +
    labs(x = "Time", y = "Plankton (log mmol N/m3)") + 
    scale_x_date(date_breaks = "6 month", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=11), axis.text.x=element_text(angle=90, hjust=1))
  
  
  listgraphs <- list(figureBest_chla, figureBest_zoo, figureBest_combined,
                     figureBest_chla_log, figureBest_zoo_log, figureBest_combined_log) 
  
  return(listgraphs)
}

bestSimulations <- function(bestID, station_int){
  count_ID <- length(bestID)
  
  readDate  <- read.csv(paste0("~/workspace/NPZ/Input data/NPZ/inputData_npzd_station", station_int, ".csv"))
  
  bestSim<-NA
  
  for(i in 1:count_ID){
    idSim <- bestID[i]
    resultsSim <-read.csv(paste0("~/workspace/NPZ/Output/Final results/NPZ/station", station_int, "_iter2/detailed_simulation_", idSim, ".csv"))
    temp_sim <-resultsSim[,c(2,12,13,17)]
    temp_sim$ID <- idSim
    temp_sim$Date <- as.Date(readDate$Date, format = "%m/%d/%Y")
    #Remove dummy years
    temp_sim<- temp_sim[-c(1:1096),]
    bestSim <- rbind(bestSim,temp_sim)
  }
  
  bestSim <- bestSim[-c(1),]
  bestSim$Date <- as.Date(bestSim$Date)
  bestSim$month <- month(bestSim$Date)
  bestSim$type <- "simulated"
  
  return(bestSim)
}

monthlyBestSimulations <- function(simulations_best_temp, chla_obs, zoo_obs){
  
  #Chl-a graph
  chla_sim <- simulations_best_temp[,c(4,6,7,8)]
  
  chla_obs_temp<-as.data.frame(matrix(nrow=nrow(chla_obs), ncol = 4))
  names(chla_obs_temp)[1] <- "chl_a"
  names(chla_obs_temp)[2] <- "Date"
  names(chla_obs_temp)[3] <- "month"
  names(chla_obs_temp)[4] <- "type"
  
  chla_obs_temp$chl_a <- chla_obs$Chlorophyll_a
  chla_obs_temp$Date <- chla_obs$Time
  chla_obs_temp$month <- month(as.Date(chla_obs_temp$Date))
  chla_obs_temp$type <- "observed"
  
  stats_chla <- count(chla_obs_temp, "month")
  stats_chla <- getStats(stats_chla)
  
  all_chla <- rbind(chla_sim, chla_obs_temp)
  all_chla$month <- as.factor(all_chla$month)
  
  graph_chla <- ggplot(data = all_chla, aes(x=month, y=chl_a, fill = type)) +
    geom_boxplot() + scale_x_discrete(labels=c("1" = paste("January (n = ", stats_chla$freq[1],")"), "2" = paste("February (n = ", stats_chla$freq[2],")"),"3" = paste("March (n = ", stats_chla$freq[3],")"), 
                                               "4" = paste("April (n = ", stats_chla$freq[4],")"),"5" = paste("May (n = ", stats_chla$freq[5],")"),
                                               "6" = paste("June (n = ", stats_chla$freq[6],")"), "7" = paste("July (n = ", stats_chla$freq[7],")"),
                                               "8" = paste("August (n = ", stats_chla$freq[8],")"), "9" = paste("September (n = ", stats_chla$freq[9],")"), 
                                               "10" = paste("October (n = ", stats_chla$freq[10],")"), "11" = paste("November (n = ", stats_chla$freq[11],")"), "12" = paste("December (n = ", stats_chla$freq[12],")"))) +
    theme(axis.text.x = element_text(angle=90))  + labs(y="Chlorophyll a (mg Chl_a/m3)")  
  
  #Zooplankton graph
  zoo_sim <- simulations_best_temp[,c(3,6,7,8)]
  
  zoo_obs_temp<-as.data.frame(matrix(nrow=nrow(zoo_obs), ncol = 4))
  names(zoo_obs_temp)[1] <- "zoo"
  names(zoo_obs_temp)[2] <- "Date"
  names(zoo_obs_temp)[3] <- "month"
  names(zoo_obs_temp)[4] <- "type"
  
  zoo_obs_temp$zoo <- zoo_obs$n_mol
  zoo_obs_temp$Date <- zoo_obs$Date
  zoo_obs_temp$month <- month(as.Date(zoo_obs_temp$Date, format = "%m/%d/%Y"))
  zoo_obs_temp$type <- "observed"
  
  stats_zoo <- count(zoo_obs_temp, "month")
  stats_zoo <- getStats(stats_zoo)
  
  all_zoo <- rbind(zoo_sim, zoo_obs_temp)
  all_zoo$month <- as.factor(all_zoo$month)
  
  graph_zoo <- ggplot(data = all_zoo, aes(x=month, y=zoo, fill = type)) +
    geom_boxplot() + scale_x_discrete(labels=c("1" = paste("January (n = ", stats_zoo$freq[1],")"), "2" = paste("February (n = ", stats_zoo$freq[2],")"),"3" = paste("March (n = ", stats_zoo$freq[3],")"), 
                                               "4" = paste("April (n = ", stats_zoo$freq[4],")"),"5" = paste("May (n = ", stats_zoo$freq[5],")"),
                                               "6" = paste("June (n = ", stats_zoo$freq[6],")"), "7" = paste("July (n = ", stats_zoo$freq[7],")"),
                                               "8" = paste("August (n = ", stats_zoo$freq[8],")"), "9" = paste("September (n = ", stats_zoo$freq[9],")"), 
                                               "10" = paste("October (n = ", stats_zoo$freq[10],")"), "11" = paste("November (n = ", stats_zoo$freq[11],")"), "12" = paste("December (n = ", stats_zoo$freq[12],")"))) +
    theme(axis.text.x = element_text(angle=90))  + labs(y="Zooplankton (m mol N/m3)")  
  
  graphs_temp <-list(graph_chla, graph_zoo)
  
  return(graphs_temp)
  
}

getStats <- function(temp_stats){
  stats_obs <- data.frame(matrix(nrow = 12, ncol = 2))
  names(stats_obs)[1] <- "month"
  names(stats_obs)[2] <- "freq"
  
  for (m in 1:12){
    stats_obs[m,1] <- m
    temp <- temp_stats[which(temp_stats$month == m),2]
    stats_obs[m,2] <- ifelse(length(temp)==0, 0, temp)
  }
  
  return(stats_obs)
  
}


