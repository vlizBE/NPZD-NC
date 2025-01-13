library(ggplot2)
library(xts)
library(reshape2)
library(lubridate)
library(stats)
library(plyr)
library(RColorBrewer)


convertZoo <- function(zoo_obs, num_station){
  #Calculate biomass from ind/m3 to mg C/m3
  
  taxa_select <- c("Calanoida", "Noctiluca", "Harpacticoida", "Appendicularia")
  
  ##############Calanoida
  #Conversion from ind/m3 to mgC/m3
  calanoida <- zoo_obs[which(zoo_obs$Taxon ==  taxa_select[1]),]
  calanoida$biomass <- calanoida$Density * 0.0006     #Median value based on Brun et al. 2016
  
  #Conver carbon to mol (mmolC/m3)
  calanoida$biomass_carbon <- calanoida$biomass / 12.0107
  
  #convert to N with C:N ratio of 5.5 to 6 from Van Nieuwerburgh et al. 2004
  #result in mmol N / m3
  calanoida$n_mol <- calanoida$biomass_carbon / 5.75
  
  
  ##############Noctiluca
  
  noctiluca <- zoo_obs[which(zoo_obs$Taxon ==  taxa_select[2]),]
  noctiluca$biomass <- noctiluca$Density * 0.0003     #value based on Loder et al 2012
  noctiluca$biomass_carbon <- noctiluca$biomass / 12.0107
  noctiluca$n_mol <- noctiluca$biomass_carbon / 3.35
  
  
  ##############Harpacticoida
  harpa <- zoo_obs[which(zoo_obs$Taxon ==  taxa_select[3]),]
  harpa$biomass <- harpa$Density * 0.001     #value based on Sautour and Castel 1995
  harpa$biomass_carbon <- harpa$biomass / 12.0107
  harpa$n_mol <- harpa$biomass_carbon / 4.5
  
  
  ##############Appendicularia
  append <- zoo_obs[which(zoo_obs$Taxon ==  taxa_select[4]),]
  append$biomass <- append$Density * 0.004     #value based on Lombard et al. 2009
  append$biomass_carbon <- append$biomass / 12.0107
  append$n_mol <- append$biomass_carbon / 4.08
  
  zoo_converted <- rbind(calanoida, noctiluca, harpa, append)
  zoo_converted$Date<-as.Date(zoo_converted$Date, format = "%Y-%m-%d")
  
  if(num_station == "offshore"){
    summary_zoo <- aggregate(n_mol ~  Date + Station, data = zoo_converted, FUN = "sum")
  }else{
    summary_zoo <- aggregate(n_mol ~  Date, data = zoo_converted, FUN = "sum")
  }
  
  summary_zoo$month <- month(summary_zoo$Date)
  
  list_results <- list(zoo_converted, summary_zoo)
  
  return(list_results)
  
}

get_validation_data <- function(){
    if (area == "bpns") {
#Get observations of Chla and zooplankton
observations <- get_abiotic_pigment(startdate, stopdate, station_code,station_region) 

observations$Station <- station_region

#Observations from LifeWatch
#Station, Date, Chlorophyll_a
chla <- observations[,c("Station", "Date", "Chla")]
chla<-chla[which(chla$Station==station_region),]

chla$Time<-as.Date(chla$Date)
chla <- chla[,c("Station", "Time", "Chla")]
chla<-subset(subset(chla, Time > startdate), Time < stopdate)
chla <<- chla[!is.na(chla$Chla),]

#Zooplankton data
taxa <- c("Calanoida", "Noctiluca", "Harpacticoida", "Appendicularia")
zoo_observations<- get_zooplankton(startdate, stopdate, taxa, station_code) #read.csv(paste0(wd, "Input data/NPZD/zooscan.csv")) # only zooscan data
zoo_observations$Station <- station_region
zoo_observations<-zoo_observations[which(zoo_observations$Station == station_region),]

zoo_observations<-subset(subset(zoo_observations, Date > startdate), Date < stopdate)

zoo_converted <- convertZoo(zoo_observations, station_region)
zoo_perTaxa <- as.data.frame(zoo_converted[[1]])
zoo_aggregated <- as.data.frame(zoo_converted[[2]])
write.csv(zoo_aggregated, paste0(wd, "Input data/zoo_converted_station", station_region, ".csv"))

#zoo_aggregated <- read.csv(paste0(wd, "Input data/zoo_converted_station", station_region, ".csv"))
#zoo_aggregated <- zoo_aggregated[,c(-1)]
zoo_aggregated$Date <- as.Date(zoo_aggregated$Date)

#Remove outlier station nearshore
zoo_aggregated <<- zoo_aggregated[which(zoo_aggregated$n_mol < 0.5),]
colnames(zoo_aggregated) <<- c("Date", "station", "zoo", "month")
        } else {
chla <<- validation_data
 colnames(chla) <<- c("Time","Chla")
zoo_aggregated <<- data.frame("Date" = NA, "zoo" = NA)
        }

# get pco2 seawater data to validate
if (!is.na(carbon_validation_folder)){
file_list <- list.files(path = paste0(wd,"Input data/",carbon_validation_folder), pattern = "*.csv", full.names = TRUE)
data_list <- lapply(file_list, function(file) {
  data <- read.csv(file)
    selected_data <- data %>%
    select(c('Date.Time','pCO2..uatm.')) # Use any_of to handle non-existent columns
    selected_data <- na.omit(selected_data)

  return(selected_data)
})
pco2w_validation <- do.call(rbind, data_list)
colnames(pco2w_validation) <- c('Date', 'pco2w')
pco2w_validation <- subset(subset(pco2w_validation, Date > startdate), Date < stopdate)         

pco2w_validation$Date <- as.Date(pco2w_validation$Date)
pco2w_validation <<-  pco2w_validation %>% group_by(Date) %>%
                                     summarise(Date = mean(Date),
                                                      pco2w = mean(pco2w)) 
    }
    
}

calculateErrors <- function(num_station, chla, zoo_aggregated, pco2w_validation, num_sim){
  
  cont <- 1
  possible_parameters <- read.csv(paste0(wd, "Output/Intermediate results/station", num_station, "_iter1/possible_parameters_iter1.csv"))
  possible_parameters <- possible_parameters[,-c(1)]
  errors_parameters <- createErrorFile(num_sim)
  
  readDate <- read.csv(paste0(wd, "Input data/inputData_npzd_station", num_station, ".csv"))
  
  

  for(i in 1:num_sim){
    
    output <- read.csv(paste0(wd, "Output/Intermediate results/station", num_station, "_iter1/detailed_simulation_", i, ".csv"))
    output$Date <- as.Date(readDate$Date)
    
    #remove dummy years (one leap year: 2012)
    #output <- output[-c(1:1096),] 
    
    if(!anyNA(output)){
      
      ################
      #calculate error
      temp_error<-calculateRMSE(output, possible_parameters[i,12],i, zoo_aggregated[,c("Date","zoo")], chla[,c("Time","Chla")], pco2w_validation, num_station)
        
      

      errors_parameters$ID_sim[cont] <- i
      errors_parameters$rmse_zoo[cont] <- temp_error[[1]]
      errors_parameters$rmse_phyto[cont] <- temp_error[[2]]
      errors_parameters$rmse_pco2w[cont] <- temp_error[[3]]
      errors_parameters$total_rmse[cont] <- temp_error[[1]] + temp_error[[2]] + temp_error[[3]]
      
      #####
      #attach set of parameters
      
      errors_parameters[cont,c(6:19)] <-possible_parameters[i,]
      
      cont <- cont + 1
      
    }
  }
  
  write.csv(errors_parameters, paste0(wd, "Output/Intermediate results/station", num_station,"_iter1/errors_parameters_iter1.csv"))
  
  return(errors_parameters)
  
}


createErrorFile <- function(numSim){
  errors_parameters <- data.frame(matrix(nrow = numSim, ncol = 19))
  names(errors_parameters)[1] <- "ID_sim"
  names(errors_parameters)[2] <- "rmse_phyto"
  names(errors_parameters)[3] <- "rmse_zoo"
  names(errors_parameters)[4] <- "rmse_pco2w"
  names(errors_parameters)[5] <- "total_rmse"
  names(errors_parameters)[6] <- "ID"
  names(errors_parameters)[7] <- "maxUptake"
  names(errors_parameters)[8] <- "ksPAR"
  names(errors_parameters)[9] <- "ksDIN"
  names(errors_parameters)[10] <- "ksP"
  names(errors_parameters)[11] <- "ksSi"
  names(errors_parameters)[12] <- "maxGrazing"
  names(errors_parameters)[13] <- "ksGrazing"
  names(errors_parameters)[14] <- "pFaeces"
  names(errors_parameters)[15] <- "excretionRate"
  names(errors_parameters)[16] <- "mortalityRate"
  names(errors_parameters)[17] <- "ChlNratio"
  names(errors_parameters)[18] <- "Tobs"
  names(errors_parameters)[19] <- "kd"
  return(errors_parameters)
}

calculateRMSE <- function(resultsSimulation1, chla_conversion, id_temp, zoo_obs, chla_obs, pco2w_validation, numStation){
  
  #get only zoo, phyto and chl_a results
  resultsSimulation1 <- resultsSimulation1[,c("phyto", "zoo", "chl_a", "pCO2w", "Date")]
  
  ######################
  #Calculate RMSE zooplankton
  if (is.na(zoo_obs[1,1])){
  rmse_sim_zoo <- 0
      } else {
  names(zoo_obs)[1]<-"Date_obs"
  names(zoo_obs)[2]<-"zoo_obs"
      
  
  temp_zoo <<-merge(resultsSimulation1, zoo_obs , by.x = "Date", by.y= "Date_obs" ,  all.x = FALSE)
  
  rmse_sim_zoo <- (sqrt(mean((temp_zoo$zoo - temp_zoo$zoo_obs)^2)))/mean(temp_zoo$zoo_obs)
  class(temp_zoo$Date)
  }
  ######################
  #Calculate error Chl_a in units of mmolN
  names(chla_obs)[1]<-"Date_obs"
  names(chla_obs)[2]<-"chla_obs"
  
  temp_chla <<- merge(resultsSimulation1, chla_obs , by.x = "Date", by.y= "Date_obs" ,  all.x = FALSE)
  
  rmse_sim_chla <- sqrt(mean((temp_chla$chl_a - temp_chla$chla_obs)^2))
  rmse_sim_phyto <- (rmse_sim_chla/chla_conversion)/(mean(temp_chla$chla_obs)/chla_conversion)
  rmse_sim_chla <- rmse_sim_chla/mean(temp_chla$chla_obs)

  ######################
    #Calculate RMSE pCO2w
   if (is.na(pco2w_validation[1,1])){
  rmse_sim_pco2w <- 0
      } else {
       pco2w_obs <<- pco2w_validation
       names(pco2w_obs)[1]<<-"Date_obs"
       names(pco2w_obs)[2]<<-"pCO2w_obs"
      
  
  temp_pco2w <<-merge(resultsSimulation1, pco2w_obs , by.x = "Date", by.y= "Date_obs" ,  all.x = FALSE)
  
  rmse_sim_pco2w <- (sqrt(mean((temp_pco2w$pCO2w - temp_pco2w$pCO2w_obs)^2))/mean(temp_pco2w$pCO2w_obs))
  }

  ######################
  # get rmse's in list
  list_results <- list(rmse_sim_zoo, rmse_sim_phyto, rmse_sim_pco2w)
    
  return(list_results)
  
}


createBestSim <- function(ID_minError, station_int, chla_obs, zoo_obs, pco2w_obs){
  
  countID <- length(ID_minError)
  
  readDate  <- read.csv(paste0(wd, "Input data/inputData_npzd_station", station_int, ".csv"))
  
  bestSim<-NA
  
  for(i in 1:countID){
    idSim <- ID_minError[i]
   
    resultsSim <-read.csv(paste0(wd, "Output/Intermediate results/station", station_int, "_iter1/detailed_simulation_", idSim, ".csv"))
    temp_sim <-resultsSim[,c("time_step", "phyto", "zoo", "chl_a", "pCO2w")]
    temp_sim$ID <- idSim
    temp_sim$Date <- as.Date(readDate$Date)
    
    bestSim <- rbind(bestSim,temp_sim)
  }
  
  bestSim <- bestSim[-c(1),]
  bestSim$Date <- as.Date(bestSim$Date)
  
  
  #####Chl_a
  bestSim_mean_chla <- aggregate( chl_a ~ Date, data = bestSim, mean )
  
  figureBest_chla <- ggplot() + 
    geom_line(data = bestSim, aes(x=as.Date(Date), y = chl_a),alpha = 0.1) +
    geom_line(data = bestSim_mean_chla, aes(x=as.Date(Date), y = chl_a),colour = "blue", linewidth = 1) +
    geom_point(data = chla_obs, aes(x = as.Date(Time), y= Chlorophyll_a), colour = "red")+
    labs(x = "Time", y = "Chl_a (mg Chl_a/m3)") + 
    scale_x_date(date_breaks = "6 month", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=11), axis.text.x=element_text(angle=90, hjust=1))
 
  ##### Zoo
  bestSim_mean_zoo <- aggregate( zoo ~ Date, data = bestSim, mean )
  
  figureBest_zoo <- ggplot() + 
    geom_line(data = bestSim, aes(x=as.Date(Date), y = zoo),alpha = 0.1) +
    geom_line(data = bestSim_mean_zoo, aes(x=as.Date(Date), y = zoo),colour = "blue", linewidth = 1) +
    geom_point(data = zoo_obs, aes(x = as.Date(Date), y= zoo), colour = "red")+
    labs(x = "Time", y = "Zooplankton (mmol N/m3)") + 
    scale_x_date(date_breaks = "6 month", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=11), axis.text.x=element_text(angle=90, hjust=1))

   ##### pCO2w
  bestSim_mean_pco2w <- aggregate( pCO2w ~ Date, data = bestSim, mean)
  
  figureBest_pco2w <- ggplot() + 
    geom_line(data = bestSim, aes(x=as.Date(Date), y = pCO2w),alpha = 0.1) +
    geom_line(data = bestSim_mean_pco2w, aes(x=as.Date(Date), y = pCO2w),colour = "blue", linewidth = 1) +
    geom_point(data = pco2w_obs, aes(x = as.Date(Date), y= pCO2w), colour = "red")+
    labs(x = "Time", y = "pCO2w µatm") + 
    scale_x_date(date_breaks = "6 month", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=11), axis.text.x=element_text(angle=90, hjust=1))
  
  #####Combined
  bestSim_mean <- aggregate( phyto ~ Date, data = bestSim, mean )
  
  figureBest_combined <- ggplot() + 
    geom_line(data = bestSim, aes(x=as.Date(Date), y = phyto),alpha = 0.1, colour = "blue") +
    geom_line(data = bestSim_mean, aes(x=as.Date(Date), y = phyto),colour = "blue", linewidth = 1) +
    geom_line(data = bestSim, aes(x=as.Date(Date), y = zoo),alpha = 0.1) +
    geom_line(data = bestSim_mean_zoo, aes(x=as.Date(Date), y = zoo),colour = "black", linewidth = 1) +
    labs(x = "Time", y = "Plankton (mmol N/m3)") + 
    scale_x_date(date_breaks = "6 month", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=11), axis.text.x=element_text(angle=90, hjust=1))
  
  ############################
  ###logarithmic scale
  
  #####Chl_a
  figureBest_chla_log <- ggplot() + 
    geom_line(data = bestSim, aes(x=as.Date(Date), y = log(chl_a)),alpha = 0.1) +
    geom_line(data = bestSim_mean_chla, aes(x=as.Date(Date), y = log(chl_a)),colour = "blue", linewidth = 1) +
    geom_point(data = chla_obs, aes(x = as.Date(Time), y= log(Chlorophyll_a)), colour = "red")+
    labs(x = "Time", y = "Chl_a (log mg Chl_a/m3)") + 
    scale_x_date(date_breaks = "6 month", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=11), axis.text.x=element_text(angle=90, hjust=1))
  
  ##### Zoo
  figureBest_zoo_log <- ggplot() + 
    geom_line(data = bestSim, aes(x=as.Date(Date), y = log(zoo)),alpha = 0.1) +
    geom_line(data = bestSim_mean_zoo, aes(x=as.Date(Date), y = log(zoo)),colour = "blue", linewidth = 1) +
    geom_point(data = zoo_obs, aes(x = as.Date(Date), y= log(zoo)), colour = "red")+
    labs(x = "Time", y = "Zooplankton (log mmol N/m3)") + 
    scale_x_date(date_breaks = "6 month", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=11), axis.text.x=element_text(angle=90, hjust=1))
  
  
  #####Combined
  figureBest_combined_log <- ggplot() + 
    geom_line(data = bestSim, aes(x=as.Date(Date), y = log(phyto)),alpha = 0.1, colour = "blue") +
    geom_line(data = bestSim_mean, aes(x=as.Date(Date), y = log(phyto)),colour = "blue", linewidth = 1) +
    geom_line(data = bestSim, aes(x=as.Date(Date), y = log(zoo)),alpha = 0.1) +
    geom_line(data = bestSim_mean_zoo, aes(x=as.Date(Date), y = log(zoo)),colour = "black", linewidth = 1) +
    labs(x = "Time", y = "Plankton (log mmol N/m3)") + 
    scale_x_date(date_breaks = "6 month", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=11), axis.text.x=element_text(angle=90, hjust=1))
  
  
  listgraphs <- list(figureBest_chla, figureBest_zoo, figureBest_combined, bestSim_mean_pco2w,
                     figureBest_chla_log, figureBest_zoo_log, figureBest_combined_log) 
  
  return(listgraphs)
}

bestSimulations <- function(bestID, station_int){
  count_ID <- length(bestID)
  
  readDate  <- read.csv(paste0(wd, "Input data/inputData_npzd_station", station_int, ".csv"))
  
  bestSim<-NA
  
  for(i in 1:count_ID){
    idSim <- bestID[i]
    resultsSim <-read.csv(paste0(wd, "Output/Intermediate results/station", station_int, "_iter1/detailed_simulation_", idSim, ".csv"))
    temp_sim <-resultsSim[,c("time_step", "phyto", "zoo", "chl_a", "pCO2w", "detritus", "detritusC")]
    temp_sim$ID <- idSim
    temp_sim$Date <- as.Date(readDate$Date)

    bestSim <- rbind(bestSim,temp_sim)
  }
  
  bestSim <- bestSim[-c(1),]
  bestSim$Date <- as.Date(bestSim$Date)
  bestSim$month <- month(bestSim$Date)
  bestSim$type <- "simulated"
  
  return(bestSim)
}

monthlyBestSimulations <- function(simulations_best_temp, chla_obs, zoo_obs, pco2w_validation){
  
  #Chl-a graph
  chla_sim <- simulations_best_temp[,c("chl_a", "Date", "month", "type")]
  
  chla_obs_temp<-as.data.frame(matrix(nrow=nrow(chla_obs), ncol = 4))
  names(chla_obs_temp)[1] <- "chl_a"
  names(chla_obs_temp)[2] <- "Date"
  names(chla_obs_temp)[3] <- "month"
  names(chla_obs_temp)[4] <- "type"
  
  chla_obs_temp$chl_a <- chla_obs$Chla
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
    if(!is.na(zoo_obs[1,1])){
  zoo_sim <- simulations_best_temp[,c("zoo","Date","month","type")]
  
  zoo_obs_temp<-as.data.frame(matrix(nrow=nrow(zoo_obs), ncol = 4))
  names(zoo_obs_temp)[1] <- "zoo"
  names(zoo_obs_temp)[2] <- "Date"
  names(zoo_obs_temp)[3] <- "month"
  names(zoo_obs_temp)[4] <- "type"
  
  zoo_obs_temp$zoo <- zoo_obs$zoo
  zoo_obs_temp$Date <- zoo_obs$Date
  zoo_obs_temp$month <- month(as.Date(zoo_obs_temp$Date))
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
  
  } 

#pco2w graph
  if(!is.na(pco2w_validation[1,1])){
  pco2w_sim <- simulations_best_temp[,c("pCO2w","Date","month","type")]
  
  pco2w_obs_temp<-as.data.frame(matrix(nrow=nrow(pco2w_validation), ncol = 4))
  names(pco2w_obs_temp)[1] <- "pCO2w"
  names(pco2w_obs_temp)[2] <- "Date"
  names(pco2w_obs_temp)[3] <- "month"
  names(pco2w_obs_temp)[4] <- "type"
  
  pco2w_obs_temp$pCO2w <- pco2w_validation$pco2w
  pco2w_obs_temp$Date <- pco2w_validation$Date
  pco2w_obs_temp$month <- month(as.Date(pco2w_obs_temp$Date))
  pco2w_obs_temp$type <- "observed"
  
  stats_pco2w <- count(pco2w_obs_temp, "month")
  stats_pco2w <- getStats(stats_pco2w)
  
  all_pco2w <- rbind(pco2w_sim, pco2w_obs_temp)
  all_pco2w$month <- as.factor(all_pco2w$month)
  
  graph_pco2w <- ggplot(data = all_pco2w, aes(x=month, y=pCO2w, fill = type)) +
    geom_boxplot() + scale_x_discrete(labels=c("1" = paste("January (n = ", stats_pco2w$freq[1],")"), "2" = paste("February (n = ", stats_pco2w$freq[2],")"),"3" = paste("March (n = ", stats_pco2w$freq[3],")"), 
                                               "4" = paste("April (n = ", stats_pco2w$freq[4],")"),"5" = paste("May (n = ", stats_pco2w$freq[5],")"),
                                               "6" = paste("June (n = ", stats_pco2w$freq[6],")"), "7" = paste("July (n = ", stats_pco2w$freq[7],")"),
                                               "8" = paste("August (n = ", stats_pco2w$freq[8],")"), "9" = paste("September (n = ", stats_pco2w$freq[9],")"), 
                                               "10" = paste("October (n = ", stats_pco2w$freq[10],")"), "11" = paste("November (n = ", stats_pco2w$freq[11],")"), "12" = paste("December (n = ", stats_pco2w$freq[12],")"))) +
    theme(axis.text.x = element_text(angle=90))  + labs(y="pCO2w (µatm)")  

  } 

    if (!is.na(zoo_obs[1,1]) & !is.na(pco2w_validation[1,1])) {
        graphs_temp <-list(graph_chla, graph_zoo, graph_pco2w)
         }
    if (!is.na(zoo_obs[1,1]) & is.na(pco2w_validation[1,1])) {
            graphs_temp <- list(graph_chla, graph_zoo)
        } 
    if (is.na(zoo_obs[1,1]) & !is.na(pco2w_validation[1,1])) {
            graphs_temp <- list(graph_chla, graph_pco2w)
        } 
    if (is.na(zoo_obs[1,1]) & is.na(pco2w_validation[1,1])) {
            graphs_temp <- graph_chla
        } 
    
  return(graphs_temp)
  
}


increaseError_springSummer <- function(simulations_error, station_int, chla_obs){
  copy_simulationsError <- simulations_error
  copy_simulationsError$total_rmse_2 <- NA
  copy_simulationsError$total_rmse_3 <- NA
  
  readDate  <- read.csv(paste0(wd, "Input data/inputData_npzd_station", station_int, ".csv"))
  
  if(station_int == "offshore"){
    chla_obs <- aggregate(chla_obs$Chla, by=list(chla_obs$Time), FUN="median",  data = chla_obs)  
    names(chla_obs)[1] <- "Time"
    names(chla_obs)[2] <- "chla_obs"
  }
  
  for(i in 1:nrow(simulations_error)){
    idSim <- simulations_error$ID_sim[i]
    resultsSim <- read.csv(paste0(wd, "Output/Intermediate results/station", station_int, "_iter1/detailed_simulation_", idSim, ".csv"))
    temp_sim <-resultsSim[,c("time_step",	"phyto", "zoo",	"chl_a")]
    temp_sim$ID <- idSim
    temp_sim$Date <- as.Date(readDate$Date)

    #extract only the simulation rows that have an observation
    #Chl_a
    only_obs_dates1 <- merge(temp_sim, chla_obs, by.x = "Date", by.y= "Time" ,  all.x = FALSE)
    
     if (area == "bpns"){
           if(station_int != "offshore"){
                  only_obs_dates1 <- only_obs_dates1[,c(1,3,5,6,8)]
           } else {
                 only_obs_dates1 <- only_obs_dates1[,c(1,3,5,6,7)]
           }
    } else {
           only_obs_dates1 <- only_obs_dates1[,c(1,3,5,6,7)]
    }
    
    only_obs_dates1$Doy <- yday(only_obs_dates1$Date)
    names(only_obs_dates1)[2]<- "phyto_sim"
    names(only_obs_dates1)[3]<- "chla_sim"
    names(only_obs_dates1)[5]<- "chla_obs"
    
    only_obs_dates1$Difference <- only_obs_dates1$chla_sim - only_obs_dates1$chla_obs
    
    #Increase error by 10 for Spring (21/03 until 21/06)
    only_obs_dates1$Difference_spring <- ifelse( only_obs_dates1$Doy >= 80  & only_obs_dates1$Doy <= 172 , only_obs_dates1$Difference * 10 , only_obs_dates1$Difference)   
    
    rmse_sim_chla_spring <- sqrt(mean((only_obs_dates1$Difference_spring)^2))
    rmse_sim_phyto2 <- rmse_sim_chla_spring/simulations_error$ChlNratio[i]
    
    copy_simulationsError$total_rmse_2[i] <- copy_simulationsError$rmse_zoo[i] + rmse_sim_phyto2
    
    #Increase error by 10 for Summer (21/06 until 21/09)
    only_obs_dates1$Difference_summer <- ifelse( only_obs_dates1$Doy > 172  & only_obs_dates1$Doy <= 264 , only_obs_dates1$Difference * 10 , only_obs_dates1$Difference)   
    
    rmse_sim_chla_summer <- sqrt(mean((only_obs_dates1$Difference_summer)^2))
    rmse_sim_phyto3 <- rmse_sim_chla_summer/simulations_error$ChlNratio[i]
    
    copy_simulationsError$total_rmse_3[i] <- copy_simulationsError$rmse_zoo[i] + rmse_sim_phyto3
    
  }
  
  return(copy_simulationsError)
}


createNewParameters_def <- function(min_error_rmse, min_q, max_q){
  newParameters <- data.frame(matrix(nrow = 13, ncol = 3))
  
  names(newParameters)[1] <- "Parameter"
  names(newParameters)[2] <- "Min"
  names(newParameters)[3] <- "Max"
  
  newParameters[1,1] <- "maxUptake"
  newParameters[2,1] <- "ksPAR"
  newParameters[3,1] <- "ksDIN"
  newParameters[4,1] <- "ksP"
  newParameters[5,1] <- "ksSi"
  newParameters[6,1] <- "maxGrazing"
  newParameters[7,1] <- "ksGrazing"
  newParameters[8,1] <- "pFaeces"
  newParameters[9,1] <- "excretionRate"
  newParameters[10,1] <- "mortalityRate"
  newParameters[11,1] <- "ChlNratio"
  newParameters[12,1] <- "Tobs"
  newParameters[13,1] <- "kd" 
  
  newParameters[1,2] <- quantile(min_error_rmse$maxUptake)[min_q]
  newParameters[1,3] <- quantile(min_error_rmse$maxUptake)[max_q]
  
  newParameters[2,2] <- trunc(quantile(min_error_rmse$ksPAR)[min_q])
  newParameters[2,3] <- trunc(quantile(min_error_rmse$ksPAR)[max_q])
  
  newParameters[3,2] <- quantile(min_error_rmse$ksDIN)[min_q]
  newParameters[3,3] <- quantile(min_error_rmse$ksDIN)[max_q]
  
  newParameters[4,2] <- quantile(min_error_rmse$ksP)[min_q]
  newParameters[4,3] <- quantile(min_error_rmse$ksP)[max_q]
  
  newParameters[5,2] <- quantile(min_error_rmse$ksSi)[min_q]
  newParameters[5,3] <- quantile(min_error_rmse$ksSi)[max_q]
  
  newParameters[6,2] <- quantile(min_error_rmse$maxGrazing)[min_q]
  newParameters[6,3] <- quantile(min_error_rmse$maxGrazing)[max_q]
  
  newParameters[7,2] <- quantile(min_error_rmse$ksGrazing)[min_q]
  newParameters[7,3] <- quantile(min_error_rmse$ksGrazing)[max_q]
  
  newParameters[8,2] <- quantile(min_error_rmse$pFaeces)[min_q]
  newParameters[8,3] <- quantile(min_error_rmse$pFaeces)[max_q]
  
  newParameters[9,2] <- quantile(min_error_rmse$excretionRate)[min_q]
  newParameters[9,3] <- quantile(min_error_rmse$excretionRate)[max_q]
  
  newParameters[10,2] <- quantile(min_error_rmse$mortalityRate)[min_q]
  newParameters[10,3] <- quantile(min_error_rmse$mortalityRate)[max_q]
  
  newParameters[11,2] <- quantile(min_error_rmse$ChlNratio)[min_q]
  newParameters[11,3] <- quantile(min_error_rmse$ChlNratio)[max_q]
  
  newParameters[12,2] <- quantile(min_error_rmse$Tobs)[min_q]
  newParameters[12,3] <- quantile(min_error_rmse$Tobs)[max_q]
  
  newParameters[13,2] <- quantile(min_error_rmse$kd)[min_q]
  newParameters[13,3] <- quantile(min_error_rmse$kd)[max_q]
  
  return(newParameters)
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
