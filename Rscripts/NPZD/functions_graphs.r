##############################
# Functions for final graphs #
##############################

# required packages
library(ggplot2)
library(xts)
library(reshape2)
library(lubridate)
library(stats)
library(dplyr)
library(viridis)
theme_set(theme_bw())

# This function searches for the best 10% of the simulations based on the RMSE. 
getBestSim <- function(id_station){
  
  #Read the file with the errors calculated for all simulations
  
  error_log <- read.csv(paste0(wd, "Output/Final results/NPZD/station", id_station,
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
  simulations_best_rmse <- simulations_best_rmse[,c("phyto", "zoo", "chl_a", "ID", "Date")]
  
  #This function calculates the mean value for each month for chlorophyll-a abundances
  bestSim_month_chla <- reduceSimulations(simulations_best_rmse[,c( "chl_a", "ID", "Date")])
  bestSim_month_chla$Station <- id_station
  
  #This function calculates the mean value for each month for zooplankton abundances
  bestSim_month_zoo <- reduceSimulations(simulations_best_rmse[,c("zoo", "ID", "Date")])
  bestSim_month_zoo$Station <- id_station
  
  #Two elements in the list, the first the best Chlorophyll-a values, and the second the zooplankton abundances
  list_bestSim <- list(bestSim_month_chla, bestSim_month_zoo)

    #Only get first element of the list, that is a dataframe with the values of Chlorophyll-a density for the best 10% best       simulations
    chla_station <- list_bestSim[[1]]

    #This data frame has the average value per date, the previous one has all the values not aggregated
    chla_station_avg <- aggregate( .~roundDate, data = chla_station, FUN = function(x) median(as.numeric(as.character(x))))
    chla_station_avg$Type <- "Chlorophyll a"
    chla_station$Type <- "Chlorophyll a"

    #Only get the second element of the list, that is a dataframe with the values of zooplankton density
    #for the best 10% best simulations
    zoo_station <- list_bestSim[[2]]
    zoo_station_avg <- aggregate( .~roundDate, data = zoo_station, FUN = function(x) median(as.numeric(as.character(x))) )
    zoo_station_avg$Type <- "Zooplankton"
    zoo_station$Type <- "Zooplankton"

    temp_station <- chla_station_avg[,c(1,3,6)]
    temp_station_2 <- zoo_station_avg[,c(1,3,6)]
    
    names(temp_station)[2] <- "Plankton"
    names(temp_station_2)[2] <- "Plankton"
    
    combine_station <- rbind(temp_station, temp_station_2)
    
    #Summary per Date
    chla_station_avg <- chla_station_avg[,c(1,3,4,5)]
    
    all_chla <<- chla_station_avg
    all_chla_complete <- chla_station
    
    names(all_chla)[2] <<- "Plankton"
    names(all_chla_complete)[3] <- "Plankton"
      
    zoo_station_avg <- zoo_station_avg[,c(1,3,4,5)]
    
    all_zoo <<- zoo_station_avg
    all_zoo_complete <- zoo_station
    
    names(all_zoo)[2] <<- "Plankton"
    names(all_zoo_complete)[3] <- "Plankton"
    
    all_output_avg <- rbind(all_chla, all_zoo)
    all_output_avg$Station <- station_code
    all_output_complete <- rbind(all_chla_complete, all_zoo_complete)
        
  return(all_output_complete)
  
}

reduceSimulations <- function(all_sim){
  
  #Create a new date to refer to the first day of the month as data will be averaged per month
  all_sim$year <- year(all_sim$Date)
  all_sim$roundDate <- floor_date(all_sim$Date, unit = "month")
  
  #Calculate average values per month mantaining the different ID of each simulation
  simulations_by_month <- aggregate( .~ID+roundDate, data = all_sim, median )
  
  #Remove the original column date that after all transformations does not contain relevant information anymore
  simulations_by_month <- subset(simulations_by_month, select = -c(Date))
  
  return(simulations_by_month)
  
}

bestSimulations <- function(bestID, station_int){
  #Extract the number of the 10% best simualtions
  count_ID <- length(bestID)
  
  #Read the input data for the correct date for each day
  readDate  <- read.csv(paste0(wd, "Input data/inputData_npzd_station", station_int, ".csv"))
  
  bestSim<-NA
  
  #Extract the details of each simulation 
  
  for(i in 1:count_ID){
    idSim <- bestID[i]
    resultsSim <-read.csv(paste0(wd, "Output/Final results/NPZD/station", station_int,
                                 "_iter2/detailed_simulation_", idSim, ".csv"))

    #Select only the columns of interest (final results, no limitation factors)
    temp_sim <-resultsSim[,c("time_step","phyto","zoo","chl_a","detritus","detritusC","pCO2w")]
    #Attach the id simulation to differentiate the 10% best simulations from each other to calculate the variance
    temp_sim$ID <- idSim
    temp_sim$Date <- as.Date(readDate$Date)

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

# Function graph Chla & Zoo
graph_phyto_zoo <- function(){
    all_zoo_2nd_scale <- all_zoo
    all_zoo_2nd_scale$Plankton <- log(all_zoo_2nd_scale$Plankton+1)*10
    all_chla_1st_scale <- all_chla
    all_chla_1st_scale$Plankton <-log(all_chla_1st_scale$Plankton+1)
    all_chla_zoo_1st_2nd_scale <- rbind(all_chla_1st_scale,all_zoo_2nd_scale)
    
    comparisonGrap_chla_zoo <<- ggplot() +
      geom_line(data = all_chla, aes(x = as.Date(roundDate), y = log(Plankton+1)), linewidth = 1) + 
      geom_line(data = all_zoo, aes(x = as.Date(roundDate), y = log(Plankton+1)*10),colour = 'grey', linewidth = 1, linetype = "dashed") + 
      scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y", expand = c(0, 0)) +
      scale_y_continuous(name = expression(paste("Chlorophyll a (log (mg Chla ", m^-3," + 1))")), sec.axis = sec_axis( trans=~./10, 
                          name=expression(paste("Zooplankton (log (mmol N ", m^-3," + 1))")))) +
      theme(text = element_text(size=10), axis.text.x=element_text(angle=90, hjust=1)) +
      theme(legend.position = "bottom",legend.direction = "horizontal",legend.background = element_rect(),legend.title = element_blank(),
            legend.text = element_text(size=12),legend.text.align = 0,axis.text.x = element_text(colour="black",size=11.5),
            axis.text.y = element_text(colour="black",size=14), axis.title.x = element_text(colour="black",size=16,face = 'bold'),
            axis.title.y = element_text(colour="black",size=16,face = 'bold'), strip.text = element_text(size = 12))+
      labs(x = "Date", y= expression(paste("Zooplankton (log (mmol N ", m^-3," + 1))")))
}
                           
# Function graph Chla
graph_phyto <- function(){
# create chla observation subset
if (area=='bpns'){                             
chla_obs_station <- subset(chla, Station==station_region)
colnames(chla_obs_station) <- c('Station', 'Time',"Chlorophyll_a")
    } else {
    chla_obs_station <- chla
    chla_obs_station$Station <- station_region
    colnames(chla_obs_station) <- c('Time',"Chlorophyll_a",'Station')
    }

# calculate min and max for uncertainty ribbon used in plot chla and zooplankton
values_min <- aggregate( Plankton~roundDate + Station +Type, data = all_output_complete, FUN = min )
values_min$Limit <-"min"
values_max <- aggregate( Plankton~roundDate + Station +Type, data = all_output_complete, FUN = max )
values_max$Limit <-"max"
values_merged <- merge(values_min,values_max, by= c("roundDate", "Station", "Type"))
values_merged <- values_merged[,c(1,2,3,4,6)]
names(values_merged)[4] <- "Min"
names(values_merged)[5] <- "Max"
values_merged_chla <- values_merged[which(values_merged$Type == "Chlorophyll a"),]
values_merged_zoo <- values_merged[which(values_merged$Type == "Zooplankton"),]

# create plot
comparisonGraphAvg_chla <<- ggplot() +
  geom_line(data = all_chla, aes(x = as.Date(roundDate), y = log(Plankton+1)), linewidth = 1) + 
  geom_ribbon(data = values_merged_chla, aes(x=as.Date(roundDate), ymin = log(Min+1), ymax = log(Max+1)), fill = 'black',  alpha = 0.15) +
  geom_point(data=chla_obs_station, aes(x=as.Date(Time), y = log(Chlorophyll_a+1), alpha = 0.25)) +
  scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y", expand = c(0, 0)) +
  theme(text = element_text(size=10), axis.text.x=element_text(angle=90, hjust=1)) +
  theme(legend.position = "none",legend.direction = "vertical",legend.background = element_rect(),legend.title = element_text(size=14,face = 'bold'),
        legend.text = element_text(size=12),legend.text.align = 0,axis.text.x = element_text(colour="black",size=11.5),
        axis.text.y = element_text(colour="black",size=14), axis.title.x = element_text(colour="black",size=16,face = 'bold'),
        axis.title.y = element_text(colour="black",size=16,face = 'bold'), strip.text = element_text(size = 12))+
  labs(x = "Date", y = expression(paste("Chlorophyll a (log (mg Chla ", m^-3," + 1))")))
}

# Function graph Zoo
graph_zoo <- function(){
#zooplankton
zoo_obs_station <- zoo_aggregated
    
# calculate min and max for uncertainty ribbon used in plot chla and zooplankton
values_min <- aggregate( Plankton~roundDate + Station +Type, data = all_output_complete, FUN = min )
values_min$Limit <-"min"
values_max <- aggregate( Plankton~roundDate + Station +Type, data = all_output_complete, FUN = max )
values_max$Limit <-"max"
values_merged <- merge(values_min,values_max, by= c("roundDate", "Station", "Type"))
values_merged <- values_merged[,c(1,2,3,4,6)]
names(values_merged)[4] <- "Min"
names(values_merged)[5] <- "Max"
values_merged_chla <- values_merged[which(values_merged$Type == "Chlorophyll a"),]
values_merged_zoo <- values_merged[which(values_merged$Type == "Zooplankton"),]

if (!is.na(zoo_obs_station[1,1])){
comparisonGraphAvg_zoo <<- ggplot() +
  geom_line(data = all_zoo, aes(x = as.Date(roundDate), y = log(Plankton+1)), linewidth = 1) + 
  geom_ribbon(data = values_merged_zoo, aes(x=as.Date(roundDate), ymin = log(Min+1), ymax = log(Max+1)),  alpha = 0.15) +
  geom_point(data=zoo_obs_station, aes(x=as.Date(Date), y = log(zoo+1))) +
  scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y", expand = c(0, 0)) +
  theme(text = element_text(size=10), axis.text.x=element_text(angle=90, hjust=1)) +
  theme(legend.position = "none",legend.direction = "vertical",legend.background = element_rect(),legend.title = element_text(size=14,face = 'bold'),
        legend.text = element_text(size=12),legend.text.align = 0,axis.text.x = element_text(colour="black",size=11.5),
        axis.text.y = element_text(colour="black",size=14), axis.title.x = element_text(colour="black",size=16,face = 'bold'),
        axis.title.y = element_text(colour="black",size=16,face = 'bold'), strip.text = element_text(size = 12))+
  labs(x = "Date", y= expression(paste("Zooplankton (log (mmol N ", m^-3," + 1))")))
    } else {
    comparisonGraphAvg_zoo <<- ggplot() +
  geom_line(data = all_zoo, aes(x = as.Date(roundDate), y = log(Plankton+1)), linewidth = 1) + 
  geom_ribbon(data = values_merged_zoo, aes(x=as.Date(roundDate), ymin = log(Min+1), ymax = log(Max+1)),  alpha = 0.15) +
  #geom_point(data=zoo_obs_station, aes(x=as.Date(Date), y = log(zoo+1))) +
  scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y", expand = c(0, 0)) +
  theme(text = element_text(size=10), axis.text.x=element_text(angle=90, hjust=1)) +
  theme(legend.position = "none",legend.direction = "vertical",legend.background = element_rect(),legend.title = element_text(size=14,face = 'bold'),
        legend.text = element_text(size=12),legend.text.align = 0,axis.text.x = element_text(colour="black",size=11.5),
        axis.text.y = element_text(colour="black",size=14), axis.title.x = element_text(colour="black",size=16,face = 'bold'),
        axis.title.y = element_text(colour="black",size=16,face = 'bold'), strip.text = element_text(size = 12))+
  labs(x = "Date", y= expression(paste("Zooplankton (log (mmol N ", m^-3," + 1))")))
    }
}
                                 
# Function graph Chla
graph_rel_contributions <- function(){
# relative contributions 
# Let's prepare our data for the graph
  relativeContributions_station <- read.csv(paste0(wd, "Output/Final results/Relative contributions/relativeContributions_month_station",station_region,".csv"))

  relativeContributions_station$Station <- station_region

  values_min <- aggregate( Contribution ~ Date + Type, data = relativeContributions_station, FUN = min )
  values_min$Limit <-"min"
  values_max <- aggregate( Contribution ~ Date + Type, data = relativeContributions_station, FUN = max )
  values_max$Limit <-"max"
  
  contributions_merged <- merge(values_min,values_max, by= c("Date", "Type"))
  
  contributions_merged <- contributions_merged[,c(1,2,3,5)]
  names(contributions_merged)[3] <- "Min"
  names(contributions_merged)[4] <- "Max"
  
  contributions_merged$Type <- factor(contributions_merged$Type, 
                                      levels = c('DIN', 
                                                 "PO4", 
                                                 "SiO4",
                                                 "PAR",
                                                 "Temperature",
                                                 "Zooplankton"))
  contributions_plot$Type <- factor(contributions_plot$Type, 
                               levels = c('DIN', 
                                          "PO4", 
                                          "SiO4",
                                          "PAR",
                                          "Temperature",
                                          "Zooplankton"))

  my_labeller <- as_labeller(c(DIN = 'DIN', PO4 = "PO[4]", SiO4 = "SiO[4]", PAR = "PAR", Temperature = "SST", Zooplankton = "Zooplankton~grazing"), default = label_parsed)
  
# Create the graph
comparisonGraph <- ggplot() +
    geom_line(data = contributions_plot, aes(x = as.Date(Date), y = percentage), linewidth = 0.8) + 
    geom_ribbon(data = contributions_merged, aes(x=as.Date(Date), ymin = Min, ymax = Max),  alpha = 0.15) +
    ylim(0,0.5) +
    facet_wrap(vars(Type), scales = "free_y", labeller = labeller(Type =  my_labeller))+
    scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=7), axis.text.x=element_text(angle=90, hjust=1)) +
    theme(legend.position = "right",legend.direction = "vertical",legend.background = element_rect(),legend.title = element_text(size=14,face = 'bold'),legend.text = element_text(size=12),legend.text.align = 0,axis.text.x = element_text(colour="black",size=12),
 axis.text.y = element_text(colour="black",size=14), axis.title.x = element_text(colour="black",size=16,face = 'bold'),
 axis.title.y = element_text(colour="black",size=16,face = 'bold'), strip.text = element_text(size = 12))+
    labs(x = "Date", y= "Averaged monthly relative contribution")

#Summary Relative contributions (boxplots)
comparisonContribution <- ggplot() +
  geom_boxplot(data = contributions_plot, aes(x = as.factor(Type), y = percentage),outlier.shape = NA) + 
  ylim(0,0.45 )+
  theme(text = element_text(size=12), axis.text.x=element_text(angle=45, hjust=1)) +
  theme(legend.position = "none",legend.direction = "vertical",legend.background = element_rect(),legend.title = element_text(size=14,face = 'bold'),
        legend.text = element_text(size=12),legend.text.align = 0,axis.text.x = element_text(colour="black",size=11.5),
        axis.text.y = element_text(colour="black",size=14), axis.title.x = element_text(colour="black",size=16,face = 'bold'),
        axis.title.y = element_text(colour="black",size=16,face = 'bold'), strip.text = element_text(size = 12))+
  labs(x = "Region", y= "Averaged monthly relative Contribution") 

}




detritus_graphs <- function(){
  
 #get average and min/max
  detritus_data <- simulations_best_rmse %>% 
    group_by(time_step) %>%
    summarise(avg_detritusc = median(detritusC),
              min_detritusc = min(detritusC),
              max_detritusc = max(detritusC),
              avg_detritus = median(detritus),
              min_detritus = min(detritus),
              max_detritus = max(detritus),
              Date = Date)
  
  
  # Create detritus plots  
  graphs_detritus <- ggplot() +
    geom_line(data = detritus_data, aes(x = as.Date(Date), y = avg_detritusc), linewidth = 1) + 
    geom_ribbon(data = detritus_data, aes(x=as.Date(Date), ymin = min_detritusc, ymax = max_detritusc), fill = 'black',  alpha = 0.15) +
    geom_line(data = detritus_data, aes(x = as.Date(Date), y = avg_detritus*50), linewidth = 1, colour = "blue") + 
    geom_ribbon(data = detritus_data, aes(x=as.Date(Date), ymin = min_detritus*50, ymax = max_detritus*50), fill = 'blue',  alpha = 0.15) +
    #geom_point(data=chla_obs_station, aes(x=as.Date(Time), y = log(Chlorophyll_a+1))) +
    scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y", expand = c(0, 0)) +
    scale_y_continuous(name = expression(paste("mg Carbon ", m^-3)), 
                       sec.axis = sec_axis( trans=~./50, 
                                            name=expression(paste("mg Nitrogen ", m^-3)))) +
    theme(text = element_text(size=10), axis.text.x=element_text(angle=90, hjust=1)) +
    theme(legend.position = "none",legend.direction = "vertical",legend.background = element_rect(),legend.title = element_text(size=14,face = 'bold'),
          legend.text = element_text(size=12),legend.text.align = 0,axis.text.x = element_text(colour="black",size=11.5),
          axis.text.y = element_text(colour="black",size=14), axis.title.x = element_text(colour="black",size=16,face = 'bold'),
          axis.title.y = element_text(colour="black",size=16,face = 'bold'), strip.text = element_text(size = 12))+
    labs(x = "Date", y = expression(paste("mg Carbon ", m^-3)))
  
  return(graphs_detritus)
  
}

pco2w_graph <- function(){
  
 #get average and min/max
  pco2w_data <- simulations_best_rmse %>% 
    group_by(time_step) %>%
    summarise(avg_pco2w = median(pCO2w),
              min_pco2w = min(pCO2w),
              max_pco2w = max(pCO2w),
              Date = Date)
  
  
  # Create pco2w plots  
  graphs_pco2w <- ggplot() +
    geom_line(data = pco2w_data, aes(x = as.Date(Date), y = avg_pco2w), linewidth = 1) + 
    geom_ribbon(data = pco2w_data, aes(x=as.Date(Date), ymin = min_pco2w, ymax = max_pco2w), fill = 'black',  alpha = 0.15) +
    geom_point(data = pco2w_validation, aes(x=as.Date(Date), y = pco2w)) +
    scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=10), axis.text.x=element_text(angle=90, hjust=1)) +
    theme(legend.position = "none",legend.direction = "vertical",legend.background = element_rect(),legend.title = element_text(size=14,face = 'bold'),
          legend.text = element_text(size=12),legend.text.align = 0,axis.text.x = element_text(colour="black",size=11.5),
          axis.text.y = element_text(colour="black",size=14), axis.title.x = element_text(colour="black",size=16,face = 'bold'),
          axis.title.y = element_text(colour="black",size=16,face = 'bold'), strip.text = element_text(size = 12))+
    labs(x = "Date", y = expression("µatm"))
  
  return(graphs_pco2w)
  
}
                                 
# create relative contributions
relative_contribution_graph <- function(type){
# Let's prepare our data for the graph
  relativeContributions_station <- read.csv(paste0(wd, "Output/Final results/Relative contributions/relativeContributions_month_station",station_region,".csv"))

  relativeContributions_station$Station <- station_region

  values_min <- aggregate( Contribution ~ Date + Type, data = relativeContributions_station, FUN = min )
  values_min$Limit <-"min"
  values_max <- aggregate( Contribution ~ Date + Type, data = relativeContributions_station, FUN = max )
  values_max$Limit <-"max"
  
  contributions_merged <- merge(values_min,values_max, by= c("Date", "Type"))
  
  contributions_merged <- contributions_merged[,c(1,2,3,5)]
  names(contributions_merged)[3] <- "Min"
  names(contributions_merged)[4] <- "Max"
  
  contributions_merged$Type <- factor(contributions_merged$Type, 
                                      levels = c('DIN', 
                                                 "PO4", 
                                                 "SiO4",
                                                 "PAR",
                                                 "Temperature",
                                                 "Zooplankton"))
  contributions_plot$Type <- factor(contributions_plot$Type, 
                               levels = c('DIN', 
                                          "PO4", 
                                          "SiO4",
                                          "PAR",
                                          "Temperature",
                                          "Zooplankton"))

  my_labeller <- as_labeller(c(DIN = 'DIN', PO4 = "PO[4]", SiO4 = "SiO[4]", PAR = "PAR", Temperature = "SST", Zooplankton = "Zooplankton~grazing"), default = label_parsed)
  
# Create the line graph
if (type == 'line'|type == 'both'){
relative_contribution_line <<- ggplot() +
    geom_line(data = contributions_plot, aes(x = as.Date(Date), y = percentage), linewidth = 0.8) + 
    geom_ribbon(data = contributions_merged, aes(x=as.Date(Date), ymin = Min, ymax = Max),  alpha = 0.15) +
    ylim(0,0.5) +
    facet_wrap(vars(Type), scales = "free_y", labeller = labeller(Type =  my_labeller))+
    scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=7), axis.text.x=element_text(angle=90, hjust=1)) +
    theme(legend.position = "right",legend.direction = "vertical",legend.background = element_rect(),legend.title = element_text(size=14,face = 'bold'),legend.text = element_text(size=12),legend.text.align = 0,axis.text.x = element_text(colour="black",size=12),
 axis.text.y = element_text(colour="black",size=14), axis.title.x = element_text(colour="black",size=16,face = 'bold'),
 axis.title.y = element_text(colour="black",size=16,face = 'bold'), strip.text = element_text(size = 12))+
    labs(x = "Date", y= "Averaged monthly relative contribution")
}
# Create the boxplot graph
if (type == 'boxplot'|type == 'both'){
relative_contribution_boxplot <<- ggplot() +
  geom_boxplot(data = contributions_plot, aes(x = as.factor(Type), y = percentage),outlier.shape = NA) + 
  ylim(0,0.45 )+
  theme(text = element_text(size=12), axis.text.x=element_text(angle=45, hjust=1)) +
  theme(legend.position = "none",legend.direction = "vertical",legend.background = element_rect(),legend.title = element_text(size=14,face = 'bold'),
        legend.text = element_text(size=12),legend.text.align = 0,axis.text.x = element_text(colour="black",size=11.5),
        axis.text.y = element_text(colour="black",size=14), axis.title.x = element_text(colour="black",size=16,face = 'bold'),
        axis.title.y = element_text(colour="black",size=16,face = 'bold'), strip.text = element_text(size = 12))+
  labs(x = "Region", y= "Averaged monthly relative Contribution") 
}
}