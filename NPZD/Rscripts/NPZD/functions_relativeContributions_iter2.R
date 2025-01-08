library(ggplot2)
library(xts)
library(reshape2)
library(lubridate)
library(dplyr)
library(ggpubr)
theme_set(theme_bw())


calculateContributions <- function(selectedSim, station){
 
  relativeContributions <- data.frame(matrix(nrow = 0, ncol = 4))
  names(relativeContributions)[1] <- "Date"
  names(relativeContributions)[2] <- "Type"
  names(relativeContributions)[3] <- "Contribution"
  names(relativeContributions)[4] <- "ID"
  
  numSim <- nrow(selectedSim)
  
  for (i in 1:numSim){
    print(paste0("ID sim: ", selectedSim$ID_sim[i]))
    #Read results simulations and remove 3 dummy years
    values <- read.csv(paste0(wd, "Output/Final results/NPZD/station", station, "_iter2/detailed_simulation_",selectedSim$ID_sim[i],".csv"))
    values <- values[-c(1:1096),]
    
    #calculate relative contributions
    contributions_all <- calculateContributions_sim(values)
    contributions_all$ID <- selectedSim$ID_sim[i]
    
    contributions_month <- reduceContributions(contributions_all)

    temp_new <- data.frame(matrix(nrow = nrow(contributions_month), ncol = 4))
    names(temp_new)[1] <- "Date"
    names(temp_new)[2] <- "Type"
    names(temp_new)[3] <- "Contribution"
    names(temp_new)[4] <- "ID"
    
    temp_new$Date <- contributions_month$roundDate
    temp_new$Type <- contributions_month$Type
    temp_new$Contribution <- contributions_month$Contribution
    temp_new$ID <- contributions_month$ID
    
    
    relativeContributions<- rbind(relativeContributions, temp_new)
  }
  
  return(relativeContributions)
}

calculateContributions_sim <- function(simulation_results){
  
  temp_cont <- data.frame(matrix(nrow = 0, ncol = 3))
  names(temp_cont)[1] <- "Date"
  names(temp_cont)[2] <- "Type"
  names(temp_cont)[3] <- "Contribution"
  
  
  for(j in 1:nrow(simulation_results)){
    #Calculate absolute limitations
    temp_limitation <- data.frame(matrix(nrow=6, ncol=2))
    
    #Limitation for DIN
    temp_limitation[1,1] <- simulation_results$din_lim[j]
    temp_limitation[1,2] <- 1-temp_limitation[1,1]
    
    #Limitation for PO4
    temp_limitation[2,1] <- simulation_results$p_lim[j]
    temp_limitation[2,2] <- 1-temp_limitation[2,1]
    
    #Limitation for SiO4
    temp_limitation[3,1] <- simulation_results$si_lim[j]
    temp_limitation[3,2] <- 1-temp_limitation[3,1]
    
    #Limitation for Temperature
    #The Temp_lim equation can give values slightly higher than 1
    temp_limitation[4,1] <- min(1,simulation_results$temp_lim[j])
    temp_limitation[4,2] <- 1-temp_limitation[4,1]
    
    #Limitation for PAR
    temp_limitation[5,1] <- simulation_results$par_lim[j]
    temp_limitation[5,2] <- 1-temp_limitation[5,1]
    
    #Limitation for Grazing zooplankton
    temp_limitation[6,1] <- simulation_results$grazing_lim[j]
    temp_limitation[6,2] <- 1-temp_limitation[6,1]
    
    #Sum all absolute limitations
    total_limitation <- sum(temp_limitation[,2])
    
    #Calculate relative contributions
    aux_cont <- data.frame(matrix(nrow = 6, ncol = 3))
    names(aux_cont)[1] <- "Date"
    names(aux_cont)[2] <- "Type"
    names(aux_cont)[3] <- "Contribution"
    
    aux_cont[,1]<-as.Date(simulation_results$Date[j])
    
    aux_cont[1,2] <- "DIN"
    aux_cont[1,3] <- temp_limitation[1,2]/total_limitation
    aux_cont[2,2] <- "PO4"
    aux_cont[2,3] <- temp_limitation[2,2]/total_limitation
    aux_cont[3,2] <- "SiO4"
    aux_cont[3,3] <- temp_limitation[3,2]/total_limitation
    aux_cont[4,2] <- "Temperature"
    aux_cont[4,3] <- temp_limitation[4,2]/total_limitation
    aux_cont[5,2] <- "PAR"
    aux_cont[5,3] <- temp_limitation[5,2]/total_limitation
    aux_cont[6,2] <- "Zooplankton"
    aux_cont[6,3] <- temp_limitation[6,2]/total_limitation
    
    temp_cont<-rbind(temp_cont, aux_cont)
  }
  
  return(temp_cont)
}

reduceContributions <- function(all_cont){
  
  all_cont$year <- year(all_cont$Date)
  all_cont$roundDate <- floor_date(all_cont$Date, unit = "month")
  
  contributions_by_month <- aggregate( .~Type+roundDate, data = all_cont, mean )
  contributions_by_month <- subset(contributions_by_month, select = -c(Date))
  
  return(contributions_by_month)
  
}

normalizeData <- function(contributions){
  data_plot <- contributions  %>%
    group_by(Date, Type) %>%
    dplyr::summarise(n = sum(Contribution)) %>%
    dplyr::mutate(percentage = n / sum(n))
  
  return(data_plot)
}

createContributions <- function(contributions_plot, num_station){

  contributions_plot$Type <- factor(contributions_plot$Type, levels=c("SiO4", "PO4", "DIN", "Temperature", "PAR", "Zooplankton"))

  graph_contributions <- ggplot(data = contributions_plot) +
    geom_area( aes(x=as.Date(Date), y = percentage, fill = Type)) + 
    scale_fill_viridis(discrete = T, direction = -1) +
    #scale_colour_grey(aesthetics = c("fill"), start = 0.8, end = 0.2)+
    scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=9), axis.text.x=element_text(angle=90, hjust=1)) +
    labs(x = "Date", y= "Relative contribution", 
         title = paste0("Relative contributions station: ", num_station)) +
    theme(plot.title = element_text(size=12, face="bold"))
  
  return(graph_contributions)
  
}

createContributions_line <- function(contributions_plot, num_station){
  
  contributions_plot$Type <- factor(contributions_plot$Type, levels=c("SiO4", "PO4", "DIN", "Temperature", "PAR", "Zooplankton"))
  
  graph_contributions <- ggplot(data = contributions_plot) +
    geom_line( aes(x=as.Date(Date), y = percentage, colour = Type)) + 
    scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=9), axis.text.x=element_text(angle=90, hjust=1)) +
    labs(x = "Date", y= "Relative contribution", 
         title = paste0("Relative contributions station: ", num_station)) +
    theme(plot.title = element_text(size=12, face="bold"))
  
  return(graph_contributions)
  
}

compareContributions_spatial <- function(){
  #Load contributions aggregated per station

  contributions_plot_130 <- read.csv("~/workspace/NPZ/Output/Final results/Relative contributions/contributionsNormalized_month_station130.csv")
  contributions_plot_130$Station <- "130"
  contributions_plot_330 <- read.csv("~/workspace/NPZ/Output/Final results/Relative contributions/contributionsNormalized_month_station330.csv")
  contributions_plot_330$Station <- "330"
  contributions_plot_offshore <- read.csv("~/workspace/NPZ/Output/Final results/Relative contributions/contributionsNormalized_month_stationoffshore.csv")
  contributions_plot_offshore$Station <- "offshore"
  
  combined_iter <- rbind(contributions_plot_130,contributions_plot_330, contributions_plot_offshore)
  combined_iter$Station <- as.factor(combined_iter$Station)
  
  relativeContributions_station130 <- read.csv("workspace/NPZD_model_SP/output/relativeContributions_smooth/relativeContributions_month_station130.csv")
  relativeContributions_station330 <- read.csv("workspace/NPZD_model_SP/output/relativeContributions_smooth/relativeContributions_month_station330.csv")
  relativeContributions_offshore <- read.csv("workspace/NPZD_model_SP/output/relativeContributions_smooth/relativeContributions_month_stationoffshore.csv")
  
  relativeContributions_station130$Station <- "130"
  relativeContributions_station330$Station <- "330"
  relativeContributions_offshore$Station <- "offshore"
  
  all_contributions <- rbind(relativeContributions_station130, relativeContributions_station330, relativeContributions_offshore)
  
  values_min <- aggregate( Contribution ~ Date + Station +Type, data = all_contributions, FUN = min )
  values_min$Limit <-"min"
  values_max <- aggregate( Contribution ~ Date + Station +Type, data = all_contributions, FUN = max )
  values_max$Limit <-"max"
  
  contributions_merged <- merge(values_min,values_max, by= c("Date", "Station", "Type"))
  
  contributions_merged <- contributions_merged[,c(1,2,3,4,6)]
  names(contributions_merged)[4] <- "Min"
  names(contributions_merged)[5] <- "Max"
  
  contributions_merged$Type <- factor(contributions_merged$Type, 
                                      levels = c('DIN', 
                                                 "PO4", 
                                                 "SiO4",
                                                 "PAR",
                                                 "Temperature",
                                                 "Zooplankton"))
  combined_iter$Type <- factor(combined_iter$Type, 
                               levels = c('DIN', 
                                          "PO4", 
                                          "SiO4",
                                          "PAR",
                                          "Temperature",
                                          "Zooplankton"))
  
  my_labeller <- as_labeller(c(DIN = 'DIN', PO4 = "PO[4]", SiO4 = "SiO[4]", PAR = "PAR", Temperature = "SST", Zooplankton = "Zooplankton~grazing"),
                             default = label_parsed)
  
  comparisonGraph <- ggplot() +
    geom_line(data = combined_iter, aes(x = as.Date(Date), y = percentage, group=Station, colour = Station), size = 0.8) + 
    geom_ribbon(data = contributions_merged, aes(x=as.Date(Date), ymin = Min, ymax = Max, group = Station, fill = Station),  alpha = 0.15) +
    ylim(0,0.5) +
    facet_wrap(vars(Type), scales = "free_y", labeller = labeller(Type =  my_labeller))+
    scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y", expand = c(0, 0)) +
    scale_colour_manual(name = "Region",labels= c("Nearshore","Midshore","Offshore"), 
                        values=c("#F8766D",'#00BA38','#619CFF'), breaks=c("130","330","offshore"))+
    scale_fill_manual(name = "Region",labels= c("Nearshore","Midshore","Offshore"), 
                        values=c("#F8766D",'#00BA38','#619CFF'), breaks=c("130","330","offshore"))+
    theme(text = element_text(size=7), axis.text.x=element_text(angle=90, hjust=1)) +
    theme(legend.position = "right",legend.direction = "vertical",legend.background = element_rect(),legend.title = element_text(size=14,face = 'bold'),
          legend.text = element_text(size=12),legend.text.align = 0,axis.text.x = element_text(colour="black",size=12),
          axis.text.y = element_text(colour="black",size=14), axis.title.x = element_text(colour="black",size=16,face = 'bold'),
          axis.title.y = element_text(colour="black",size=16,face = 'bold'), strip.text = element_text(size = 12))+
    labs(x = "Date", y= "Relative contribution")
  
  return(comparisonGraph)
  
}
