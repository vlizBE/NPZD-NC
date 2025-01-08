library(mgcv)
library(lubridate)
library(ggplot2)
library(stats)
library(xts)
library(plyr)

##################
##Functions for GAM Models
# functions to load packages that are needed for NPZ model (Steven Pint)

install_needed_packages_GAM <- function(){ 
  not_installed <- setdiff(c("mgcv",
                             "lubridate",
                             "ggplot2",
                             "stats",
                             "xts",
                             "lwdataexplorer",
                             "plyr",
                             "dplyr",
                             "devtools")
                           , rownames(installed.packages()))
  
  if (length(not_installed) != 0) {
    if ("lwdataexplorer" %in% not_installed){
      install.packages(not_installed[which(not_installed != "lwdataexplorer")])
      devtools::install_github("lifewatch/lwdataexplorer")
      } else {
    install.packages(not_installed)
      }
  } else {
    print("All packages that are needed are installed")
  }
}


test_model <-function(y_values, day, year, k_day, k_year){
  
  model_1 <- gam(y_values ~ s(day, k = k_day) + s(year, k = k_year), method = "REML")
  
  #print(summary(model_1))
  #plot(model_1, pages=1, all.terms = TRUE, residuals = TRUE)
  #print(gam.check(model_1))
  
  #print(paste("AIC:", AIC(model_1)))
  
  return(model_1)
}

test_model_gamma <-function(y_values, day, year, k_day, k_year, g1){
  
  model_1 <- gam(y_values ~ s(day, k = k_day) + s(year, k = k_year), method = "REML", gamma = g1)
  
  #print(summary(model_1))
  #plot(model_1, pages=1, all.terms = TRUE, residuals = TRUE)
  #print(gam.check(model_1))
  
  print(paste("AIC:", AIC(model_1)))
  
  return(model_1)
}

comparisonPoint <- function(day, year, y_value, model_base) {
  
  compare<- data.frame(day, year, y_value)
  
  names(compare)[1] <- "day"
  names(compare)[2] <- "year"
  names(compare)[3] <- "observation"
  
  compare <- compare[which(!is.na(compare$observation)),]
  
  pred<-predict.gam(model_base, newdata = compare[,c(1,2)], type = "response")
  
  compare$pred<-pred
  
  ##################################
  #error measurement
  
  rmse_value <- sqrt(mean((compare$pred - compare$observation)^2))
  print(paste("RMSE model: ",rmse_value))
  
  mean_obs <- mean(compare$observation)
  num1 <- sum((compare$pred-mean_obs)^2)
  den1 <- sum((compare$observation-mean_obs)^2)
  coefficient_determination <- num1 / den1
  print(paste("Coefficient of determination: ", coefficient_determination))
  
  ###########################
  #Create graph
  graph <- ggplot(data=compare) +
    geom_point(aes(x=observation, y=pred)) + 
    geom_abline(slope = 1, intercept = 0, color = "blue") +
    labs(x = "Observed value", y = "Predicted value") 
  
  return(graph)
}

comparisonTimeSeries <- function (predicted, observed){
  
  graph_DIN<-ggplot() + 
    geom_line(aes(x=as.Date(Date) , y=DIN), data = predicted) +
    geom_point(aes(x=as.Date(Date), y=DIN), data = observed, colour = "red", size = 1) +
    scale_x_date(date_breaks = "6 month", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=11), axis.text.x=element_text(angle=90, hjust=1)) +
    labs(x = "Date", y= "DIN (mmol N/m3)") 
  
  graph_po4<-ggplot() + 
    geom_line(aes(x=as.Date(Date) , y=po4), data = predicted) +
    geom_point(aes(x=as.Date(Date), y=PO4), data = observed, colour = "red", size = 1) +
    scale_x_date(date_breaks = "6 month", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=11), axis.text.x=element_text(angle=90, hjust=1)) +
    labs(x = "Date", y= "PO4 (mmol P/m3)") 
  
  graph_sio4<-ggplot() + 
    geom_line(aes(x=as.Date(Date) , y=sio4), data = predicted) +
    geom_point(aes(x=as.Date(Date), y=SiO4), data = observed, colour = "red", size = 1) +
    scale_x_date(date_breaks = "6 month", date_labels =  "%b %Y", expand = c(0, 0)) +
    theme(text = element_text(size=11), axis.text.x=element_text(angle=90, hjust=1)) +
    labs(x = "Date", y= "SiO4 (mmol Si/m3)") 
  
  listgraphs <- list(graph_DIN, graph_po4, graph_sio4)
  
  return(listgraphs)
  
}

compareGraphs <- function(data_obs, data_sim, number_of_days){
  
  data_obs$type <- "Observed"
  data_sim$type <- "GAM"
  
  data_obs2 <- as.data.frame(matrix(nrow=nrow(data_obs), ncol = 6))
  names(data_obs2)[1] <- "Date"
  names(data_obs2)[2] <- "PO4"
  names(data_obs2)[3] <- "SiO4"
  names(data_obs2)[4] <- "month"
  names(data_obs2)[5] <- "DIN"
  names(data_obs2)[6] <- "type"
  
  data_obs2$Date <- data_obs$Date
  data_obs2$PO4 <- data_obs$PO4
  data_obs2$SiO4 <- data_obs$SiO4
  data_obs2$month <- data_obs$month
  data_obs2$DIN <- data_obs$DIN
  data_obs2$type <- data_obs$type
  
  data_obs2$month <- as.factor(data_obs2$month)
  data_obs2$Date <- as.Date(data_obs2$Date)
  
  data_sim2 <- as.data.frame(matrix(nrow=number_of_days, ncol = 6))
  names(data_sim2)[1] <- "Date"
  names(data_sim2)[2] <- "PO4"
  names(data_sim2)[3] <- "SiO4"
  names(data_sim2)[4] <- "month"
  names(data_sim2)[5] <- "DIN"
  names(data_sim2)[6] <- "type"
  
  data_sim2$Date <- data_sim$Date
  data_sim2$PO4 <- data_sim$po4
  data_sim2$SiO4 <- data_sim$sio4
  data_sim2$month <- as.factor(data_sim$month)
  data_sim2$DIN <- data_sim$DIN
  data_sim2$type <- data_sim$type
  
  all_data_graph <- rbind(data_obs2, data_sim2)
  all_data_graph$month <- factor(all_data_graph$month, levels=c(1:12))
  
  #comparison graph DIN
  
  #Calculate the number of observations per month 
  data_DIN <- all_data_graph[!is.na(all_data_graph$DIN),]
  din_obs<-data_DIN[which(data_DIN$type == "Observed"),]
  stats_din <- count(din_obs, "month")
  stats_din <- getStats(stats_din)
  print(paste("total DIN observations: ",sum(count(din_obs, "month")$freq)))
  
  graphDIN <- ggplot(data = all_data_graph, aes(x=month, y=DIN, fill = type)) +
    geom_boxplot() + scale_x_discrete(labels=c("1" = paste("January (n = ", stats_din$freq[1],")"), "2" = paste("February (n = ", stats_din$freq[2],")"),"3" = paste("March (n = ", stats_din$freq[3],")"), 
                                               "4" = paste("April (n = ", stats_din$freq[4],")"),"5" = paste("May (n = ", stats_din$freq[5],")"),
                                               "6" = paste("June (n = ", stats_din$freq[6],")"), "7" = paste("July (n = ", stats_din$freq[7],")"),
                                               "8" = paste("August (n = ", stats_din$freq[8],")"), "9" = paste("September (n = ", stats_din$freq[9],")"), 
                                               "10" = paste("October (n = ", stats_din$freq[10],")"), "11" = paste("November (n = ", stats_din$freq[11],")"), "12" = paste("December (n = ", stats_din$freq[12],")"))) +
    theme(axis.text.x = element_text(angle=90))  + labs(y="DIN (m_moles_N/m3)")
  
  #comparison graph for PO4
  po4_obs<-data_obs2[!is.na(data_obs2$PO4),c(2,4)]
  stats_po4 <- count(po4_obs, "month")
  stats_po4 <- getStats(stats_po4)
  print(paste("total PO4 observations: ",sum(count(po4_obs, "month")$freq)))
  
  graphP <- ggplot(data = all_data_graph, aes(x=month, y=PO4, fill = type)) +
    geom_boxplot() + scale_x_discrete(labels=c("1" = paste("January (n = ", stats_po4$freq[1],")"), "2" = paste("February (n = ", stats_po4$freq[2],")"),"3" = paste("March (n = ", stats_po4$freq[3],")"), 
                                               "4" = paste("April (n = ", stats_po4$freq[4],")"),"5" = paste("May (n = ", stats_po4$freq[5],")"),
                                               "6" = paste("June (n = ", stats_po4$freq[6],")"), "7" = paste("July (n = ", stats_po4$freq[7],")"),
                                               "8" = paste("August (n = ", stats_po4$freq[8],")"), "9" = paste("September (n = ", stats_po4$freq[9],")"), 
                                               "10" = paste("October (n = ", stats_po4$freq[10],")"), "11" = paste("November (n = ", stats_po4$freq[11],")"), "12" = paste("December (n = ", stats_po4$freq[12],")"))) +
    theme(axis.text.x = element_text(angle=90)) + labs(y="PO4 (m_moles_P/m3)")
  
  
  #Comparison graph for SiO4
  si_obs<-data_obs2[!is.na(data_obs2$SiO4),c(3,4)]
  stats_si <- count(si_obs, "month")
  stats_si <-getStats(stats_si)
  print(paste("total SiO4 observations: ",sum(count(si_obs, "month")$freq)))
  
  graphSi <- ggplot(data = all_data_graph, aes(x=month, y=SiO4, fill = type)) +
    geom_boxplot() + scale_x_discrete(labels=c("1" = paste("January (n = ", stats_si$freq[1],")"), "2" = paste("February (n = ", stats_si$freq[2],")"),"3" = paste("March (n = ", stats_si$freq[3],")"), 
                                               "4" = paste("April (n = ", stats_si$freq[4],")"),"5" = paste("May (n = ", stats_si$freq[5],")"),
                                               "6" = paste("June (n = ", stats_si$freq[6],")"), "7" = paste("July (n = ", stats_si$freq[7],")"),
                                               "8" = paste("August (n = ", stats_si$freq[8],")"), "9" = paste("September (n = ", stats_si$freq[9],")"), 
                                               "10" = paste("October (n = ", stats_si$freq[10],")"), "11" = paste("November (n = ", stats_si$freq[11],")"), "12" = paste("December (n = ", stats_si$freq[12],")"))) +
    theme(axis.text.x = element_text(angle=90)) + labs(y="SiO4 (m_moles_Si/m3)")
  
  listgraphs <- list(graphDIN, graphP, graphSi) 
  
  return(listgraphs)
  
}

getStats <- function(temp_stats){
  stats_obs <- data.frame(matrix(nrow = 12, ncol = 2))
  names(stats_obs)[1] <- "month"
  names(stats_obs)[2] <- "freq"

  #in case there are no observations in a given month, this for makes sure that there is 0 registered in the dataframe
  #this has to be made manually as the count function in R does not incorporate this
  for (m in 1:12){
    stats_obs[m,1] <- m
    temp <- temp_stats[which(temp_stats$month == m),2]
    stats_obs[m,2] <- ifelse(length(temp)==0, 0, temp)
  }
  
  return(stats_obs)

}
