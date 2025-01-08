
create_parameters <- function(numSim, parameters_values, id_param){
  
  parameters <- data.frame(matrix(nrow = numSim, ncol = 14))
  
  names(parameters)[1] <- "ID"
  names(parameters)[2] <- "maxUptake"
  names(parameters)[3] <- "ksPAR"
  names(parameters)[4] <- "ksDIN"
  names(parameters)[5] <- "ksP"
  names(parameters)[6] <- "ksSi"
  names(parameters)[7] <- "maxGrazing"
  names(parameters)[8] <- "ksGrazing"
  names(parameters)[9] <- "pFaeces"
  names(parameters)[10] <- "excretionRate"
  names(parameters)[11] <- "mortalityRate"
  names(parameters)[12] <- "ChlNratio"
  names(parameters)[13] <- "Tobs"
  names(parameters)[14] <- "kd"
  
  parameters[1] <- c(1:numSim)                        #ID set of parameters
  parameters[2] <- runif(numSim, min=parameters_values$Min[1], max=parameters_values$Max[1])        # MaxUptake: /day
  parameters[3] <- floor(runif(numSim, min=parameters_values$Min[2], max=parameters_values$Max[2]))      # ksPAR: ?Einst/m2/s
  parameters[4] <- runif(numSim, min=parameters_values$Min[3], max=parameters_values$Max[3])       # ksDIN: mmolN/m3
  parameters[5] <- runif(numSim, min=parameters_values$Min[4], max=parameters_values$Max[4])        # ksP
  parameters[6] <- runif(numSim, min=parameters_values$Min[5], max=parameters_values$Max[5])       # ksSi: Lancelot et al. (2005)
  parameters[7] <- runif(numSim, min=parameters_values$Min[6], max=parameters_values$Max[6])        # maxGrazing: /day
  parameters[8] <- runif(numSim, min=parameters_values$Min[7], max=parameters_values$Max[7])        # ksGrazing: mmolN/m3
  parameters[9] <- runif(numSim, min=parameters_values$Min[8], max=parameters_values$Max[8])        # pFaeces: %
  parameters[10] <- runif(numSim, min=parameters_values$Min[9], max=parameters_values$Max[9])       # excretionRate: /day
  parameters[11] <- runif(numSim, min=parameters_values$Min[10], max=parameters_values$Max[10])       # mortalityRate: /(mmolN/m3)/day
  parameters[12] <- runif(numSim, min=parameters_values$Min[11], max=parameters_values$Max[11])      # ChlNratio: mgChl/mmolN
  parameters[13] <- runif(numSim, min=parameters_values$Min[12], max=parameters_values$Max[12])      # Tobs and Lancelot et al. 2005
  parameters[14] <- runif(numSim, min=parameters_values$Min[13], max=parameters_values$Max[13])     #Kd based on quantiles 25% and 75%
  
  return(parameters)
}

npzd_run <- function(time_steps, parameters1, parameters2, int_station, inputData){
  
  times <- time_steps

  #simulation output
  simulation_output <- data.frame(matrix(nrow = times, ncol = 18))
  
  names(simulation_output)[1] <- "time_step"
  names(simulation_output)[2] <- "par_lim"
  names(simulation_output)[3] <- "din_lim"
  names(simulation_output)[4] <- "p_lim"
  names(simulation_output)[5] <- "si_lim"
  names(simulation_output)[6] <- "temp_lim"
  names(simulation_output)[7] <- "nuptake"
  names(simulation_output)[8] <- "grazing_lim"
  names(simulation_output)[9] <- "dphyto"
  names(simulation_output)[10] <- "dzoo"
  names(simulation_output)[11] <- "phyto"
  names(simulation_output)[12] <- "zoo"
  names(simulation_output)[13] <- "din_sim"
  names(simulation_output)[14] <- "din_gam"
  names(simulation_output)[15] <- "ddin"
  names(simulation_output)[16] <- "chl_a"
  names(simulation_output)[17] <- "detritus"  
  names(simulation_output)[18] <- "Date"
  
  simulation_output$Date <- as.Date(inputData$Date)
  
  #initial values, taken from book
  PHYTO     <- 1.0
  ZOO       <- 0.1
  DIN       <- 5.0
  DETRITUS  <- 5.0
  
  #validation exponential growth
  flagExp <- 0
  
  
  for(j in 1:times){
    
    #set time step in output
    simulation_output[j,"time_step"]<-j
    
    day_step <- yday(simulation_output$Date[j])
    
    #Day 172 of the year corresponds to 21/06 
    if(day_step<=172){
      parameters <- parameters1
    }else {
      parameters<- parameters2
    }
    
    #Set limitation factors and effects
    par_0 <- 0.5*(540+440*sin(2*pi*j/365-1.4))
    PAR <- par_0*exp(-parameters$kd*3)    #Nutrient data is sampled at 3m depth
    PARlim <- PAR/(PAR+parameters$ksPAR)
    simulation_output[j,"par_lim"]<-PARlim
    
    #Replace zero and negative values based on the minimum level of detection (NH4, NO2, NOx)
    minDIN <- runif(1, min=0.0135, max=0.027)
    DIN <- max(minDIN,inputData$DIN[j])
    DINlim <- DIN /(DIN+parameters$ksDIN)
    simulation_output[j,"din_lim"]<-DINlim
    
    P <- inputData$po4[j]
    Plim <- P/(P + parameters$ksP)
    simulation_output[j,"p_lim"]<-Plim
    
    Si <- inputData$sio4[j]
    Silim <- Si / (Si + parameters$ksSi)
    simulation_output[j,"si_lim"]<-Silim
    
    Temp <- inputData$Temp[j]
    theta <- 1.185 - 0.00729*Temp
    Templim <- theta^(Temp - parameters$Tobs)
    simulation_output[j,"temp_lim"]<-Templim
    
    #Calculate values of the variables
    
    Nuptake        <- parameters$maxUptake * PARlim * DINlim * Templim * Plim * Silim * PHYTO
    simulation_output[j,"nuptake"]<-Nuptake
    
    PHYTOlim <-PHYTO/(PHYTO+parameters$ksGrazing)
    Grazing        <- parameters$maxGrazing * PHYTOlim * ZOO
    simulation_output[j,"grazing_lim"]<- PHYTOlim           #limiting factor grazing
    
    Faeces         <- parameters$pFaeces * Grazing
    Excretion      <- parameters$excretionRate * ZOO
    Mortality      <- parameters$mortalityRate * ZOO * ZOO
    Mineralisation <- 0.1 * DETRITUS
    
    dPHYTO    <- Nuptake - Grazing
    simulation_output[j,"dphyto"]<-dPHYTO
    
    dZOO      <- Grazing - Faeces - Excretion - Mortality
    simulation_output[j,"dzoo"]<-dZOO
    
    dDIN      <- Mineralisation + Excretion - Nuptake
    simulation_output[j,"ddin"]<-dDIN
    
    dDETRITUS <- Mortality + Faeces -  Mineralisation
    
    PHYTO     <- PHYTO + dPHYTO
    simulation_output[j,"phyto"]<-PHYTO
    
    ZOO       <- ZOO + dZOO
    simulation_output[j,"zoo"]<-ZOO
    
    DIN       <- DIN + dDIN
    simulation_output[j,"din_sim"]<-DIN
    simulation_output[j,"din_gam"]<-max(minDIN,inputData$DIN[j])
    
    DETRITUS  <- DETRITUS + dDETRITUS
    simulation_output[j,"detritus"] <- DETRITUS
    
    Chlorophyll    <- parameters$ChlNratio * PHYTO
    simulation_output[j,"chl_a"]<-Chlorophyll
    
    if(PHYTO > 10000 | PHYTO < -10000){
      flagExp <- j
      break
    }
  }
  
  #If phyto grows exponentially, the simulation will have an error code
  if(flagExp>0){
    simulation_output[flagExp:times,11] <- 99999
  }
  
  id_simulation<-parameters[1,"ID"]
  nameFile <- paste("~/workspace/NPZ/Output/Final results/NPZ/station", int_station, "_iter2/test/detailed_simulation_", id_simulation, ".csv", sep = "")
  
  write.csv(simulation_output, file = nameFile)
  
  return(NA)
  
}