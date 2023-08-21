
createErrorFile <- function(numSim){
  errors_parameters <- data.frame(matrix(nrow = numSim, ncol = 19))
  names(errors_parameters)[1] <- "ID_sim"
  names(errors_parameters)[2] <- "rmse_phyto"
  names(errors_parameters)[3] <- "rmse_zoo"
  names(errors_parameters)[4] <- "total_rmse"
  names(errors_parameters)[5] <- "Rsquared"
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

create_parameters <- function(numSim, numStation){
  
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
  
  parameters_values <- read.csv(paste0("~/workspace/NPZ/Input data/NPZ/parameters_min_max_", numStation, ".csv"))
  
  parameters$ID <- c(1:numSim)                        #ID set of parameters
  parameters$maxUptake <- runif(numSim, min=parameters_values$Min[1], max=parameters_values$Max[1])        # MaxUptake: /day
  parameters$ksPAR <- floor(runif(numSim, min=parameters_values$Min[2], max=parameters_values$Max[2]))      # ksPAR: ?Einst/m2/s
  parameters$ksDIN <- runif(numSim, min=parameters_values$Min[3], max=parameters_values$Max[3])       # ksDIN: mmolN/m3
  parameters$ksP <- runif(numSim, min=parameters_values$Min[4], max=parameters_values$Max[4])        # ksP
  parameters$ksSi <- runif(numSim, min=parameters_values$Min[5], max=parameters_values$Max[5])       # ksSi: Lancelot et al. (2005)
  parameters$maxGrazing <- runif(numSim, min=parameters_values$Min[6], max=parameters_values$Max[6])        # maxGrazing: /day
  parameters$ksGrazing <- runif(numSim, min=parameters_values$Min[7], max=parameters_values$Max[7])        # ksGrazing: mmolN/m3
  parameters$pFaeces <- runif(numSim, min=parameters_values$Min[8], max=parameters_values$Max[8])        # pFaeces: %
  parameters$excretionRate <- runif(numSim, min=parameters_values$Min[9], max=parameters_values$Max[9])       # excretionRate: /day
  parameters$mortalityRate <- runif(numSim, min=parameters_values$Min[10], max=parameters_values$Max[10])       # mortalityRate: /(mmolN/m3)/day
  parameters$ChlNratio <- runif(numSim, min=parameters_values$Min[11], max=parameters_values$Max[11])      # ChlNratio: mgChl/mmolN
  parameters$Tobs <- runif(numSim, min=parameters_values$Min[12], max=parameters_values$Max[12])      # Tobs and Lancelot et al. 2005
  parameters$kd <- runif(numSim, min=parameters_values$Min[13], max=parameters_values$Max[13])     #Kd based on quantiles 25% and 75%
  
  return(parameters)
}

npzd_run <- function(time_steps, parameters, numStation, inputData){
  
  times <- time_steps
  
  #simulation output
  simulation_output <- data.frame(matrix(nrow = times, ncol = 17))
  
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
  
  #initial values, taken from book
  PHYTO     <- 1.0
  ZOO       <- 0.1
  DIN       <- 5.0
  DETRITUS  <- 5.0
  
  #validation exponential growth
  flagExp <- 0
  
  for(j in 1:times){
    
    #set time step in output
    simulation_output[j,1]<-j
    
    #Set limitation factors and effects
    par_0 <- 0.5*(540+440*sin(2*pi*j/365-1.4))
    PAR <- par_0*exp(-parameters[14]*3)    #Nutrient data is sampled at 3m depth
    PARlim <- PAR/(PAR+parameters[3])
    simulation_output[j,2]<-PARlim
    
    #Replace zero and negative values based on the minimum level of detection (NH4, NO2, NOx)
    minDIN <- runif(1, min=0.0135, max=0.027)
    DIN <- max(minDIN,inputData$DIN[j])
    DINlim <- DIN /(DIN+parameters[4])
    simulation_output[j,3]<-DINlim
    
    P <- inputData$po4[j]
    Plim <- P/(P + parameters[5])
    simulation_output[j,4]<-Plim
    
    Si <- inputData$sio4[j]
    Silim <- Si / (Si + parameters[6])
    simulation_output[j,5]<-Silim
    
    Temp <- inputData$Temp[j]
    theta <- 1.185 - 0.00729*Temp
    Templim <- theta^(Temp - parameters[13])
    simulation_output[j,6]<-Templim
    
    #Calculate values of the variables
    
    Nuptake        <- parameters[2] * PARlim * DINlim * Templim * Plim * Silim * PHYTO
    simulation_output[j,7]<-Nuptake
    
    PHYTOlim <-PHYTO/(PHYTO+parameters[8])
    Grazing        <- parameters[7] * PHYTOlim * ZOO
    simulation_output[j,8]<- PHYTOlim           #limiting factor grazing
    
    Faeces         <- parameters[9] * Grazing
    Excretion      <- parameters[10] * ZOO
    Mortality      <- parameters[11] * ZOO * ZOO
    Mineralisation <- 0.1 * DETRITUS
    
    dPHYTO    <- Nuptake - Grazing
    simulation_output[j,9]<-dPHYTO
    
    dZOO      <- Grazing - Faeces - Excretion - Mortality
    simulation_output[j,10]<-dZOO
    
    dDIN      <- Mineralisation + Excretion - Nuptake
    simulation_output[j,15]<-dDIN
    
    dDETRITUS <- Mortality + Faeces -  Mineralisation
    
    PHYTO     <- PHYTO + dPHYTO
    simulation_output[j,11]<-PHYTO
    
    ZOO       <- ZOO + dZOO
    simulation_output[j,12]<-ZOO
    
    DIN       <- DIN + dDIN
    simulation_output[j,13]<-DIN
    simulation_output[j,14]<-max(minDIN,inputData$DIN[j])
    
    DETRITUS  <- DETRITUS + dDETRITUS
    simulation_output[j,17] <- DETRITUS
    
    Chlorophyll    <- parameters[12] * PHYTO
    simulation_output[j,16]<-Chlorophyll
    
    if(PHYTO > 10000 | PHYTO < -10000){
      flagExp <- j
      break
    }
  }
  
  #If phyto grows exponentially, the simulation will have an error code
  if(flagExp>0){
    simulation_output[flagExp:times,11] <- 99999
  }
  
  id_simulation<-parameters[1,1]
  nameFile <- paste("~/workspace/NPZ/Output/Intermediate results/station", numStation, "_iter1/detailed_simulation_", id_simulation, ".csv", sep = "")
  
  write.csv(simulation_output, file = nameFile)
  
  return(parameters)
  
}
