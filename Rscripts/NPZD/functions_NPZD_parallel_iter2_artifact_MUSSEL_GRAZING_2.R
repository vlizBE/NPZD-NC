
create_parameters <- function(numSim, parameters_values, id_param){
  
  parameters <- data.frame(matrix(nrow = numSim, ncol = 14))
  
  names(parameters)[1] <- paste0("ID_", id_param)
  names(parameters)[2] <- paste0("maxUptake_", id_param)
  names(parameters)[3] <- paste0("ksPAR_", id_param)
  names(parameters)[4] <- paste0("ksDIN_", id_param)
  names(parameters)[5] <- paste0("ksP_", id_param)
  names(parameters)[6] <- paste0("ksSi_", id_param)
  names(parameters)[7] <- paste0("maxGrazing_", id_param)
  names(parameters)[8] <- paste0("ksGrazing_", id_param)
  names(parameters)[9] <- paste0("pFaeces_", id_param)
  names(parameters)[10] <- paste0("excretionRate_", id_param)
  names(parameters)[11] <- paste0("mortalityRate_", id_param)
  names(parameters)[12] <- paste0("ChlNratio_", id_param)
  names(parameters)[13] <- paste0("Tobs_", id_param)
  names(parameters)[14] <- paste0("kd_" , id_param)
  
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
  
  simulation_output$Date <- as.Date(inputData$Date, format= "%m/%d/%Y")
  
  #initial values, taken from book
  PHYTO     <- 1.0
  ZOO       <- 0.1
  DIN       <- 5.0
  DETRITUS  <- 5.0
  
  #create average parameters 
  parameters_avg <- average_param(parameters1, parameters2)

  #validation exponential growth
  flagExp <- 0
  
  for(j in 1:times){
    
    #set time step in output
    simulation_output[j,1]<-j
    
    day_step <- yday(simulation_output$Date[j])
    
    #Use the parameters according to the date of the year using the solstice as reference periods
    #Day 172 corresponds to 21/06
    #Day 354 corresponds to 21/12
    
    #Reference periods for summer solstice
    #Day 157 corresponds to 7/06
    #Day 187 corresponds to 5/07
    
    #Reference periods for winter solstice
    #Day 340 corresponds to 5/12
    #Day 5 corresponds to 5/01
    
    if(day_step<157 & day_step>5){
      parameters <- parameters1
      #print(paste0("day: ", day_step, " parameters1"))
    }else {
      if(day_step>=157 & day_step<=187){
        parameters <- parameters_avg
       # print(paste0("day: ", day_step, " parameters average"))
      }else{
        if(day_step>187 & day_step<340){
          parameters <- parameters2
        #  print(paste0("day: ", day_step, " parameters2"))
        }else{
          parameters <- parameters_avg
         # print(paste0("day: ", day_step, " parameters average"))
        }
      }
    }
    
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
    
    ###addition mussel grazing ###
   if (times >1095) {
     if(PHYTO > (8/parameters[12])){    #8ug chla/L = 8mg/m3 --> 8mmolN/m3 //  mgChl/mmolN
      Grazing <- parameters[7] * PHYTOlim * ZOO
    }else {
      Grazing <- (parameters[7] * PHYTOlim * ZOO)+2 #40g chla/h + mosselen 24h actief
    }
     } else {Grazing <- parameters[7] * PHYTOlim * ZOO}
   
     ###
    
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
  nameFile <- paste("~/workspace/NPZ/Output/Final results/NPZ/station", int_station, "_iter2/detailed_simulation_", id_simulation, ".csv", sep = "")
  
  write.csv(simulation_output, file = nameFile)
  
  return(simulation_output)
  
}

average_param <- function(parameters1, parameters2){
  
  parameters <- data.frame(matrix(nrow = 3, ncol = 14))
  id_param <- "avg"
  
  names(parameters)[1] <- paste0("ID_", id_param)
  names(parameters)[2] <- paste0("maxUptake_", id_param)
  names(parameters)[3] <- paste0("ksPAR_", id_param)
  names(parameters)[4] <- paste0("ksDIN_", id_param)
  names(parameters)[5] <- paste0("ksP_", id_param)
  names(parameters)[6] <- paste0("ksSi_", id_param)
  names(parameters)[7] <- paste0("maxGrazing_", id_param)
  names(parameters)[8] <- paste0("ksGrazing_", id_param)
  names(parameters)[9] <- paste0("pFaeces_", id_param)
  names(parameters)[10] <- paste0("excretionRate_", id_param)
  names(parameters)[11] <- paste0("mortalityRate_", id_param)
  names(parameters)[12] <- paste0("ChlNratio_", id_param)
  names(parameters)[13] <- paste0("Tobs_", id_param)
  names(parameters)[14] <- paste0("kd_" , id_param)
  
  
  parameters[1,] <- parameters1
  parameters[2,] <- parameters2

  parameters[3,] <- apply(parameters[c(1,2),], 2, FUN = mean)
  
  return(parameters[3,])
  
}