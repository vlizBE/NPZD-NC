
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

npzd_run <- function(time_steps, parameters1, parameters2, int_station, inputData, station_lat, station_lon, sample_depth){

  times <- time_steps

  #simulation output
  simulation_output <- data.frame(matrix(nrow = times, ncol = 27))
  
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
  #names(simulation_output)[18] <- "grazingC_lim"  
  names(simulation_output)[18] <- "dphytoc"  
  names(simulation_output)[19] <- "dzooc"  
  names(simulation_output)[20] <- "ddic" 
  names(simulation_output)[21] <- "zooc"  
  names(simulation_output)[22] <- "dic_sim"  
  #names(simulation_output)[24] <- "dic_gam"  
  names(simulation_output)[23] <- "detritusC"  
  names(simulation_output)[24] <- "cuptake"  
  names(simulation_output)[25] <- "phytoc"
  names(simulation_output)[26] <- "pCO2w"
  names(simulation_output)[27] <- "Cflux"

  simulation_output$Date <- as.Date(inputData$Date)
  
  #initial values, taken from book
  redfield <- 6.625
  PHYTO     <- 1.0
  PHYTOC     <- PHYTO*redfield
  ZOO       <- 0.1
  ZOOC      <- ZOO * 4.5 #qhet ?=? CNzooconst = 4.5
  DIN       <- 5.0
  DIC       <- 2200 # 2200 umol/l = mmol/m³
  DETRITUS  <- 5.0
  DETRITUSC  <- 5.0 * redfield # niet de juiste conversie factor (Redfield)
  #resDOC    <- DIC * 0.03 # DOC-DIC ratio only 2-4%  (https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2021GL092904#:~:text=The%20concentration%20of%20dissolved%20organic,sinks%20than%20that%20for%20DIC.)
  #resdoc kan eruit -> eventueel wel
  ChlCratio  <-  0.055 ## origially 0.1
  
  #create average parameters 
  parameters_avg <- average_param(parameters1, parameters2)

  #validation exponential growth
  flagExp <- 0
  
  for(j in 1:times){
    
    #set time step in output
    simulation_output[j,"time_step"]<-j
    
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
    par_0 <- mean(na.omit(solrad::OpenRadiation(seq(j, (j+1)-1/24,1/24), Lat = station_lat, Lon = station_lon, SLon = station_lon, DS=0, Elevation = 0))) # solrad package, previously -> 0.5*(540+440*sin(2*pi*j/365-1.4))
    PAR <- par_0*exp(-parameters$kd*sample_depth)    #Nutrient data is sampled at ?m depth
    PARlim <- PAR/(PAR+parameters$ksPAR)
    simulation_output[j,"par_lim"]<-PARlim
    
    #Replace zero and negative values based on the minimum level of detection (NH4, NO2, NOx)
    minDIN <- runif(1, min=0.0135, max=0.027)
    DIN <- max(minDIN,inputData$DIN[j]) #DIN OP BASIS VAN MODEL ZETTEN?? BEKIJKEN OF DIT GEDAAN MODEL BEGRENZEN? MOMENTEEL MODEL TREKT WAARNEMINGEN 
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
    # temp info for carbon
    Tref <- 288.15
    TempK <- Temp + 273.15                 # temperature °C to Kelvin 
    Tf <- exp(-4500 * ((1/TempK) - (1/Tref))) # T in Kelvin - slope of Arrhenius relation: 4500K
    # 1 temp gebruiken?
    
    # CARBON FLUX CALCULATIONS check hoe te implementeren
    # Calculate pCO2w from simulated DIC using seacarb package
    # flag 9:  var 1 pH and var 2 DIC (this can be in situ pH and modelled DIC)
    # S =  in situ salinity
    # T =  in situ temperature
    # Patm = 1 default setting
    # long & lat can be filled in using coord of the station, but no difference observed
    # P Hydrostatic pressure in bar (surface = 0)
    # Pt Concentration of total phosphate in mol/kg; set to 0 if NA: can use PO4 concentration from input
    # Sit Concentration of total silicate in mol/kg; set to 0 if NA: can use SiO4 concentration from input
    carb_df <- seacarb::carb(flag = 9, var1 = carbon_inputData$ph[j], var2=DIC*10^-6/1.025, S=carbon_inputData$sal[j], T=Temp, Patm=1, P=0, Pt=P, Sit=Si,
    k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential",
    warn="y", eos="eos80", long=1.e20, lat=1.e20)
   
    FC <- FC_buoy_validation(u10=carbon_inputData$wind[j],Temp=Temp,Sal=carbon_inputData$sal[j],pCO2W=carb_df$pCO2,pCO2A=carbon_inputData$pco2atm[j])

    simulation_output[j,"pCO2w"] <- carb_df$pCO2
    simulation_output[j,"Cflux"] <- FC

    #marelac package bekijken
    #Calculate values of the variables
    
    Nuptake        <- parameters$maxUptake * PARlim * DINlim * Templim * Plim * Silim * PHYTO
    simulation_output[j,"nuptake"] <- Nuptake
    
    CNratio <- PHYTOC/PHYTO
    Rphot <- ((1/CNratio) - 0.043) / (0.171 - 0.043)  #CNratio -> cellullar molar? N:C=16:106 fixed or variable?
    umaxc <- 2.3 * Rphot * Templim#Tf
    # parameter for C specific photosynthesis (uc): 0.1-6.4 d^-1 -> best 2.3
    Cuptake         <- umaxc * (1 - exp((-ChlCratio * 0.3 * PAR) / umaxc))
    
    #chla specific photosynthetic efficiency (alpha): 0.3 mmol C (mg chla)^-1 m² W^-1 d^-1
    #thetac is chl:C ratio phyto
    #PAR is included in the equation instead of I (Schartau) 
    simulation_output[j,"cuptake"] <- Cuptake
    
    PHYTOlim <-PHYTO/(PHYTO+parameters$ksGrazing)
    
    Grazing        <- parameters$maxGrazing * PHYTOlim * ZOO
    GrazingC        <- (parameters$maxGrazing * PHYTOlim * ZOO) / (1/CNratio) # #oorspronklijk was het delen door NCratio
    simulation_output[j,"grazing_lim"]<- PHYTOlim           #limiting factor grazing
    #simulation_output[j,"grazingC_lim"]<- PHYTOClim           #limiting factor grazing  
    
    Faeces         <- parameters$pFaeces * Grazing
    FaecesC         <- parameters$pFaeces * GrazingC # correct? or just use the CNratio?
    
    Excretion      <- parameters$excretionRate * ZOO
    
    ExcretionCmin  = 0.01 *Templim* ZOOC   #Tf      # respiration rate ZOO: 0.01 d^-1 ------- from Tom's model
    homeostasis = 0 #GrazingC * (CNratio - 4.5)  #!!in Tom's script it's CNratioPhy - CNratioHet     # instantaneous loss of excess carbon
    if (homeostasis >= 0){  
      respZOOExcessC = homeostasis
    } else {
      respZOOExcessC  = 0.00
    }
    ExcretionC <- ExcretionCmin + respZOOExcessC
    #ExcretionC <- parameters$excretionRate * ZOO * 4.5
    
    Mortality      <- parameters$mortalityRate * ZOO * ZOO
    MortalityC      <- Mortality * 4.5  #parameters$mortalityRate * ZOOC * ZOOC # CNratio zoo = 4.5
    
    Mineralisation <- 0.1 * DETRITUS
    MineralisationC <- 0.004 * Templim * DETRITUSC # Tf # remineralisation of C (d^-1): best value 0.004 mean 0.018 std 0.019
    
    dPHYTO    <- Nuptake - Grazing
    simulation_output[j,"dphyto"]<-dPHYTO
    
    dPHYTOC    <- Cuptake - GrazingC
    simulation_output[j,"dphytoc"]<-dPHYTOC
    
    dZOO      <- Grazing - Faeces - Excretion - Mortality
    simulation_output[j,"dzoo"]<-dZOO
    
    dZOOC     <- GrazingC - FaecesC - ExcretionC - MortalityC
    simulation_output[j,"dzooc"]<-dZOOC
    
    dDIN      <- Mineralisation + Excretion - Nuptake
    simulation_output[j,"ddin"]<-dDIN
    
    #dDIC    <- - (Cuptake - Calcif)*PHYTOC  + respZOOC + MineralisationC + FC # geprobeert ovet te nemen van Tom's model
    #resDOC <- 0.250 * (1 - 0.64) * PHYTOC + 0.001 * ZOOC + MineralisationC - 0.372 * Tf * resDOC
    # phyto loss carbon rate Yc = 0.250; polysaccharide fraction of total DOC exdates fpcho = 0.64
    # loss rate hetrotrophs Yhet = 0.001 d^-1; remineralisation of resDOC pc = 0.372 d^-1
    
    Rnc <- 1 - exp(-1000 *(abs((1/CNratio)-0.171)-((1/CNratio)-0.171))^2) 
    # slope parameter for DIN uptake regulation (sigmanc): 1000 mmol N²/ mmol C²
    
    Vnc <- umaxc * 0.171 * Rnc * DINlim
    # max cellular CN ratio qmax 0.171
    
    rphy <- 0.01 + 2.3 * Vnc
    # maintanance respiration rate (rc) 0.01 d^-1; biosynthetic cost (biosyncost) 2.3 mmol C/mmol N
    pc <- 0.372
    dDIC <- (rphy - Cuptake) * PHYTOC + ExcretionC + MineralisationC + FC  # respiration_phy, carbon_phtosynthesis, carbon_phyto, respiration_zoo, carbon_zoo, polysacchariden, temperature dependent,residual carbon, carbon flux
    #remineralisation of resDOC pc = 0.372 d^-1
    #replace  pc * Tf * resDOC by just MineralisationC similar to dDIC?
    
    simulation_output[j,"ddic"]<-dDIC
    
    dDETRITUS <- Mortality + Faeces -  Mineralisation
    
    dDETRITUSC <- MortalityC + FaecesC -  MineralisationC
    
    PHYTO     <- PHYTO + dPHYTO
    simulation_output[j,"phyto"]<-PHYTO
    
    PHYTOC    <- PHYTOC + dPHYTOC
    simulation_output[j,"phytoc"]<-PHYTOC
    
    ZOO       <- ZOO + dZOO
    simulation_output[j,"zoo"]<-ZOO
    
    ZOOC      <- ZOOC + dZOOC
    simulation_output[j,"zooc"]<-ZOOC
    
    DIN       <- DIN + dDIN
    simulation_output[j,"din_sim"]<-DIN
    simulation_output[j,"din_gam"]<-max(minDIN,inputData$DIN[j])
    
    DIC       <- DIC + dDIC                         
    simulation_output[j,"dic_sim"]<-DIC
    #simulation_output[j,"dic_gam"]<-max(minDIC,inputData$DIC[j])
    
    DETRITUS  <- DETRITUS + dDETRITUS
    simulation_output[j,"detritus"] <- DETRITUS
    
    DETRITUSC  <- DETRITUSC + dDETRITUSC
    simulation_output[j,"detritusC"] <- DETRITUSC  
    
    Chlorophyll    <- parameters$ChlNratio * PHYTO
    simulation_output[j,"chl_a"]<-Chlorophyll

          
    if(PHYTO > 10000 | PHYTO < -10000){
      flagExp <- j
      break
    }
    
    #if(ZOOC > 10000 | ZOOC < 0){
    #  flagExp <- j
    #  break
    #}
    
    
  }
  
  #If phyto grows exponentially, the simulation will have an error code
  if(flagExp>0){
    simulation_output[flagExp:times,11] <- 99999
  }
  
  id_simulation<-parameters[1,"ID"]
  nameFile <- paste(wd,"Output/Final results/NPZD/station", int_station, "_iter2/detailed_simulation_", id_simulation, ".csv", sep = "")
  
  write.csv(simulation_output, file = nameFile)
  
  return(NA)
  
}

average_param <- function(parameters1, parameters2){
  
  parameters <- data.frame(matrix(nrow = 3, ncol = 14))
  id_param <- "avg"
  
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
  
  
  parameters[1,] <- parameters1
  parameters[2,] <- parameters2

  parameters[3,] <- apply(parameters[c(1,2),], 2, FUN = mean)
  
  return(parameters[3,])
  
}