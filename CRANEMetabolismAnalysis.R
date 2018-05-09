####################################################
### CRANE Data Analysis for Biogeochemical Responses####
### Created by Nyssa Silbiger                   ###
### Created on 5/03/2016                       ###
### Edited on 5/8/2018                               ###
### Edited by NJS                               ###
###################################################

# Clear workspace --------------------------------
rm(list=ls())

# Add Libraries ---------------------------------
library('plyr')
library('oce')
library('lme4')
library('lmerTest')
library('reshape2')
library('seacarb')
library('MuMIn')
library('effects')
library('RColorBrewer')

# Functions-----------------------------------------
#easy errorbar barplots
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

NECCalc<-function(HeaderTA,TankTA,ResidenceTime,SurfaceArea, TankVolume=5678,SWDenstiy=1.023, HeaderN=0, TankN=0, HeaderP=0, TankP=0){
 
   #NEC<-0.5*(HeaderTA-TankTA)*TankVolume*SWDenstiy/(ResidenceTime*SurfaceArea)/1000
   deltaN<-HeaderN-TankN # account for uptake of N
   deltaP<-HeaderP-TankP #account for updakte of P
   deltaTA<-(HeaderTA-TankTA)
   NEC<-0.5*(deltaTA-deltaN-2*deltaP)*(TankVolume/SurfaceArea)*(SWDenstiy/ResidenceTime)/1000
   
   
   return(NEC)
  #NEC calc calculated the net ecosystem calcification of a flow through mesocosm system at a given time point
  #this uses residence time (lagrangian) rather than change in TA over time (Eulerian-- this is what I used for the 
  #biogeochemistry paper)
  #NEC is in mmol cm^-2 hr-1  (or umol g-1 hr-1 if using one of the other measurements)
  
  #HeaderTA is TA from the header in umol/kg
  #TankTAis TA from the tank in umol/kg
  #Residence time is the residence time in hours
  #Surface area is SA of the substrate in cm2
  #TankVolume is the volume in cm3 = (default = 5678)
  #SWDensity is density of seawater in g/cm3 (default =1.023)
# divide by 1000 to make umol g-1 h-1
   
}  

NCPCalc<-function(HeaderDIC,TankDIC,ResidenceTime,SurfaceArea, TankVolume=5678,SWDenstiy=1.023, NEC){
  
  deltaDIC<-(TankVolume*SWDenstiy*(HeaderDIC-TankDIC)/(ResidenceTime*SurfaceArea))/1000
  #substract calcification rate
  NCP<-deltaDIC-NEC
  return(NCP)
  #NEC calc calculated the net ecosystem calcification of a flow through mesocosm system at a given time point
  #this uses residence time (lagrangian) rather than change in TA over time (Eulerian-- this is what I used for the 
  #biogeochemistry paper)
  #NCP is in umol C cm^-2 hr-1  (or umol g-1 hr-1 if using one of the other measurements)
  
  #HeaderDIC is DIC from the header in umol/kg
  #TankDICis DIC from the tank in umol/kg
  #Residence time is the residence time in hours
  #Surface area is SA of the substrate in cm2
  #TankVolume is the volume in cm3 = (default = 5678)
  #SWDensity is density of seawater in kg/cm3 (default =1.023)
  #NEC is the calcification rate at that time 
  
}  

##Calculate Nutrient uptake
NutCalc<-function(HeaderN,TankN,ResidenceTime,SurfaceArea, TankVolume=5678,SWDenstiy=1.023){
  
  #NEC<-0.5*(HeaderTA-TankTA)*TankVolume*SWDenstiy/(ResidenceTime*SurfaceArea)/1000
  
  NutRate<-(HeaderN-TankN)*(TankVolume/SurfaceArea)*(SWDenstiy/ResidenceTime)/1000
  
  
  return(NutRate)
  #Nutcalc calculated the nutrient uptake of a flow through mesocosm system at a given time point
  #this uses residence time (lagrangian) r
  #Nutrate is in umol cm^-2 hr-1  (or umol g-1 hr-1 if using one of the other measurements)
  
  
}  

# Load Data-----------------------------------------
#Chem Data
ChemData<-read.csv('Data/AllChemData_noCorrect_Nuts.csv') #read in 1st 12 columns only
ChemData<-ChemData[,-c(24,25)]

#Biology Data

#Coral
Coral <- read.csv("data/CoralSets_Rprocessed.csv",header=TRUE)

#Rubble
Rubble <- read.csv("data/Rubble_Rprocessed.csv",header=TRUE)

#Algae
Algae <- read.csv("data/Algae_Rprocessed.csv",header=TRUE)

#Sand
Sand<- read.csv("data/Sand_Rprocessed.csv",header=TRUE)


## TA Normalization----------------------------------------------------------------------
# 

ChemData$HeaderTA.norm<-ChemData$HeaderTA # I switched to nutrient normalization step to the NEC function see line 31

ChemData$TankTA.norm<-ChemData$TankTA 

#calculuate all the carbonate parameters with seacarb---------------------------------------
#Tanks
TankCO2<-carb(flag=8, ChemData$TankpH, ChemData$TankTA/1000000, S=ChemData$TankSalinity, T=ChemData$TankTemp, Patm=1, P=0, Pt=0, Sit=0,
     k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")
#TA is divided by 1000 because all calculations are in mol/kg in the seacarb package

#convert CO2, HCO3, CO3, DIC, and Alk back to micromol for easier interpretation
TankCO2[,c("CO2","HCO3","CO3","DIC","ALK")]<-TankCO2[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

#add these data to the Dataframe
ChemData[,c("TankCO2","TankHCO3","TankCO3","TankDIC","TankOmegaArag","TankOmegaCalcite","TankpCO2","TankfCO2")]<-
 TankCO2[,c("CO2","HCO3","CO3","DIC","OmegaAragonite","OmegaCalcite","pCO2","fCO2")]

#headers
HeaderCO2<-carbb(flag=8, ChemData$HeaderpH, ChemData$HeaderTA/1000000, S=ChemData$HeaderSalinity, T=ChemData$HeaderTemp, Patm=1, P=0, Pt=0, Sit=0,
              k1k2="x", kf="x", ks="d", pHscale="T", b="u74", gas="potential")
HeaderCO2[,c("CO2","HCO3","CO3","DIC","ALK")]<-HeaderCO2[,c("CO2","HCO3","CO3","DIC","ALK")]*1000000

ChemData[,c("HeaderCO2","HeaderHCO3","HeaderCO3","HeaderDIC","HeaderOmegaArag","HeaderOmegaCalcite","HeaderpCO2","HeaderfCO2")]<-
  HeaderCO2[,c("CO2","HCO3","CO3","DIC","OmegaAragonite","OmegaCalcite","pCO2","fCO2")]

## convert flow to residence time
#Flow is now in ml per min
TankVol<-5678 #hole was at 6 quarts... this is the volume of each tank in ml

ChemData$ResTime<-(1/60)*(1/ChemData$Flow)*TankVol


#Add the algae bits to the final biomass (pieces of algae fell out of the vexarcages, but stayed in aquaria
#during experiment this makes sure that the total biomass is used in the calculation )
Algae$FinalSA<-Algae$FinalSA+Algae$BitsSA #surface area
Algae$FinalWW<-Algae$FinalWW+Algae$BitsWW # wet weight
Algae$FinalVol<-Algae$FinalVol+Algae$BitsVol # volume
Algae$AFDW<-Algae$AFDW+Algae$AFDWbits # organic biomass
Algae$DW<-Algae$DW+Algae$DWbits # dry weight


#calculate the average residence time by aquarium for each experiment (1-36 is exp 1 and 37-72 is exp 2).
ResTime.mean <- ddply(ChemData, c("Aquarium"), summarise,
                          ResTime.mean = mean(ResTime, na.rm = T)
                        
                           )

#calculate the average flow rate by aquarium for each experiment (1-36 is exp 1 and 37-72 is exp 2).
Flow.mean <- ddply(ChemData, c("Aquarium"), summarise,
                      Flow.mean = mean(Flow, na.rm = T)
                      
)


## Sum up all the biological data by aquarium for each experiment
Coral.Exp1Summary <- ddply(Coral, c("Aq_Ex1"), summarise,
                          SA = sum(SA, na.rm = T), #the SA for coral and algae are in mm2 while Sand is in cm2
                          AFDW = sum(AFDW, na.rm = T),
                          DW = sum(DW, na.rm = T),
                          Volume = sum(Volume, na.rm = T)
)

Coral.Exp2Summary <- ddply(Coral, c("Aq_Ex2"), summarise,
                           SA = sum(SA, na.rm = T),
                           AFDW = sum(AFDW, na.rm = T),
                           DW = sum(DW, na.rm = T),
                           Volume = sum(Volume, na.rm = T)
)


Rubble.Exp1Summary <- ddply(Rubble, c("Aq_Ex1"), summarise,
                          SA = sum(SA, na.rm = T),
                           AFDW = sum(AFDW, na.rm = T),
                           DW = sum(DW, na.rm = T),
                           Volume = sum(Volume, na.rm = T)
)

Rubble.Exp2Summary <- ddply(Rubble, c("Aq_Ex2"), summarise,
                            SA = sum(SA, na.rm = T),
                            AFDW = sum(AFDW, na.rm = T),
                            DW = sum(DW, na.rm = T),
                            Volume = sum(Volume, na.rm = T)
)

Algae.Exp1Summary <- ddply(Algae, c("Aq_Ex1"), summarise,
                            SA = sum(FinalSA, na.rm = T),
                            AFDW = sum(AFDW, na.rm = T),
                            DW = sum(DW, na.rm = T),
                            Volume = sum(FinalVol, na.rm = T)
)

Algae.Exp2Summary <- ddply(Algae, c("Aq_Ex2"), summarise,
                           SA = sum(FinalSA, na.rm = T),
                           AFDW = sum(AFDW, na.rm = T),
                           DW = sum(DW, na.rm = T),
                           Volume = sum(FinalVol, na.rm = T)
)

Sand.Exp1Summary <- ddply(Sand, c("Aq_Ex1"), summarise,
                           SA = sum(SA, na.rm = T),
                           AFDW = sum(AFDW, na.rm = T),
                           DW = sum(DW, na.rm = T),
                           Volume = sum(Vol, na.rm = T)
)

Sand.Exp2Summary <- ddply(Sand, c("Aq_Ex2"), summarise,
                          SA = sum(SA, na.rm = T),
                          AFDW = sum(AFDW, na.rm = T),
                          DW = sum(DW, na.rm = T),
                          Volume = sum(Vol, na.rm = T)
)
#join all the biology together
biology<-rbind(Rubble.Exp1Summary,Coral.Exp1Summary, Algae.Exp1Summary,Sand.Exp1Summary) #exp 1

biology2<-rbind(Rubble.Exp2Summary,Coral.Exp2Summary, Algae.Exp2Summary,Sand.Exp2Summary) #exp2

#this sums up all 4 parts for each aquarium for total biomass etc per aquarium
Exp2biology <- ddply(biology2, c("Aq_Ex2"), summarise,
                            SA = sum(SA, na.rm = T),
                            AFDW = sum(AFDW, na.rm = T),
                            DW = sum(DW, na.rm = T),
                            Volume = sum(Volume, na.rm = T)
)

colnames(biology)[1]<-'Aquarium'
colnames(Exp2biology)[1]<-'Aquarium'

biology<-rbind(biology,Exp2biology)

#add the mean residence times
AllData<-merge(ChemData, ResTime.mean, all.x=TRUE)

#add the mean Flow rates times
AllData<-merge(AllData, Flow.mean, all.x=TRUE)

#Make one huge dataset with all the biology and chem data together
AllData<-merge(AllData,biology, by='Aquarium', all.x=TRUE)
#format the dates
AllData$DateTime<-as.POSIXct(paste(AllData$Date, AllData$Time), format="%m/%d/%Y %H:%M")

##sort all the data by time, substrate, nutrient level, and experiment
AllData<-AllData[order(AllData$Experiment, AllData$DateTime, AllData$Substrate, AllData$NutLevel),]

# change the order of the factors for nutrient level so that it is ambient then medium then high

AllData$NutLevel <- factor(AllData$NutLevel, levels = c("Ambient", "Med", "High"))


#add a column for "Day" or "Night"
time<-unique(AllData$Time)
AllData$DayNight<-ifelse(AllData$Time==time[1]|AllData$Time==time[2]|AllData$Time==time[3], 'Day', 'Night')

#Calculate NEC---------------------------------------------

Nuts<-unique(AllData$NutLevel)[c(1,3,2)] #puts the order for loops as ambient, then med, then high
sub<-unique(AllData$Substrate)

#sort data by aquarium so that I can calculate NEC easier

AllData<-AllData[order(AllData$Aquarium),]


#NEC using lagrangian method with AFDW normalization--THIS IS WHAT I AM USING
AllData$NEC.AFDW<-NECCalc(HeaderTA = AllData$HeaderTA.norm, 
                    TankTA = AllData$TankTA.norm, 
                    ResidenceTime = AllData$ResTime.mean, 
                    SurfaceArea = AllData$AFDW,
                    HeaderN = AllData$HeaderN, TankN = AllData$TankN, 
                    HeaderP = AllData$HeaderP, TankP = AllData$TankP)

  
  #NEC using lagrangian method with dry weight normalization
  AllData$NEC.DW<-NECCalc(HeaderTA = AllData$HeaderTA.norm, 
                            TankTA = AllData$TankTA.norm, 
                            ResidenceTime = AllData$ResTime.mean, 
                            SurfaceArea = AllData$DW,
                          HeaderN = AllData$HeaderN, TankN = AllData$TankN, 
                          HeaderP = AllData$HeaderP, TankP = AllData$TankP)
  
  #NEC using lagrangian method with volume normalization
  AllData$NEC.Vol<-NECCalc(HeaderTA = AllData$HeaderTA.norm, 
                            TankTA = AllData$TankTA.norm, 
                            ResidenceTime = AllData$ResTime.mean, 
                            SurfaceArea = AllData$Volume,
                           HeaderN = AllData$HeaderN, TankN = AllData$TankN, 
                           HeaderP = AllData$HeaderP, TankP = AllData$TankP)  
  
  #NEC using lagrangian method with SA normalization
  AllData$NEC.SA<-NECCalc(HeaderTA = AllData$HeaderTA.norm, 
                           TankTA = AllData$TankTA.norm, 
                           ResidenceTime = AllData$ResTime.mean, 
                           SurfaceArea = AllData$SA,
                          HeaderN = AllData$HeaderN, TankN = AllData$TankN, 
                          HeaderP = AllData$HeaderP, TankP = AllData$TankP) 

 
  ###NCP calcs----------------------------------------------
  ## DIC
  AllData$TankDIC.norm<-AllData$TankDIC
  AllData$HeaderDIC.norm<-AllData$HeaderDIC
  
 
  #NCP using lagrangian method with AFDW normalization
  AllData$NCP.AFDW<-NCPCalc(HeaderDIC = AllData$HeaderDIC.norm, 
                            TankDIC = AllData$TankDIC.norm, 
                            ResidenceTime = AllData$ResTime.mean, 
                            SurfaceArea = AllData$AFDW,
                            NEC = AllData$NEC.AFDW)
  
  
  #NCP using lagrangian method with dry weight normalization
  AllData$NCP.DW<-NCPCalc(HeaderDIC = AllData$HeaderDIC.norm, 
                            TankDIC = AllData$TankDIC.norm, 
                            ResidenceTime = AllData$ResTime.mean, 
                            SurfaceArea = AllData$DW,
                            NEC = AllData$NEC.DW)
  
  #NCP using lagrangian method with volume normalization
  AllData$NCP.Vol<-NCPCalc(HeaderDIC = AllData$HeaderDIC.norm, 
                            TankDIC = AllData$TankDIC.norm, 
                            ResidenceTime = AllData$ResTime.mean, 
                            SurfaceArea = AllData$Vol,
                            NEC = AllData$NEC.Vol) 
  
  #NEC using lagrangian method with SA normalization
  AllData$NCP.SA<-NCPCalc(HeaderDIC = AllData$HeaderDIC.norm, 
                           TankDIC = AllData$TankDIC.norm, 
                           ResidenceTime = AllData$ResTime.mean, 
                           SurfaceArea = AllData$SA,
                           NEC = AllData$NEC.SA) 
  
  
  
   #### calculate means------------------

  #calcification mean by substrate, treatment, and sampling time--- this just takes the means over the replicate aquaria
NEC.mean <- ddply(AllData, c("Substrate","NutLevel","DateTime"), summarise,
                     Mean.AFDW = mean(NEC.AFDW, na.rm = T),
                     N=sum(!is.na(NEC.AFDW)),
                     SE.AFDW= sd(NEC.AFDW, na.rm = T)/sqrt(N),
                     Mean.SA = mean(NEC.SA, na.rm = T),
                     SE.SA= sd(NEC.SA, na.rm = T)/sqrt(N),
                     Mean.DW = mean(NEC.DW, na.rm = T),
                     SE.DW= sd(NEC.DW, na.rm = T)/sqrt(N),
                     Mean.Vol = mean(NEC.Vol, na.rm = T),
                     SE.Vol= sd(NEC.Vol, na.rm = T)/sqrt(N)
)
#production
  NCP.mean <- ddply(AllData, c("Substrate","NutLevel","DateTime"), summarise,
                    Mean.AFDW = mean(NCP.AFDW, na.rm = T),
                    N=sum(!is.na(NCP.AFDW)),
                    SE.AFDW= sd(NCP.AFDW, na.rm = T)/sqrt(N),
                    Mean.SA = mean(NCP.SA, na.rm = T),
                    SE.SA= sd(NCP.SA, na.rm = T)/sqrt(N),
                    Mean.DW = mean(NCP.DW, na.rm = T),
                    SE.DW= sd(NCP.DW, na.rm = T)/sqrt(N),
                    Mean.Vol = mean(NCP.Vol, na.rm = T),
                    SE.Vol= sd(NCP.Vol, na.rm = T)/sqrt(N)
  )
  
  

##----------------------------------------------------------------------
#calculate daily average per substrate and nutrient
#add a day night column
time<-unique(NEC.mean$DateTime)
NEC.mean$DayNight<-ifelse(NEC.mean$DateTime==time[1]|NEC.mean$DateTime==time[2]|NEC.mean$DateTime==time[3]|
                            NEC.mean$DateTime==time[7]|NEC.mean$DateTime==time[8]|
                            NEC.mean$DateTime==time[9]|NEC.mean$DateTime==time[10]|
                            NEC.mean$DateTime==time[14], 'Day', 'Night')
NCP.mean$DayNight<-NEC.mean$DayNight

#summarise across all the different times
NEC.mean.Net <- ddply(NEC.mean, c("Substrate","NutLevel"), summarise,
                  Mean.AFDW2 = mean(Mean.AFDW, na.rm = T),
                  N2=sum(!is.na(Mean.AFDW)),
                  SE.AFDW2= sd(Mean.AFDW, na.rm = T)/sqrt(N2),
                  Mean.SA2 = mean(Mean.SA, na.rm = T),
                  SE.SA2= sd(Mean.SA, na.rm = T)/sqrt(N2),
                  Mean.DW2 = mean(Mean.DW, na.rm = T),
                  SE.DW2= sd(Mean.DW, na.rm = T)/sqrt(N2),
                  Mean.Vol2 = mean(Mean.Vol, na.rm = T),
                  SE.Vol2= sd(Mean.Vol, na.rm = T)/sqrt(N2)
)

NCP.mean.Net <- ddply(NCP.mean, c("Substrate","NutLevel"), summarise,
                      Mean.AFDW2 = mean(Mean.AFDW, na.rm = T),
                      N2=sum(!is.na(Mean.AFDW)),
                      SE.AFDW2= sd(Mean.AFDW, na.rm = T)/sqrt(N2),
                      Mean.SA2 = mean(Mean.SA, na.rm = T),
                      SE.SA2= sd(Mean.SA, na.rm = T)/sqrt(N2),
                      Mean.DW2 = mean(Mean.DW, na.rm = T),
                      SE.DW2= sd(Mean.DW, na.rm = T)/sqrt(N2),
                      Mean.Vol2 = mean(Mean.Vol, na.rm = T),
                      SE.Vol2= sd(Mean.Vol, na.rm = T)/sqrt(N2)
)



#calculate averages for day and night separately
NEC.mean.DayNight <- ddply(NEC.mean, c("Substrate","NutLevel", "DayNight"), summarise,
                      Mean.AFDW2 = mean(Mean.AFDW, na.rm = T),
                      N2=sum(!is.na(Mean.AFDW)),
                      SE.AFDW2= sd(Mean.AFDW, na.rm = T)/sqrt(N2),
                      Mean.SA2 = mean(Mean.SA, na.rm = T),
                      SE.SA2= sd(Mean.SA, na.rm = T)/sqrt(N2),
                      Mean.DW2 = mean(Mean.DW, na.rm = T),
                      SE.DW2= sd(Mean.DW, na.rm = T)/sqrt(N2),
                      Mean.Vol2 = mean(Mean.Vol, na.rm = T),
                      SE.Vol2= sd(Mean.Vol, na.rm = T)/sqrt(N2)
)

NCP.mean.DayNight <- ddply(NCP.mean, c("Substrate","NutLevel", "DayNight"), summarise,
                      Mean.AFDW2 = mean(Mean.AFDW, na.rm = T),
                      N2=sum(!is.na(Mean.AFDW)),
                      SE.AFDW2= sd(Mean.AFDW, na.rm = T)/sqrt(N2),
                      Mean.SA2 = mean(Mean.SA, na.rm = T),
                      SE.SA2= sd(Mean.SA, na.rm = T)/sqrt(N2),
                      Mean.DW2 = mean(Mean.DW, na.rm = T),
                      SE.DW2= sd(Mean.DW, na.rm = T)/sqrt(N2),
                      Mean.Vol2 = mean(Mean.Vol, na.rm = T),
                      SE.Vol2= sd(Mean.Vol, na.rm = T)/sqrt(N2))

#PR
NCP.mean.PRbyTank <- ddply(AllData, c("Substrate","NutLevel", "Tank", "Aquarium"), summarise,
                          
                           Mean.AFDW2.p = mean(NCP.AFDW[DayNight=='Day'], na.rm=T),
                           Mean.AFDW2.r =mean(abs(NCP.AFDW[DayNight=='Night']), na.rm = T),       
                           Mean.AFDW2.NCP = mean(NCP.AFDW, na.rm=T),
                           GPP=Mean.AFDW2.p+Mean.AFDW2.r,
                           PR.AFDW2=GPP/Mean.AFDW2.r,
                           Mean.AFDW2.Day = mean(NEC.AFDW[DayNight=='Day'], na.rm=T), #day NEC
                           Mean.AFDW2.Night = mean(NEC.AFDW[DayNight=='Night'], na.rm=T), #Night NEC
                           Mean.AFDW2.NEC = mean(NEC.AFDW, na.rm=T) #Night NEC
                           
)



NCP.mean.PR<-ddply(NCP.mean.PRbyTank, c("Substrate","NutLevel"), summarise,
                Mean.AFDW2=mean(PR.AFDW2, na.rm=T),
                Mean.GPP=mean(GPP, na.rm=T),
                N2=3,
                SE.AFDW2=sd(PR.AFDW2)/sqrt(N2),
                SE.GPP=sd(GPP)/sqrt(N2)
)




## Stats-------------------------------------
#I took the average across times then ran GLMs for each inidividual substrate using black tank as a random effect

#NEC net

#Take the average of everything across aquariums so its easier to code in the model.  I need to drop the factors for that to work
AllData2<-AllData
levels(AllData2$Substrate)<-c(1:5)
AllData2$Substrate<-as.numeric(droplevels(AllData2$Substrate))
levels(AllData2$NutLevel)<-c(1:3)
AllData2$NutLevel<-as.numeric(droplevels(AllData2$NutLevel))

Mean.Time<-aggregate(AllData2, list(AllData2$Aquarium), mean)
#now add the factors back for the GLM
Mean.Time$NutLevel<-as.factor(Mean.Time$NutLevel)
levels(Mean.Time$NutLevel)<-c('Ambient','Med','High')
Mean.Time$Substrate<-as.factor(Mean.Time$Substrate)
levels(Mean.Time$Substrate)<-c('Algae','Coral','Mixed','Rubble','Sand')


#add mean respiration rates by aquarium to all data so that I can add it to the photosynthesis
#to get GPP across each time point

AllData <-merge(AllData, NCP.mean.PRbyTank[,c(4,6)], all.x=TRUE, by.x = "Aquarium")
AllData$GPP<-AllData$NCP.AFDW+AllData$Mean.AFDW2.r


#Algae
model.NECNet.algae<-lmer(NEC.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Algae' ,]) 
model.NECDay.algae<-lmer(NEC.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Algae' & AllData$DayNight=='Day',]) 
model.NECNight.algae<-lmer(NEC.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Algae' & AllData$DayNight=='Night',]) 
model.NCPNet.algae<-lmer(NCP.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Algae' ,]) 
model.R.algae<-lmer(NCP.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Algae' & AllData$DayNight=='Night',]) 
model.GCP.algae<-lmer(GPP~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Algae' & AllData$DayNight=='Day',]) 
model.PR.algae<-lmer(PR.AFDW2~NutLevel +(1|Tank), data=NCP.mean.PRbyTank[NCP.mean.PRbyTank$Substrate=='Algae',])

#Coral
model.NECNet.Coral<-lmer(NEC.AFDW~NutLevel +(1|Tank) +(1|DateTime), data=AllData[AllData$Substrate=='Coral' ,]) 
model.NECDay.Coral<-lmer(NEC.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Coral' & AllData$DayNight=='Day',]) 
model.NECNight.Coral<-lmer(NEC.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Coral' & AllData$DayNight=='Night',]) 
model.NCPNet.Coral<-lmer(NCP.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Coral' ,]) 
model.R.Coral<-lmer(NCP.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Coral' & AllData$DayNight=='Night',]) 
model.GCP.Coral<-lmer(GPP~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Coral' & AllData$DayNight=='Day',]) 
model.PR.Coral<-lmer(PR.AFDW2~NutLevel +(1|Tank), data=NCP.mean.PRbyTank[NCP.mean.PRbyTank$Substrate=='Coral',])


#Rubble
model.NECNet.Rubble<-lmer(NEC.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Rubble' ,]) 
model.NECDay.Rubble<-lmer(NEC.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Rubble' & AllData$DayNight=='Day',]) 
model.NECNight.Rubble<-lmer(NEC.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Rubble' & AllData$DayNight=='Night',]) 
model.NCPNet.Rubble<-lmer(NCP.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Rubble' ,]) 
model.R.Rubble<-lmer(NCP.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Rubble' & AllData$DayNight=='Night',]) 
model.GCP.Rubble<-lmer(GPP~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Rubble' & AllData$DayNight=='Day',]) 
model.PR.Rubble<-lmer(PR.AFDW2~NutLevel +(1|Tank), data=NCP.mean.PRbyTank[NCP.mean.PRbyTank$Substrate=='Rubble',])

#Sand
model.NECNet.Sand<-lmer(NEC.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Sand' ,]) 
model.NECDay.Sand<-lmer(NEC.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Sand' & AllData$DayNight=='Day',]) 
model.NECNight.Sand<-lmer(NEC.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Sand' & AllData$DayNight=='Night',]) 
model.NCPNet.Sand<-lmer(NCP.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Sand' ,]) 
model.R.Sand<-lmer(NCP.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Sand' & AllData$DayNight=='Night',]) 
model.GCP.Sand<-lmer(GPP~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Sand' & AllData$DayNight=='Day',]) 
model.PR.Sand<-lmer(PR.AFDW2~NutLevel +(1|Tank), data=NCP.mean.PRbyTank[NCP.mean.PRbyTank$Substrate=='Sand',])

#Mixed
model.NECNet.Mixed<-lmer(NEC.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Mixed' ,]) 
model.NECDay.Mixed<-lmer(NEC.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Mixed' & AllData$DayNight=='Day',]) 
model.NECNight.Mixed<-lmer(NEC.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Mixed' & AllData$DayNight=='Night',]) 
model.NCPNet.Mixed<-lmer(NCP.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Mixed' ,]) 
model.R.Mixed<-lmer(NCP.AFDW~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Mixed' & AllData$DayNight=='Night',]) 
model.GCP.Mixed<-lmer(GPP~NutLevel +(1|Tank)+(1|DateTime), data=AllData[AllData$Substrate=='Mixed' & AllData$DayNight=='Day',]) 
model.PR.Mixed<-lmer(PR.AFDW2~NutLevel +(1|Tank), data=NCP.mean.PRbyTank[NCP.mean.PRbyTank$Substrate=='Mixed',])

#omega vs NEC------------------- varying intercepts for aquarium within black tank and varying slopes (and intercept) by nutrient level

model.NECOmega.Algae<-lmer(NEC.AFDW~ TankOmegaArag*NutLevel+ (1|Tank/Aquarium), data=AllData[AllData$Substrate=='Algae',])

model.NECOmega.Coral<-lmer(NEC.AFDW~ TankOmegaArag*NutLevel+ (1|Tank/Aquarium), data=AllData[AllData$Substrate=='Coral',])

model.NECOmega.Sand<-lmer(NEC.AFDW~ TankOmegaArag*NutLevel+ (1|Tank/Aquarium), data=AllData[AllData$Substrate=='Sand',])

model.NECOmega.Rubble<-lmer(NEC.AFDW~ TankOmegaArag*NutLevel+ (1|Tank/Aquarium), data=AllData[AllData$Substrate=='Rubble',])

model.NECOmega.Mixed<-lmer(NEC.AFDW~ TankOmegaArag*NutLevel+(1|Tank/Aquarium), data=AllData[AllData$Substrate=='Mixed',])

# models for pH vs NCP
#Ranef to account for repeated measures of aquarium within tank and time  
model.pH.NCP<-lmer(TankpH~NCP.AFDW*NutLevel*Substrate+ (1|Tank/Aquarium) , data=AllData)
 


#### Change color for pub quality figures
mypalette<-brewer.pal(9,"Blues") # light pink to purple for ambient to high nutrients
mypalette<-mypalette[c(3,6,9)] # pull out 3 contrasting colors for ambient, med, and high

#add the colors to AllData
AllData$colors<-NA
  for (i in 1:length(Nuts)){
AllData$colors[AllData$NutLevel ==Nuts[i]]<-mypalette[i] 
  }

###Looking at feedbacks--------------------------------------------
# Plot NCP versus pH
pdf(file = 'Plots/MSplots/NCPvspH.pdf', height = 6, width = 6, useDingbats = FALSE)
par(mfrow=c(1,1), pty='s')
#plot(AllData$NCP.AFDW, AllData$TankpH, col=AllData$NutLevel, xlab='NCP', ylab='pH')
#legend('topleft',legend=c('Ambient',"Medium","High"), col=unique(AllData$NutLevel), pch=19, bty = 'n')
plot(AllData$NCP.AFDW[AllData$Substrate=='Coral'], AllData$TankpH[AllData$Substrate=='Coral'], col=AllData$colors[AllData$Substrate=='Coral'], xlab='NCP', ylab='pH', pch=0)
points(AllData$NCP.AFDW[AllData$Substrate=='Algae'], AllData$TankpH[AllData$Substrate=='Algae'], col=AllData$colors[AllData$Substrate=='Algae'],  pch=1)
points(AllData$NCP.AFDW[AllData$Substrate=='Rubble'], AllData$TankpH[AllData$Substrate=='Rubble'], col=AllData$colors[AllData$Substrate=='Rubble'],  pch=3)
points(AllData$NCP.AFDW[AllData$Substrate=='Sand'], AllData$TankpH[AllData$Substrate=='Sand'], col=AllData$colors[AllData$Substrate=='Sand'],  pch=2)
points(AllData$NCP.AFDW[AllData$Substrate=='Mixed'], AllData$TankpH[AllData$Substrate=='Mixed'], col=AllData$colors[AllData$Substrate=='Mixed'],  pch=19, cex=0.2)
legend('topleft', legend=c('Coral','Algae','Rubble','Sand','Mixed'),pch=c(0,1,3,2,19), bty='n')
legend('bottomright',legend=c('Ambient',"Medium","High"), col=unique(AllData$colors), pch=19, bty = 'n')
abline(v=0, lty=2, col='grey')
dev.off()


## NEC vs Omega plot-------------------------------------
##plot for the NEC versus aragonite saturation state relationships by substrate and treatment

pdf("plots/MSplots/NECvsOmega.pdf", width=6, height=8.5, useDingbats = FALSE )
par(mfrow=c(3,2))
#cols <- c(unique(NEC.mean$NutLevel))
cols <- mypalette
y<-AllData$NEC.AFDW
yse<-NEC.mean$SE.AFDW
for (i in 1:length(sub)){
  plot(NA, xlab=expression(paste(Omega)[arag]),ylim=c(min(y[AllData$Substrate==sub[i]]), max(y[AllData$Substrate==sub[i]])), 
       ylab=expression(paste("NCC ",mu,"mol CaCO"[3], " g"^{-1}," hr"^{-1})), main = sub[i], xlim=c(1.5,6), cex.lab=1.5, cex.main = 2)
  
  abline(h=0, lty=2)
  for (j in 1:length(Nuts)){
    par(new = TRUE)
    plot(AllData$TankOmegaArag[AllData$Substrate==sub[i] & AllData$NutLevel==Nuts[j]], yaxt='n',
         y[AllData$Substrate==sub[i] & AllData$NutLevel==Nuts[j]], col = cols[j],
         pch=19, type="p", xaxt='n', xlim=c(1.5,6), ylab='', xlab='',ylim=c(min(y[AllData$Substrate==sub[i]]), max(y[AllData$Substrate==sub[i]])))
   #adjustcolor( cols[j], alpha.f = 0.5)
    
    model<-lm(y[AllData$Substrate==sub[i] & AllData$NutLevel==Nuts[j]]~AllData$TankOmegaArag[AllData$Substrate==sub[i] & AllData$NutLevel==Nuts[j]])
    #if the p>0.05 make it a dashed line (or don't include)
    p<-anova(model)$`Pr(>F)`[1]
    if(p<=0.05){
      lines(AllData$TankOmegaArag[AllData$Substrate==sub[i] & AllData$NutLevel==Nuts[j]], model$fitted.values, col=cols[j], lwd=3, lty=1)
    }
      }
}
#empty plot to put the legend
plot(NA, ylim=c(0,1), xlim=c(0,1), axes=FALSE,  ylab="", xlab="")
legend('center', horiz = FALSE, legend=unique(NEC.mean$NutLevel), col=mypalette, pch=19, bty = 'n')

dev.off()

## pH plots----------------------------------------------------------------
#pH means across substrate, treatment, and time
deltapHMeans <- ddply(AllData, c("Substrate","NutLevel", "DayNight"), summarise,
                           pHMean = mean(TankpH-HeaderpH, na.rm = T),
                           N2=sum(!is.na(TankpH)),
                           pHSE= sd(TankpH-HeaderpH, na.rm = T)/sqrt(N2)#,
                          # pH.NEC=mean(pHNEC.treat),
                          # pH.NECSE= sd(pHNEC.treat, na.rm = T)/sqrt(N2),
                           #pH.NCP=mean(pHNCP.treat),
                           #pH.NCPSE= sd(pHNCP.treat, na.rm = T)/sqrt(N2)
)



#absolute change in pH
deltapHMeans.net <- ddply(AllData, c("Substrate","NutLevel"), summarise,
                      pHMean = mean(abs(TankpH-HeaderpH), na.rm = T),
                      N2=sum(!is.na(TankpH)),
                      pHSE= sd(abs(TankpH-HeaderpH), na.rm = T)/sqrt(N2)
                      
)


 #pH across time by substrate
  
  deltapHMeans.time <- ddply(AllData, c("Substrate","NutLevel", "DateTime"), summarise,
                        pHMean = mean(TankpH-HeaderpH, na.rm = T),
                        N2=sum(!is.na(TankpH)),
                        pHSE= sd(TankpH-HeaderpH, na.rm = T)/sqrt(N2)
                      
  )
  
  #Delta pH across time
  pdf(file = 'Plots/MSplots/DeltapHTime.pdf', height = 8, width = 6, useDingbats = FALSE)
  par(mfrow=c(3,2))
  rm(y)
  rm(yse)
  y<-deltapHMeans.time$pHMean
  yse<-deltapHMeans.time$pHSE
  for (i in 1:length(sub)){
    plot(NA, xaxt='n', xlab="Time",ylim=c(min(y), max(y)+0.1), ylab=expression(paste(Delta,"pH")), main = sub[i])
    
    abline(h=0)
    par(new = TRUE)
    #cols <- unique(NCP.mean$NutLevel)
    #
    for (j in 1:length(Nuts)){
      par(new = TRUE)
        
      plot(as.numeric(deltapHMeans.time$DateTime [deltapHMeans.time$Substrate==sub[i] & deltapHMeans.time$NutLevel==Nuts[j]]),
           y[deltapHMeans.time$Substrate==sub[i] & deltapHMeans.time$NutLevel==Nuts[j]], col = cols[j],
           pch=19, type="b", xaxt='n', ylab='', xlab='',ylim=c(min(y), max(y)+.1), yaxt='n')
      
      arrows(unique(deltapHMeans.time $DateTime), y[deltapHMeans.time $Substrate==sub[i] & deltapHMeans.time $NutLevel==Nuts[j]]
             + yse[deltapHMeans.time $Substrate==sub[i] & deltapHMeans.time $NutLevel==Nuts[j]], 
             unique(deltapHMeans.time $DateTime), y[deltapHMeans.time $Substrate==sub[i] & deltapHMeans.time$NutLevel==Nuts[j]]
             - yse[deltapHMeans.time $Substrate==sub[i] & deltapHMeans.time $NutLevel==Nuts[j]], 
             angle=90, code=3, length = 0.05, col = cols[j])
      start<-ifelse(i<=4,c(1),c(8))
      stops<-ifelse(i<=4,c(7),c(14))
      
    }
    axis(1, at=unique(deltapHMeans.time$DateTime)[start:stops], labels=c('10:00',"14:00","18:00","22:00","02:00","06:00","10:00"))
    
    #shaded area for night
    a<-ifelse(i<=4,3,10) #because mixed has different dates than the rest of the substrats
    b<-ifelse(i<=4,6,13)
    
    rect(unique(deltapHMeans.time $DateTime)[a]+3600,min(y),unique(deltapHMeans.time$DateTime)[b]+3600,max(y)+.1,col = rgb(0.5,0.5,0.5,1/4), border = NA)
  }
  legend('bottomleft', legend=unique(deltapHMeans.time$NutLevel), col=mypalette, pch=19, bty = 'n')
  
  dev.off()
 
  #calculate the mean nutrient conditions
  meanNuts<-ddply(AllData, 'NutLevel', summarise,
        meanNN=mean(HeaderN, na.rm=T),
        SENN=sd(HeaderN, na.rm=T)/sqrt(7),
        meanP=mean(HeaderP, na.rm=T),
        SEP=sd(HeaderP, na.rm=T)/sqrt(7))
  pdf(file = 'Plots/MSplots/NutrientsColor.pdf', height = 5, width = 8, useDingbats = FALSE)      
  par(mfrow=c(1,2)) #average nitrate for experiment
  x<-barplot(meanNuts$meanNN, main=expression(paste('NO'[3]^{'-'},'+ NO'[2]^{'-'})), 
             ylab=expression(paste(mu,'M')), ylim=c(0,10), col = mypalette)
  errorbars(x,meanNuts$meanNN,0,meanNuts$SENN)
  axis(1, at=x, labels=c("Ambient","Medium","High"))
  
  #average phosphate for experiment
  x<-barplot(meanNuts$meanP, main=expression(paste('PO'[4]^{'3-'})), 
             ylab='', ylim=c(0,3.0), col= mypalette)
  errorbars(x,meanNuts$meanP,0,meanNuts$SEP)
  axis(1, at=x, labels=c("Ambient","Medium","High"))
  dev.off()



 ##### Mixed effects models for the change in BW for algae, rubble, coral 
 #coral
 # NOTE THE pc DATA IS PROPORTIONAL, NOT PERCENT SO MULTUPLY EVERYTHING BY 100
 Coral$pcDeltaBW<-Coral$pcDeltaBW*100
 Algae$pcDeltaWW<-Algae$pcDeltaWW*100
 Rubble$pcDeltaBW<-Rubble$pcDeltaBW*100
 
 CoralBW.mod<-lmer(pcDeltaBW ~ Nuts+Species + (1|Tank), data= Coral) # I looked for an interaction, but found none
 #CoralBW.mod<-lmer(pcDeltaBW ~ Nuts + (1|Tank/Species), data= Coral) # this puts species nested within tank as a random effect
 
 # order the nutrient levels
  Coral$Nuts<-ordered(Coral$Nuts, levels =c('Ambient','Medium','High') )
 # calculate means and SE for plotting
 Coral.mean<-ddply(Coral, c("Species","Nuts"), summarise,
       pcBW.mean = mean(pcDeltaBW, na.rm=TRUE),
       n = sum(!is.na(pcDeltaBW)),
       pcBW.se = sd(pcDeltaBW, na.rm=TRUE)/sqrt(n),
       CoralGrowth.mean = mean((1000/50) *DeltaBW/AFDW, na.rm=TRUE), #mg/ gram AFDW rate per day
       CoralGrowth.SE = sd((1000/50) *DeltaBW/AFDW, na.rm=TRUE)/sqrt(n),
       CoralGrowth.meanSA = mean((1000/50) *DeltaBW/SA, na.rm=TRUE), #gram/ SA AFDW rate
       CoralGrowthSA.SE = sd((1000/50) *DeltaBW/SA, na.rm=TRUE)/sqrt(n)
       )
 
 #algae
 AlgaeBW.mod<-lmer(pcDeltaWW ~ Nuts + (1|Tank), data= Algae) # I looked for an interaction, but found none
 Algae$Nuts<-ordered(Algae$Nuts, levels =c('Ambient','Medium','High') )
 Algae.mean<-ddply(Algae, c("Nuts"), summarise,
                   pcBW.mean = mean(pcDeltaWW, na.rm=TRUE),
                   n = sum(!is.na(pcDeltaWW)),
                   pcBW.se = sd(pcDeltaWW, na.rm=TRUE)/sqrt(n),
                   AlgaeGrowth.mean = mean((1000/50) *DeltaWW/(AFDW+AFDWbits), na.rm=TRUE), #mg/ gram AFDW rate per day
                   AlgaeGrowth.SE = sd((1000/50) *DeltaWW/(AFDW+AFDWbits), na.rm=TRUE)/sqrt(n),
                   AlgaeGrowth.meanSA = mean((1000/50) *DeltaWW/FinalSA, na.rm=TRUE), #gram/ SA AFDW rate
                   AlgaeGrowthSA.SE = sd((1000/50) *DeltaWW/FinalSA, na.rm=TRUE)/sqrt(n))
 
 # rubble
RubbleBW.mod<-lmer(pcDeltaBW ~ Nuts + (1|Tank), data= Rubble) # I looked for an interaction, but found none
Rubble$Nuts<-ordered(Rubble$Nuts, levels =c('Ambient','Medium','High') )
Rubble.mean<-ddply(Rubble, c("Nuts"), summarise,
                  pcBW.mean = mean(pcDeltaBW, na.rm=TRUE),
                  n = sum(!is.na(pcDeltaBW)),
                  pcBW.se = sd(pcDeltaBW, na.rm=TRUE)/sqrt(n),
                  RubbleGrowth.mean = mean((1000/50) *DeltaBW/AFDW, na.rm=TRUE), #mg/ gram AFDW rate per day
                  RubbleGrowth.SE = sd((1000/50) *DeltaBW/AFDW, na.rm=TRUE)/sqrt(n),
                  RubbleGrowth.meanSA = mean((1000/50) *DeltaBW/SA, na.rm=TRUE), #gram/ SA AFDW rate
                  RubbleGrowthSA.SE = sd((1000/50) *DeltaBW/SA, na.rm=TRUE)/sqrt(n))


pdf('plots/MSplots/BuoyantWeight_dot.pdf', width = 4, height = 8.5, useDingbats = FALSE)
par(mfrow=c(3,1))

plot(c(1:3), Coral.mean$CoralGrowth.mean[1:3], pch=19, ylim=c(10, 35), xlab="", xaxt='n', ylab = expression(paste('Coral Growth (mg g ', AFDW^-1, day^-1,')')), col = mypalette, cex=1.5)
lines(c(1:3), Coral.mean$CoralGrowth.mean[1:3], col = 'black', type = "c")
arrows(c(1:3), Coral.mean$CoralGrowth.mean + Coral.mean$CoralGrowth.SE,  # x , mean + SE
       c(1:3), Coral.mean$CoralGrowth.mean - Coral.mean$CoralGrowth.SE,  # x, mean - SE
       angle=90, code=3, length = 0.05)
points(c(1:3), Coral.mean$CoralGrowth.mean[4:6], pch=15, ylim=c(0, 5), col = mypalette, cex=1.5)
lines(c(1:3), Coral.mean$CoralGrowth.mean[4:6], col = 'black', type = "c")
legend('topright',legend = c(expression(italic('M. capitata')), expression(italic('P. compressa'))), bty='n',
       pch = c(19,15), col = 'black')

  
# algae
plot(1:3,Algae.mean$AlgaeGrowth.mean,   ylim = c(240,340), xlab="",xaxt='n',ylab = expression(paste('Algal Growth (mg g ', AFDW^-1, day^-1,')')), pch = 19, col = mypalette, cex=1.5)
lines(c(1:3), Algae.mean$AlgaeGrowth.mean, col = 'black', type = "c")
arrows(1:3, Algae.mean$AlgaeGrowth.mean + Algae.mean$AlgaeGrowth.SE,  # x , mean + SE
       1:3, Algae.mean$AlgaeGrowth.mean - Algae.mean$AlgaeGrowth.SE,  # x, mean - SE
       angle=90, code=3, length = 0.05)
       

#Rubble
plot(Rubble.mean$RubbleGrowth.mean, xaxt='n',  ylab = expression(paste('Rubble Growth (mg g ', AFDW^-1, day^-1,')')),xlab="",
     pch = 19, col = mypalette, cex=1.5, ylim=c(-3,2))
arrows(1:3, Rubble.mean$RubbleGrowth.mean + Rubble.mean$RubbleGrowth.SE,  # x , mean + SE
       1:3, Rubble.mean$RubbleGrowth.mean - Rubble.mean$RubbleGrowth.SE,  # x, mean - SE
       angle=90, code=3, length = 0.05)
abline(h=0, lty=2)
lines(c(1:3), Rubble.mean$RubbleGrowth.mean, col = 'black', type = "c")
text(x = 1:3, par("usr")[3] ,  labels = c("Ambient","Medium","High"), pos = 1, xpd = TRUE)

dev.off()


### NEC and NCP plots and stats file #####
substrate<-c('Coral','Algae','Rubble','Sand','Mixed')

# create a list of all the models
mods.NEC<-list(model.NECNet.Coral,model.NECNet.algae, model.NECNet.Rubble,model.NECNet.Sand,model.NECNet.Mixed)
mods.NECD<-list(model.NECDay.Coral, model.NECDay.algae, model.NECDay.Rubble, model.NECDay.Sand,model.NECDay.Mixed)            
mods.NECN<-list(model.NECNight.Coral, model.NECNight.algae, model.NECNight.Rubble, model.NECNight.Sand,model.NECNight.Mixed)            

#  create a dataframe of the stats to report in a the paper
NEC.net.stats<-data.frame(matrix(NA,nrow=5,ncol=7, dimnames = list(c(1:5),c('Substrate','SS','MS','NumDF','DenDF','F','P'))))
NEC.day.stats<-data.frame(matrix(NA,nrow=5,ncol=7, dimnames = list(c(1:5),c('Substrate','SS','MS','NumDF','DenDF','F','P'))))
NEC.night.stats<-data.frame(matrix(NA,nrow=5,ncol=7, dimnames = list(c(1:5),c('Substrate','SS','MS','NumDF','DenDF','F','P'))))

# Make a function so I can minimize errors
Nutplot.NEC<-function(species,i, Main = TRUE, YLAB = 'NEC'){
  
  #Day
  y<-summary(mods.NECD[[i]])$coefficients[,1:2]
  a<-anova(mods.NECD[[i]])
  ef<-as.data.frame(effect("NutLevel", mods.NECD[[i]]))
  
  SE.upper<-ef$fit+NEC.mean.DayNight$SE.AFDW2[NEC.mean.DayNight$Substrate==species & NEC.mean.DayNight$DayNight=='Day']
  SE.lower<-ef$fit-NEC.mean.DayNight$SE.AFDW2[NEC.mean.DayNight$Substrate==species& NEC.mean.DayNight$DayNight=='Day']
  
  # make plotting y axis beter
  ef2<-as.data.frame(effect("NutLevel", mods.NECN[[i]]))
  SE.low<-ef2$fit-NEC.mean.DayNight$SE.AFDW2[NEC.mean.DayNight$Substrate==species& NEC.mean.DayNight$DayNight=='Night']
  
  #plot
  plot(1:3-0.2,ef$fit, pch = 2,  cex = 3, main=ifelse(Main ==TRUE, species,NA),
       ylim=c(c(ifelse(min(SE.low)>0,0,floor(min(SE.low))),ceiling(max(SE.upper)))),
       yaxp=c(c(ifelse(min(SE.low)>0,0,floor(min(SE.low))),ceiling(max(SE.upper))), 5),
       ylab=ifelse(YLAB=='NEC', expression(paste("NEC ",mu,"mol CaCO"[3], " g"^{-1}," hr"^{-1})), 
                   expression(paste("NCP ",mu,"mol C g"^{-1}," hr"^{-1}))),
       cex.main=2, cex.axis=2, cex.lab=2, xlab = "", col=mypalette, xaxt='n', xlim=c(0.25,3.75))
  
  lines(1:3-0.2,ef$fit, col = 'black', type = 'c' )
  arrows(1:3-0.2, SE.upper,1:3-0.2, SE.lower, length=0.05, angle=90, code=3,
         col=mypalette, lwd=2)
  
  # add star or NS for significance
  if(a$`Pr(>F)`<=0.055){text(3+0.6,ef$fit[3],'*', cex=3)}   #  legend("top",'*', cex=2, bty="n")
  if(a$`Pr(>F)`>=0.055){text(3+0.6,ef$fit[3],'NS')}
 
  #calculate 95%CI using effects
  y<-summary(mods.NEC[[i]])$coefficients[,1:2]
  a<-anova(mods.NEC[[i]])
  #calculate 95%CI using effects
  ef<-as.data.frame(effect("NutLevel", mods.NEC[[i]]))
  
  SE.upper<-ef$fit+NEC.mean.Net$SE.AFDW2[NEC.mean.Net$Substrate==species]
  SE.lower<-ef$fit-NEC.mean.Net$SE.AFDW2[NEC.mean.Net$Substrate==species]
  
  points(1:3,ef$fit, main=ifelse(Main ==TRUE, species,NA), pch=19,
         col=mypalette, xaxt='n', cex=3)
  lines(1:3,ef$fit, col = 'black', type = 'c' )
  
  arrows(1:3, SE.upper,1:3, SE.lower, length=0.05, angle=90, code=3,
         col=mypalette, lwd=2)
  
  abline(h=0, lty=2, col='grey')
 
  #add a star to the graph if it is statistically significant
  if(a$`Pr(>F)`<=0.055){text(3+0.6,ef$fit[3],'*', cex=3)}   #  legend("top",'*', cex=2, bty="n")
  if(a$`Pr(>F)`>=0.055){text(3+0.6,ef$fit[3],'NS')}
  
  # night
  y<-summary(mods.NECN[[i]])$coefficients[,1:2]
  a<-anova(mods.NECN[[i]])
  ef<-as.data.frame(effect("NutLevel", mods.NECN[[i]]))
  
  SE.upper<-ef$fit+NEC.mean.DayNight$SE.AFDW2[NEC.mean.DayNight$Substrate==species& NEC.mean.DayNight$DayNight=='Night']
  SE.lower<-ef$fit-NEC.mean.DayNight$SE.AFDW2[NEC.mean.DayNight$Substrate==species& NEC.mean.DayNight$DayNight=='Night']
  
  points(1:3+0.2,ef$fit, cex=3, pch = 17, col = mypalette)
  lines(1:3+0.2,ef$fit, col = 'black', type = 'c' )
  arrows(1:3+0.2, SE.upper,1:3+0.2, SE.lower, length=0.05, angle=90, code=3,
         col=mypalette, lwd=2)
  if(a$`Pr(>F)`<=0.055){text(3+0.6,ef$fit[3],'*', cex=3)}   #  legend("top",'*', cex=2, bty="n")
  if(a$`Pr(>F)`>=0.055){text(3+0.6,ef$fit[3],'NS')}
  
  axis(1, at=1:3, labels = c("Ambient","Medium","High"), las = 2, cex.lab=1.5,cex.axis=1.5)
  
}

pdf("plots/MSplots/MeanRates_NEC.pdf", width=18, height=5, useDingbats = FALSE)
j<-2
par(bg=NA) 
par(pty="m")
##  plot of all metabolic rates by substrate
par(mfrow=c(1,5))
par(mar=c(6.2,6.1,4.1,1.5))
par(lwd = 2) 
for(i in 1:5){
  Nutplot.NEC(substrate[i],i,'TRUE','NEC')
  # save the stats
  #day
  a<-anova(mods.NECD[[i]])
  NEC.day.stats[i,1:7]<-t(matrix(c(substrate[i],a[1:6])))
  #night
  a<-anova(mods.NECN[[i]])
  NEC.night.stats[i,1:7]<-t(matrix(c(substrate[i],a[1:6])))
  #net
  a<-anova(mods.NEC[[i]])
  NEC.net.stats[i,1:7]<-t(matrix(c(substrate[i],a[1:6])))
  
}
legend('topright', legend = c('Day NCC','Night NCC','Net Diel NCC'), pch = c(2,17,19), col = mypalette, bty='n', cex=1.5)
dev.off()
###############

## NCP
#  create a dataframe of the stats to report in a the paper
NCP.net.stats<-data.frame(matrix(NA,nrow=5,ncol=7, dimnames = list(c(1:5),c('Substrate','SS','MS','NumDF','DenDF','F','P'))))
NCP.day.stats<-data.frame(matrix(NA,nrow=5,ncol=7, dimnames = list(c(1:5),c('Substrate','SS','MS','NumDF','DenDF','F','P'))))
NCP.night.stats<-data.frame(matrix(NA,nrow=5,ncol=7, dimnames = list(c(1:5),c('Substrate','SS','MS','NumDF','DenDF','F','P'))))

Nutplot.NCP<-function(species,i, Main = TRUE, YLAB = 'NCP'){
  
  #Day
  y<-summary(mods.NECD[[i]])$coefficients[,1:2]
  a<-anova(mods.NECD[[i]])
  ef<-as.data.frame(effect("NutLevel", mods.NECD[[i]]))
  
  SE.upper<-ef$fit+NCP.mean.DayNight$SE.AFDW2[NCP.mean.DayNight$Substrate==species & NCP.mean.DayNight$DayNight=='Day']
  SE.lower<-ef$fit-NCP.mean.DayNight$SE.AFDW2[NCP.mean.DayNight$Substrate==species& NCP.mean.DayNight$DayNight=='Day']
  
  # make plotting y axis beter
  ef2<-as.data.frame(effect("NutLevel", mods.NECN[[i]]))
  SE.low<-ef2$fit-NCP.mean.DayNight$SE.AFDW2[NCP.mean.DayNight$Substrate==species& NCP.mean.DayNight$DayNight=='Night']
  
  #plot
  plot(1:3-0.2,ef$fit, pch = 2,  cex = 3, main=ifelse(Main ==TRUE, species,NA),
       ylim=c(c(ifelse(min(SE.low)>0,0,floor(min(SE.low))),ceiling(max(SE.upper)))),
       yaxp=c(c(ifelse(min(SE.low)>0,0,floor(min(SE.low))),ceiling(max(SE.upper))), 5),
       ylab=ifelse(YLAB=='NEC', expression(paste("NEC ",mu,"mol CaCO"[3], " g"^{-1}," hr"^{-1})), 
                   expression(paste("NCP ",mu,"mol C g"^{-1}," hr"^{-1}))),
       cex.main=2, cex.axis=2, cex.lab=2, xlab = "", col=mypalette, xaxt='n', xlim=c(0.25,3.75))
       
  lines(1:3-0.2,ef$fit, col = 'black', type = 'c' )
  arrows(1:3-0.2, SE.upper,1:3-0.2, SE.lower, length=0.05, angle=90, code=3,
         col=mypalette, lwd=2)
  
  # add star or NS for significance
  if(a$`Pr(>F)`<=0.055){text(3+0.6,ef$fit[3],'*', cex=3)}   #  legend("top",'*', cex=2, bty="n")
  if(a$`Pr(>F)`>=0.055){text(3+0.6,ef$fit[3],'NS')}
  
  #calculate 95%CI using effects
  y<-summary(mods.NEC[[i]])$coefficients[,1:2]
  a<-anova(mods.NEC[[i]])
  #calculate 95%CI using effects
  ef<-as.data.frame(effect("NutLevel", mods.NEC[[i]]))
  
  SE.upper<-ef$fit+NCP.mean.Net$SE.AFDW2[NCP.mean.Net$Substrate==species]
  SE.lower<-ef$fit-NCP.mean.Net$SE.AFDW2[NCP.mean.Net$Substrate==species]
  
  points(1:3,ef$fit, main=ifelse(Main ==TRUE, species,NA), pch=19,
         col=mypalette, xaxt='n', cex=3)
  lines(1:3,ef$fit, col = 'black', type = 'c' )
  
  arrows(1:3, SE.upper,1:3, SE.lower, length=0.05, angle=90, code=3,
         col=mypalette, lwd=2)
  
  abline(h=0, lty=2, col='grey')
 
  #add a star to the graph if it is statistically significant
  if(a$`Pr(>F)`<=0.055){text(3+0.6,ef$fit[3],'*', cex=3)}   #  legend("top",'*', cex=2, bty="n")
  if(a$`Pr(>F)`>=0.055){text(3+0.6,ef$fit[3],'NS')}
  
  # night
  y<-summary(mods.NECN[[i]])$coefficients[,1:2]
  a<-anova(mods.NECN[[i]])
  ef<-as.data.frame(effect("NutLevel", mods.NECN[[i]]))
  
  SE.upper<-ef$fit+NCP.mean.DayNight$SE.AFDW2[NCP.mean.DayNight$Substrate==species& NCP.mean.DayNight$DayNight=='Night']
  SE.lower<-ef$fit-NCP.mean.DayNight$SE.AFDW2[NCP.mean.DayNight$Substrate==species& NCP.mean.DayNight$DayNight=='Night']
  
  points(1:3+0.2,ef$fit, cex=3, pch = 17, col = mypalette)
  lines(1:3+0.2,ef$fit, col = 'black', type = 'c' )
  arrows(1:3+0.2, SE.upper,1:3+0.2, SE.lower, length=0.05, angle=90, code=3,
         col=mypalette, lwd=2)
  if(a$`Pr(>F)`<=0.055){text(3+0.6,ef$fit[3],'*', cex=3)}   #  legend("top",'*', cex=2, bty="n")
  if(a$`Pr(>F)`>=0.055){text(3+0.6,ef$fit[3],'NS')}
  
  axis(1, at=1:3, labels = c("Ambient","Medium","High"), las = 2, cex.lab=1.5,cex.axis=1.5)
  
}

pdf("plots/MSplots/MeanRates_NCP.pdf", width=18, height=5, useDingbats = FALSE)
j<-2
par(bg=NA) 
par(pty="m")
##  plot of all metabolic rates by substrate
par(mfrow=c(1,5))
par(mar=c(6.2,6.1,4.1,1.5))
par(lwd = 2) 

# create a list of all the NCP models
mods.NEC<-list(model.NCPNet.Coral,model.NCPNet.algae, model.NCPNet.Rubble,model.NCPNet.Sand,model.NCPNet.Mixed)
mods.NECD<-list(model.GCP.Coral, model.GCP.algae, model.GCP.Rubble, model.GCP.Sand,model.GCP.Mixed)            
mods.NECN<-list(model.R.Coral, model.R.algae, model.R.Rubble, model.R.Sand,model.R.Mixed)            

for(i in 1:5){
  Nutplot.NCP(substrate[i],i, 'TRUE', 'NCP')
  a<-anova(mods.NECD[[i]])
  NCP.day.stats[i,1:7]<-t(matrix(c(substrate[i],a[1:6])))
  #night
  a<-anova(mods.NECN[[i]])
  NCP.night.stats[i,1:7]<-t(matrix(c(substrate[i],a[1:6])))
  #net
  a<-anova(mods.NEC[[i]])
  NCP.net.stats[i,1:7]<-t(matrix(c(substrate[i],a[1:6])))
}
legend('topright', legend = c('GCP','R','NCP'), pch = c(2,17,19), col = mypalette, bty='n', cex=1.5)
dev.off()

# write the stats to a csv
# NEC and NCP as a function of nutrients models
write.csv(NEC.day.stats, 'stats/NEC.day.stats.csv')
write.csv(NEC.night.stats, 'stats/NEC.night.stats.csv')
write.csv(NEC.net.stats, 'stats/NEC.net.stats.csv')
write.csv(NCP.day.stats, 'stats/NCP.day.stats.csv')
write.csv(NCP.night.stats, 'stats/NCP.night.stats.csv')
write.csv(NCP.net.stats, 'stats/NCP.net.stats.csv')

# pH as a function of species
write.csv(data.frame(summary(model.pH.NCP)$coefficients), 'stats/pH.NCP.stats.csv')

# NEC, omega, nutrient ancova
write.csv(data.frame(summary(model.NECOmega.Algae)$coefficients),'stats/algae.omega.stats.csv')
write.csv(data.frame(summary(model.NECOmega.Coral)$coefficients),'stats/coral.omega.stats.csv')
write.csv(data.frame(summary(model.NECOmega.Rubble)$coefficients),'stats/rubble.omega.stats.csv')
write.csv(data.frame(summary(model.NECOmega.Sand)$coefficients),'stats/sand.omega.stats.csv')
write.csv(data.frame(summary(model.NECOmega.Mixed)$coefficients),'stats/mixed.omega.stats.csv')


#### make plot of temperature and light for supplement
# read in the data and convert the date
TempData<-read.csv('data/Temperature/20151224_Crane_Tank_1_A.csv',stringsAsFactors=FALSE, header=TRUE)
TempData$Date.Time<-as.POSIXct(TempData$Date.Time, format="%m/%d/%Y %H:%M")
#convert all the columns to numeric since excel put in some funky things

# bring in the LiCOR data from Mike to convert LUX to uMol m-2 s-1
G<-read.csv('data/light/PAR_HIMB_MFox.csv')
G$Date.Time<-as.POSIXct(paste(G$Date, G$Time), format="%m/%d/%y %H:%M:%S")

#merge the files so that I can calibrate
TotalLight<-merge(TempData,G, by = 'Date.Time')

#plot Ave quantum yeild by light intensity
plot(TotalLight$IntensityB,TotalLight$Avg.Quantum..uE.sec.m2.)

#fit an exponetial function following Long et al. 2012
mod <- nls(Avg.Quantum..uE.sec.m2. ~ A1*exp(-IntensityB/t1) +y0, 
           data = TotalLight, start = list(A1 = -100, t1 = 100, y0 = 100))

points(TotalLight$IntensityB, predict(mod, list(IntensityB = TotalLight$IntensityB)), col = 'red', pch =19)

# calibrate the light intensities to PAR based on this fit
m<-summary(mod)$coefficients[,1] # pull out the coefficients from the model
TempData[,5:7]<-m[1]*exp(-TempData[,5:7]/m[2]) +m[3] #convert data to PAR

# plot the temperature and light across time for each tank
pdf('plots/MSplots/TempLight.pdf', width = 8, height = 8, useDingbats = FALSE)
par(mfrow=c(2,1))
#Temp
# Time-series
plot(TempData$Date.Time, TempData[,2], col = 'black', type = 'l', xlab = 'Date', ylab = expression(paste('Temperature (',~degree~C, ')')), ylim = c(24,30))
lines(TempData$Date.Time, TempData[,3], col = 'blue', type = 'l')
lines(TempData$Date.Time, TempData[,4], col = 'green', type = 'l')
legend('topright', legend = c('Incubation Tank A','Incubation Tank B','Incubation Tank C'), lty=1, col=c('black','blue','green'), bty='n')
# average
#means
M<-colMeans(TempData[,c(2:4)])
#SDs
S<-apply(TempData[,c(2:4)], 2, sd)
#x<-barplot(M, names.arg = c('Tank A','Tank B','Tank C'), col = c('black','blue','green'), ylim = c(0,30), ylab = expression(paste('Mean Temperature (',~degree~C, ')')))
#arrows(x, M+S, x,M-S, length = 0.05, angle = 90, code = 3)

# light

plot(TempData$Date.Time, TempData[,5], col = 'black', type = 'l', xlab = 'Date', 
     ylab = expression(paste('PAR (',mu,'mol m'^{-2},'s'^{-1},')')))
lines(TempData$Date.Time, TempData[,6], col = 'blue', type = 'l')
lines(TempData$Date.Time, TempData[,7], col = 'green', type = 'l')

#means
M<-apply(TempData[,c(5:7)], 2, max) # max light intensity
#SDs
S<-apply(TempData[,c(5:7)], 2, sd)
#x<-barplot(M, names.arg = c('Tank A','Tank B','Tank C'), col = c('black','blue','green'), ylim = c(0,80000), ylab = 'Max Light Intensity (Lux)')
#arrows(x, M+S, x,M-S, length = 0.05, angle = 90, code = 3)
dev.off()

##Calculate % organic matter from the sediment samples.
Sand$OM<-100*Sand$AFDW/Sand$DW
# calculate means
aggregate(Sand$OM~Sand$Nuts, FUN = mean)
