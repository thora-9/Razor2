require(Kendall)
require(hydroTSM)
require(lubridate)
require(tibble)
require(dplyr)
source('juli_ET.R')
source('juli_ET_v2.R')
source('stage1.R')
source('euler.R')
source('SCS_curve_v2.R')
source('CN_calendar.R')
source('Kc_calendar.R')
source('Percolation.R')
source('Percolation_Rice.R')



########################## Long term rainfall
input_file = "Tirumungulam_rainfall_modified_v2.csv"
input = read.csv(input_file,header=TRUE)
in_date = as.Date(input[,1],"%Y-%m-%d")
x = zoo(input[,2],in_date)
t = window(x,start=as.Date(input[1,1]))

########################## 2013 field collected rainfall
input_file = "Hourly_rain.csv"
input2 = read.csv(input_file,header=TRUE)
date_seq = seq(as.POSIXct("2013-09-26 00:00"),by="hour",length.out = 5324)

t2 = zoo(input2[4],date_seq)

t_daily = subdaily2daily(t2,FUN=sum)*1000
t_daily = t_daily[1:163]
#t_daily[1:length(t_daily)]=10

#############################
##Box 1: Runoff generation
#############################
LU_details = tribble(
                    ~LU1, ~LU2, ~LU3, ~Per_Area,
                    'Fallow','Fallow','Fallow',30,
                    'Rice', 'Cotton', 'Fallow', 50,
                    'KH_Millet', 'Fallow','Fallow',0,
                    'Juliflora','Fallow','Fallow',20)

LU_num = 3                    
              
#For now assuming constant crop depletion factor

LU_Parameters = tribble(
              ~LU,~plant_month,~LI,~LD,~LM,~LL,~KCI,~KCM,~KCL,~RD,~CN,
              'Rice',10,20,30,30,25,1,1.20,0.9,300,9999,
              'KH_Millet',10,20,35,40,30,0.3,1.2,0.5,300,70,
              'Fallow',0,0,0,0,0,0,0,0,100,58,
              'Juliflora',1,365,0,0,0,1.2,0,0,300,80,
              'Cotton',2,30,50,60,55,0.35,1.15,0.6,300,66)
#Rice CN is 9999 (basically implying that the SCS method is invalid for it)
CN_cal = CN_calendar(LU_details,LU_num,LU_Parameters)
Kc_cal = Kc_calendar(LU_details,LU_num,LU_Parameters)

AMC = 0

#This basically converts the orignal curve number formula to mm
#S = (25400/CN)-254

#Input the PET rate from Madurai for every month
PET = c(4.1,4.4,5.7,5.2,5.3,4.7,4.4,4.6,4.6,3.9,3.5,3.7)

#############################
##Box 2: Soil Moisture box
#############################
#Used to calculate the readily available moisture in the soil
rho1 = 0.3
#The max percolation rate calculated using Saturated hydraulic conductivity
max_percolation = (5*10^-5)*(24*60*60) #mm/day; from Gowing paper
Rice_max_percolation = 7 #mm/day; from Gowing paper

#############################
##Box 3: Groundwater box
#############################
AQ1_max = 600
AQ1_ini = 500

#############################
##Soil Parameters
#############################
RD = 500
Soil_WP = 0.1725
Soil_FC = 0.2825
Soil_sat = 0.415
max_percolation =  4.32 #mm/day; from Gowing paper
TAW1 = Soil_FC - Soil_WP
RAW1 = rho1*TAW1

pad_BH = 70
soil_paddy1 = vector()
soil_paddy1[1] = 0

#############################
##Variables
#############################
Q1f = vector()
ET1 = vector()
ET2 = vector()
Qw = vector()
Qr = vector()
Qx = vector()
S1 = vector()
S2 = vector()
S3 = vector()
Qf = vector()
Qu = vector()
Sc1 = vector()
Sc2 = vector()
runoff = matrix(nrow = length(t_daily),ncol = nrow(LU_details))
DP1 = matrix(nrow = length(t_daily),ncol = nrow(LU_details))
ET1 = matrix(nrow = length(t_daily),ncol = nrow(LU_details))
SM1 = matrix(nrow = length(t_daily),ncol = nrow(LU_details))
IF1 = matrix(nrow = length(t_daily),ncol = nrow(LU_details))
AQ1 = vector()
runoff_vol = matrix(nrow = length(t_daily),ncol = nrow(LU_details))


#Tank variables
inflow_f1 = vector()
inflow_s1 = vector()
t1_inflow = vector()
t1_precip = vector()
t1_area = vector()
t1_area0 = 0
t1_vol = vector()
t1_vol0 = 0
t1_stage = vector()
t1_area = vector()
t1_spill = vector()
t1_sluice = vector()
t1_GW = vector()
t1_ET = vector()
t1_all = data.frame(matrix(ncol = 8, nrow = 0))
t1_const = as.data.frame(cbind(5e6,3.595,30,276405))
colnames(t1_const) = c("max_catch","weir_height","spill_len","max_volume")#Units c(m2,meter,meter,m3)
#Com Constants
HRU_Areas1 = t(LU_details$Per_Area)*t1_const$max_catch/100
  

#############################
#Initialize
#############################
samay = length(t_daily)
i = 1

for(i in 1:samay) {
  
  cur_date = t_daily[i]
  doy = yday(cur_date)
  cur_month = month(index(cur_date))#months[i]
  month_name = month.abb[cur_month]
  cur_P = coredata(t_daily[i]) #100#test[i]
  
  #Runoff Generation
  #Estimate the antecendant conditions of the catchment
  if(i > 5){
    #rain_5 = sum(t_daily[(i-5):(i-1)]) * 1000
    rain_5 = 500
  } else {rain_5 = 30}
  
  j = 1
  
  for (j in 1:ncol(CN_cal)){
    cur_CN = CN_cal[doy,j]
    cur_pars = filter(LU_Parameters,CN==cur_CN)
    cur_LU = cur_pars$LU
    cur_kc = Kc_cal[doy,j]
    if(cur_LU != 'Rice'){
      runoff[i,j] = SCS_curve(cur_CN,cur_P,rain_5)
      runoff_vol[i,j] = (LU_details$Per_Area[j] * t1_const$max_catch) * runoff[i,j] * (1/1000)
      if (i > 1) {
      #Keep an eye on if I should be calculating Percolation and ET based on the modified SM vs exact SM from previous timestep  
      temp_SM = SM1[i-1,j] + cur_P - runoff[i,j]
      #Remove the interflow volume right away (if applicable)
      if (temp_SM>cur_pars$RD*Soil_sat){
        IF1[i,j] = temp_SM - cur_pars$RD*Soil_sat
        temp_SM = temp_SM - IF1[i,j]
      } else {
        IF1[i,j] = 0
      }
      #Basically call the function Percolation that uses one of the methods described in Raven
      DP1[i,j] = Percolation(temp_SM,max_percolation,cur_pars)
      #Estimate the evapotranspiration by first estimating the Water stress factor Ks
      TAW1 = (cur_pars$RD * Soil_FC)-(cur_pars$RD * Soil_WP)
      RAW1 = rho1 * TAW1
      #Adjusted soil moisture is basically the total SM - WP
      AdSM = temp_SM - (cur_pars$RD * Soil_WP)
      if (AdSM < (TAW1-RAW1)){
        #Keep an eye on whether to use SM1[i-1,j] or temp_SM
        Ks = AdSM / (TAW1-RAW1)
        ET1[i,j] = Ks * cur_kc * PET[cur_month]
        if (ET1[i,j] < 0) {ET1[i,j] = 0}
      } else {
        ET1[i,j] = 1*cur_kc*PET[cur_month]
      }
      #The way it is setup right now the SM will never end the day at saturation point
      SM1[i,j] = temp_SM - DP1[i,j] - ET1[i,j]
      } else {
        SM1[i,j] = cur_P - runoff[i,j]
        if (SM1[i,j]>cur_pars$RD*Soil_sat){
           DP1[i,j] = max_percolation
           IF1[i,j] = SM1[i,j] - cur_pars$RD * Soil_sat - max_percolation
           ET1[i,j] = 1 * cur_kc * PET[cur_month]
           SM1[i,j] = SM1[i,j] - DP1[i,j] - IF1[i,j] - ET1[i,j]
        } else if (SM1[i,j] < cur_pars$RD * Soil_sat){
           DP1[i,j] = Percolation(SM1[i,j],max_percolation,cur_pars)
           IF1[i,j] = 0
           #Need to add a value to scale the AET based on the soil moisture level 
           #Ignore the values generated at i=1
           ET1[i,j] = 1 * cur_kc * PET[cur_month]
           SM1[i,j] = SM1[i,j] - DP1[i,j] - IF1[i,j] - ET1[i,j]
           if (SM1[i,j] <= 0){SM1[i,j] = 0}
      } else {
        (print('Error due to Deep percolation calculation'))
        } 
      } 
    } else if (cur_LU=='Rice'){
      temp_SM = SM1[i-1,j] + cur_P
      Rice_SM_sat = 0.415*cur_pars$RD
      Rice_run_thresh = Rice_SM_sat + pad_BH
      #Calculate the runoff from the rice fields here. Runoff only starts when the total SM value on the rice fields is greater than 
      #the bund height + soil saturation water content
      if (temp_SM > Rice_run_thresh){
        runoff[i,j] = temp_SM - Rice_run_thresh
        runoff_vol[i,j] = (LU_details$Per_Area[j] * t1_const$max_catch) * runoff[i,j] * (1/1000)
        temp_SM = temp_SM - runoff[i,j]
      } else if (temp_SM < Rice_run_thresh){
        runoff[i,j] = 0
        runoff_vol[i,j] = (LU_details$Per_Area[j] * t1_const$max_catch) * runoff[i,j] * (1/1000)
      } else {print('Error due to Rice runoff calculation')}
      #Start calculating deep percolation from the rice fields here using the method in Gowing paper
      DP1[i,j] = Percolation_Rice(temp_SM,Rice_max_percolation,cur_pars)
      #Estimate the evapotranspiration by first estimating the Water stress factor Ks
      TAW1 = (cur_pars$RD * Soil_FC)-(cur_pars$RD * Soil_WP)
      RAW1 = rho1 * TAW1
      #Adjusted soil moisture is basically the total SM - WP
      if (temp_SM > Rice_SM_sat){
        AdSM = Rice_SM_sat -(cur_pars$RD * Soil_WP)
      } else if (temp_SM < Rice_SM_sat) {
        AdSM = temp_SM -(cur_pars$RD * Soil_WP)
      } else {print('Error due to Soil Moisture--ET calculation')}
      if (AdSM < (TAW1-RAW1)){
        #Keep an eye on whether to use SM1[i-1,j] or temp_SM
        Ks = AdSM / (TAW1-RAW1)
        ET1[i,j] = Ks * cur_kc * PET[cur_month]
        if (ET1[i,j] < 0) {ET1[i,j] = 0}
      } else {
        ET1[i,j] = 1*cur_kc*PET[cur_month]
      }
      #Interflow will always be zero
      IF1[i,j] = 0
      #The way it is setup right now the SM will never end the day at saturation point
      SM1[i,j] = temp_SM - DP1[i,j] - ET1[i,j]
      }
  }
  #Estimate the deep percolation from all the HRU's
  if (i > 1) {AQ1[i] = AQ1[i-1] + sum(DP1[i,])
  } else {AQ1[1] = AQ1_ini + sum(DP1[i,])}
  
  #Tank1 Start
  
  #Get the tank area from the previous timestep
  if (i == 1) {
    t1_area[1] = 0
    t1_vol[1] = 0
    t1_temp_area = 0
    t1_temp_vol = 0
  } else {
    t1_temp_area = t1_area[i-1]
    t1_temp_vol = t1_vol[i-1]
  } 
  #Estimate the area for each of the HRU's (special HRU is Fallow which converts to tank)
  HRU_Areas1_temp = HRU_Areas1
  HRU_Areas1_temp[1] = HRU_Areas1_temp[1] - t1_temp_area
  t1_inflow_temp = (IF1[i,] + runoff[i,]) * HRU_Areas1_temp * (1/1000)
  t1_inflow[i] = sum(t1_inflow_temp)
  
  #Estimate the water added to the tank by direct precipitation
  t1_precip[i] = cur_P * t1_temp_area *(1/1000)
  
  #Update the temp_tank volume to include the inputs calculated above
  t1_temp_vol = t1_temp_vol + t1_precip[i] + t1_inflow[i]
  
  #Update the t1_temp area and then subsequently the fallow area
  #Stage-Volume relationship from Mike data
  t1_temp_stage = (t1_temp_vol/22914)^(1/1.9461)
  #Stage Area 
  t1_temp_area=42942*(t1_temp_stage)^1.0993
  #HRU_Areas update
  HRU_Areas1_temp[1] = HRU_Areas1[1] - t1_temp_area
  
  #ET from tank
  t1_ET[i] = t1_temp_area * PET[cur_month] * (1/1000) #m3/day
  
  ####Compare the tank capacity and current volume of water in tank.
  vol_diff1 = t1_temp_vol - t1_const$max_volume
  if (vol_diff1 >= 0){
    t1_spill[i] = vol_diff1
  } else{ t1_spill[i] = 0 }
  
  ####Sluice Outflow
  Qo1a = (( t1_temp_stage - 0.785 ) * 5.1903 ) * 86.4 #86.4 Converst L/s to m3/d
  Qo1b = ((( t1_temp_stage - 1.185 ) * 9.6768 ) + (( t1_temp_stage - 1.185 ) * 4.9196 )) * 86.4 #86.4 Converst L/s to m3/d
  
  if ( Qo1a < 0 ){ Qo1a = 0 }
  if ( Qo1b < 0 ) { Qo1b = 0 }
  
  if ( t1_temp_stage < 0.785 ){
    t1_sluice[i] = 0 } else if (0.785 < t1_temp_stage & t1_temp_stage < 1.185) {
      t1_sluice[i] = Qo1a} else if( t1_temp_stage > 1.185 ){
        t1_sluice[i] = ( Qo1a + Qo1b )}
  
  ###Spillage from upstream tank
  t1_spill_add = 0
  
  ####GW exchange- Mike paper
  GW_loss1 = 8.6 * t1_temp_stage - 6.5 #mm/day
  if ( GW_loss1 < 0) { GW_loss1 = 0}
  t1_GW[i] = ( GW_loss1/1000 ) * t1_temp_area #m3/day
  
  #Total Storage change
  t1_vol[i] = t1_temp_vol - ( t1_ET[i] + t1_sluice[i] + t1_spill[i] + t1_GW[i] )
  
  #Stage-Volume relationship from Mike data
  t1_stage[i] = (t1_vol[i]/22914) ^ (1/1.9461)
  #Stage Area 
  t1_area[i] = 42942 * (t1_stage[i]) ^ 1.0993
  
  #cur_all=cbind(t1_stage[j],t1_area[j],t1_vol[j],inflow1[j],coredata(inflow_vol),t1_ET[j],t1_sluice[j],t1_spill[j],t1_GW[j])
  #t1_all=rbind(t1_all,cur_all)
  
  
  i=i+1
 }
  


plot(t1_stage[12:163],type = 'l',ylim=c(0,4))
lines(t1_field_stage$T1_stage.m.[12:174],col = 'purple')

plot(AQ1,type = 'l')
  
  #Estimating the Soil Moisture Balance for each HRU
  
inflow1 = IF1 + runoff

inflow2 = apply(inflow1, 1, mean)

plot(inflow2, type = 'l')
  
t1_field_stage = read.csv('tank1_fieldstage.csv',stringsAsFactors = FALSE)


plot(t1_stage,type = 'l')
lines(t1_field_stage$T1_stage.m.,col = 'purple')

apply(input2,2,sum)

