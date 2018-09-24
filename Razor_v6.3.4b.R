require(Kendall)
require(hydroTSM)
require(lubridate)
require(tibble)
require(dplyr)
require(hydroGOF)
#source('juli_ET.R')
source('juli_ET_v2.R')
#source('stage1.R')
#source('euler.R')
source('SCS_curve_v2.R')
source('CN_calendar.R')
source('Kc_calendar.R')
source('Percolation.R')
source('Percolation_Rice.R')
source('baseflow.R')
source('lateralflow.R')
source('GW_yield.R')



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

########################### Field Data
t1_field_stage = read.csv('tank1_fieldstage.csv',stringsAsFactors = FALSE)

#############################
##Box 1: Runoff generation
#############################

#Catchment 1
LU1_details = tribble(
                    ~LU1, ~LU2, ~LU3, ~Per_Area, ~GW_Irr, ~SW_Irr, ~Pumps,
                    'Fallow','Fallow','Fallow',30, 'N', 'N',0,
                    'Rice', 'Cotton', 'Fallow', 5, 'N', 'N',0,
                    'Rice', 'Cotton', 'Fallow', 25, 'Y', 'N', 50,
                    'KH_Millet', 'Fallow','Fallow',30, 'N', 'N',0,
                    'Juliflora','Fallow','Fallow',10, 'N', 'N',0)

#The number of LU each HRU can have
LU1_num = 3

#Catchment 2
#The command area Per Area is different from the rest of the catchment area; So the HRU with SW_irr 
#as Y will need to add up separately to 100
LU2_details = tribble(
                    ~LU1, ~LU2, ~LU3, ~Per_Area, ~GW_Irr, ~SW_Irr, ~Pumps,
                    'Fallow','Fallow','Fallow',30, 'N', 'N',0,
                    'Rice', 'Cotton', 'Fallow', 60, 'Y', 'N',30,
                    'KH_Millet', 'Fallow','Fallow',0, 'N', 'N',0,
                    'Juliflora','Fallow','Fallow',10, 'N', 'N',0,
                    'Rice', 'Cotton', 'Fallow', 60, 'Y', 'Y',10,
                    'Juliflora','Fallow','Fallow',40, 'N', 'Y',0)

LU2_num = 3
              
#For now assuming constant crop depletion factor

LU_Parameters = tribble(
              ~LU,~plant_month,~LI,~LD,~LM,~LL,~KCI,~KCM,~KCL,~RD,~CN,
              'Rice',10,20,30,30,25,1,1.20,0.9,300,9999,
              'KH_Millet',10,20,35,40,30,0.3,1.2,0.5,300,70,
              'Fallow',0,0,0,0,0,0,0,0,100,58,
              'Juliflora',1,365,0,0,0,1.2,0,0,300,80,
              'Cotton',2,30,50,60,55,0.35,1.15,0.6,300,66)

#Rice CN is 9999 (basically implying that the SCS method is invalid for it)
CN1_cal = CN_calendar(LU1_details,LU1_num,LU_Parameters)
Kc1_cal = Kc_calendar(LU1_details,LU1_num,LU_Parameters)

CN2_cal = CN_calendar(LU2_details,LU2_num,LU_Parameters)
Kc2_cal = Kc_calendar(LU2_details,LU2_num,LU_Parameters)

AMC = 0

#This basically converts the orignal curve number formula to mm
#S = (25400/CN)-254

#Input the PET rate from Madurai for every month; Estimated using the monthly PET values from VDSA database
PET = c(4,4.7,5.2,5,5.16,5.43,5.16,4.81,4.7,4,3.6,3.7)

#############################
##Box 2: Soil Moisture box
#############################
#Used to calculate the readily available moisture in the soil
rho = 0.3
#The max percolation rate calculated using Saturated hydraulic conductivity
max_percolation = (5*10^-5)*(24*60*60) #mm/day; from Gowing paper
Rice_max_percolation = 7 #mm/day; from Gowing paper

#############################
##Box 3: Groundwater box
#############################
AQ1_max = 400
AQ1_ini = 250
LF1_thresh = 0.3 #percent of max aquifer storage
LF1_max = 2 #mm/day
LF_coeff = 1 #default 
BF1_thresh = 0.7 #percent of max aquifer storage
BF1_max = 3 #mm/day
BF_coeff = 1 #default 

AQ2_max = 400
AQ2_ini = 250
LF2_thresh = 0.3 #percent of max aquifer storage
LF2_max = 2 #mm/day
LF_coeff = 1 #default 
BF2_thresh = 0.7 #percent of max aquifer storage
BF2_max = 3 #mm/day
BF_coeff = 1 #default 


#############################
##Soil Parameters
#############################
RD = 500
Soil_WP = 0.1725
Soil_FC = 0.2825
Soil_sat = 0.415
max_percolation =  4.32 #mm/day; from Gowing paper
TAW1 = Soil_FC - Soil_WP
RAW1 = rho*TAW1

TAW2 = Soil_FC - Soil_WP
RAW2 = rho*TAW2

pad_BH = 70
soil_paddy1 = vector()
soil_paddy1[1] = 0

soil_paddy2 = vector()
soil_paddy2[1] = 0

#Groundwater Irrigation Parameters

#Used to control the point at which the well yields start declining more rapidly with declining aquifer water levels
wl_thresh = 0.5
well_max = 20 #m3/hour; basically the maximum volume that wells can pump
target_SM = 1.3*Soil_FC
electricity = 8

#############################
##Variables
#############################
Q1f = vector()
ET1 = vector()
S1 = vector()
S2 = vector()
S3 = vector()
Qf = vector()
Qu = vector()
RF1 = matrix(nrow = length(t_daily),ncol = nrow(LU1_details))
DP1 = matrix(nrow = length(t_daily),ncol = nrow(LU1_details))
ET1 = matrix(nrow = length(t_daily),ncol = nrow(LU1_details))
SM1 = matrix(nrow = length(t_daily),ncol = nrow(LU1_details))
IF1 = matrix(nrow = length(t_daily),ncol = nrow(LU1_details))
AQ1 = vector()
BF1 = vector()
LF1 = vector()
GW1_used = matrix(0,nrow = length(t_daily),ncol = nrow(LU1_details))
GW2_used = matrix(0,nrow = length(t_daily),ncol = nrow(LU1_details))
RF1_vol = matrix(nrow = length(t_daily),ncol = nrow(LU1_details))
mod_loss = 0

#############################
Q2f = vector()
ET2 = vector()
S1 = vector()
S2 = vector()
S3 = vector()
Qf = vector()
Qu = vector()
Sc2 = vector()
Sc1 = vector()
RF2 = matrix(nrow = length(t_daily),ncol = nrow(LU2_details))
DP2 = matrix(nrow = length(t_daily),ncol = nrow(LU2_details))
ET2 = matrix(nrow = length(t_daily),ncol = nrow(LU2_details))
SM2 = matrix(nrow = length(t_daily),ncol = nrow(LU2_details))
IF2 = matrix(nrow = length(t_daily),ncol = nrow(LU2_details))
AQ2 = vector()
BF2 = vector()
LF2 = vector()
RF2_vol = matrix(nrow = length(t_daily),ncol = nrow(LU2_details))
mod_loss = 0

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
t1_BF = vector()
t1_all = data.frame(matrix(ncol = 8, nrow = 0))
t1_const = as.data.frame(cbind(5e6,3.595,30,276405,27e4))
colnames(t1_const) = c("max_catch","weir_height","spill_len","max_volume","command")#Units c(m2,meter,meter,m3)

inflow_f2 = vector()
inflow_s2 = vector()
t2_inflow = vector()
t2_precip = vector()
t2_area = vector()
t2_area0 = 0
t2_vol = vector()
t2_vol0 = 0
t2_stage = vector()
t2_area = vector()
t2_spill = vector()
t2_sluice = vector()
t2_GW = vector()
t2_ET = vector()
t2_BF = vector()
t2_const = as.data.frame(cbind(16.2e6,3.595,30,407513,45e4))
colnames(t2_const) = c("max_catch","weir_height","spill_len","max_volume","command")#Units c(m2,meter,meter,m3)


#Com Constants
#Divide by 100 to convert the Percent_Area 
HRU_Areas1 = t(LU1_details$Per_Area)*t1_const$max_catch/100

LU2_com = subset(LU2_details,LU2_details$SW_Irr=='Y')
HRU_Areas2a = t(LU2_com$Per_Area)*t1_const$command/100

LU2_catch = subset(LU2_details,LU2_details$SW_Irr=='N')
HRU_Areas2b = t(LU2_catch$Per_Area)*t2_const$max_catch/100

HRU_Areas2 = c(HRU_Areas2b, HRU_Areas2a)


#############################
#Initialize
#############################
samay = length(t_daily)
i = 1
#For comparing potential vs actual ET
PET_max = vector()
TP_max = vector()

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
  
  j1 = 2
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  
  for (j1 in 1:ncol(CN1_cal)){
    cur_CN = CN1_cal[doy,j1]
    cur_pars = filter(LU_Parameters,CN==cur_CN)
    cur_LU = cur_pars$LU
    cur_kc = Kc1_cal[doy,j1]
    GW_irrigation = LU1_details$GW_Irr[j1]
    if(cur_LU != 'Rice'){
      RF1[i,j1] = SCS_curve(cur_CN,cur_P,rain_5)
      RF1_vol[i,j1] = (LU1_details$Per_Area[j1] * t1_const$max_catch) * RF1[i,j1] * (1/1000)
      if (i > 1) {
      #Keep an eye on if I should be calculating Percolation and ET based on the modified SM vs exact SM from previous timestep  
      temp_SM = SM1[i-1,j1] + cur_P - RF1[i,j1]
      #Remove the interflow volume right away (if applicable)
      if (temp_SM>cur_pars$RD*Soil_sat){
        #Purely to ensure that the change in LU doesn't lead to interflow when there is no rainfall due to changing root depth
        if (cur_P == 0 & (CN1_cal[doy,j1] != CN1_cal[doy-1,j1])) {
          mod_loss_temp = (temp_SM - cur_pars$RD*Soil_sat)
          mod_loss = mod_loss + (temp_SM - cur_pars$RD*Soil_sat)
          IF1[i,j1] = 0 
          temp_SM = temp_SM - IF1[i,j1] - mod_loss_temp
          } else {
        IF1[i,j1] = temp_SM - cur_pars$RD*Soil_sat
        temp_SM = temp_SM - IF1[i,j1]
          }
      } else {
        IF1[i,j1] = 0
      }
      #Basically call the function Percolation that uses one of the methods described in Raven
      DP1[i,j1] = Percolation(temp_SM,max_percolation,cur_pars)
      
      #Estimate the evapotranspiration by first estimating the Water stress factor Ks
      TAW1 = (cur_pars$RD * Soil_FC)-(cur_pars$RD * Soil_WP)
      RAW1 = rho * TAW1
      #Adj1usted soil moisture is basically the total SM - WP
      AdSM = temp_SM - (cur_pars$RD * Soil_WP)
      #GW Irrigation Component
      if (cur_LU != 'Fallow') {
        if (GW_irrigation == 'Y'){
          if (temp_SM < target_SM*cur_pars$RD){
            #Estimate the difference between target and current SM levels (in mm)
            #then multiply that by the total Area of that HRU to estimate the volume needed
            GW_volume_req = (target_SM*cur_pars$RD - temp_SM)*HRU_Areas1[j1]*(1/1000)
            GW_available = GW_yield(AQ1[i-1], AQ1_max, well_max, wl_thresh)*LU1_details$Pumps[j1]*electricity
            GW_diff = GW_volume_req - GW_available
            if (GW_diff > 0) {
              temp_SM = temp_SM + (GW_available / HRU_Areas1[j1])*1000
              GW1_used[i,j1] = (GW_available / HRU_Areas1[j1])*1000
            } else if (GW_diff <= 0) {
              temp_SM = target_SM*cur_pars$RD
              GW1_used[i,j1] = (GW_volume_req / HRU_Areas1[j1])*1000
            }
          }
        }
      }
      AdSM = temp_SM - (cur_pars$RD * Soil_WP)
      if (AdSM < (TAW1-RAW1)){
        #Keep an eye on whether to use SM1[i-1,j] or temp_SM
        Ks = AdSM / (TAW1-RAW1)
        ET1[i,j1] = Ks * cur_kc * PET[cur_month]
        if (ET1[i,j1] < 0) {ET1[i,j1] = 0}
      } else {
        ET1[i,j1] = 1*cur_kc*PET[cur_month]
      }
      #The way it is setup right now the SM will never end the day at saturation point
      SM1[i,j1] = temp_SM - DP1[i,j1] - ET1[i,j1]
      } else {
        SM1[i,j1] = (Soil_WP*cur_pars$RD) + cur_P - RF1[i,j1]
        if (SM1[i,j1]>cur_pars$RD*Soil_sat){
           DP1[i,j1] = max_percolation
           IF1[i,j1] = SM1[i,j1] - cur_pars$RD * Soil_sat - max_percolation
           ET1[i,j1] = 1 * cur_kc * PET[cur_month]
             SM1[i,j1] = SM1[i,j1] - DP1[i,j1] - IF1[i,j1] - ET1[i,j1]
        } else if (SM1[i,j1] < cur_pars$RD * Soil_sat){
           DP1[i,j1] = Percolation(SM1[i,j1],max_percolation,cur_pars)
           IF1[i,j1] = 0
           #Need to add a value to scale the AET based on the soil moisture level 
           #Ignore the values generated at i=1
           ET1[i,j1] = 1 * cur_kc * PET[cur_month]
           SM1[i,j1] = SM1[i,j1] - DP1[i,j1] - IF1[i,j1] - ET1[i,j1]
           if (SM1[i,j1] <= 0){SM1[i,j1] = 0}
      } else {
        (print('Error due to Deep percolation calculation'))
        } 
      } 
    } else if (cur_LU=='Rice'){
      temp_SM = SM1[i-1,j1] + cur_P
      Rice_SM_sat = 0.415*cur_pars$RD
      Rice_run_thresh = Rice_SM_sat + pad_BH
      #Calculate the runoff from the rice fields here. Runoff only starts when the total SM value on the rice fields is greater than 
      #the bund height + soil saturation water content
      if (temp_SM > Rice_run_thresh){
        RF1[i,j1] = temp_SM - Rice_run_thresh
        RF1_vol[i,j1] = (LU1_details$Per_Area[j1] * t1_const$max_catch) * RF1[i,j1] * (1/1000)
        temp_SM = temp_SM - RF1[i,j1]
      } else if (temp_SM < Rice_run_thresh){
        RF1[i,j1] = 0
        RF1_vol[i,j1] = (LU1_details$Per_Area[j1] * t1_const$max_catch) * RF1[i,j1] * (1/1000)
      } else {print('Error due to Rice runoff calculation')}
      #Start calculating deep percolation from the rice fields here using the method in Gowing paper
      DP1[i,j1] = Percolation_Rice(temp_SM,Rice_max_percolation,cur_pars)
      #Estimate the evapotranspiration by first estimating the Water stress factor Ks
      TAW1 = (cur_pars$RD * Soil_FC)-(cur_pars$RD * Soil_WP)
      RAW1 = rho * TAW1
      #Adj1usted soil moisture is basically the total SM - WP
      #GW Irrigation Component
      if (cur_LU != 'Fallow') {
        if (GW_irrigation == 'Y'){
          if (temp_SM < target_SM*cur_pars$RD){
            #Estimate the difference between target and current SM levels (in mm)
            #then multiply that by the total Area of that HRU to estimate the volume needed
            GW_volume_req = (target_SM*cur_pars$RD - temp_SM)*HRU_Areas1[j1]*(1/1000)
            #GW yield is in m3/hr, so convert that to per day based on the electricity constant, times the number of pumps
            GW_available = GW_yield(AQ1[i-1], AQ1_max, well_max, wl_thresh)*LU1_details$Pumps[j1]*electricity
            GW_diff = GW_volume_req - GW_available
            if (GW_diff > 0) {
              temp_SM = temp_SM + (GW_available / HRU_Areas1[j1])*1000
              GW1_used[i,j1] = (GW_available / HRU_Areas1[j1])*1000
            } else if (GW_diff <= 0) {
              temp_SM = target_SM*cur_pars$RD
              GW1_used[i,j1] = (GW_volume_req / HRU_Areas1[j1])*1000
            }
          }
        }
      }
      if (temp_SM > Rice_SM_sat){
        AdSM = Rice_SM_sat -(cur_pars$RD * Soil_WP)
      } else if (temp_SM < Rice_SM_sat) {
        AdSM = temp_SM -(cur_pars$RD * Soil_WP)
      } else {print('Error due to Soil Moisture--ET calculation')}
      if (AdSM < (TAW1-RAW1)){
        #Keep an eye on whether to use SM1[i-1,j] or temp_SM
        Ks = AdSM / (TAW1-RAW1)
        ET1[i,j1] = Ks * cur_kc * PET[cur_month]
        if (ET1[i,j1] < 0) {ET1[i,j1] = 0}
      } else {
        ET1[i,j1] = 1*cur_kc*PET[cur_month]
      }
      #Interflow will always be zero
      IF1[i,j1] = 0
      #The way it is setup right now the SM will never end the day at saturation point
      SM1[i,j1] = temp_SM - DP1[i,j1] - ET1[i,j1]
      }
  }
  ############################################################  ############################################################
  ############################################################  ############################################################
  
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
  t1_inflow_temp = (IF1[i,] + RF1[i,]) * HRU_Areas1_temp * (1/1000)
  t1_inflow[i] = sum(t1_inflow_temp)
  
  
  # Use the tank area to estimate the percolation occuring from different HRU's in the catchment area.
  if (i > 1) {
    #Basically calculate the volumetric input into the aquifer bucket based on the size of the HRU's
    #then divide it by the total catchment of the tank (including tank area) to estimate the change in mm
    AQ1_inputs = (DP1[i,] * HRU_Areas1_temp) /(t1_const$max_catch)
    AQ1[i] = AQ1[i-1] + sum(AQ1_inputs)
    LF1[i] = lateralflow (AQ1[i], AQ1_max, LF1_thresh, LF1_max, LF_coeff)
    BF1[i] = baseflow (AQ1[i], AQ1_max, BF1_thresh, BF1_max, BF_coeff)
    GW1_Irr = sum(GW1_used[i,])
    AQ1[i] = AQ1[i] - BF1[i] - LF1[i] - GW1_Irr
  } else {
    AQ1[i] = AQ1_ini + sum(DP1[i,])
    BF1[i] = 0 
    LF1[i] = 0
  }

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
  
  if (Qo1a < 0){Qo1a = 0}
  if (Qo1b < 0) {Qo1b = 0}
  
  if ( t1_temp_stage < 0.785 ){
    t1_sluice[i] = 0 } else if (0.785 < t1_temp_stage & t1_temp_stage < 1.185) {
      t1_sluice[i] = Qo1a} else if( t1_temp_stage > 1.185 ){
        t1_sluice[i] = ( Qo1a + Qo1b )}
  
  ###Spillage from upstream tank
  t1_spill_add = 0
  
  ####Net GW exchange- Mike paper
  GW_loss1 = 8.6 * t1_temp_stage - 6.5 #mm/day
  if (GW_loss1 < 0) {GW_loss1 = 0}
  t1_GW[i] = (GW_loss1/1000) * t1_temp_area #m3/day
  
  ####Baseflow into tank
  t1_BF[i] = (BF1[i]/1000) * t1_temp_area
  
  #Total Storage change
  t1_vol[i] = t1_temp_vol - ( t1_ET[i] + t1_sluice[i] + t1_spill[i] + t1_GW[i] )
  
  #Stage-Volume relationship from Mike data
  t1_stage[i] = (t1_vol[i]/22914) ^ (1/1.9461)
  #Stage Area 
  t1_area[i] = 42942 * (t1_stage[i]) ^ 1.0993
  
  #cur_all=cbind(t1_stage[j],t1_area[j],t1_vol[j],inflow1[j],coredata(inflow_vol),t1_ET[j],t1_sluice[j],t1_spill[j],t1_GW[j])
  #t1_all=rbind(t1_all,cur_all)
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  j2=5
  for (j2 in 1:ncol(CN2_cal)){
    cur_CN = CN2_cal[doy,j2]
    cur_pars = filter(LU_Parameters,CN==cur_CN)
    cur_LU = cur_pars$LU
    cur_kc = Kc2_cal[doy,j2]
    SW2_irrigation = LU2_details$SW_Irr[j2]
    GW2_irrigation = LU2_details$GW_Irr[j2]
    if(cur_LU != 'Rice'){
    #Add the sluice from the upstream tank here for areas that are irrigated by Surface water   
      if (SW2_irrigation == 'Y') {
        #Split the sluice calculated above based on the percent command area of this HRU
        HRU_per_area = LU2_details$Per_Area[j2]/100
        #Convert the sluice volume to depth by dividing by the area of the command HRU
        sluice_depth = t1_sluice[i]*100*HRU_per_area / HRU_Areas2 [j2]
        rain_sluice = cur_P + sluice_depth
      } else {rain_sluice = cur_P}
      
      RF2[i,j2] = SCS_curve(cur_CN, rain_sluice, rain_5)
      RF2_vol[i,j2] = (LU2_details$Per_Area[j2] * t2_const$max_catch) * RF2[i,j2] * (1/1000)
      
      if (i > 1) {
        #Keep an eye on if I should be calculating Percolation and ET based on the modified SM vs exact SM from previous timestep  
        temp_SM = SM2[i-1,j2] + rain_sluice - RF2[i,j2] 
        #Remove the interflow volume right away (if applicable)
        if (temp_SM>cur_pars$RD*Soil_sat){
          #Purely to ensure that the change in LU doesn't lead to interflow when there is no rainfall due to changing root depth
          if (cur_P == 0 & (CN2_cal[doy,j2] != CN2_cal[doy-1,j2])) {
            mod_loss_temp = (temp_SM - cur_pars$RD*Soil_sat)
            mod_loss = mod_loss + (temp_SM - cur_pars$RD*Soil_sat)
            IF2[i,j2] = 0 
            temp_SM = temp_SM - IF2[i,j2] - mod_loss_temp
          } else {
            IF2[i,j2] = temp_SM - cur_pars$RD*Soil_sat
            temp_SM = temp_SM - IF2[i,j2]
          }
        } else {
          IF2[i,j2] = 0
        }
        #Basically call the function Percolation that uses one of the methods described in Raven
        DP2[i,j2] = Percolation(temp_SM,max_percolation,cur_pars)
        
        #Estimate the evapotranspiration by first estimating the Water stress factor Ks
        TAW2 = (cur_pars$RD * Soil_FC)-(cur_pars$RD * Soil_WP)
        RAW2 = rho * TAW2
        #Adj2usted soil moisture is basically the total SM - WP
        AdSM = temp_SM - (cur_pars$RD * Soil_WP)
        if (cur_LU != 'Fallow') {
          if (GW2_irrigation == 'Y'){
            if (temp_SM < target_SM*cur_pars$RD){
              #Estimate the difference between target and current SM levels (in mm)
              #then multiply that by the total Area of that HRU to estimate the volume needed
              GW_volume_req = (target_SM*cur_pars$RD - temp_SM)*HRU_Areas2[j2]*(1/1000)
              GW_available = GW_yield(AQ2[i-1], AQ2_max, well_max, wl_thresh)*LU2_details$Pumps[j2]*electricity
              GW_diff = GW_volume_req - GW_available
              if (GW_diff > 0) {
                temp_SM = temp_SM + (GW_available / HRU_Areas2[j2])*1000
                GW2_used[i,j2] = (GW_available / HRU_Areas2[j2])*1000
              } else if (GW_diff <= 0) {
                temp_SM = target_SM*cur_pars$RD
                GW2_used[i,j2] = (GW_volume_req / HRU_Areas2[j2])*1000
              }
            }
          }
        }
        if (AdSM < (TAW2-RAW2)){
          #Keep an eye on whether to use SM1[i-1,j] or temp_SM
          Ks = AdSM / (TAW2-RAW2)
          ET2[i,j2] = Ks * cur_kc * PET[cur_month]
          if (ET2[i,j2] < 0) {ET2[i,j2] = 0}
        } else {
          ET2[i,j2] = 1*cur_kc*PET[cur_month]
        }
        #The way it is setup right now the SM will never end the day at saturation point
        SM2[i,j2] = temp_SM - DP2[i,j2] - ET2[i,j2]
      } else {
        SM2[i,j2] = (Soil_WP*cur_pars$RD)+ rain_sluice - RF2[i,j2]
        if (SM2[i,j2]>cur_pars$RD*Soil_sat){
          DP2[i,j2] = max_percolation
          IF2[i,j2] = SM2[i,j2] - cur_pars$RD * Soil_sat - max_percolation
          ET2[i,j2] = 1 * cur_kc * PET[cur_month]
          SM2[i,j2] = SM2[i,j2] - DP2[i,j2] - IF2[i,j2] - ET2[i,j2]
        } else if (SM2[i,j2] < cur_pars$RD * Soil_sat){
          DP2[i,j2] = Percolation(SM2[i,j2],max_percolation,cur_pars)
          IF2[i,j2] = 0
          #Need to add a value to scale the AET based on the soil moisture level 
          #Ignore the values generated at i=1
          ET2[i,j2] = 1 * cur_kc * PET[cur_month]
          SM2[i,j2] = SM2[i,j2] - DP2[i,j2] - IF2[i,j2] - ET2[i,j2]
          if (SM2[i,j2] <= 0){SM2[i,j2] = 0}
        } else {
          (print('Error due to Deep percolation calculation'))
        } 
      } 
    } else if (cur_LU=='Rice'){
      
      #Add the sluice from the upstream tank here for areas that are irrigated by Surface water   
      if (SW2_irrigation == 'Y') {
        #Split the sluice calculated above based on the percent command area of this HRU
        HRU_per_area = LU2_details$Per_Area[j2]/100
        #Convert the sluice volume to depth by dividing by the area of the command HRU; convert meter to mm
        sluice_depth = t1_sluice[i]*100*HRU_per_area / HRU_Areas2 [j2]
        rain_sluice = cur_P + sluice_depth
      } else {rain_sluice = cur_P}
      
      temp_SM = SM2[i-1,j2] + rain_sluice
      Rice_SM_sat = 0.415*cur_pars$RD
      Rice_run_thresh = Rice_SM_sat + pad_BH
      #Calculate the runoff from the rice fields here. Runoff only starts when the total SM value on the rice fields is greater than 
      #the bund height + soil saturation water content
      if (temp_SM > Rice_run_thresh){
        RF2[i,j2] = temp_SM - Rice_run_thresh
        RF2_vol[i,j2] = (LU2_details$Per_Area[j2] * t2_const$max_catch) * RF2[i,j2] * (1/1000)
        temp_SM = temp_SM - RF2[i,j2]
      } else if (temp_SM < Rice_run_thresh){
        RF2[i,j2] = 0
        RF2_vol[i,j2] = (LU2_details$Per_Area[j2] * t2_const$max_catch) * RF2[i,j2] * (1/1000)
      } else {print('Error due to Rice runoff calculation')}
      #Start calculating deep percolation from the rice fields here using the method in Gowing paper
      DP2[i,j2] = Percolation_Rice(temp_SM,Rice_max_percolation,cur_pars)
      #Estimate the evapotranspiration by first estimating the Water stress factor Ks
      TAW2 = (cur_pars$RD * Soil_FC)-(cur_pars$RD * Soil_WP)
      RAW2 = rho * TAW2
      #Adj2usted soil moisture is basically the total SM - WP
      if (cur_LU != 'Fallow') {
        if (GW2_irrigation == 'Y'){
          if (temp_SM < target_SM*cur_pars$RD){
            #Estimate the difference between target and current SM levels (in mm)
            #then multiply that by the total Area of that HRU to estimate the volume needed
            GW_volume_req = (target_SM*cur_pars$RD - temp_SM)*HRU_Areas2[j2]*(1/1000)
            GW_available = GW_yield(AQ2[i-1], AQ2_max, well_max, wl_thresh)*LU2_details$Pumps[j2]*electricity
            GW_diff = GW_volume_req - GW_available
            if (GW_diff > 0) {
              temp_SM = temp_SM + (GW_available / HRU_Areas2[j2])*1000
              GW2_used[i,j2] = (GW_available / HRU_Areas2[j2])*1000
            } else if (GW_diff <= 0) {
              temp_SM = target_SM*cur_pars$RD
              GW2_used[i,j2] = (GW_volume_req / HRU_Areas2[j2])*1000
            }
          }
        }
      }
      if (temp_SM > Rice_SM_sat){
        AdSM = Rice_SM_sat -(cur_pars$RD * Soil_WP)
      } else if (temp_SM < Rice_SM_sat) {
        AdSM = temp_SM -(cur_pars$RD * Soil_WP)
      } else {print('Error due to Soil Moisture--ET calculation')}
      if (AdSM < (TAW2-RAW2)){
        #Keep an eye on whether to use SM1[i-1,j] or temp_SM
        Ks = AdSM / (TAW2-RAW2)
        ET2[i,j2] = Ks * cur_kc * PET[cur_month]
        if (ET2[i,j2] < 0) {ET2[i,j2] = 0}
      } else {
        ET2[i,j2] = 1*cur_kc*PET[cur_month]
      }
      #Interflow will always be zero
      IF2[i,j2] = 0
      #The way it is setup right now the SM will never end the day at saturation point
      SM2[i,j2] = temp_SM - DP2[i,j2] - ET2[i,j2]
    }
  }
  
  HRU_Areas2_temp = HRU_Areas2
  if (i > 1) {
    #Basically calculate the volumetric input into the aquifer bucket based on the size of the HRU's
    #then divide it by the total catchment of the tank (including tank area) to estimate the change in mm
    AQ2_inputs = (DP2[i,] * HRU_Areas2_temp) /(t2_const$max_catch)
    AQ2[i] = AQ2[i-1] + sum(AQ2_inputs) + 1000*(t1_GW[i]/(t2_const$max_catch))
    LF2[i] = lateralflow (AQ2[i], AQ2_max, LF2_thresh, LF2_max, LF_coeff)
    BF2[i] = baseflow (AQ2[i], AQ2_max, BF2_thresh, BF2_max, BF_coeff)
    GW2_Irr = sum(GW2_used[i,])
    #GW1_Irr = sum(GW1_used[i,])
    AQ2[i] = AQ2[i] - BF2[i] - LF2[i] - GW2_Irr
  } else {
    AQ2[i] = AQ1_ini + sum(DP2[i,])
    BF2[i] = 0 
    LF2[i] = 0
  }
  
  
  PET_max[i] = Kc2_cal[doy,5]*PET[cur_month]
  TP_max[i] = Kc2_cal[i,5]
  
  i=i+1
}


par(mfrow = c(2,2))  

plot(t_daily,type = 'l')

plot(t1_sluice,type = 'l')

plot(SM1[,2],type = 'l',ylim = c(0,200))
lines(SM1[,3],col='green')
lines(SM2[,5])

plot(ET1[,3],type = 'l',col='green',ylim = c(0,7))
lines(ET1[,2])
lines(ET2[,5])
lines(PET_max,col='red')

plot(t1_stage[12:110],type = 'l',ylim=c(0,4))
lines(t1_field_stage$T1_stage.m.[12:110],col = 'purple')

NSE.default(t1_stage[12:110],t1_field_stage$T1_stage.m.[12:110])
dev.off()
plot(AQ1,col='green',type = 'l',ylim=c(0,400))
lines(AQ2)

plot(AQ2,type = 'l',ylim=c(0,400))

plot(GW1_used[,3],type = 'l')

AQ1old = AQ1


i=82
j2=5

