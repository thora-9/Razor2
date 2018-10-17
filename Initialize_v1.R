#v6.3.4d had continuous irrigation based on the most they could pump 
#v6.3.4c/b probably had the SM deficit based irrigation

require(Kendall)
require(hydroTSM)
require(lubridate)
require(tibble)
require(dplyr)
require(hydroGOF)
require(ggplot2)
require(ggpubr)
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
source('Crop_WD.R')
source('beta_func.R')



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

t2 = zoo(input2[7],date_seq)

t_daily = subdaily2daily(t2,FUN=sum)*1000
t_daily = t_daily[2:163]
#write.csv(t_daily,'Rainfall_daily_tank2_mm.csv')
#t_daily[1:length(t_daily)]=10

########################### Field Data
t1_field_stage = read.csv('tank1_fieldstage.csv',stringsAsFactors = FALSE)

#############################
##Field Data
#############################
field_all = read.csv('tank_stage_rain.csv',stringsAsFactors = FALSE)
rainfall = zoo(field_all$Rainfall,field_all$Date)
t_daily = rainfall

#############################
##Box 1: Runoff generation
#############################

#Catchment 1
LU1_details = tribble(
  ~LU1, ~LU2, ~LU3, ~Per_Area, ~GW_Irr, ~SW_Irr, ~Pumps,
  'Fallow','Fallow','Fallow',20, 'N', 'N',0,
  'Rice', 'Cotton', 'Fallow', 5, 'N', 'N',0,
  'Rice', 'Cotton', 'Fallow', 65, 'Y', 'N', 40,
  'KH_Millet', 'Fallow','Fallow',0, 'N', 'N',0,
  'Juliflora','Fallow','Fallow',10, 'N', 'N',0)

#The number of LU each HRU can have
LU1_num = 3

#Catchment 3
#The command area Per Area is different from the rest of the catchment area; So the HRU with SW_irr 
#as Y will need to add up separately to 100
LU2_details = tribble(
  ~LU1, ~LU2, ~LU3, ~Per_Area, ~GW_Irr, ~SW_Irr, ~Pumps,
  'Fallow','Fallow','Fallow',20, 'N', 'N',0,
  'Rice', 'Cotton', 'Fallow', 60, 'Y', 'N',40,
  'KH_Millet', 'Fallow','Fallow',0, 'N', 'N',0,
  'Juliflora','Fallow','Fallow',20, 'N', 'N',0,
  'Rice', 'Cotton', 'Fallow', 80, 'Y', 'Y',10,
  'Juliflora','Fallow','Fallow',20, 'N', 'Y',0)

LU2_num = 3

#Catchment 3
#The command area Per Area is different from the rest of the catchment area; So the HRU with SW_irr 
#as Y will need to add up separately to 100
LU3_details = tribble(
  ~LU1, ~LU2, ~LU3, ~Per_Area, ~GW_Irr, ~SW_Irr, ~Pumps,
  'Fallow','Fallow','Fallow',20, 'N', 'N',0,
  'Rice', 'Cotton', 'Fallow', 30, 'Y', 'N',40,
  'KH_Millet', 'Fallow','Fallow',30, 'N', 'N',0,
  'Juliflora','Fallow','Fallow',20, 'N', 'N',0,
  'Rice', 'Cotton', 'Fallow', 50, 'Y', 'Y',10,
  'Juliflora','Fallow','Fallow',50, 'N', 'Y',0)

LU3_num = 3

#Catchment 4
#The command area Per Area is different from the rest of the catchment area; So the HRU with SW_irr 
#as Y will need to add up separately to 100
LU4_details = tribble(
  ~LU1, ~LU2, ~LU3, ~Per_Area, ~GW_Irr, ~SW_Irr, ~Pumps,
  'Fallow','Fallow','Fallow',20, 'N', 'N',0,
  'Rice', 'Cotton', 'Fallow', 30, 'Y', 'N',40,
  'KH_Millet', 'Fallow','Fallow',30, 'N', 'N',0,
  'Juliflora','Fallow','Fallow',20, 'N', 'N',0,
  'Rice', 'Cotton', 'Fallow', 50, 'Y', 'Y',10,
  'Juliflora','Fallow','Fallow',50, 'N', 'Y',0)

LU4_num = 3

#For now assuming constant crop depletion factor

LU_Parameters = tribble(
  ~LU,~plant_month,~LI,~LD,~LM,~LL,~KCI,~KCM,~KCL,~RD,~CN,
  'Rice',10,30,35,30,30,1,1.20,0.9,300,9999,
  'KH_Millet',10,20,35,40,30,0.3,1.2,0.5,500,70,
  'Fallow',0,0,0,0,0,0,0,0,100,58,
  'Juliflora',1,365,0,0,0,1.2,0,0,500,80,
  'Cotton',3,30,50,60,55,0.35,1.15,0.6,500,66)

#Rice CN is 9999 (basically implying that the SCS method is invalid for it)
CN1_cal = CN_calendar(LU1_details,LU1_num,LU_Parameters)
Kc1_cal = Kc_calendar(LU1_details,LU1_num,LU_Parameters)

CN2_cal = CN_calendar(LU2_details,LU2_num,LU_Parameters)
Kc2_cal = Kc_calendar(LU2_details,LU2_num,LU_Parameters)

CN3_cal = CN_calendar(LU3_details,LU3_num,LU_Parameters)
Kc3_cal = Kc_calendar(LU3_details,LU3_num,LU_Parameters)

CN4_cal = CN_calendar(LU4_details,LU4_num,LU_Parameters)
Kc4_cal = Kc_calendar(LU4_details,LU4_num,LU_Parameters)

AMC = 0

#This basically converts the orignal curve number formula to mm
#S = (25400/CN)-254

#Input the daily PET rate from Madurai for every month; Estimated using the monthly PET values from VDSA database
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
porosity = 0.01

AQ2_max = 400
AQ2_ini = 250
LF2_thresh = 0.3 #percent of max aquifer storage
LF2_max = 2 #mm/day
LF_coeff = 1 #default 
BF2_thresh = 0.7 #percent of max aquifer storage
BF2_max = 3 #mm/day
BF_coeff = 1 #default 

AQ3_max = 400
AQ3_ini = 250
LF3_thresh = 0.3 #percent of max aquifer storage
LF3_max = 2 #mm/day
LF_coeff = 1 #default 
BF3_thresh = 0.7 #percent of max aquifer storage
BF3_max = 3 #mm/day
BF_coeff = 1 #default 

AQ4_max = 400
AQ4_ini = 250
LF4_thresh = 0.3 #percent of max aquifer storage
LF4_max = 2 #mm/day
LF_coeff = 1 #default 
BF4_thresh = 0.7 #percent of max aquifer storage
BF4_max = 3 #mm/day
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

soil_paddy3 = vector()
soil_paddy3[1] = 0

soil_paddy4 = vector()
soil_paddy4[1] = 0

#############################
##Juliflora specific parameters
#############################
#c corresponds to the e-folding parameter (m^-1) which desribes the decline of relative root density with depth z (Vervoort 2009)
c = 2
LAI = 2.3
PET_J = 10 #(mm/day)
AW_critical = 0.84*Soil_FC
light_par = 0.4

#############################
#Groundwater Irrigation Parameters
#############################

#Used to control the point at which the well yields start declining more rapidly with declining aquifer water levels
wl_thresh = 0.1
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
GW1_change=vector()
RF1_vol = matrix(nrow = length(t_daily),ncol = nrow(LU1_details))
mod_loss = 0
ET_GW = vector()
WT_Z = vector()

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
Spill2 = vector()
RF2_vol = matrix(nrow = length(t_daily),ncol = nrow(LU2_details))
mod_loss = 0
GW2_used = matrix(0,nrow = length(t_daily),ncol = nrow(LU2_details))
GW2_change=vector()

#############################
Q3f = vector()
ET3 = vector()
S1 = vector()
S2 = vector()
S3 = vector()
Qf = vector()
Qu = vector()
Sc3 = vector()
Sc1 = vector()
RF3 = matrix(nrow = length(t_daily),ncol = nrow(LU3_details))
DP3 = matrix(nrow = length(t_daily),ncol = nrow(LU3_details))
ET3 = matrix(nrow = length(t_daily),ncol = nrow(LU3_details))
SM3 = matrix(nrow = length(t_daily),ncol = nrow(LU3_details))
IF3 = matrix(nrow = length(t_daily),ncol = nrow(LU3_details))
AQ3 = vector()
BF3 = vector()
LF3 = vector()
Spill3 = vector()
RF3_vol = matrix(nrow = length(t_daily),ncol = nrow(LU3_details))
mod_loss = 0
GW3_used = matrix(0,nrow = length(t_daily),ncol = nrow(LU3_details))
GW3_change=vector()

#############################
Q4f = vector()
ET4 = vector()
S1 = vector()
S2 = vector()
S4 = vector()
Qf = vector()
Qu = vector()
Sc4 = vector()
Sc1 = vector()
RF4 = matrix(nrow = length(t_daily),ncol = nrow(LU4_details))
DP4 = matrix(nrow = length(t_daily),ncol = nrow(LU4_details))
ET4 = matrix(nrow = length(t_daily),ncol = nrow(LU4_details))
SM4 = matrix(nrow = length(t_daily),ncol = nrow(LU4_details))
IF4 = matrix(nrow = length(t_daily),ncol = nrow(LU4_details))
AQ4 = vector()
BF4 = vector()
LF4 = vector()
Spill4 = vector()
RF4_vol = matrix(nrow = length(t_daily),ncol = nrow(LU4_details))
mod_loss = 0
GW4_used = matrix(0,nrow = length(t_daily),ncol = nrow(LU4_details))
GW4_change=vector()


#Tank variables

#Tank 1
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
#Catchment calculated from (Van Meter, 2016) by subtracting previous tanks catchment
t1_const = as.data.frame(cbind(5e6,3.1842,30,276405,27e4))
colnames(t1_const) = c("max_catch","weir_height","spill_len","max_volume","command")#Units c(m2,meter,meter,m3)

#Tank 2
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
#Catchment calculated from (Van Meter, 2016) by subtracting previous tanks catchment
t2_const = as.data.frame(cbind(11.2e6,3.423,30,407513,45e4))
colnames(t2_const) = c("max_catch","weir_height","spill_len","max_volume","command")#Units c(m2,meter,meter,m3)

#Tank 3
inflow_f3 = vector()
inflow_s3 = vector()
t3_inflow = vector()
t3_precip = vector()
t3_area = vector()
t3_area0 = 0
t3_vol = vector()
t3_vol0 = 0
t3_stage = vector()
t3_area = vector()
t3_spill = vector()
t3_sluice = vector()
t3_GW = vector()
t3_ET = vector()
t3_BF = vector()
#Catchment calculated from (Van Meter, 2016) by subtracting previous tanks catchment
t3_const = as.data.frame(cbind(6.3e6,4.0,30,217633,19e4))
colnames(t3_const) = c("max_catch","weir_height","spill_len","max_volume","command")#Units c(m2,meter,meter,m3)

#Tank 4
inflow_f4 = vector()
inflow_s4 = vector()
t4_inflow = vector()
t4_precip = vector()
t4_area = vector()
t4_area0 = 0
t4_vol = vector()
t4_vol0 = 0
t4_stage = vector()
t4_area = vector()
t4_spill = vector()
t4_sluice = vector()
t4_GW = vector()
t4_ET = vector()
t4_BF = vector()
#Catchment calculated from (Van Meter, 2016) by subtracting previous tanks catchment
t4_const = as.data.frame(cbind(5.9e6,3.3,30,139270,24e4))
colnames(t4_const) = c("max_catch","weir_height","spill_len","max_volume","command")#Units c(m2,meter,meter,m3)




#Com Constants
#Divide by 100 to convert the Percent_Area 
HRU_Areas1 = t(LU1_details$Per_Area)*t1_const$max_catch/100
#############################
LU2_com = subset(LU2_details,LU2_details$SW_Irr=='Y')
HRU_Areas2a = t(LU2_com$Per_Area)*t1_const$command/100

LU2_catch = subset(LU2_details,LU2_details$SW_Irr=='N')
HRU_Areas2b = t(LU2_catch$Per_Area)*(t2_const$max_catch-t1_const$command)/100

HRU_Areas2 = c(HRU_Areas2b, HRU_Areas2a)
#############################
LU3_com = subset(LU3_details,LU3_details$SW_Irr=='Y')
HRU_Areas3a = t(LU3_com$Per_Area)*t2_const$command/100

LU3_catch = subset(LU3_details,LU3_details$SW_Irr=='N')
HRU_Areas3b = t(LU3_catch$Per_Area)*(t3_const$max_catch-t2_const$command)/100

HRU_Areas3 = c(HRU_Areas3b, HRU_Areas3a)
#############################
LU4_com = subset(LU4_details,LU4_details$SW_Irr=='Y')
HRU_Areas4a = t(LU4_com$Per_Area)*t2_const$command/100

LU4_catch = subset(LU4_details,LU4_details$SW_Irr=='N')
HRU_Areas4b = t(LU4_catch$Per_Area)*(t4_const$max_catch-t2_const$command)/100

HRU_Areas4 = c(HRU_Areas4b, HRU_Areas4a)

#############################
#Initialize
#############################
samay = length(t_daily)
i = 40
#For comparing potential vs actual ET
PET_max = vector()
TP_max = vector()

source('RWH_Model_v1.1.R')

par(mfrow =  c(2,2))  

plot(as.Date(index(t_daily)),t1_stage,type='l')
lines(as.Date(index(t_daily)),field_all$Tank.1,col='red')

plot(as.Date(index(t_daily)),t2_stage,type='l')
lines(as.Date(index(t_daily)),field_all$Tank.2,col='red')

plot(as.Date(index(t_daily)),t3_stage,type='l',ylim = c(0,4.5))
lines(as.Date(index(t_daily)),field_all$Tank.3,col='red')

plot(as.Date(index(t_daily)),t4_stage,type='l',ylim = c(0,4.5))
lines(as.Date(index(t_daily)),field_all$Tank.4,col='red')

dev.off()

sum(t1_spill)
sum(t2_spill)

lines(t2_stage,col='red')

source('Plot_v1.R')

NSE.default(t1_stage,field_all$Tank.1)
NSE.default(t2_stage,field_all$Tank.2)
NSE.default(t3_stage,field_all$Tank.3)
NSE.default(t4_stage,field_all$Tank.4)


dev.off()
