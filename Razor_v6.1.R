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
#t_daily[1:length(t_daily)]=10

#############################
##Box 1: Runoff generation
#############################
LU_details = tribble(
                    ~LU1, ~LU2, ~LU3, ~Per_Area,
                    'Rice', 'Cotton', 'Fallow', 50,
                    'KH_Millet', 'Fallow','Fallow',20,
                    'Fallow','Fallow','Fallow',20,
                    'Juliflora','Fallow','Fallow',10)

LU_num = 3                    
              
#For now assuming constant crop depletion factor

LU_Parameters = tribble(
              ~LU,~plant_month,~LI,~LD,~LM,~LL,~KCI,~KCM,~KCL,~RD,~CN,
              'Rice',10,20,30,30,25,1,1.20,0.9,300,9999,
              'KH_Millet',10,20,35,40,30,0.3,1.2,0.5,300,70,
              'Fallow',0,0,0,0,0,0,0,0,100,58,
              'Juliflora',1,365,0,0,0,1.2,0,0,300,80,
              'Cotton',2,30,50,60,55,0.35,1.15,0.6,300,66)

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

#############################
##Box 3: Groundwater box
#############################
maxS3 = 600
ini_S3 = 500

#############################
##Soil Parameters
#############################
RD = 500
Soil_WP = 0.1725
Soil_FC = 0.2825
Soil_sat = 0.415
max_percolation = 4.32 #mm/day; from Gowing paper
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
runoff_vol = matrix(nrow = length(t_daily),ncol = nrow(LU_details))

#Tank variables
inflow_f1 = vector()
inflow_s1 = vector()
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


#############################
#Initialize
#############################
samay = length(t_daily)
i = 1

for(i in 1:10) {
  
  cur_date = t_daily[i]
  doy = yday(cur_date)
  cur_month = month(index(cur_date))#months[i]
  month_name = month.abb[cur_month]
  cur_P = 100#test[i]
  
  #Runoff Generation
  #Estimate the antecendant conditions of the catchment
  if(i > 5){
    rain_5 = sum(t_daily[(i-5):(i-1)]) * 1000
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
      temp_SM = SM1[i-1,j] + cur_P - runoff[i,j]
      #Basically call the function Percolation that uses one of the methods described in Raven
      DP1[i,j] = Percolation(temp_SM,max_percolation,cur_pars)
      #Estimate the evapotranspiration by first estimating the Water stress factor Ks
      TAW1 = (cur_pars$RD * Soil_FC)-(cur_pars$RD * Soil_WP)
      RAW1 = rho1 * TAW1
      if (SM1[i-1,j] < (TAW1-RAW1)){
        Ks = SM1[i-1,j] / (TAW1-RAW1)
        ET1[i,j] = Ks * cur_kc * PET[cur_month]
      } else {
        ET1[i,j] = 1*cur_kc*PET[cur_month]
        }
      } else {
        SM1[i,j] = cur_P - runoff[i,j]
        if (SM1[i,j]>cur_pars$RD*Soil_sat){
           DP1[i,j] = max_percolation
           IF1[i,j] = SM1[i,j] - cur_pars$RD * Soil_sat - max_percolation
           ET1[i,j] = 1 * cur_kc * PET[cur_month]
           SM1[i,j] = SM1[i,j] - DP1[i,j] - IF1[i,j] - ET1[i,j]
        } else if (SM1[i,j] < cur_pars$RD * Soil_sat){
           DP1[i,j]=Percolation(SM1[i,j],max_percolation,cur_pars)
           IF1[i,j] = 0
           ET1[i,j] = 1 * cur_kc * PET[cur_month]
           SM1[i,j] = SM1[i,j] - DP1[i,j] - IF1[i,j] - ET1[i,j]
      } else {
        (print('Error due to Deep percolation calculation'))
        } 
      } 
    } 
  }
}
    else if (cur_LU=='Rice'){
      if (i==1){
        soil_paddy1[i]=0
        
      }
      }
      
    }
  
  
  #Estimating the Soil Moisture Balance for each HRU
  
  
  
  if (cur_P==0){}
  #Bucket 1  
  if (i>1){
    cur_S1=S1[i-1]
  } else {cur_S1=ini_S1}
  
  Q1f[i]=(cur_S1-maxS1)
  if (Q1f[i]<0){Q1f[i]=0}
  Qw[i]=(cur_S1/tw)*h #units of tw are days, so reduce that to timestep of FE
  ET1[i]=cur_PET*(cur_S1/maxS1)*0
  
  S1[i]=cur_S1+(cur_P-Q1f[i]-Qw[i]-ET1[i])
  #Bucker 2  
  if (i>1){
    cur_S2=S2[i-1]
  } else {cur_S2=ini_S2}
  
  ET2[i]=cur_PET*(cur_S2/Se)
  S2[i]=cur_S2+(Qw[i])#-ET2[i])
  if (S2[i]>theta_fc){
    Qr[i]=S2[i]-theta_fc
  } else {Qr[i]=0}
  
  #Bucket 3
  if (i>1){
    cur_S3=S3[i-1]
  } else {cur_S3=ini_S3}
  
  Qx[i]=ET2[i]
  
  S3[i]=cur_S3+(Qr[i]-Qx[i])
  
  #i=i+1
  #print(i)
  
  #Start the tank water balance every after every (1/h) steps  
  if (i%%(1/h)==0){
    j=i*h
    inflow_f1[j]=sum(Q1f[(i-(1/h)):i])
    inflow_s1[j]=0#sum(Q2u[(i-(1/h)):i])*0
    inflow1[j]=inflow_f1[j]+inflow_s1[j]
    
    if(j==1){
      t1_area_cur=t1_area0
    } else {t1_area_cur=t1_area[j-1]}
    
    if(j==1){
      t1_vol_cur=t1_vol0
    } else {t1_vol_cur=t1_vol[j-1]}
    #Rainfall directly on the tank in that timestep
    rain_t1=t1_area_cur*t_daily[j]*(1/1000) #in meters
    
    inflow_vol=rain_t1+((inflow1[j])*(t1_const[1]-t1_area_cur)*(1/1000))
    t1_vol[j]=t1_vol_cur+rain_t1+inflow_vol#((inflow1[j])*(t1_const[1]-t1_area_cur)*(1/1000))
    
    #Stage-Volume relationship from Mike data
    t1_stage[j]=(t1_vol[j]/22914)^(1/1.9461)
    #Stage Area 
    t1_area[j]=42942*(t1_stage[j])^1.0993
    
    # ####Spillway taken from Sakti paper modified to take absolute value of difference
    # if(tank1_const[2]-stage1<0){
    #   form_spill1=1.7*tank1_const[3]*(abs(tank1_const[2]-stage1))^(3/2)*(24*60*60)
    # } else{form_spill1=0}
    
    ####ET from tanks
    t1_ET[j]=t1_area[j]*juli_ET(cur_month)*(1/1000) #m3/day
    
    ####Compare the tank capacity and current volume of water in tank.
    vol_diff1=t1_vol[j]-tank1_const[4]
    if (vol_diff1>=0){
      t1_spill[j]=vol_diff1
    } else{t1_spill[j]=0}
    
    
    ####Sluice Outflow
    Qo1a = ((t1_stage[j]-.785)*5.1903)*86.4 #86.4 Converst L/s to m3/d
    Qo1b = (((t1_stage[j]-1.185)*9.6768)+((t1_stage[j]-1.185)*4.9196))*86.4 #86.4 Converst L/s to m3/d
    
    if (Qo1a<0){Qo1a = 0}
    if (Qo1b<0) {Qo1b = 0}
    
    if (t1_stage[j]<0.785){
      t1_sluice[j] = 0} else if(0.785<t1_stage[j] & t1_stage[j]<1.185) {
        t1_sluice[j] = Qo1a} else if(t1_stage[j]>1.185){
          t1_sluice[j] = (Qo1a + Qo1b)}
    
    ###Spillage from upstream tank
    spill_add1=0
    
    ####GW exchange- Mike paper
    GW_loss1=8.6*t1_stage[j]-6.5 #mm/day
    t1_GW[j]=(GW_loss1/1000)*t1_area[j]#m3/day
    
    #Total Storage change
    t1_vol[j]=t1_vol[j]-(t1_ET[j]+t1_sluice[j]+t1_spill[j]+t1_GW[j])
    
    #Stage-Volume relationship from Mike data
    t1_stage[j]=(t1_vol[j]/22914)^(1/1.9461)
    #Stage Area 
    t1_area[j]=42942*(t1_stage[j])^1.0993
    
    cur_all=cbind(t1_stage[j],t1_area[j],t1_vol[j],inflow1[j],coredata(inflow_vol),t1_ET[j],t1_sluice[j],t1_spill[j],t1_GW[j])
    t1_all=rbind(t1_all,cur_all)
  }
  
}




