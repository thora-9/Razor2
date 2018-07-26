###################################
#Using curve number to get the maximum soil water retention
CN=80 #from Mike
#This basically converts the orignal curve number formula to mm
S=(25400/CN)-254
##################################
i=2
Qr=vector()

for (i in 2:222){
  cur_month=month(t_daily[i])
  cur_rain=coredata(t_daily[i])#*1000  #in mm
  #Making sure the rain >=0
  if(cur_rain<0){
    cur_rain=0
  }
  if(i>5){
    rain_5=sum(t_daily[(i-5):(i-1)])#*1000
    if(rain_5<(0.96*25.4)){
      S_AMC=2.281*S
    } else if(rain_5>(1.6*25.4)){
      S_AMC=0.427*S
    } else{S_AMC=S}
  }
  if(i<5){
    S_AMC=S
  }
  
  ################################TANK 1################################
  
  if(cur_rain>0 & cur_rain>(0.2*S_AMC))
  { Qr[i]=(cur_rain-0.2*S_AMC)^2/(cur_rain+0.8*S_AMC)
    } else {Qr[i]=0}
}
