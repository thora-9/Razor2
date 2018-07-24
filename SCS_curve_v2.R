SCS_curve <- function(cur_CN,cur_P,rain_5){
  #Convert the CN to retention parameter using emperical equation
  S=(25400/cur_CN)-254
  #Change the retnetion parameter based on the antecedent conditions
  #AMC 1
  if(rain_5<(0.96*25.4)){
    S_AMC=2.281*S
  } else if(rain_5>(1.6*25.4)){ #AMC 3
    S_AMC=0.427*S
  } else{S_AMC=S} #AMC 2
  
  #Estimate the daily runoff based on the retention parameter value
  if(cur_P>0 & cur_P>(0.2*S_AMC)){
     Qr=(cur_P-0.2*S_AMC)^2/(cur_P+0.8*S_AMC)
  } else {Qr=0}
  
  return(Qr)
  }