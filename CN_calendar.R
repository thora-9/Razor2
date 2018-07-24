CN_calendar <- function (LU_details, LU_num, LU_Parameters){
  cal_out=vector()
  LU_out=vector()
  # a=4
  # b=1
  for (a in 1:nrow(LU_details)){
    for (b in 1:LU_num){
    cur_LU=LU_details[[a,b]]
    cur_pars=filter(LU_Parameters,LU==cur_LU)
    plant_month=cur_pars$plant_month
    if(plant_month>0){
    grow_start=as.Date(paste(plant_month,"-15",sep = ""),'%m-%d')
    doy=yday(grow_start)
    grow_end=doy+cur_pars$LI+cur_pars$LD+cur_pars$LM+cur_pars$LL
    #Yearly calendar
    #Get the curve number for Fallow
    fall_CN=LU_Parameters %>% filter(LU=='Fallow') %>% select(CN) %>% as.integer()
    cal1=rep(0,doy) 
    cal2=rep(cur_pars$CN,abs(doy-grow_end))
    if((365-grow_end)>0){
      cal3=rep(0,365-grow_end)
    } else {
      tmp1=grow_end-365
      cal1[1:tmp1]=cur_pars$CN
      cal2=cal2[1:(length(cal2)-tmp1)]
      cal3=NULL}
    cal_all=c(cal1,cal2,cal3)
    } else {cal_all=rep(0,365)}
    #Transpose individual calendars
    cal_out=cbind(cal_out,cal_all)
    colnames(cal_out)[b]=cur_LU
    }
  LU_temp=apply(cal_out,1,sum)
  cal_out=vector()
  LU_out=cbind(LU_out,LU_temp)
  colnames(LU_out)[a]=paste('HRU',a,sep = '_')
  }
  fall_CN=LU_Parameters %>% filter(LU=='Fallow') %>% select(CN) %>% as.integer()
  LU_out[LU_out==0]<-fall_CN
  return(LU_out)
  
}
