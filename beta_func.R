beta_func<- function (S,AW,WP,Sat,RD){
  RD=cur_pars$RD
  AW_mm = AW*RD
  WP_mm = WP*RD
  Sat_mm = Sat*RD
  if (S < WP_mm){
    beta = 0
  } else if (S > WP_mm && S < AW_mm){
    beta = (S - WP_mm)/(AW_mm - WP_mm)
  } else if (S > AW_mm) {
    beta = 1
  }
  
  return (beta)
}