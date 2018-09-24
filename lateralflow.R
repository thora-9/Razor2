#Input: 1) the current aquifer storage
# 2) the maximum possible aquifer storage
# 3) the threshold at which baseflow starts
# 4) the max baseflow
# 4) Power law coefficient 
# This relationship has been taken from Raven
#Output: the baseflow amount 
lateralflow <- function (AQ_cur, AQ_max, lf_thresh, lf_max, lf_coeff){
  #Threshold driven relationship obtained from Raven manual
  lf_thresh_stor = lf_thresh * AQ_max
  lf_thresh_sat = lf_thresh_stor / AQ_max
  cur_lf = lf_max * (((AQ_cur / AQ_max) - lf_thresh_sat)/(1 - lf_thresh_sat)) ^ lf_coeff
  #Conditions to prevent the baseflow from exceeding max rate or dropping below 0
  if (cur_lf > lf_max) {
    cur_lf = lf_max
  } else if (cur_lf < 0){
    cur_lf = 0
  }
  return(cur_lf)
}
