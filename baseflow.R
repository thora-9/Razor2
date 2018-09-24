#Input: 1) the current aquifer storage
# 2) the maximum possible aquifer storage
# 3) the threshold at which baseflow starts
# 4) the max baseflow
# 4) Power law coefficient 
# This relationship has been taken from Raven
#Output: the baseflow amount 
baseflow <- function (AQ_cur, AQ_max, bf_thresh, bf_max, bf_coeff){
  #Threshold driven relationship obtained from Raven manual
  bf_thresh_stor = bf_thresh * AQ_max
  bf_thresh_sat = bf_thresh_stor / AQ_max
  cur_bf = bf_max * (((AQ_cur / AQ_max) - bf_thresh_sat)/(1 - bf_thresh_sat)) ^ bf_coeff
  #Conditions to prevent the baseflow from exceeding max rate or dropping below 0
  if (cur_bf > bf_max) {
    cur_bf = bf_max
  } else if (cur_bf < 0){
    cur_bf = 0
  }
  return(cur_bf)
}

