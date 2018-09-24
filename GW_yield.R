#Input: 1) the current soil moisture level
# 2) the maximum possible deep percolation rate
# 3) the current set of parameters associated with the land-use
#Output: the deep percolation amount 
GW_yield <- function (AQ, AQ_max, well_max, wl_thresh){
  # AQ = 400
  # AQ_max = 600
  # well_max = 27 #m3/hour
  # wl_thresh = 0.8 # 80%
  # max_DP=max_percolation
  AQ_ratio = AQ/AQ_max
  if (AQ_ratio < wl_thresh){
    cur_yield = well_max*exp(-4*(1-(AQ_ratio/(wl_thresh))))
  } else {
    cur_yield = well_max*exp(-1*(1-(AQ/AQ_max)))
  }
   #Limit the percolation to max_DP at most
  if (cur_yield > well_max){cur_yield = well_max}
  if (cur_yield < 0){cur_yield = 0}
  return(cur_yield) 
  }
