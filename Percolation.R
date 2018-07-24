#Input: 1) the current soil moisture level
# 2) the maximum possible deep percolation rate
# 3) the current set of parameters associated with the land-use
#Output: the deep percolation amount 
Percolation <- function (SM, max_DP,cur_pars){
  # SM=SM1[i,j]
  # max_DP=max_percolation
  cur_RD=cur_pars$RD
  WP=0.1725*cur_RD
  FC=0.2825*cur_RD
  sat=0.415*cur_RD
  if (SM > FC){
  cur_dp=max_DP*(SM-FC)/(sat-FC)
  return(cur_dp)
  } else {return(0)}
  }
  
  