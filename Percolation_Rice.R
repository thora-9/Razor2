#Input: 1) the current soil moisture level
# 2) the maximum possible deep percolation rate
# 3) the current set of parameters associated with the land-use
#Output: the deep percolation amount 
Percolation_Rice <- function (SM, max_DP,cur_pars){
  #SM=200#temp_SM#SM1[i,j]
  #max_DP=Rice_max_percolation #max_percolation
  cur_RD=cur_pars$RD
  #These number are obtained from the Gowing paper 
  WP=0.1725*cur_RD
  FC=0.2825*cur_RD
  sat=0.415*cur_RD

  if (SM > FC){
    cur_dp = max_DP * (SM - FC)/(sat - FC)
    #Limit the percolation to max_DP at most
    if (cur_dp > max_DP){cur_dp = max_DP}
    return(cur_dp)
  } else {return(0)}
}
