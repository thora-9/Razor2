juli_ET<-function(mth){
  grp1=c(6,7,8,9)
  grp2=c(9,10,11,12)
  grp3=c(1,2)
  grp4=c(3,4)
  if (mth %in% grp1){
    cur_jET=4.57
  } else if (mth %in% grp2){
    cur_jET=3.9
  } else if (mth %in% grp3){
    cur_jET=5.31
  } else if (mth %in% grp4){
    cur_jET=6.8
  }
  return(cur_jET)
}

# juli_ET<-function(mth){
#   grp1=c(1)
#   grp2=c(2,3,4,5,6)
#   grp3=c(7,8,9,10)
#   grp4=c(11,12)
#   if (mth %in% grp1){
#     cur_jET=9.1
#   } else if (mth %in% grp2){
#     cur_jET=6.8
#   } else if (mth %in% grp3){
#     cur_jET=8.1
#   } else if (mth %in% grp4){
#     cur_jET=9.1
#   }
#   return(cur_jET)
# }