Crop_WD<- function (Crop, LU_Parameters,cur_kc){
  if(Crop=='Rice'){
    if (cur_kc == LU_Parameters$KCI){
      #Use the nursery water requirement + Main field preparation requirement from TNAU webpage
      cwd = (40+200)/(LU_Parameters$LI)
    } else if (cur_kc == LU_Parameters$KCM){
      #Use the nursery water requirement + Main field preparation requirement from TNAU webpage
      cwd = (458+417)/(LU_Parameters$LD+LU_Parameters$LM)
    } else if (cur_kc == LU_Parameters$KCL){
      #Use the nursery water requirement + Main field preparation requirement from TNAU webpage
      cwd = (125)/(LU_Parameters$LL)
    }
  } else if (Crop == 'Cotton'){
    #Just use the stock all season water requirement from TNAU website and divide by growth period
    cwd = 600/(165)
  } 
  
  return(cwd) 
}
