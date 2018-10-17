for(i in 1:samay) {
  
  cur_date = t_daily[i]
  doy = yday(cur_date)
  cur_month = month(index(cur_date))#months[i]
  month_name = month.abb[cur_month]
  cur_P = coredata(t_daily[i]) #100#test[i]
  
  #Runoff Generation
  #Estimate the antecendant conditions of the catchment
  if(i > 5){
    #rain_5 = sum(t_daily[(i-5):(i-1)]) * 1000
    rain_5 = 500
  } else {rain_5 = 30}
  
  j1 = 3
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  
  for (j1 in 1:ncol(CN1_cal)){
    cur_CN = CN1_cal[doy,j1]
    cur_pars = filter(LU_Parameters,CN==cur_CN)
    cur_LU = cur_pars$LU
    cur_kc = Kc1_cal[doy,j1]
    GW_irrigation = LU1_details$GW_Irr[j1]
    if(cur_LU != 'Rice'){
      RF1[i,j1] = SCS_curve(cur_CN,cur_P,rain_5)
      RF1_vol[i,j1] = (LU1_details$Per_Area[j1] * t1_const$max_catch) * RF1[i,j1] * (1/1000)
      if (i > 1) {
        #Keep an eye on if I should be calculating Percolation and ET based on the modified SM vs exact SM from previous timestep  
        temp_SM = SM1[i-1,j1] + cur_P - RF1[i,j1]
        #GW Irrigation Component
        #if (cur_LU != 'Fallow') {
        if (GW_irrigation == 'Y' && cur_LU != 'Fallow'){
          #Estimate the difference between target and current SM levels (in mm)
          #then multiply that by the total Area of that HRU to estimate the volume needed
          GW_available = GW_yield(AQ1[i-1], AQ1_max, well_max, wl_thresh)*LU1_details$Pumps[j1]*electricity
          cwd = Crop_WD(cur_LU,cur_pars,cur_kc)
          if (cur_P < cwd) {
            GW_req = (cwd - cur_P)*(1/1000)*(HRU_Areas1[j1])
            if (GW_req > GW_available){
              temp_SM = temp_SM + (GW_available / HRU_Areas2[j2])*1000
              GW1_used[i,j1] = (GW_available / t1_const$max_catch)*1000
              GW1_change[i] = (GW_available / HRU_Areas1[j1])*1000
            } else if (GW_req < GW_available) {
              temp_SM = temp_SM + (GW_req / HRU_Areas2[j2])*1000
              GW1_used[i,j1] = (GW_req / t1_const$max_catch)*1000
              GW1_change[i] = (GW_req / HRU_Areas1[j1])*1000
            } 
          } else if (cur_P > cwd) {
            GW_req = 0
            temp_SM = temp_SM 
            GW1_used[i,j1] = 0
            GW1_change[i] = 0
            
          }
        }
        #}
        #Remove the interflow volume right away (if applicable)
        if (temp_SM>cur_pars$RD*Soil_sat){
          #Purely to ensure that the change in LU doesn't lead to interflow when there is no rainfall due to changing root depth
          if (cur_P == 0 & (CN1_cal[doy,j1] != CN1_cal[doy-1,j1])) {
            mod_loss_temp = (temp_SM - cur_pars$RD*Soil_sat)
            mod_loss = mod_loss + (temp_SM - cur_pars$RD*Soil_sat)
            IF1[i,j1] = 0 
            temp_SM = temp_SM - IF1[i,j1] - mod_loss_temp
          } else {
            IF1[i,j1] = temp_SM - cur_pars$RD*Soil_sat
            temp_SM = temp_SM - IF1[i,j1]
          }
        } else {
          IF1[i,j1] = 0
        }
        #Basically call the function Percolation that uses one of the methods described in Raven
        DP1[i,j1] = Percolation(temp_SM,max_percolation,cur_pars)
        
        #Estimate the evapotranspiration by first estimating the Water stress factor Ks
        TAW1 = (cur_pars$RD * Soil_FC)-(cur_pars$RD * Soil_WP)
        RAW1 = rho * TAW1
        #Adj1usted soil moisture is basically the total SM - WP
        AdSM = temp_SM - (cur_pars$RD * Soil_WP)
        if (AdSM < (TAW1-RAW1)){
          #Keep an eye on whether to use SM1[i-1,j] or temp_SM
          Ks = AdSM / (TAW1-RAW1)
          ET1[i,j1] = Ks * cur_kc * PET[cur_month]
          if (ET1[i,j1] < 0) {ET1[i,j1] = 0}
        } else {
          ET1[i,j1] = 1*cur_kc*PET[cur_month]
        }
        #The way it is setup right now the SM will never end the day at saturation point
        SM1[i,j1] = temp_SM - DP1[i,j1] - ET1[i,j1]
      } else {
        SM1[i,j1] = (Soil_WP*cur_pars$RD) + cur_P - RF1[i,j1]
        if (SM1[i,j1]>cur_pars$RD*Soil_sat){
          DP1[i,j1] = max_percolation
          IF1[i,j1] = SM1[i,j1] - cur_pars$RD * Soil_sat - max_percolation
          ET1[i,j1] = 1 * cur_kc * PET[cur_month]
          SM1[i,j1] = SM1[i,j1] - DP1[i,j1] - IF1[i,j1] - ET1[i,j1]
        } else if (SM1[i,j1] < cur_pars$RD * Soil_sat){
          DP1[i,j1] = Percolation(SM1[i,j1],max_percolation,cur_pars)
          IF1[i,j1] = 0
          #Need to add a value to scale the AET based on the soil moisture level 
          #Ignore the values generated at i=1
          ET1[i,j1] = 1 * cur_kc * PET[cur_month]
          SM1[i,j1] = SM1[i,j1] - DP1[i,j1] - IF1[i,j1] - ET1[i,j1]
          if (SM1[i,j1] <= 0){SM1[i,j1] = 0}
        } else {
          (print('Error due to Deep percolation calculation'))
        } 
      } 
    } else if (cur_LU=='Rice'){
      temp_SM = SM1[i-1,j1] + cur_P
      Rice_SM_sat = 0.415*cur_pars$RD
      Rice_run_thresh = Rice_SM_sat + pad_BH
      #GW irrigation input
      #if (cur_LU != 'Fallow') {
      if (GW_irrigation == 'Y' && cur_LU != 'Fallow'){
        #Estimate the difference between target and current SM levels (in mm)
        #then multiply that by the total Area of that HRU to estimate the volume needed
        GW_available = GW_yield(AQ1[i-1], AQ1_max, well_max, wl_thresh)*LU1_details$Pumps[j1]*electricity
        cwd = Crop_WD(cur_LU,cur_pars,cur_kc)
        if (cur_P < cwd) {
          GW_req = (cwd - cur_P)*(1/1000)*(HRU_Areas1[j1])
          if (GW_req > GW_available){
            temp_SM = temp_SM + (GW_available / HRU_Areas1[j1])*1000
            GW1_used[i,j1] = (GW_available / t1_const$max_catch)*1000
            GW1_change[i] = (GW_available / HRU_Areas1[j1])*1000
          } else if (GW_req < GW_available) {
            temp_SM = temp_SM + (GW_req / HRU_Areas1[j1])*1000
            GW1_used[i,j1] = (GW_req / t1_const$max_catch)*1000
            GW1_change[i] = (GW_req / HRU_Areas1[j1])*1000
          } 
        } else if (cur_P > cwd) {
          GW_req = 0
          temp_SM = temp_SM 
          GW1_used[i,j1] = 0
          GW1_change[i] = 0
          
        }
      }
      #}
      #Calculate the runoff from the rice fields here. Runoff only starts when the total SM value on the rice fields is greater than 
      #the bund height + soil saturation water content
      if (temp_SM > Rice_run_thresh){
        RF1[i,j1] = temp_SM - Rice_run_thresh
        RF1_vol[i,j1] = (LU1_details$Per_Area[j1] * t1_const$max_catch) * RF1[i,j1] * (1/1000)
        temp_SM = temp_SM - RF1[i,j1]
      } else if (temp_SM < Rice_run_thresh){
        RF1[i,j1] = 0
        RF1_vol[i,j1] = (LU1_details$Per_Area[j1] * t1_const$max_catch) * RF1[i,j1] * (1/1000)
      } else {print('Error due to Rice runoff calculation')}
      #Start calculating deep percolation from the rice fields here using the method in Gowing paper
      DP1[i,j1] = Percolation_Rice(temp_SM,Rice_max_percolation,cur_pars)
      #Estimate the evapotranspiration by first estimating the Water stress factor Ks
      TAW1 = (cur_pars$RD * Soil_FC)-(cur_pars$RD * Soil_WP)
      RAW1 = rho * TAW1
      #Adj1usted soil moisture is basically the total SM - WP
      if (temp_SM > Rice_SM_sat){
        AdSM = Rice_SM_sat -(cur_pars$RD * Soil_WP)
      } else if (temp_SM < Rice_SM_sat) {
        AdSM = temp_SM -(cur_pars$RD * Soil_WP)
      } else {print('Error due to Soil Moisture--ET calculation')}
      if (AdSM < (TAW1-RAW1)){
        #Keep an eye on whether to use SM1[i-1,j] or temp_SM
        Ks = AdSM / (TAW1-RAW1)
        ET1[i,j1] = Ks * cur_kc * PET[cur_month]
        if (ET1[i,j1] < 0) {ET1[i,j1] = 0}
      } else {
        ET1[i,j1] = 1*cur_kc*PET[cur_month]
      }
      #Interflow will always be zero
      IF1[i,j1] = 0
      #The way it is setup right now the SM will never end the day at saturation point
      SM1[i,j1] = temp_SM - DP1[i,j1] - ET1[i,j1]
    }
  }
  ############################################################  ############################################################
  ############################################################  ############################################################
  
  #Tank1 Start
  #
  #Get the tank area from the previous timestep
  if (i == 1) {
    t1_area[1] = 0
    t1_vol[1] = 0
    t1_temp_area = 0
    t1_temp_vol = 0
  } else {
    t1_temp_area = t1_area[i-1]
    t1_temp_vol = t1_vol[i-1]
  } 
  #Estimate the area for each of the HRU's (special HRU is Fallow which converts to tank)
  HRU_Areas1_temp = HRU_Areas1
  HRU_Areas1_temp[1] = HRU_Areas1_temp[1] - t1_temp_area
  t1_inflow_temp = (IF1[i,] + RF1[i,]) * HRU_Areas1_temp * (1/1000)
  t1_inflow[i] = sum(t1_inflow_temp)
  
  
  # Use the tank area to estimate the percolation occuring from different HRU's in the catchment area.
  if (i > 1) {
    #Basically calculate the volumetric input into the aquifer bucket based on the size of the HRU's
    #then divide it by the total catchment of the tank (including tank area) to estimate the change in mm
    AQ1_inputs = (DP1[i,] * HRU_Areas1_temp) /(t1_const$max_catch)
    AQ1[i] = AQ1[i-1] + sum(AQ1_inputs)
    LF1[i] = lateralflow (AQ1[i], AQ1_max, LF1_thresh, LF1_max, LF_coeff)
    BF1[i] = baseflow (AQ1[i], AQ1_max, BF1_thresh, BF1_max, BF_coeff)
    GW1_Irr = sum(GW1_used[i,])
    AQ1[i] = AQ1[i] - BF1[i] - LF1[i] - GW1_Irr
  } else {
    AQ1[i] = AQ1_ini + sum(DP1[i,])
    BF1[i] = 0 
    LF1[i] = 0
  }
  
  #Estimate the water added to the tank by direct precipitation
  t1_precip[i] = cur_P * t1_temp_area *(1/1000)
  
  
  
  #Update the temp_tank volume to include the inputs calculated above
  t1_temp_vol = t1_temp_vol + t1_precip[i] + t1_inflow[i]
  
  #Update the t1_temp area and then subsequently the fallow area
  #Stage-Volume relationship from Mike data
  t1_temp_stage = (t1_temp_vol/22914)^(1/1.9461)
  #Stage Area 
  t1_temp_area=42942*(t1_temp_stage)^1.0993
  #HRU_Areas update
  HRU_Areas1_temp[1] = HRU_Areas1[1] - t1_temp_area
  
  #ET from tank
  t1_ET[i] = t1_temp_area * PET[cur_month] * (1/1000) #m3/day
  
  ####Compare the tank capacity and current volume of water in tank.
  vol_diff1 = t1_temp_vol - t1_const$max_volume
  if (vol_diff1 >= 0){
    t1_spill[i] = vol_diff1
  } else{ t1_spill[i] = 0 }
  
  ####Sluice Outflow
  Qo1a = (( t1_temp_stage - 0.785 ) * 5.1903 ) * 86.4 #86.4 Converst L/s to m3/d
  Qo1b = ((( t1_temp_stage - 1.185 ) * 9.6768 ) + (( t1_temp_stage - 1.185 ) * 4.9196 )) * 86.4 #86.4 Converst L/s to m3/d
  
  if (Qo1a < 0){Qo1a = 0}
  if (Qo1b < 0) {Qo1b = 0}
  
  if ( t1_temp_stage < 0.785 ){
    t1_sluice[i] = 0 } else if (0.785 < t1_temp_stage & t1_temp_stage < 1.185) {
      t1_sluice[i] = Qo1a} else if( t1_temp_stage > 1.185 ){
        t1_sluice[i] = ( Qo1a + Qo1b )}
  
  ###Spillage from upstream tank
  t1_spill_add = 0
  
  ####Net GW exchange- Mike paper
  GW_loss1 = 8.6 * t1_temp_stage - 6.5 #mm/day
  if (GW_loss1 < 0) {GW_loss1 = 0}
  t1_GW[i] = (GW_loss1/1000) * t1_temp_area #m3/day
  
  ####Baseflow into tank
  t1_BF[i] = (BF1[i]/1000) * t1_temp_area
  
  #Total Storage change
  t1_vol[i] = t1_temp_vol - ( t1_ET[i] + t1_sluice[i] + t1_spill[i] + t1_GW[i] )
  
  #Stage-Volume relationship from Mike data
  t1_stage[i] = (t1_vol[i]/22914) ^ (1/1.9461)
  #Stage Area 
  t1_area[i] = 42942 * (t1_stage[i]) ^ 1.0993
  #
  #cur_all=cbind(t1_stage[j],t1_area[j],t1_vol[j],inflow1[j],coredata(inflow_vol),t1_ET[j],t1_sluice[j],t1_spill[j],t1_GW[j])
  #t1_all=rbind(t1_all,cur_all)
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  j2=5
  for (j2 in 1:ncol(CN2_cal)){
    cur_CN = CN2_cal[doy,j2]
    cur_pars = filter(LU_Parameters,CN==cur_CN)
    cur_LU = cur_pars$LU
    cur_kc = Kc2_cal[doy,j2]
    SW2_irrigation = LU2_details$SW_Irr[j2]
    GW2_irrigation = LU2_details$GW_Irr[j2]
    if(cur_LU != 'Rice'){
      #Add the sluice from the upstream tank here for areas that are irrigated by Surface water   
      if (SW2_irrigation == 'Y') {
        #Split the sluice calculated above based on the percent command area of this HRU
        HRU_per_area = LU2_details$Per_Area[j2]/100
        #Convert the sluice volume to depth by dividing by the area of the command HRU
        sluice_depth = t1_sluice[i]*1000*HRU_per_area / HRU_Areas2 [j2]
        rain_sluice = cur_P + sluice_depth
      } else {rain_sluice = cur_P}
      
      RF2[i,j2] = SCS_curve(cur_CN, rain_sluice, rain_5)
      RF2_vol[i,j2] = (LU2_details$Per_Area[j2] * t2_const$max_catch) * RF2[i,j2] * (1/1000)
      
      if (i > 1) {
        #Keep an eye on if I should be calculating Percolation and ET based on the modified SM vs exact SM from previous timestep  
        temp_SM = SM2[i-1,j2] + rain_sluice - RF2[i,j2] 
        #GW_irrigation
        #if (cur_LU != 'Fallow') {
        if (GW2_irrigation == 'Y' && cur_LU != 'Fallow'){
          #Estimate the difference between target and current SM levels (in mm)
          #then multiply that by the total Area of that HRU to estimate the volume needed
          GW_available = GW_yield(AQ2[i-1], AQ2_max, well_max, wl_thresh)*LU2_details$Pumps[j2]*electricity
          cwd = Crop_WD(cur_LU,cur_pars,cur_kc)
          if (rain_sluice < cwd) {
            GW_req = (cwd - rain_sluice )*(1/1000)*(HRU_Areas2[j2])
            if (GW_req > GW_available){
              temp_SM = temp_SM + (GW_available / HRU_Areas2[j2])*1000
              GW2_used[i,j2] = (GW_available / t2_const$max_catch)*1000
              GW2_change[i] = (GW_available / HRU_Areas2[j2])*1000
            } else if (GW_req < GW_available) {
              temp_SM = temp_SM + (GW_req / HRU_Areas2[j2])*1000
              GW2_used[i,j2] = (GW_req / t2_const$max_catch)*1000
              GW2_change[i] = (GW_req / HRU_Areas2[j2])*1000
            } 
          } else if (rain_sluice  > cwd) {
            GW_req = 0
            temp_SM = temp_SM 
            GW2_used[i,j2] = 0
            GW2_change[i] = 0
          }
        } 
        #}
        #Remove the interflow volume right away (if applicable)
        if (temp_SM>cur_pars$RD*Soil_sat){
          #Purely to ensure that the change in LU doesn't lead to interflow when there is no rainfall due to changing root depth
          if (cur_P == 0 & (CN2_cal[doy,j2] != CN2_cal[doy-1,j2])) {
            mod_loss_temp = (temp_SM - cur_pars$RD*Soil_sat)
            mod_loss = mod_loss + (temp_SM - cur_pars$RD*Soil_sat)
            IF2[i,j2] = 0 
            temp_SM = temp_SM - IF2[i,j2] - mod_loss_temp
          } else {
            IF2[i,j2] = temp_SM - cur_pars$RD*Soil_sat
            temp_SM = temp_SM - IF2[i,j2]
          }
        } else {
          IF2[i,j2] = 0
        }
        #Basically call the function Percolation that uses one of the methods described in Raven
        DP2[i,j2] = Percolation(temp_SM,max_percolation,cur_pars)
        
        #Estimate the evapotranspiration by first estimating the Water stress factor Ks
        TAW2 = (cur_pars$RD * Soil_FC)-(cur_pars$RD * Soil_WP)
        RAW2 = rho * TAW2
        #Adj2usted soil moisture is basically the total SM - WP
        AdSM = temp_SM - (cur_pars$RD * Soil_WP)
        if (AdSM < (TAW2-RAW2)){
          #Keep an eye on whether to use SM1[i-1,j] or temp_SM
          Ks = AdSM / (TAW2-RAW2)
          ET2[i,j2] = Ks * cur_kc * PET[cur_month]
          if (ET2[i,j2] < 0) {ET2[i,j2] = 0}
        } else {
          ET2[i,j2] = 1*cur_kc*PET[cur_month]
        }
        #The way it is setup right now the SM will never end the day at saturation point
        SM2[i,j2] = temp_SM - DP2[i,j2] - ET2[i,j2]
      } else {
        SM2[i,j2] = (Soil_WP*cur_pars$RD)+ rain_sluice - RF2[i,j2]
        if (SM2[i,j2]>cur_pars$RD*Soil_sat){
          DP2[i,j2] = max_percolation
          IF2[i,j2] = SM2[i,j2] - cur_pars$RD * Soil_sat - max_percolation
          ET2[i,j2] = 1 * cur_kc * PET[cur_month]
          SM2[i,j2] = SM2[i,j2] - DP2[i,j2] - IF2[i,j2] - ET2[i,j2]
        } else if (SM2[i,j2] < cur_pars$RD * Soil_sat){
          DP2[i,j2] = Percolation(SM2[i,j2],max_percolation,cur_pars)
          IF2[i,j2] = 0
          #Need to add a value to scale the AET based on the soil moisture level 
          #Ignore the values generated at i=1
          ET2[i,j2] = 1 * cur_kc * PET[cur_month]
          SM2[i,j2] = SM2[i,j2] - DP2[i,j2] - IF2[i,j2] - ET2[i,j2]
          if (SM2[i,j2] <= 0){SM2[i,j2] = 0}
        } else {
          (print('Error due to Deep percolation calculation'))
        } 
      } 
    } else if (cur_LU=='Rice'){
      
      #Add the sluice from the upstream tank here for areas that are irrigated by Surface water   
      if (SW2_irrigation == 'Y') {
        #Split the sluice calculated above based on the percent command area of this HRU
        HRU_per_area = LU2_details$Per_Area[j2]/100
        #Convert the sluice volume to depth by dividing by the area of the command HRU; convert meter to mm
        sluice_depth = t1_sluice[i]*1000*HRU_per_area / HRU_Areas2 [j2]
        rain_sluice = cur_P + sluice_depth
      } else {rain_sluice = cur_P}
      
      temp_SM = SM2[i-1,j2] + rain_sluice
      Rice_SM_sat = 0.415*cur_pars$RD
      Rice_run_thresh = Rice_SM_sat + pad_BH
      #GW_irrigation
      #if (cur_LU != 'Fallow') {
      if (GW2_irrigation == 'Y' && cur_LU != 'Fallow'){
        #Estimate the difference between target and current SM levels (in mm)
        #then multiply that by the total Area of that HRU to estimate the volume needed
        GW_available = GW_yield(AQ2[i-1], AQ2_max, well_max, wl_thresh)*LU2_details$Pumps[j2]*electricity
        cwd = Crop_WD(cur_LU,cur_pars,cur_kc)
        if (rain_sluice < cwd) {
          GW_req = (cwd - rain_sluice )*(1/1000)*(HRU_Areas2[j2])
          if (GW_req > GW_available){
            temp_SM = temp_SM + (GW_available / HRU_Areas2[j2])*1000
            GW2_used[i,j2] = (GW_available / t2_const$max_catch)*1000
            GW2_change[i] = (GW_available / HRU_Areas2[j2])*1000
          } else if (GW_req < GW_available) {
            temp_SM = temp_SM + (GW_req / HRU_Areas2[j2])*1000
            GW2_used[i,j2] = (GW_req / t2_const$max_catch)*1000
            GW2_change[i] = (GW_req / HRU_Areas2[j2])*1000
          } 
        } else if (rain_sluice  > cwd) {
          GW_req = 0
          temp_SM = temp_SM 
          GW2_used[i,j2] = 0
          GW2_change[i] = 0
        }
      } 
      # }
      #Calculate the runoff from the rice fields here. Runoff only starts when the total SM value on the rice fields is greater than 
      #the bund height + soil saturation water content
      if (temp_SM > Rice_run_thresh){
        RF2[i,j2] = temp_SM - Rice_run_thresh
        RF2_vol[i,j2] = (LU2_details$Per_Area[j2] * t2_const$max_catch) * RF2[i,j2] * (1/1000)
        temp_SM = temp_SM - RF2[i,j2]
      } else if (temp_SM < Rice_run_thresh){
        RF2[i,j2] = 0
        RF2_vol[i,j2] = (LU2_details$Per_Area[j2] * t2_const$max_catch) * RF2[i,j2] * (1/1000)
      } else {print('Error due to Rice runoff calculation')}
      #Start calculating deep percolation from the rice fields here using the method in Gowing paper
      DP2[i,j2] = Percolation_Rice(temp_SM,Rice_max_percolation,cur_pars)
      #Estimate the evapotranspiration by first estimating the Water stress factor Ks
      TAW2 = (cur_pars$RD * Soil_FC)-(cur_pars$RD * Soil_WP)
      RAW2 = rho * TAW2
      #Adj2usted soil moisture is basically the total SM - WP
      if (temp_SM > Rice_SM_sat){
        AdSM = Rice_SM_sat -(cur_pars$RD * Soil_WP)
      } else if (temp_SM < Rice_SM_sat) {
        AdSM = temp_SM -(cur_pars$RD * Soil_WP)
      } else {print('Error due to Soil Moisture--ET calculation')}
      if (AdSM < (TAW2-RAW2)){
        #Keep an eye on whether to use SM1[i-1,j] or temp_SM
        Ks = AdSM / (TAW2-RAW2)
        ET2[i,j2] = Ks * cur_kc * PET[cur_month]
        if (ET2[i,j2] < 0) {ET2[i,j2] = 0}
      } else {
        ET2[i,j2] = 1*cur_kc*PET[cur_month]
      }
      #Interflow will always be zero
      IF2[i,j2] = 0
      #The way it is setup right now the SM will never end the day at saturation point
      SM2[i,j2] = temp_SM - DP2[i,j2] - ET2[i,j2]
    }
  }
  
  HRU_Areas2_temp = HRU_Areas2
  if (i > 1) {
    #Basically calculate the volumetric input into the aquifer bucket based on the size of the HRU's
    #then divide it by the total catchment of the tank (including tank area) to estimate the change in mm
    AQ2_inputs = (DP2[i,] * HRU_Areas2_temp) /(t2_const$max_catch)
    AQ2[i] = AQ2[i-1] + sum(AQ2_inputs) + 1000*(t1_GW[i]/(t2_const$max_catch))
    LF2[i] = lateralflow (AQ2[i], AQ2_max, LF2_thresh, LF2_max, LF_coeff)
    BF2[i] = baseflow (AQ2[i], AQ2_max, BF2_thresh, BF2_max, BF_coeff)
    #Assume that 30% of the tank spillage recharges the aquifer (conveyance loss)
    Spill2[i] = 0.3*t1_spill[i] / (t2_const$max_catch)
    GW2_Irr = sum(GW2_used[i,])
    #GW1_Irr = sum(GW1_used[i,])
    AQ2[i] = AQ2[i] - BF2[i] - LF2[i] - GW2_Irr + Spill2[i]
  } else {
    AQ2[i] = AQ1_ini + sum(DP2[i,])
    BF2[i] = 0 
    LF2[i] = 0
  }
  
  
  PET_max[i] = Kc2_cal[doy,5]*PET[cur_month]
  TP_max[i] = Kc2_cal[doy,5]
  
  i=i+1
}