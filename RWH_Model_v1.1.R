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
        if (cur_LU != 'Juliflora'){
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
        } else if (cur_LU == 'Juliflora'){
          #Use the Vervoort and Van Der Zee method here; adapted for aquifer bucket instead of actual water table
          #WT_Z is the pseudo-water table
            WT_Z[i] = (AQ1[i-1])/(porosity*1000)+(cur_pars$RD/1000)#(AQ1_max- AQ1[i-1])/porosity
          #These parameters are the same as the paper
            Rz_1 = (100/c)*(exp(-cur_pars$RD/((100/c)*1000)) - exp(-WT_Z[i]/(100/c)))
            Rz_2 = 50
            Rc_Z = (WT_Z[i])/(WT_Z[i]-(cur_pars$RD/1000))*(Rz_1/Rz_2)
            fr = (100/c)*(1 - exp(-cur_pars$RD/((100/c)*1000)))
            #test_SM = Soil_sat*cur_pars$RD
            cur_beta = beta_func (temp_SM, AW_critical, Soil_WP, Soil_sat, cur_pars$RD)
            ET_soil = fr*cur_beta*(1-exp(-light_par*LAI))*PET_J
            temp_min = min((1-fr*cur_beta),Rc_Z)
            ET_GW[i] = temp_min*(1-exp(-light_par*LAI))*PET_J
            ET1[i,j1] = ET_soil + ET_GW[i]
        }
        #The way it is setup right now the SM will never end the day at saturation point
        SM1[i,j1] = temp_SM - DP1[i,j1] - ET_soil#ET1[i,j1]
      } else {
        ET_GW[i] = 0
        SM1[i,j1] = (Soil_WP*cur_pars$RD) + cur_P - RF1[i,j1]
        if (SM1[i,j1]>cur_pars$RD*Soil_sat){
          DP1[i,j1] = max_percolation
          IF1[i,j1] = SM1[i,j1] - cur_pars$RD * Soil_sat - max_percolation
          ET1[i,j1] = 0 * cur_kc * PET[cur_month]
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
    ET_AQ1[i] = ET_GW[i]*(HRU_Areas1_temp[5])/(t1_const$max_catch)
    AQ1[i] = AQ1[i] - BF1[i] - LF1[i] - GW1_Irr - ET_AQ1[i]
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
  
  ############################################################  ############################################################
  ############################################################  ############################################################
  
  #Tank2 Start
  #
  #Get the tank area from the previous timestep
  if (i == 1) {
    t2_area[1] = 0
    t2_vol[1] = 0
    t2_temp_area = 0
    t2_temp_vol = 0
  } else {
    t2_temp_area = t2_area[i-1]
    t2_temp_vol = t2_vol[i-1]
  } 
  #Estimate the area for each of the HRU's (special HRU is Fallow which converts to tank)
  HRU_Areas2_temp = HRU_Areas2
  HRU_Areas2_temp[1] = HRU_Areas2_temp[1] - t2_temp_area
  t2_inflow_temp = (IF2[i,] + RF2[i,]) * HRU_Areas2_temp * (1/1000)
  t2_inflow[i] = sum(t2_inflow_temp)
  
  
  #Estimate the water added to the tank by direct precipitation
  t2_precip[i] = cur_P * t2_temp_area *(1/1000)
  
  #Update the temp_tank volume to include the inputs calculated above
  t2_temp_vol = t2_temp_vol + t2_precip[i] + t2_inflow[i]
  
  #Update the t2_temp area and then subsequently the fallow area
  #Stage-Volume relationship from Mike data
  t2_temp_stage = (t2_temp_vol/4852.2)^(1/3.60)
  #Stage Area 
  t2_temp_area=14674*(t2_temp_stage)^2.88
  #HRU_Areas update
  HRU_Areas2_temp[1] = HRU_Areas2[1] - t2_temp_area
  
  #ET from tank
  t2_ET[i] = t2_temp_area * PET[cur_month] * (1/1000) #m3/day
  
  ####Compare the tank capacity and current volume of water in tank.
  vol_diff2 = t2_temp_vol - t2_const$max_volume
  if (vol_diff2 >= 0){
    t2_spill[i] = vol_diff2
  } else{ t2_spill[i] = 0 }
  
  ####Sluice Outflow
  Qo2a = (( t2_temp_stage - 0.785 ) * 16.98 ) * 86.4 #86.4 Converst L/s to m3/d
  Qo2b = ((( t2_temp_stage - 1.185 ) * 2.35 ) + (( t2_temp_stage - 1.185 ) * 4.9196 )) * 86.4 #86.4 Converst L/s to m3/d
  
  if (Qo2a < 0){Qo2a = 0}
  if (Qo2b < 0) {Qo2b = 0}
  
  if ( t2_temp_stage < 0.785 ){
    t2_sluice[i] = 0 } else if (0.785 < t2_temp_stage & t2_temp_stage < 1.185) {
      t2_sluice[i] = Qo2a} else if( t2_temp_stage > 1.185 ){
        t2_sluice[i] = ( Qo2a + Qo2b )}
  
  ###Spillage from upstream tank
  t2_spill_add = 0
  
  ####Net GW exchange- Mike paper
  GW_loss2 = 8.6 * t2_temp_stage - 6.9 #mm/day; (Van Meter et.al., 2016)
  if (GW_loss2 < 0) {GW_loss2 = 0}
  t2_GW[i] = (GW_loss2/1000) * t2_temp_area #m3/day
  
  ####Baseflow into tank
  t2_BF[i] = (BF2[i]/1000) * t2_temp_area
  
  #Total Storage change
  t2_vol[i] = t2_temp_vol - ( t2_ET[i] + t2_sluice[i] + t2_spill[i] + t2_GW[i] )
  
  #Stage-Volume relationship from Mike data
  t2_stage[i] = (t2_vol[i]/4852.2)^(1/3.60)
  #Stage Area 
  t2_area[i] = 14674 * (t2_stage[i]) ^ 2.88
  #
  #cur_all=cbind(t1_stage[j],t1_area[j],t1_vol[j],inflow1[j],coredata(inflow_vol),t1_ET[j],t1_sluice[j],t1_spill[j],t1_GW[j])
  #t1_all=rbind(t1_all,cur_all)
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  j3=5
  for (j3 in 1:ncol(CN3_cal)){
    cur_CN = CN3_cal[doy,j3]
    cur_pars = filter(LU_Parameters,CN==cur_CN)
    cur_LU = cur_pars$LU
    cur_kc = Kc3_cal[doy,j3]
    SW3_irrigation = LU3_details$SW_Irr[j3]
    GW3_irrigation = LU3_details$GW_Irr[j3]
    if(cur_LU != 'Rice'){
      #Add the sluice from the upstream tank here for areas that are irrigated by Surface water   
      if (SW3_irrigation == 'Y') {
        #Split the sluice calculated above based on the percent command area of this HRU
        HRU_per_area = LU3_details$Per_Area[j3]/100
        #Convert the sluice volume to depth by dividing by the area of the command HRU
        sluice_depth = t2_sluice[i]*1000*HRU_per_area / HRU_Areas3 [j3]
        rain_sluice = cur_P + sluice_depth
      } else {rain_sluice = cur_P}
      
      RF3[i,j3] = SCS_curve(cur_CN, rain_sluice, rain_5)
      RF3_vol[i,j3] = (LU3_details$Per_Area[j3] * t3_const$max_catch) * RF3[i,j3] * (1/1000)
      
      if (i > 1) {
        #Keep an eye on if I should be calculating Percolation and ET based on the modified SM vs exact SM from previous timestep  
        temp_SM = SM3[i-1,j3] + rain_sluice - RF3[i,j3] 
        #GW_irrigation
        #if (cur_LU != 'Fallow') {
        if (GW3_irrigation == 'Y' && cur_LU != 'Fallow'){
          #Estimate the difference between target and current SM levels (in mm)
          #then multiply that by the total Area of that HRU to estimate the volume needed
          GW_available = GW_yield(AQ3[i-1], AQ3_max, well_max, wl_thresh)*LU3_details$Pumps[j3]*electricity
          cwd = Crop_WD(cur_LU,cur_pars,cur_kc)
          if (rain_sluice < cwd) {
            GW_req = (cwd - rain_sluice )*(1/1000)*(HRU_Areas3[j3])
            if (GW_req > GW_available){
              temp_SM = temp_SM + (GW_available / HRU_Areas3[j3])*1000
              GW3_used[i,j3] = (GW_available / t3_const$max_catch)*1000
              GW3_change[i] = (GW_available / HRU_Areas3[j3])*1000
            } else if (GW_req < GW_available) {
              temp_SM = temp_SM + (GW_req / HRU_Areas3[j3])*1000
              GW3_used[i,j3] = (GW_req / t3_const$max_catch)*1000
              GW3_change[i] = (GW_req / HRU_Areas3[j3])*1000
            } 
          } else if (rain_sluice  > cwd) {
            GW_req = 0
            temp_SM = temp_SM 
            GW3_used[i,j3] = 0
            GW3_change[i] = 0
          }
        } 
        #}
        #Remove the interflow volume right away (if applicable)
        if (temp_SM>cur_pars$RD*Soil_sat){
          #Purely to ensure that the change in LU doesn't lead to interflow when there is no rainfall due to changing root depth
          if (cur_P == 0 & (CN3_cal[doy,j3] != CN3_cal[doy-1,j3])) {
            mod_loss_temp = (temp_SM - cur_pars$RD*Soil_sat)
            mod_loss = mod_loss + (temp_SM - cur_pars$RD*Soil_sat)
            IF3[i,j3] = 0 
            temp_SM = temp_SM - IF3[i,j3] - mod_loss_temp
          } else {
            IF3[i,j3] = temp_SM - cur_pars$RD*Soil_sat
            temp_SM = temp_SM - IF3[i,j3]
          }
        } else {
          IF3[i,j3] = 0
        }
        #Basically call the function Percolation that uses one of the methods described in Raven
        DP3[i,j3] = Percolation(temp_SM,max_percolation,cur_pars)
        
        #Estimate the evapotranspiration by first estimating the Water stress factor Ks
        TAW3 = (cur_pars$RD * Soil_FC)-(cur_pars$RD * Soil_WP)
        RAW3 = rho * TAW3
        #Adj3usted soil moisture is basically the total SM - WP
        AdSM = temp_SM - (cur_pars$RD * Soil_WP)
        if (AdSM < (TAW3-RAW3)){
          #Keep an eye on whether to use SM1[i-1,j] or temp_SM
          Ks = AdSM / (TAW3-RAW3)
          ET3[i,j3] = Ks * cur_kc * PET[cur_month]
          if (ET3[i,j3] < 0) {ET3[i,j3] = 0}
        } else {
          ET3[i,j3] = 1*cur_kc*PET[cur_month]
        }
        #The way it is setup right now the SM will never end the day at saturation point
        SM3[i,j3] = temp_SM - DP3[i,j3] - ET3[i,j3]
      } else {
        SM3[i,j3] = (Soil_WP*cur_pars$RD)+ rain_sluice - RF3[i,j3]
        if (SM3[i,j3]>cur_pars$RD*Soil_sat){
          DP3[i,j3] = max_percolation
          IF3[i,j3] = SM3[i,j3] - cur_pars$RD * Soil_sat - max_percolation
          ET3[i,j3] = 1 * cur_kc * PET[cur_month]
          SM3[i,j3] = SM3[i,j3] - DP3[i,j3] - IF3[i,j3] - ET3[i,j3]
        } else if (SM3[i,j3] < cur_pars$RD * Soil_sat){
          DP3[i,j3] = Percolation(SM3[i,j3],max_percolation,cur_pars)
          IF3[i,j3] = 0
          #Need to add a value to scale the AET based on the soil moisture level 
          #Ignore the values generated at i=1
          ET3[i,j3] = 1 * cur_kc * PET[cur_month]
          SM3[i,j3] = SM3[i,j3] - DP3[i,j3] - IF3[i,j3] - ET3[i,j3]
          if (SM3[i,j3] <= 0){SM3[i,j3] = 0}
        } else {
          (print('Error due to Deep percolation calculation'))
        } 
      } 
    } else if (cur_LU=='Rice'){
      
      #Add the sluice from the upstream tank here for areas that are irrigated by Surface water   
      if (SW3_irrigation == 'Y') {
        #Split the sluice calculated above based on the percent command area of this HRU
        HRU_per_area = LU3_details$Per_Area[j3]/100
        #Convert the sluice volume to depth by dividing by the area of the command HRU; convert meter to mm
        sluice_depth = t2_sluice[i]*1000*HRU_per_area / HRU_Areas3 [j3]
        rain_sluice = cur_P + sluice_depth
      } else {rain_sluice = cur_P}
      
      temp_SM = SM3[i-1,j3] + rain_sluice
      Rice_SM_sat = 0.415*cur_pars$RD
      Rice_run_thresh = Rice_SM_sat + pad_BH
      #GW_irrigation
      #if (cur_LU != 'Fallow') {
      if (GW3_irrigation == 'Y' && cur_LU != 'Fallow'){
        #Estimate the difference between target and current SM levels (in mm)
        #then multiply that by the total Area of that HRU to estimate the volume needed
        GW_available = GW_yield(AQ3[i-1], AQ3_max, well_max, wl_thresh)*LU3_details$Pumps[j3]*electricity
        cwd = Crop_WD(cur_LU,cur_pars,cur_kc)
        if (rain_sluice < cwd) {
          GW_req = (cwd - rain_sluice )*(1/1000)*(HRU_Areas3[j3])
          if (GW_req > GW_available){
            temp_SM = temp_SM + (GW_available / HRU_Areas3[j3])*1000
            GW3_used[i,j3] = (GW_available / t3_const$max_catch)*1000
            GW3_change[i] = (GW_available / HRU_Areas3[j3])*1000
          } else if (GW_req < GW_available) {
            temp_SM = temp_SM + (GW_req / HRU_Areas3[j3])*1000
            GW3_used[i,j3] = (GW_req / t3_const$max_catch)*1000
            GW3_change[i] = (GW_req / HRU_Areas3[j3])*1000
          } 
        } else if (rain_sluice  > cwd) {
          GW_req = 0
          temp_SM = temp_SM 
          GW3_used[i,j3] = 0
          GW3_change[i] = 0
        }
      } 
      # }
      #Calculate the runoff from the rice fields here. Runoff only starts when the total SM value on the rice fields is greater than 
      #the bund height + soil saturation water content
      if (temp_SM > Rice_run_thresh){
        RF3[i,j3] = temp_SM - Rice_run_thresh
        RF3_vol[i,j3] = (LU3_details$Per_Area[j3] * t3_const$max_catch) * RF3[i,j3] * (1/1000)
        temp_SM = temp_SM - RF3[i,j3]
      } else if (temp_SM < Rice_run_thresh){
        RF3[i,j3] = 0
        RF3_vol[i,j3] = (LU3_details$Per_Area[j3] * t3_const$max_catch) * RF3[i,j3] * (1/1000)
      } else {print('Error due to Rice runoff calculation')}
      #Start calculating deep percolation from the rice fields here using the method in Gowing paper
      DP3[i,j3] = Percolation_Rice(temp_SM,Rice_max_percolation,cur_pars)
      #Estimate the evapotranspiration by first estimating the Water stress factor Ks
      TAW3 = (cur_pars$RD * Soil_FC)-(cur_pars$RD * Soil_WP)
      RAW3 = rho * TAW3
      #Adj3usted soil moisture is basically the total SM - WP
      if (temp_SM > Rice_SM_sat){
        AdSM = Rice_SM_sat -(cur_pars$RD * Soil_WP)
      } else if (temp_SM < Rice_SM_sat) {
        AdSM = temp_SM -(cur_pars$RD * Soil_WP)
      } else {print('Error due to Soil Moisture--ET calculation')}
      if (AdSM < (TAW3-RAW3)){
        #Keep an eye on whether to use SM1[i-1,j] or temp_SM
        Ks = AdSM / (TAW3-RAW3)
        ET3[i,j3] = Ks * cur_kc * PET[cur_month]
        if (ET3[i,j3] < 0) {ET3[i,j3] = 0}
      } else {
        ET3[i,j3] = 1*cur_kc*PET[cur_month]
      }
      #Interflow will always be zero
      IF3[i,j3] = 0
      #The way it is setup right now the SM will never end the day at saturation point
      SM3[i,j3] = temp_SM - DP3[i,j3] - ET3[i,j3]
    }
  }
  
  HRU_Areas3_temp = HRU_Areas3
  if (i > 1) {
    #Basically calculate the volumetric input into the aquifer bucket based on the size of the HRU's
    #then divide it by the total catchment of the tank (including tank area) to estimate the change in mm
    AQ3_inputs = (DP3[i,] * HRU_Areas3_temp) /(t3_const$max_catch)
    AQ3[i] = AQ3[i-1] + sum(AQ3_inputs) + 1000*(t1_GW[i]/(t3_const$max_catch))
    LF3[i] = lateralflow (AQ3[i], AQ3_max, LF3_thresh, LF3_max, LF_coeff)
    BF3[i] = baseflow (AQ3[i], AQ3_max, BF3_thresh, BF3_max, BF_coeff)
    #Assume that 30% of the tank spillage recharges the aquifer (conveyance loss)
    Spill3[i] = 0.3*t1_spill[i] / (t3_const$max_catch)
    GW3_Irr = sum(GW3_used[i,])
    #GW1_Irr = sum(GW1_used[i,])
    AQ3[i] = AQ3[i] - BF3[i] - LF3[i] - GW3_Irr + Spill3[i]
  } else {
    AQ3[i] = AQ1_ini + sum(DP3[i,])
    BF3[i] = 0 
    LF3[i] = 0
  }
  
  ############################################################  ############################################################
  ############################################################  ############################################################
  
  #Tank3 Start
  #
  #Get the tank area from the previous timestep
  if (i == 1) {
    t3_area[1] = 0
    t3_vol[1] = 0
    t3_temp_area = 0
    t3_temp_vol = 0
  } else {
    t3_temp_area = t3_area[i-1]
    t3_temp_vol = t3_vol[i-1]
  } 
  #Estimate the area for each of the HRU's (special HRU is Fallow which converts to tank)
  HRU_Areas3_temp = HRU_Areas3
  HRU_Areas3_temp[1] = HRU_Areas3_temp[1] - t3_temp_area
  t3_inflow_temp = (IF3[i,] + RF3[i,]) * HRU_Areas3_temp * (1/1000)
  t3_inflow[i] = sum(t3_inflow_temp)
  
  
  #Estimate the water added to the tank by direct precipitation
  t3_precip[i] = cur_P * t3_temp_area *(1/1000)
  
  #Update the temp_tank volume to include the inputs calculated above
  t3_temp_vol = t3_temp_vol + t3_precip[i] + t3_inflow[i]
  
  #Update the t3_temp area and then subsequently the fallow area
  #Stage-Volume relationship from Mike data
  t3_temp_stage = (t3_temp_vol/3584.7)^(1/2.9)
  #Stage Area 
  t3_temp_area=11676*(t3_temp_stage)^1.79
  #HRU_Areas update
  HRU_Areas3_temp[1] = HRU_Areas3[1] - t3_temp_area
  
  #ET from tank
  t3_ET[i] = t3_temp_area * PET[cur_month] * (1/1000) #m3/day
  
  ####Compare the tank capacity and current volume of water in tank.
  vol_diff3 = t3_temp_vol - t3_const$max_volume
  if (vol_diff3 >= 0){
    t3_spill[i] = vol_diff3
  } else{ t3_spill[i] = 0 }
  
  ####Sluice Outflow
  #Does the L/s to m3/day actually necessary? Based on thesis appendix diagram seems necessary
  Qo3a = (( t3_temp_stage - 0.785 ) * 6.49 ) * 86.4 #86.4 Converst L/s to m3/d
  Qo3b = 0#((( t3_temp_stage - 1.185 ) * 2.35 ) + (( t3_temp_stage - 1.185 ) * 4.9196 )) * 86.4 #86.4 Converst L/s to m3/d
  
  if (Qo3a < 0){Qo3a = 0}
  if (Qo3b < 0) {Qo3b = 0}
  
  if ( t3_temp_stage < 0.785 ){
    t3_sluice[i] = 0 
    } else if (0.785 < t3_temp_stage) {
      t3_sluice[i] = Qo3a}
  
  ###Spillage from upstream tank
  t3_spill_add = 0
  
  ####Net GW exchange- Mike paper
  GW_loss3 = 8.6 * t3_temp_stage - 6.9 #mm/day; (Van Meter et.al., 2016)
  if (GW_loss3 < 0) {GW_loss3 = 0}
  t3_GW[i] = (GW_loss3/1000) * t3_temp_area #m3/day
  
  ####Baseflow into tank
  t3_BF[i] = (BF3[i]/1000) * t3_temp_area
  
  #Total Storage change
  t3_vol[i] = t3_temp_vol - ( t3_ET[i] + t3_sluice[i] + t3_spill[i] + t3_GW[i] )
  
  #Stage-Volume relationship from Mike data
  t3_stage[i] = (t3_vol[i]/3584.7)^(1/2.9)
  #Stage Area 
  t3_area[i] = 11676 * (t3_stage[i]) ^ 1.79
  #
  #cur_all=cbind(t1_stage[j],t1_area[j],t1_vol[j],inflow1[j],coredata(inflow_vol),t1_ET[j],t1_sluice[j],t1_spill[j],t1_GW[j])
  #t1_all=rbind(t1_all,cur_all)
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  j4=5
  for (j4 in 1:ncol(CN4_cal)){
    cur_CN = CN4_cal[doy,j4]
    cur_pars = filter(LU_Parameters,CN==cur_CN)
    cur_LU = cur_pars$LU
    cur_kc = Kc4_cal[doy,j4]
    SW4_irrigation = LU4_details$SW_Irr[j4]
    GW4_irrigation = LU4_details$GW_Irr[j4]
    if(cur_LU != 'Rice'){
      #Add the sluice from the upstream tank here for areas that are irrigated by Surface water   
      if (SW4_irrigation == 'Y') {
        #Split the sluice calculated above based on the percent command area of this HRU
        HRU_per_area = LU4_details$Per_Area[j4]/100
        #Convert the sluice volume to depth by dividing by the area of the command HRU
        sluice_depth = t3_sluice[i]*1000*HRU_per_area / HRU_Areas4 [j4]
        rain_sluice = cur_P + sluice_depth
      } else {rain_sluice = cur_P}
      
      RF4[i,j4] = SCS_curve(cur_CN, rain_sluice, rain_5)
      RF4_vol[i,j4] = (LU4_details$Per_Area[j4] * t4_const$max_catch) * RF4[i,j4] * (1/1000)
      
      if (i > 1) {
        #Keep an eye on if I should be calculating Percolation and ET based on the modified SM vs exact SM from previous timestep  
        temp_SM = SM4[i-1,j4] + rain_sluice - RF4[i,j4] 
        #GW_irrigation
        #if (cur_LU != 'Fallow') {
        if (GW4_irrigation == 'Y' && cur_LU != 'Fallow'){
          #Estimate the difference between target and current SM levels (in mm)
          #then multiply that by the total Area of that HRU to estimate the volume needed
          GW_available = GW_yield(AQ4[i-1], AQ4_max, well_max, wl_thresh)*LU4_details$Pumps[j4]*electricity
          cwd = Crop_WD(cur_LU,cur_pars,cur_kc)
          if (rain_sluice < cwd) {
            GW_req = (cwd - rain_sluice )*(1/1000)*(HRU_Areas4[j4])
            if (GW_req > GW_available){
              temp_SM = temp_SM + (GW_available / HRU_Areas4[j4])*1000
              GW4_used[i,j4] = (GW_available / t4_const$max_catch)*1000
              GW4_change[i] = (GW_available / HRU_Areas4[j4])*1000
            } else if (GW_req < GW_available) {
              temp_SM = temp_SM + (GW_req / HRU_Areas4[j4])*1000
              GW4_used[i,j4] = (GW_req / t4_const$max_catch)*1000
              GW4_change[i] = (GW_req / HRU_Areas4[j4])*1000
            } 
          } else if (rain_sluice  > cwd) {
            GW_req = 0
            temp_SM = temp_SM 
            GW4_used[i,j4] = 0
            GW4_change[i] = 0
          }
        } 
        #}
        #Remove the interflow volume right away (if applicable)
        if (temp_SM>cur_pars$RD*Soil_sat){
          #Purely to ensure that the change in LU doesn't lead to interflow when there is no rainfall due to changing root depth
          if (cur_P == 0 & (CN4_cal[doy,j4] != CN4_cal[doy-1,j4])) {
            mod_loss_temp = (temp_SM - cur_pars$RD*Soil_sat)
            mod_loss = mod_loss + (temp_SM - cur_pars$RD*Soil_sat)
            IF4[i,j4] = 0 
            temp_SM = temp_SM - IF4[i,j4] - mod_loss_temp
          } else {
            IF4[i,j4] = temp_SM - cur_pars$RD*Soil_sat
            temp_SM = temp_SM - IF4[i,j4]
          }
        } else {
          IF4[i,j4] = 0
        }
        #Basically call the function Percolation that uses one of the methods described in Raven
        DP4[i,j4] = Percolation(temp_SM,max_percolation,cur_pars)
        
        #Estimate the evapotranspiration by first estimating the Water stress factor Ks
        TAW4 = (cur_pars$RD * Soil_FC)-(cur_pars$RD * Soil_WP)
        RAW4 = rho * TAW4
        #Adj4usted soil moisture is basically the total SM - WP
        AdSM = temp_SM - (cur_pars$RD * Soil_WP)
        if (AdSM < (TAW4-RAW4)){
          #Keep an eye on whether to use SM1[i-1,j] or temp_SM
          Ks = AdSM / (TAW4-RAW4)
          ET4[i,j4] = Ks * cur_kc * PET[cur_month]
          if (ET4[i,j4] < 0) {ET4[i,j4] = 0}
        } else {
          ET4[i,j4] = 1*cur_kc*PET[cur_month]
        }
        #The way it is setup right now the SM will never end the day at saturation point
        SM4[i,j4] = temp_SM - DP4[i,j4] - ET4[i,j4]
      } else {
        SM4[i,j4] = (Soil_WP*cur_pars$RD)+ rain_sluice - RF4[i,j4]
        if (SM4[i,j4]>cur_pars$RD*Soil_sat){
          DP4[i,j4] = max_percolation
          IF4[i,j4] = SM4[i,j4] - cur_pars$RD * Soil_sat - max_percolation
          ET4[i,j4] = 1 * cur_kc * PET[cur_month]
          SM4[i,j4] = SM4[i,j4] - DP4[i,j4] - IF4[i,j4] - ET4[i,j4]
        } else if (SM4[i,j4] < cur_pars$RD * Soil_sat){
          DP4[i,j4] = Percolation(SM4[i,j4],max_percolation,cur_pars)
          IF4[i,j4] = 0
          #Need to add a value to scale the AET based on the soil moisture level 
          #Ignore the values generated at i=1
          ET4[i,j4] = 1 * cur_kc * PET[cur_month]
          SM4[i,j4] = SM4[i,j4] - DP4[i,j4] - IF4[i,j4] - ET4[i,j4]
          if (SM4[i,j4] <= 0){SM4[i,j4] = 0}
        } else {
          (print('Error due to Deep percolation calculation'))
        } 
      } 
    } else if (cur_LU=='Rice'){
      
      #Add the sluice from the upstream tank here for areas that are irrigated by Surface water   
      if (SW4_irrigation == 'Y') {
        #Split the sluice calculated above based on the percent command area of this HRU
        HRU_per_area = LU4_details$Per_Area[j4]/100
        #Convert the sluice volume to depth by dividing by the area of the command HRU; convert meter to mm
        sluice_depth = t3_sluice[i]*1000*HRU_per_area / HRU_Areas4 [j4]
        rain_sluice = cur_P + sluice_depth
      } else {rain_sluice = cur_P}
      
      temp_SM = SM4[i-1,j4] + rain_sluice
      Rice_SM_sat = 0.415*cur_pars$RD
      Rice_run_thresh = Rice_SM_sat + pad_BH
      #GW_irrigation
      #if (cur_LU != 'Fallow') {
      if (GW4_irrigation == 'Y' && cur_LU != 'Fallow'){
        #Estimate the difference between target and current SM levels (in mm)
        #then multiply that by the total Area of that HRU to estimate the volume needed
        GW_available = GW_yield(AQ4[i-1], AQ4_max, well_max, wl_thresh)*LU4_details$Pumps[j4]*electricity
        cwd = Crop_WD(cur_LU,cur_pars,cur_kc)
        if (rain_sluice < cwd) {
          GW_req = (cwd - rain_sluice )*(1/1000)*(HRU_Areas4[j4])
          if (GW_req > GW_available){
            temp_SM = temp_SM + (GW_available / HRU_Areas4[j4])*1000
            GW4_used[i,j4] = (GW_available / t4_const$max_catch)*1000
            GW4_change[i] = (GW_available / HRU_Areas4[j4])*1000
          } else if (GW_req < GW_available) {
            temp_SM = temp_SM + (GW_req / HRU_Areas4[j4])*1000
            GW4_used[i,j4] = (GW_req / t4_const$max_catch)*1000
            GW4_change[i] = (GW_req / HRU_Areas4[j4])*1000
          } 
        } else if (rain_sluice  > cwd) {
          GW_req = 0
          temp_SM = temp_SM 
          GW4_used[i,j4] = 0
          GW4_change[i] = 0
        }
      } 
      # }
      #Calculate the runoff from the rice fields here. Runoff only starts when the total SM value on the rice fields is greater than 
      #the bund height + soil saturation water content
      if (temp_SM > Rice_run_thresh){
        RF4[i,j4] = temp_SM - Rice_run_thresh
        RF4_vol[i,j4] = (LU4_details$Per_Area[j4] * t4_const$max_catch) * RF4[i,j4] * (1/1000)
        temp_SM = temp_SM - RF4[i,j4]
      } else if (temp_SM < Rice_run_thresh){
        RF4[i,j4] = 0
        RF4_vol[i,j4] = (LU4_details$Per_Area[j4] * t4_const$max_catch) * RF4[i,j4] * (1/1000)
      } else {print('Error due to Rice runoff calculation')}
      #Start calculating deep percolation from the rice fields here using the method in Gowing paper
      DP4[i,j4] = Percolation_Rice(temp_SM,Rice_max_percolation,cur_pars)
      #Estimate the evapotranspiration by first estimating the Water stress factor Ks
      TAW4 = (cur_pars$RD * Soil_FC)-(cur_pars$RD * Soil_WP)
      RAW4 = rho * TAW4
      #Adj4usted soil moisture is basically the total SM - WP
      if (temp_SM > Rice_SM_sat){
        AdSM = Rice_SM_sat -(cur_pars$RD * Soil_WP)
      } else if (temp_SM < Rice_SM_sat) {
        AdSM = temp_SM -(cur_pars$RD * Soil_WP)
      } else {print('Error due to Soil Moisture--ET calculation')}
      if (AdSM < (TAW4-RAW4)){
        #Keep an eye on whether to use SM1[i-1,j] or temp_SM
        Ks = AdSM / (TAW4-RAW4)
        ET4[i,j4] = Ks * cur_kc * PET[cur_month]
        if (ET4[i,j4] < 0) {ET4[i,j4] = 0}
      } else {
        ET4[i,j4] = 1*cur_kc*PET[cur_month]
      }
      #Interflow will always be zero
      IF4[i,j4] = 0
      #The way it is setup right now the SM will never end the day at saturation point
      SM4[i,j4] = temp_SM - DP4[i,j4] - ET4[i,j4]
    }
  }
  
  HRU_Areas4_temp = HRU_Areas4
  if (i > 1) {
    #Basically calculate the volumetric input into the aquifer bucket based on the size of the HRU's
    #then divide it by the total catchment of the tank (including tank area) to estimate the change in mm
    AQ4_inputs = (DP4[i,] * HRU_Areas4_temp) /(t4_const$max_catch)
    AQ4[i] = AQ4[i-1] + sum(AQ4_inputs) + 1000*(t1_GW[i]/(t4_const$max_catch))
    LF4[i] = lateralflow (AQ4[i], AQ4_max, LF4_thresh, LF4_max, LF_coeff)
    BF4[i] = baseflow (AQ4[i], AQ4_max, BF4_thresh, BF4_max, BF_coeff)
    #Assume that 30% of the tank spillage recharges the aquifer (conveyance loss)
    Spill4[i] = 0.3*t1_spill[i] / (t4_const$max_catch)
    GW4_Irr = sum(GW4_used[i,])
    #GW1_Irr = sum(GW1_used[i,])
    AQ4[i] = AQ4[i] - BF4[i] - LF4[i] - GW4_Irr + Spill4[i]
  } else {
    AQ4[i] = AQ1_ini + sum(DP4[i,])
    BF4[i] = 0 
    LF4[i] = 0
  }
  
  ############################################################  ############################################################
  ############################################################  ############################################################
  
  #Tank4 Start
  #
  #Get the tank area from the previous timestep
  if (i == 1) {
    t4_area[1] = 0
    t4_vol[1] = 0
    t4_temp_area = 0
    t4_temp_vol = 0
  } else {
    t4_temp_area = t4_area[i-1]
    t4_temp_vol = t4_vol[i-1]
  } 
  #Estimate the area for each of the HRU's (special HRU is Fallow which converts to tank)
  HRU_Areas4_temp = HRU_Areas4
  HRU_Areas4_temp[1] = HRU_Areas4_temp[1] - t4_temp_area
  t4_inflow_temp = (IF4[i,] + RF4[i,]) * HRU_Areas4_temp * (1/1000)
  t4_inflow[i] = sum(t4_inflow_temp)
  
  
  #Estimate the water added to the tank by direct precipitation
  t4_precip[i] = cur_P * t4_temp_area *(1/1000)
  
  #Update the temp_tank volume to include the inputs calculated above
  t4_temp_vol = t4_temp_vol + t4_precip[i] + t4_inflow[i]
  
  #Update the t4_temp area and then subsequently the fallow area
  #Stage-Volume relationship from Mike data
  t4_temp_stage = (t4_temp_vol/1390.1)^(1/3.76)
  #Stage Area 
  t4_temp_area=4227.8*(t4_temp_stage)^3.21
  #HRU_Areas update
  HRU_Areas4_temp[1] = HRU_Areas4[1] - t4_temp_area
  
  #ET from tank
  t4_ET[i] = t4_temp_area * PET[cur_month] * (1/1000) #m4/day
  
  ####Compare the tank capacity and current volume of water in tank.
  vol_diff4 = t4_temp_vol - t4_const$max_volume
  if (vol_diff4 >= 0){
    t4_spill[i] = vol_diff4
  } else{ t4_spill[i] = 0 }
  
  ####Sluice Outflow
  #Does the L/s to m4/day actually necessary? Based on thesis appendix diagram seems necessary
  Qo4a = (( t4_temp_stage - 0.785 ) * 6.49 ) * 86.4 #86.4 Converst L/s to m3/d
  Qo4b = 0#((( t4_temp_stage - 1.185 ) * 2.35 ) + (( t3_temp_stage - 1.185 ) * 4.9196 )) * 86.4 #86.4 Converst L/s to m3/d
  
  if (Qo4a < 0){Qo4a = 0}
  if (Qo4b < 0) {Qo4b = 0}
  
  if ( t4_temp_stage < 0.785 ){
    t4_sluice[i] = 0 
  } else if (0.785 < t4_temp_stage) {
    t4_sluice[i] = Qo4a}
  
  ###Spillage from upstream tank
  t4_spill_add = 0
  
  ####Net GW exchange- Mike paper
  GW_loss4 = 8.6 * t4_temp_stage - 6.9 #mm/day; (Van Meter et.al., 2016)
  if (GW_loss4 < 0) {GW_loss4 = 0}
  t4_GW[i] = (GW_loss4/1000) * t4_temp_area #m3/day
  
  ####Baseflow into tank
  t4_BF[i] = (BF4[i]/1000) * t4_temp_area
  
  #Total Storage change
  t4_vol[i] = t4_temp_vol - ( t4_ET[i] + t4_sluice[i] + t4_spill[i] + t4_GW[i] )
  
  #Stage-Volume relationship from Mike data
  t4_stage[i] = (t4_vol[i]/1390.1)^(1/3.76)
  #Stage Area 
  t4_area[i] = 4227.8 * (t4_stage[i]) ^ 3.21
  #
  #cur_all=cbind(t1_stage[j],t1_area[j],t1_vol[j],inflow1[j],coredata(inflow_vol),t1_ET[j],t1_sluice[j],t1_spill[j],t1_GW[j])
  #t1_all=rbind(t1_all,cur_all)
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  ############################################################  ############################################################
  
  PET_max[i] = Kc2_cal[doy,5]*PET[cur_month]
  TP_max[i] = Kc2_cal[doy,5]
  
  i=i+1
}