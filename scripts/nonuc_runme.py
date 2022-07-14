''' ==================================================================================================
This code creates the SPARC11b.rxn file to add SOM precursor so that it runs prep.exe and creates
saprc14_rev1.mod and SAPRC14_rev1.f files

written by Ali Akherati October 2018 Colorado State University
================================================================================================== '''

import numpy as np
import pandas as pd
import pygsheets as pyg
import time
import re
import os
import sys
import datetime as dt
import math
from scipy.optimize import curve_fit
from radiation import rad
from PBL_script import pbl,pbl2
#from BG_sizedist import SD_smps, SD_smps_fims, SD_nanosmps
from BG_sizedist import SD_smps, SD_nanosmps
from functions import g2,g3,closest,h
from matplotlib import pyplot as plt 
from OH_concentrations import OH_profile

# Clean ipython
# ====================================================================================================
#from IPython import get_ipython
#get_ipython().magic('reset -sf') 

startTime = time.time()

# ====================================================================================================
# Run the python script to run the actual model
# ====================================================================================================
#maxjobs = 100 # how many jobs you can run. For ozone, max is 16
#queue = 'defaultfaculty.q@@students'
queue = 'defaultfaculty.q'

# INPUTS
# ====================================================================================================
days = [#'04252016',
#       '05012016',
       '05162016',
       '05192016'
#       '09042016',
#       '09132016',
#       '09202016'
        ]

# ===== ===============================================================================================
for day  in days:
  emiss_sheetname = day  # Google sheet name
  emiss_scheme = 3       # 1, 2, 3, 4
  OH_proxy = 2           # 1, 2, 3
  NAME = '%s_E%s_OH%s_RHtest'%(day,emiss_scheme,OH_proxy) # Name for output files
  use_fims = False       # Use the vertical profile from FIMS data? True/False

  # nucleation parameters
  # ====================================================================================================
  fion = 8.0       # Ion recombination coefficient [cm-3 s-]
  cstar_nuc = 0.0 #1e-9 # Volatility bound for organics to contribute to nucleation [ug m-3]
  org_nuc = 1      # switch for organic nucleation [0 or 1]
  inorg_nuc = 1    # switch for inorganic nucleation [0 or 1]
  # ====================================================================================================

  #           BNZ,   TOL,   XYL,    ISP,   TRP
  #DLVP_v = '1.97 1.77 2.05 2.25 1.97 1.97 1.97' 
  DLVP_v = '2.075 1.966 1.755 2.245 1.571 1.571 1.571' 
  
  aadt = 30     # microphysics timestep [s]
  nintern = 1   # frequency of printing values or number of internal microphysics steps
  # adt = aadt*nintern # timestep for writing output [seconds]  
  nlayers = 20  # Number of vertical layers in model 
  FMT = '(%2i(F6.2,2x))'%nlayers  # This is a formatt string - this could be eliminated 
  
  COAG = 1    # [0 or 1] the switch for On/Off coagulation
  VWL  = [1]  # [0 or 1] the switch for On/Off vapor wall loss
  PWL  = [0]  # [0 or 1] the switch for On/Off particle wall loss
  
  ke      = 0.0        # loss rate constant
  kw0     = 0.0        # loss rate constant
  
  alpha   = 1.0        # accommodation coefficient 
  storg   = 0.025      # [N/m] surface tension
  Dbk     = [1.0E-10]  # particle-phase diffusion coefficient [m2/s] (could try 3.4e-15 from Charles) 
  kc      = 0.0        # first-order loss rate of species in the particle phase [1/s]
  boxvol  = 10000000.0 # teflon [cm3] - [CalTech 24 m3, CSU 10 m3, CMU ?? m3]
  Kz      = 2000       # Vertical diffusion coefficient 
  stptemp = 273.15     # STP temperature 
  stppres = 101325.0   # STP pressure 
  
  # Set the parameters for the dat chosen 
  # ====================================================================================================
  
  if day == '04252016':
    endtime = [10.0]        # Length of model run
    upper_pres = 70000.0    # Pressure of upper box
    LR_K = 9.8              # Lapse rate - currently not used 
    NOx = 1.7               # NOx concentration for OH proxy #3
    orgfrac_bg = 0.8        # Organic fraction of background aerosols  
    nh3_ppt = [5000.]       # NH3 concentration [ppt]
    OH_lower = dt.datetime(2016,4,25,9,0)   # This is essentially the start time 
#    OH_sheetname1 = 'sgpbeflux1longC1.c1.20160427.000000.custom'  
#    OH_sheetname2 = 'sgpbeflux1longC1.c1.20160428.000000.custom'
    #SMPS_sheetname = 'HiScaleSMPSb_SGP_20160424_R0'
    SMPS_sheetname = 'SMPS_merged_20160425'
    SMPS_lower = dt.datetime(2016,4,25,9,0)  # time at which to get the background size dist. 
    PBL_time = dt.datetime(2016,4,24,14,0)     # average time of HI-SCALE flight for PBL height
  elif day == '05012016':
    endtime = [10.0]
    upper_pres = 70000.0
    LR_K = 9.8 
    NOx = 1.65
    orgfrac_bg = 0.8 
    nh3_ppt = [5000.] 
    OH_lower = dt.datetime(2016,5,1,10,20)
#    OH_sheetname1 = 'sgpbeflux1longC1.c1.20160428.000000.custom'
#    OH_sheetname2 = 'sgpbeflux1longC1.c1.20160429.000000.custom'
    #SMPS_sheetname = 'HiScaleSMPSb_SGP_20160511_R0'
    SMPS_sheetname = 'HiScaleSMPSb_SGP_20160424_R0'
    SMPS_lower = dt.datetime(2016,5,1,10,20)
    PBL_time = dt.datetime(2016,4,30,14,30)
  elif day == '05162016':
    endtime = [10.0]
    upper_pres = 70000.0
    LR_K = 9.8 
    NOx = 1.65
    orgfrac_bg = 0.8 
    nh3_ppt = [300.] 
    OH_lower = dt.datetime(2016,5,16,9,0)
    OH_sheetname1 = 'sgpbeflux1longC1.c1.20160516.000000.custom'
    OH_sheetname2 = 'sgpbeflux1longC1.c1.20160517.000000.custom'
    #SMPS_sheetname = 'HiScaleSMPSb_SGP_20160424_R0'
    SMPS_sheetname = 'SMPS_merged_20160516'
    SMPS_lower = dt.datetime(2016,5,16,9,0)
    PBL_time = dt.datetime(2016,5,15,14,30)
    aadt = 10.0
    #OH_scale = 1.1342047483015671 #2
    OH_scale = 0.920395683297855*1.8908160768981273 #1
    #OH_scale = 1.0
  elif day == '05192016':
    endtime = [10.0]
    upper_pres = 70000.0
    LR_K = 9.8
    NOx = 1.4 
    orgfrac_bg = 0.7
    nh3_ppt = [300.]
    OH_lower = dt.datetime(2016, 5, 19,9,0)
    OH_sheetname1 = 'sgpbeflux1longC1.c1.20160519.000000.custom'
    OH_sheetname2 = 'sgpbeflux1longC1.c1.20160520.000000.custom'
    #SMPS_sheetname = 'HiScaleSMPSb_SGP_20160511_R0'
    SMPS_sheetname = 'SMPS_merged_20160519'
    SMPS_lower = dt.datetime(2016, 5, 19,9,0)
    PBL_time = dt.datetime(2016,5,18,12,30)
    aadt = 10.0
    #OH_scale = 0.6391390454033439 #2
    OH_scale = 0.06257969473152183*5.206133676195754 #1
    #OH_scale = 1.0
  elif day == '09042016':
    endtime = [10.0]
    upper_pres = 70000.0
    LR_K = 9.8
    NOx = 1.4
    orgfrac_bg = 0.7
    nh3_ppt = [5000.]
    OH_lower = dt.datetime(2016, 9, 4,9,0)
#    OH_sheetname1 = 'sgpbeflux1longC1.c1.20160514.000000.custom'
#    OH_sheetname2 = 'sgpbeflux1longC1.c1.20160515.000000.custom' 
#    SMPS_sheetname = 'SMPS_merged_20160904'
    SMPS_sheetname = 'HiScaleSMPSb_SGP_20160827_R1'
    SMPS_lower = dt.datetime(2016, 9, 4,9,0)
    PBL_time = dt.datetime(2016,9,3,12,30)
  elif day == '09132016':
    endtime = [10.0]
    upper_pres = 70000.0
    LR_K = 9.8 
    NOx = 1.8
    orgfrac_bg = 0.7 
    nh3_ppt = [5000.] 
    OH_lower = dt.datetime(2016, 9, 13,9,0)
#    OH_sheetname1 = 'sgpbeflux1longC1.c1.20160911.000000.custom'
#    OH_sheetname2 = 'sgpbeflux1longC1.c1.20160912.000000.custom'
    SMPS_sheetname = 'HiScaleSMPSb_SGP_20160827_R1'
    SMPS_lower = dt.datetime(2016, 9, 13,9,0)
    PBL_time = dt.datetime(2016,9,12,11,30)
  elif day == '09202016':
    endtime = [10.0]
    upper_pres = 70000.0
    LR_K = 9.8 
    NOx = 1.7
    orgfrac_bg = 0.9 
    nh3_ppt = [5000.] 
    OH_lower = dt.datetime(2016,9,20,9,0)
#    OH_sheetname1 = 'sgpbeflux1longC1.c1.20160917.000000.custom'
#    OH_sheetname2 = 'sgpbeflux1longC1.c1.20160918.000000.custom'
    SMPS_sheetname = 'HiScaleSMPSb_SGP_20160827_R1'
    SMPS_lower = dt.datetime(2016,9,20,9,0)
    PBL_time = dt.datetime(2016,9,19,13,30)
  else:
    sys.exit('Invalid day or no day selected')
  
  density_bg = orgfrac_bg*1400.0+(1-orgfrac_bg)*1770 # [kg/m3]
  
  
  # directory information
  # ====================================================================================================
  script_directory = os.popen('pwd').read()[:-1]
  src_directory = script_directory[:-(len('scripts'))] + 'src'
  run_directory = script_directory[:-(len('scripts'))] + 'runs'
  
  gc = pyg.authorize(service_file='../inputs/SOMTOMAS-Sam-GoogleSheetAccess.json')
  
  # Reading radiation Google sheet for OH calculations - added by SamO - 7/2/2020
  #=====================================================================================================
  # This requires two days worth of data since raw data are in UTC, so there can be carry-over from day to day. 
  if OH_proxy == 1 or OH_proxy ==3:
    DWSW = rad(OH_sheetname1,OH_sheetname2,endtime,gc,OH_lower,aadt)
    
    # Writing data to an input file to be read in box.f
    f3 = open('../inputs/%s_SW_Rad'%day,'w')
    for i in range(len(DWSW)):
      if DWSW[i]<0.0 or math.isnan(DWSW[i]) == True:
        DWSW[i] = 0.0
      f3.write('%s\n'%str(DWSW[i]))
    f3.close()
    
  #sys.exit('You suck!!!')
  
  # Pressure - Calculate the pressure levels based on surface pressure and number of layers
  #=====================================================================================================
  
  sfcpdata = '../scripts/MERRA2_P_data'    # Pressure data file
  Date = []                                # Date list
  fid3 = open(sfcpdata,'r')
  sfcpres_tmp = []
  for line in fid3.readlines():
    spl_line = line.split(' ')
    Date.append(dt.datetime.strptime(str(spl_line[0]+'-'+spl_line[1]),'%Y-%m-%d-%H:%M:%S' )-dt.timedelta(hours=5))
    sfcpres_tmp.append(spl_line [2]) #[2:-2:2])
  sfcpres_tmp = np.array(sfcpres_tmp)
 
  upper = OH_lower + dt.timedelta(hours = endtime[0])
  Date = np.array(Date)
  sfclow = np.where(Date>OH_lower)[0][0]
  sfcup = np.where(Date>upper)[0][0]
  P_cut = sfcpres_tmp[sfclow:sfcup]
  P_cut = np.array(P_cut)
  P_cut = P_cut.astype(np.float)
  print('P_cut:',P_cut)
  sfcpres = P_cut[0]
  
  #sys.exit('Sucker')
  
  x1 = np.linspace(0,len(P_cut)-1,int(endtime[0]*60*60/aadt)+1)
  x2 = np.arange(0,len(P_cut))
  P_cut = np.interp(x1,x2,P_cut)
  #sys.exit('I hate you')
  Dp_pa = (sfcpres - upper_pres)/nlayers  # Define pressure differential  
  pres_lay,temp_lay,pres_lay_edge,Z_lay = np.zeros([nlayers]),np.zeros([nlayers]),np.zeros([nlayers+1]),np.zeros([nlayers])
  pres_lay[0] = sfcpres - 0.5*Dp_pa 
  pres_lay_edge[0] = sfcpres 
  
  for i in range(nlayers-1):
    pres_lay[i+1] = pres_lay[i] - Dp_pa
  for i in range(nlayers):
    pres_lay_edge[i+1] = pres_lay_edge[i] - Dp_pa
  
  f6 = open('../inputs/%s_SfcPres'%day,'w')
  for i in range(len(P_cut)):
    f6.write('%s\n'%str(P_cut[i]))
  f6.close()
  
  #sys.exit('Sucker')
  
  # PBL Height - Get the PBL height from MERRA2 data, which is used for vertical diffusion profile
  #=====================================================================================================
  
  #PBL, RL_top, ML_top = pbl(OH_lower,endtime,PBL_time,aadt)
  PBL = pbl2(OH_lower,endtime,PBL_time,aadt)
  print('Initial PBL height =',PBL[0])
  f4 = open('../inputs/%s_PBL_Height'%day,'w')
  for i in range(len(PBL)):
    f4.write('%s\n'%str(PBL[i]))
  f4.close()
  #sys.exit('You suck')

  # Temperature - Get temperature data for the given timeframe at the right pressure levels
  #=====================================================================================================
  #sys.exit('I hate you for what you have done')
  tempdata = '../scripts/MERRA2_T_72layer_data'
  qvdata = '../scripts/MERRA2_QV_72layer_data'
  etamid_data = '../scripts/MERRA2_ETA_Mid'
  etaedge_data = '../scripts/MERRA2_ETA_Edge'
  
  Date = []
  Temp_tmp = []
  QV_tmp = []
  pres_lev = []
  eta_edge = []
  eta_lay = pres_lay/sfcpres
  eta_lay_edge = pres_lay_edge/sfcpres
  
  fid = open(tempdata,'r')
  fid2 = open(etamid_data,'r')
  fid3 = open(etaedge_data,'r')
  fid4 = open(qvdata,'r')
  
  for line in fid2.readlines():
    spl_line = line.split(' ')
    pres_lev.append(float(spl_line[0]))
  pres_lev = np.flip(np.array(pres_lev))
  
  for line in fid3.readlines():
    spl_line = line.split(' ')
    eta_edge.append(float(spl_line[0]))
  
  eta_edge = np.flip(np.array(eta_edge))
  
  for line in fid.readlines():
    spl_line = line.split(' ')
    Date.append(dt.datetime.strptime(str(spl_line[0]+'-'+spl_line[1]),'%Y-%m-%d-%H:%M:%S' )-dt.timedelta(hours=5))
    Temp_tmp.append(spl_line[2:-2:2])
  Temp_tmp = np.array(Temp_tmp)

  for line in fid4.readlines():
    spl_line = line.split(' ')
    QV_tmp.append(spl_line[2:-2:2])
  QV_tmp = np.array(QV_tmp)

  indecies = []
  for i in range(nlayers):
    indecies.append(closest(pres_lev,eta_lay[i]))
  
  Date = np.array(Date)
  low = np.where(Date>SMPS_lower)[0][0]
  up = np.where(Date>upper)[0][0]
  T_cut = Temp_tmp[low:up,indecies]
  T_cut = np.array(T_cut)
  T_cut = T_cut.astype(np.float)
  QV_cut = QV_tmp[low:up,indecies]
  QV_cut = np.array(QV_cut)
  QV_cut = QV_cut.astype(float)

  x1 = np.linspace(0,len(T_cut[:,0])-1,int(endtime[0]*60*60/aadt)+1)
  x2 = np.arange(0,len(T_cut[:,0]))
  TempK = np.zeros([nlayers,len(x1)])
  QV = np.zeros([nlayers,len(x1)])
  for i in range(nlayers):
    TempK[i] = np.interp(x1,x2,T_cut[:,i])
    QV[i] = np.interp(x1,x2,QV_cut[:,i])
 
  temp_lay[:] = TempK[:,0]

  f5 = open('../inputs/%s_Temperatures'%day,'w')
  for i in range(len(TempK[0])):
    writer = ''
    for j in range(nlayers):
      writer = writer + '%6.2f  '%TempK[j,i]
    f5.write('%s \n'%writer)
  f5.close()
  print('temp_lay =',temp_lay)

  f7 = open('../inputs/%s_QV'%day,'w')
  for i in range(len(QV[0])):
    writer = ''
    for j in range(nlayers):
      writer = writer + '%11.9f  '%QV[j,i]
    f7.write('%s \n'%writer)
  f7.close()


  Z_lay[0] = (287.0*temp_lay[0]/9.8 * np.log(pres_lay_edge[0]/pres_lay_edge[1]))*0.5
  
  for i in range(nlayers-1):
    Z_lay[i+1] = Z_lay[i] + (287.0*temp_lay[i]/9.8 * np.log(pres_lay_edge[i]/pres_lay_edge[i+1]))*0.5 + (287.0*temp_lay[i+1]/9.8 * np.log(pres_lay_edge[i+1]/pres_lay_edge[i+2]))*0.5
  
  print('Z_lay =',Z_lay)
  #sys.exit('Just a temporary stop')
  
  #=====================================================================================================
  if OH_proxy == 2:
    OH_profile(OH_lower,endtime,nlayers,aadt,eta_lay,day)
  
  # Background size distribution 
  #=====================================================================================================
  if use_fims==True:
    if day=='04272016':
      No_bg1,Dpm_bg1,sigma_bg1,No_bg2,Dpm_bg2,sigma_bg2,No_bg3,Dpm_bg3,sigma_bg3 = SD_smps(SMPS_sheetname,SMPS_lower,nlayers,endtime,gc,Z_lay,temp_lay,stppres,stptemp,day,pres_lay)
    else: 
      No_bg1,Dpm_bg1,sigma_bg1,No_bg2,Dpm_bg2,sigma_bg2,No_bg3,Dpm_bg3,sigma_bg3  = SD_smps_fims(SMPS_sheetname,SMPS_lower,nlayers,endtime,gc,Z_lay,temp_lay,stppres,stptemp,day,pres_lay,ML_top)
  
  elif use_fims==False:
    #if OH_lower.month==9:
    No_bg1,Dpm_bg1,sigma_bg1,No_bg2,Dpm_bg2,sigma_bg2,No_bg3,Dpm_bg3,sigma_bg3  = SD_smps(SMPS_sheetname,SMPS_lower,nlayers,endtime,gc,Z_lay,temp_lay,stppres,stptemp,day,pres_lay)
    #elif OH_lower.month==4 or OH_lower.month==5:
     # No_bg1,Dpm_bg1,sigma_bg1,No_bg2,Dpm_bg2,sigma_bg2,No_bg3,Dpm_bg3,sigma_bg3  = SD_nanosmps(SMPS_sheetname,SMPS_lower,nlayers,endtime,gc,Z_lay,temp_lay,stppres,stptemp,day,pres_lay)

  
  # ====================================================================================================
  # Regime and Parameter names - Changed by SamO 
  # ====================================================================================================
  regime = str(cstar_nuc)
  params = temp_lay[0]
  
         
  # --- read input sheet --------------------------------------------------------------  
  No1     = 0.0 #df_input.loc[df_input['variables']=='No_1', 'Manish_10_1.5e6']# [# cm-3] background 1st number conc.
  Dpm1    = 80e-3 #df_input.loc[df_input['variables']=='Dpm_1', 'Manish_10_1.5e6'].iloc[0]*1.0e-3 # [microns] background 1st median diameter
  sigma1  = 1.8 #df_input.loc[df_input['variables']=='sigma_1', 'Manish_10_1.5e6'].iloc[0]# 1st background sigma 
  No2     = 0.0 #df_input.loc[df_input['variables']=='No_2', 'Manish_10_1.5e6'].iloc[0]# [# cm-3] background 2nd number conc.
  Dpm2    = 80e-3 #df_input.loc[df_input['variables']=='Dpm_2', 'Manish_10_1.5e6'].iloc[0]*1.0e-3# [microns] background 2nd median diameter
  sigma2  = 1.8 #df_input.loc[df_input['variables']=='sigma_2', 'Manish_10_1.5e6'].iloc[0]# 2nd background sigma 
  
  # the following values are for OH concentration equation
  lights_on = 0.0 #df_input.loc[df_input['variables']=='lights on', 'Manish'].iloc[0] # [s] of the day
  a_oh      = [0.0]#, 3e6, 5e6, 1e7, 5e7, 1e8] #df_input.loc[df_input['variables']=='a_oh', 'Manish'].iloc[0]
  ax_oh     = 0.0 #df_input.loc[df_input['variables']=='ax_oh', 'Manish'].iloc[0]
  b_oh      = 0.0 #df_input.loc[df_input['variables']=='b_oh', 'Manish'].iloc[0]
  bx_oh     = 0.0 #df_input.loc[df_input['variables']=='bx_oh', 'Manish'].iloc[0]
  OH0       = np.array(a_oh)*np.exp(-1.*ax_oh*0.) + b_oh*np.exp(-1.*bx_oh*0.)
  #OH_scale  = 1.0 #2.2  # scalar for OH concentration 
  
  
  nsomprec    = 7 # number of som precursors for parameterizations or SOM grids ***** 1 for now
  somprecname = 'BNZSOMG TOLSOMG XYLSOMG ISPSOMG TRPSOMG IVOSOMG SVOSOMG' # som precursor parameterization's name
  dlvp        = DLVP_v
  
  seed_dens = 1.4 #df_input.loc[df_input['variables']=='POA density', 'AGU2019_F007'].iloc[0]*1.0e+3
  poa_1stname = 'SVO_C10'
  poa_1st_lenname = len(poa_1stname)
  
  #===================================================================================
#  gsh_emission_ML = gc.open('ML_PBL_SOM-TOMAS.2')
#  df_emiss_ML = gsh_emission_ML.worksheet_by_title(title='%s'%emiss_sheetname).get_as_df()
  
#  gsh_emission_ML = gc.open('PTRMS_0.0LOD_SOM-TOMAS.2')
#  df_emiss_ML = gsh_emission_ML.worksheet_by_title(title='%s'%emiss_sheetname).get_as_df()
  
#  gsh_emission_ML = gc.open('Lower_PTRMS_LOD_SOM-TOMAS.2')
#  df_emiss_ML = gsh_emission_ML.worksheet_by_title(title='%s'%emiss_sheetname).get_as_df()
  gsh_emission_ML = gc.open('No_Nuc_Days')
  df_emiss_ML = gsh_emission_ML.worksheet_by_title(title='%s'%emiss_sheetname).get_as_df()

  c5_scale_gc = gc.open('%s_C5s_scale'%day)
  df_c5_scale = np.array(c5_scale_gc.worksheet_by_title(title='Sheet1').get_as_df())
  
  c8_scale_gc = gc.open('%s_C8andC9_scale'%day)
  df_c8_scale = np.array(c8_scale_gc.worksheet_by_title(title='Sheet1').get_as_df())
  
  so2_scale_gc = gc.open('%s_SO2_scale'%day)
  df_so2_scale = np.array(so2_scale_gc.worksheet_by_title(title='Sheet1').get_as_df())
  #so2_scale_gc = gc.open('%s_CIMSso2_scale'%day)
  #df_so2_scale = np.array(so2_scale_gc.worksheet_by_title(title='Sheet1').get_as_df())
  
#  rh_vert_gc = gc.open('%s_RH'%day)
#  df_rh = np.array(rh_vert_gc.worksheet_by_title(title='Sheet1').get_as_df())[:,1]

  #===================================================================================
   
  lin_alt = np.array(df_so2_scale[:,1])
  lin_c5_scale = np.array(df_c5_scale[:,2])
  lin_c8_scale = np.array(df_c8_scale[:,2])
  lin_so2_scale = np.array(df_so2_scale[:,2])
  c5_scale = np.empty(nlayers)
  c8_scale = np.empty(nlayers)
  so2_scale = np.empty(nlayers)
  lin_alt_scale = np.empty(nlayers)
#  rh_lay = np.empty(nlayers)
  
  for i in range(nlayers):
    indx = closest(lin_alt,Z_lay[i])
    c5_scale[i] = np.nanmean(lin_c5_scale[indx-10:indx+10])
    c8_scale[i] = np.nanmean(lin_c8_scale[indx-10:indx+10])
    so2_scale[i] = np.nanmedian(lin_so2_scale[indx-20:indx+20])
#    rh_lay[i] = np.nanmean(df_rh[indx-5:indx+5])
    #rh_lay[i] = np.nanmean(df_rh[0:10])
#    if rh_lay[i]>0.9:
#      rh_lay[i]=0.9
#  print('rh_lay=',rh_lay)
  lin_alt_scale = lin_alt
  #sys.exit('Sucker')
  
  print('so2_scale =',so2_scale)
  print('c5_scale =',c5_scale)
  print('c8_scale =',c8_scale)
  print('PBL[0] =',PBL[0])
  #print('RL_top =',RL_top)
  #print('ML_top =',ML_top)
  
  iorg  = 456  # SamO changed from 893
  ibins = 40   # number of bins
  idiag = 2    # number of diagnostic species including sulfate and ammonia (2)
  
  ML_top = PBL[0]
  
  index  = np.where(df_emiss_ML.iloc[:,0]!='')[0].shape[0] 
  spname = df_emiss_ML['species'].iloc[:index]
  g_ippm_ML = df_emiss_ML['gas_id'].iloc[:index]
  p_frac = df_emiss_ML['pfrac_id'].iloc[:index] 
 
  nspemiss = len(spname)
  emiss_spname =  ''
  emiss_ippm   = ['' for i in range(nlayers)]
  seed_frac    =  ''
  
  if emiss_scheme==1:
    for i in range(len(spname)):
      emiss_spname = emiss_spname + '%15s'%(spname.iloc[i])
      seed_frac =    seed_frac + '%15.5E'%(p_frac.iloc[i])
      for j in range(nlayers):
        if Z_lay[j]<PBL[0]:
          if spname.iloc[i] == 'SO2':
            emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*np.nanmean(lin_so2_scale[np.where(lin_alt_scale<ML_top)[0][0]:np.where(lin_alt_scale>ML_top)[0][0]]))
          elif spname.iloc[i]=='BENZENE' or spname.iloc[i]=='TOLUENE' or spname.iloc[i]=='XYLENE' or spname.iloc[i]=='ISOPRENE':
            emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*np.nanmean(lin_c5_scale[np.where(lin_alt_scale<ML_top)[0][0]:np.where(lin_alt_scale>ML_top)[0][0]]))
          elif spname.iloc[i]=='TRIMBNZ' or spname.iloc[i]=='TERPENE' or spname.iloc[i]=='SESQTRP' or spname.iloc[i]=='IVOC' or spname.iloc[i]=='SVOC':
            emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*np.nanmean(lin_c8_scale[np.where(lin_alt_scale<ML_top)[0][0]:np.where(lin_alt_scale>ML_top)[0][0]]))
#        elif Z_lay[j]>PBL[0] and Z_lay[j]<RL_top:
#          if np.where(lin_alt_scale>ML_top)[0][0]!=np.where(lin_alt_scale>RL_top)[0][0]:
#            if spname.iloc[i] == 'SO2':
#              emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*np.nanmean(lin_so2_scale[np.where(lin_alt_scale>ML_top)[0][0]:np.where(lin_alt_scale>RL_top)[0][0]]))
#            elif spname.iloc[i]=='BENZENE' or spname.iloc[i]=='TOLUENE' or spname.iloc[i]=='XYLENE' or spname.iloc[i]=='ISOPRENE':
#                emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*np.nanmean(lin_c5_scale[np.where(lin_alt_scale>ML_top)[0][0]:np.where(lin_alt_scale>RL_top)[0][0]]))
#            elif spname.iloc[i]=='TRIMBNZ' or spname.iloc[i]=='TERPENE' or spname.iloc[i]=='SESQTRP' or spname.iloc[i]=='IVOC' or spname.iloc[i]=='SVOC':
#                emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*np.nanmean(lin_c8_scale[np.where(lin_alt_scale>ML_top)[0][0]:np.where(lin_alt_scale>RL_top)[0][0]]))
        else:
          if spname.iloc[i] == 'SO2':
            emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*lin_so2_scale[np.where(lin_alt_scale>ML_top)[0][0]])
          elif spname.iloc[i]=='BENZENE' or spname.iloc[i]=='TOLUENE' or spname.iloc[i]=='XYLENE' or spname.iloc[i]=='ISOPRENE':
            emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*lin_c5_scale[np.where(lin_alt_scale>ML_top)[0][0]])
          elif spname.iloc[i]=='TRIMBNZ' or spname.iloc[i]=='TERPENE' or spname.iloc[i]=='SESQTRP' or spname.iloc[i]=='IVOC' or spname.iloc[i]=='SVOC':
            emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*lin_c8_scale[np.where(lin_alt_scale>ML_top)[0][0]])
#        elif Z_lay[j]>RL_top:
#          if spname.iloc[i] == 'SO2':
#            emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*np.nanmean(lin_so2_scale[np.where(lin_alt_scale>RL_top)[0][0]:]))
#          elif spname.iloc[i]=='BENZENE' or spname.iloc[i]=='TOLUENE' or spname.iloc[i]=='XYLENE' or spname.iloc[i]=='ISOPRENE':
#            emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*np.nanmean(lin_c5_scale[np.where(lin_alt_scale>RL_top)[0][0]:]))
#          elif spname.iloc[i]=='TRIMBNZ' or spname.iloc[i]=='TERPENE' or spname.iloc[i]=='SESQTRP' or spname.iloc[i]=='IVOC' or spname.iloc[i]=='SVOC':
#            emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*np.nanmean(lin_c8_scale[np.where(lin_alt_scale>RL_top)[0][0]:]))
#        if Z_lay[j]>3000.0:
#            emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*0.01)
  
  if emiss_scheme==2:
    for i in range(len(spname)):
      emiss_spname = emiss_spname + '%15s'%(spname.iloc[i])
      seed_frac =    seed_frac + '%15.5E'%(p_frac.iloc[i])
      for j in range(nlayers):
        emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i])

  if emiss_scheme==3:
    for i in range(len(spname)):
      emiss_spname = emiss_spname + '%15s'%(spname.iloc[i])
      seed_frac =    seed_frac + '%15.5E'%(p_frac.iloc[i])
      for j in range(nlayers):
        if spname.iloc[i] == 'SO2':
          emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*so2_scale[j])
        elif spname.iloc[i]=='BENZENE' or spname.iloc[i]=='TOLUENE' or spname.iloc[i]=='XYLENE' or spname.iloc[i]=='ISOPRENE':
          emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*c5_scale[j])
        elif spname.iloc[i]=='TRIMBNZ' or spname.iloc[i]=='TERPENE' or spname.iloc[i]=='SESQTRP' or spname.iloc[i]=='IVOC' or spname.iloc[i]=='SVOC':
          emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*c8_scale[j])
  
  
  if emiss_scheme==4:
    for i in range(len(spname)):
      emiss_spname = emiss_spname + '%15s'%(spname.iloc[i])
      seed_frac =    seed_frac + '%15.5E'%(p_frac.iloc[i])
      for j in range(nlayers):
        if Z_lay[j]<ML_top:
#          if i == 0:
          emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i])
#          else:
#            emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i])
#        elif Z_lay[j]>ML_top and Z_lay[j]<RL_top:
        else:
#          if np.where(lin_alt_scale>ML_top)[0][0]!=np.where(lin_alt_scale>RL_top)[0][0]:
#            if i == 0:
          emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*0.5)
#            else:
#              emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*4.0)
#          else:
#            if i == 0:
#              emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*4.0)
#            else:
#              emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*4.0)
#        elif Z_lay[j]>RL_top:
#          if i == 0:
#            emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*0.5)
#          else:
#            emiss_ippm[j] = emiss_ippm[j] + '%15.5E'%(g_ippm_ML.iloc[i]*0.5)
  
  
  # Dilution
  #firesize = [1e2, 1e-4, 1e-2, 1e0] # [km2]
  #p_dilt1 = [0.01335, -8381.0, -106.5, -0.9539]
  #p_dilt2 = [-0.09625, 4.799e+04, 608.2, 4.61]
  #p_dilt3 = [0.3197, 2.992e+04, 312.7, 8.166]
  #p_dilt4 = [1.0, 1.0, 1.0, 1.0]
  firesize = [1e0] # [km2]
  p_dilt1  = [-0.9539]
  p_dilt2  = [4.61]
  p_dilt3  = [8.166]
  p_dilt4  = [1.0]
  
  # Running 
  # ====================================================================================================
  print('General information for the following simulations')  
  print('======================================================')
  print('Simulation: %s'%NAME)            
  print('inbins = %s'%ibins)
  print('COAG = %s'%COAG)
  print('boxvol = %s'%boxvol, '[cm3]')
  print('alpha = %s'%alpha)
  print('surface pressure = %s'%sfcpres, '[Pa]')
  print('1st layer temperature = %s'%temp_lay[0], '[K]')
  print('nh3_ppt = %s'%nh3_ppt)
  # ====================================================================================================
  
  ctr = 0
  for vwl in VWL:
   for pwl in PWL:
    for OH_conc in a_oh:
     for tend in endtime:
      for db in Dbk:
       for NH3 in nh3_ppt:
        ctr+=1
        rname = '%s_Cstar%s_T%6.2f_NH3%4.2e_nlay%i_hr%4.2e_ricco%i_dunn%i_bg10'%(NAME,regime,params,NH3,nlayers,tend,org_nuc,inorg_nuc)
        if pwl==1:
         rname = 'PWL_'+rname
         #if vwl==1:
         #    rname = 'VWL_'+rname
         #rname = precursor+'_'+regime+'_'+tend+'_'+no+'_'+dp+'_'+oh+'_'+ppm
         
        print('%i - Simulation'%ctr)
        print('---------------')
        print('VWL = %s'%vwl)
        print('PWL = %s'%pwl)
        print('Db = %s [m2/s]'%db)
        print('Kc = %s [1/s]'%kc)
        print('ke = %s [1/s]'%ke)
        print('OH = %s [molec/cm3]'%OH_conc)
        print('No1 = %s [#/cm3]'%No1)
        print('Dpm1 = %s [um]'%Dpm1)
        print('sigma1 = %s'%sigma1)
        print('No2 = %s [#/cm3]'%No2)
        print('Dpm2 = %s [um]'%Dpm2)
        print('sigma2 = %s'%sigma2)
        print('time = %s [h]'%tend)
        print('runname = %s'%rname)
        
        if os.path.exists('../runs/%s'%(rname)):
         os.system('rm -r ../runs/%s'%(rname))
          
        os.system('mkdir ../runs/%s'%(rname))
        os.system('cp ../src/box.exe ../runs/%s/'%(rname))
        os.system('cp ../src/saprc14_rev1.mod ../runs/%s/saprc14_rev1.mod'%(rname))
          
        f1 = open('../runs/%s/input'%(rname),'w')
        f1.write('%s\n'%rname)
        f1.write('%f\n'%aadt)
        f1.write('%i\n'%nintern)
        f1.write('%1i\n'%COAG)
        f1.write('%1i\n'%vwl)
        f1.write('%1i\n'%pwl)
        f1.write('%f\n'%ke)
        f1.write('%f\n'%kw0)
        f1.write('%f\n'%lights_on)
        f1.write('%i\n'%OH_proxy)
        f1.write('%e\n'%OH_conc) # ***
        f1.write('%e\n'%ax_oh) # ***
        f1.write('%e\n'%b_oh) # ***
        f1.write('%e\n'%bx_oh) # ***
        f1.write('%e\n'%OH_scale)
        f1.write('%e\n'%NOx) 
        f1.write('%15.4f\n'%fion)
        f1.write('%8.1e\n'%cstar_nuc)
        f1.write('%i\n'%org_nuc)
        f1.write('%i\n'%inorg_nuc)
        f1.write('%8.1e\n'%NH3)
        f1.write('%10.1f\n'%boxvol)
        f1.write('%5.2f\n'%tend) # ***
        f1.write('%8.5f\n'%alpha)
        f1.write('%f\n'%kc)
        f1.write('%8.5f\n'%storg)
        f1.write('%2i\n'%nlayers)
        f1.write('%10.3f\n'%upper_pres)
        f1.write('%s\n'%FMT)
        f1.write('%10.3f\n'%Kz)
        #f1.write('%6.2f\n'%temp_lay[0])
        f1.write('%6.2f\n'%LR_K)
        #for i in range(nlayers):
        #  f1.write('%8.5f\n'%rh_lay[i])
        f1.write('%8.5f\n'%No1) # ***
        f1.write('%8.5f\n'%Dpm1) # ***
        f1.write('%8.5g\n'%sigma1)
        f1.write('%8.5f\n'%No2) # ***
        f1.write('%8.5f\n'%Dpm2) # ***
        f1.write('%8.5g\n'%sigma2)
        f1.write('%02i\n'%nsomprec)
        f1.write('%s\n'%somprecname)
        f1.write('%s\n'%dlvp)
        f1.write('%s\n'%poa_1stname)
        f1.write('%s\n'%poa_1st_lenname)
        f1.write('%03i\n'%nspemiss)
        f1.write('%s\n'%emiss_spname)
        for i in range(nlayers):
          f1.write('%s\n'%emiss_ippm[i])
        f1.write('%s\n'%seed_frac)
        f1.write('%s\n'%seed_dens)
        f1.write('%s\n'%No_bg1)
        f1.write('%s\n'%Dpm_bg1) 
        f1.write('%s\n'%sigma_bg1)
        f1.write('%s\n'%No_bg2)
        f1.write('%s\n'%Dpm_bg2)
        f1.write('%s\n'%sigma_bg2)
        f1.write('%s\n'%No_bg3)
        f1.write('%s\n'%Dpm_bg3)
        f1.write('%s\n'%sigma_bg3)
        f1.write('%8.5f\n'%density_bg)
        f1.write('%8.6f\n'%orgfrac_bg)
        f1.write('%20.4f\n'%p_dilt1[0])
        f1.write('%20.4f\n'%p_dilt2[0])
        f1.write('%20.4f\n'%p_dilt3[0])
        f1.write('%20.4f\n'%p_dilt4[0])
        f1.close()
        
        # info file
        # ---------------------------------------------------------
        if os.path.exists('../outputs/%s.input'%(rname)):
         print('YES - %s.input exists!'%(rname))
         os.system('rm ../outputs/%s*'%(rname))
        f2 = open('../outputs/%s.input'%rname,'w')
        f2.write('filename    = %s\n'%rname)
        #f2.write('parameterization: %s_%s\n'%(precursor,regime))
        f2.write('aadt        = %s\n'%aadt)
        f2.write('nintern     = %s\n'%nintern)
        f2.write('ibins       = %s\n'%ibins)
        f2.write('COAG        = %s\n'%COAG)
        f2.write('VWL         = %s\n'%vwl)
        f2.write('PWL         = %s\n'%pwl)
        f2.write('ke          = %s [1/s]\n'%ke)
        f2.write('kw0         = %s [1/s]\n'%kw0)
        f2.write('LightsOn    = %s [s]\n'%lights_on)
        f2.write('OH_proxy    = %s \n'%OH_proxy)
        f2.write('OH eqn. = a_oh*exp(-1.*ax_oh*t) + b_oh*exp(-1.*bx_oh*t)\n') # ***
        f2.write('a_oh        = %s\n'%OH_conc) # ***
        f2.write('ax_oh       = %s\n'%ax_oh) # ***
        f2.write('b_oh        = %s\n'%b_oh) # ***
        f2.write('bx_oh       = %s\n'%bx_oh) # ***
        f2.write('OH_scale    = %s\n'%OH_scale)
        f2.write('NOx         = %s\n'%NOx)
        f2.write('fion        = %s\n'%fion)
        f2.write('cstar_nuc   = %s\n'%cstar_nuc)
        f2.write('org_nuc     = %s\n'%org_nuc)
        f2.write('iorg_nuc    = %s\n'%inorg_nuc)
        f2.write('nh3_ppt     = %s\n'%NH3)
        f2.write('boxvol      = %s [cm3]\n'%boxvol)
        f2.write('endtime     = %s [hours]\n'%tend) # ***
        f2.write('alpha       = %s\n'%alpha)
        f2.write('Db          = %s [m2/s]\n'%db)
        f2.write('Kc          = %s [1/s]\n'%kc)
        f2.write('storg       = %s [N/m]\n'%storg)
        f2.write('nlayers     = %s\n'%nlayers)
        f2.write('upper pressure= %s [Pa]\n'%upper_pres)
        f2.write('Format string= %s\n'%FMT)
        f2.write('lev-1 temperature = %s [K]\n'%temp_lay[0])
        f2.write('LR_K        = %s [K/km]\n'%LR_K)
        #for i in range(nlayers):
       #   f2.write('RH lay    = %s\n'%rh_lay[i])
        f2.write('No1         = %s [# cm-3]\n'%No1) # ***
        f2.write('Dpm1        = %s [um]\n'%Dpm1) # ***
        f2.write('Sigma1      = %s\n'%sigma1)
        f2.write('No2         = %s [# cm-3]\n'%No2) # ***
        f2.write('Dpm2        = %s [um]\n'%Dpm2) # ***
        f2.write('Sigma2      = %s\n'%sigma2)
        f2.write('nsomprec    = %s\n'%nsomprec)
        f2.write('somprecname = %s\n'%somprecname)
        f2.write('dlvp        = %s\n'%dlvp)
        f2.write('poa_1stname = %s\n'%poa_1stname)
        f2.write('poa_1st_len = %s\n'%poa_1st_lenname)
        f2.write('nspemiss    = %s\n'%nspemiss)
        f2.write('emiss_spname= %s\n'%emiss_spname)
        for i in range(nlayers):
          f2.write('emiss_ippm  = %s\n'%emiss_ippm[i])
        f2.write('seed_frac   = %s\n'%seed_frac)
        f2.write('seed_dens   = %s\n'%seed_dens)
        f2.write('No_bg1      = %s\n'%No_bg1)
        f2.write('Dpm_bg1     = %s\n'%Dpm_bg1) 
        f2.write('sigma_bg1   = %s\n'%sigma_bg1)
        f2.write('No_bg2      = %s\n'%No_bg2)
        f2.write('Dpm_bg2     = %s\n'%Dpm_bg2)
        f2.write('sigma_bg2   = %s\n'%sigma_bg2)
        f2.write('No_bg3      = %s\n'%No_bg3)
        f2.write('Dpm_bg3     = %s\n'%Dpm_bg3)
        f2.write('sigma_bg3   = %s\n'%sigma_bg3)
        f2.write('density_bg  = %s\n'%density_bg)
        f2.write('orgfrac_bg  = %s\n'%orgfrac_bg)
        f2.write('p_dilt1     = %s\n'%p_dilt1[0])
        f2.write('p_dilt2     = %s\n'%p_dilt2[0])
        f2.write('p_dilt3     = %s\n'%p_dilt3[0])
        f2.write('p_dilt4     = %s\n'%p_dilt4[0])
        
        f2.close()
        
        #sys.exit('Sucker')
        # ---------------------------------------------------------
        # ************************************************
        #c=os.popen('qstat').read()
        #b=re.findall('alia', c)
        #njobs=len(b)
        #print('b=',b)
        #print('njobs=',len(b))
        
        #while njobs >= maxjobs:
        #   os.system('sleep 5')
        #   c=os.popen('ps -u alia').read()
        #   b=re.findall('box.exe', c)
        #   njobs=len(b)
        # ************************************************
        
        #os.system('sleep 5')
  
        df4 = open('../inputs/raw_runme.sh', mode='r') # reading the default header file of TOMAS
        runme_line = []
        for i in df4.readlines():
         runme_line.append(i.strip('\n'))
         df4.close()
                 
        ind1 = np.where(np.array(runme_line)=='cd run_direct')[0][0]
        runme_line[ind1] = 'cd %s/%s'%(run_directory, rname)
        
        ind2 = np.where(np.array(runme_line)=="echo 'rname'")[0][0]
        runme_line[ind2] = 'echo %s'%(rname)
        
        ind3 = np.where(np.array(runme_line)=='./box.exe <input> /dev/null')[0][0]
        #runme_line[ind3] = './box.exe <input> /dev/null'
        runme_line[ind3] = './box.exe <input> %s.out'%(rname)
        
        #if os.path.exists('runme.sh'):
        #   os.system('rm runme.sh')
        
        runme_out = open('../runs/%s/runme.sh'%(rname), mode='w')
        for i in runme_line:
         runme_out.write(i)
         runme_out.write('\n')
        runme_out.close()
  
        os.system('cd ../runs/%s/; qsub -cwd -V -pe MPI 1 -q %s ./runme.sh'%(rname, queue))
  
        #os.system('sleep 1')
          
  #keep the compute node from logging out before each box.exe has finished
  #test=2
  #while test>=1:
  #    c=os.popen('ps -u alia').read()
  #    b=re.findall('box.exe', c)
  #    test=len(b)
  #    os.system('sleep 2')



