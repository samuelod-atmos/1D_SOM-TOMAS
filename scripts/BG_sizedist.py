import numpy as np
import datetime as dt
from functions import g1,g2,g3,h,h2,closest,smoothdist
from scipy.optimize import curve_fit
import pygsheets as pyg
import sys
from matplotlib import pyplot


gc = pyg.authorize(service_file='../inputs/SOMTOMAS-Sam-GoogleSheetAccess.json')

############################################################################################################

def SD_smps(SMPS_sheetname,SMPS_lower,nlayers,endtime,gc,Z_lay,temp_lay,stppres,stptemp,day1,pres_lay):

  if day1[0:2] == '09' or day1 == '05012016':
    gsh_smps = gc.open('HiScaleSMPS')
    df_1 = gsh_smps.worksheet_by_title(title='%s'%SMPS_sheetname).get_as_df()
    
    temper = df_1.to_numpy()
    time3 = temper[:,0]
    Bins = temper[:110,-1]
    flag = temper[:,-2]
    sizedist = temper[:,1:-2].astype(float)
    
    year = int(SMPS_sheetname[17:-7])
    month = int(SMPS_sheetname[-7:-5])
    day = int(SMPS_sheetname[-5:-3])
    
    startT = dt.datetime(year,month,day,0,0) - dt.timedelta(hours=5)
    date = []
    
    for i in range(len(time3)):
      date.append(startT + dt.timedelta(seconds=time3[i]))
      if flag[i]>0.:
        sizedist[i].fill('NaN')
    SMPS_upper = SMPS_lower + dt.timedelta(hours=int(endtime[0]))
    
    init_date,init_sd = ([] for i in range(2))
    
    for i in range(len(date)):   
      if date[i] >= SMPS_lower and date[i] <= SMPS_lower+dt.timedelta(minutes=30): 
        init_date.append(date[i])
        init_sd.append((sizedist[i]))
  #    if date[i] >= SMPS_upper and date[i] <= SMPS_upper+dt.timedelta(minutes=15): 
  #      final_date.append(date[i])
  #      final_sd.append((sizedist[i]))
  
  else:
    gsh_smps = gc.open('SMPS_merged')
#    df_1 = gsh_smps.worksheet_by_title(title='%s'%SMPS_sheetname).get_as_df()
    year = SMPS_lower.year
    month = SMPS_lower.month
    day = SMPS_lower.day
  
    df_1 = gsh_smps.worksheet_by_title(title='SMPS_merged_%s'%(str(year)+str(month).zfill(2)+str(day).zfill(2))).get_as_df()
    
    temper = df_1.to_numpy()
    time3 = temper[:,0].astype(float)
  #  print('time3 =',time3)
  #  sys.exit('You suck')
  #  flag = temper[:,-2]
    sizedist = temper[:,1:].astype(float)
    
    Bins = []
    bin_data = 'SMPS_merged_dp.txt'
    bin_fid = open(bin_data, 'r')
    for line in bin_fid.readlines():
      #print(line)
      spl_line = line.split('	')
    #print(type(spl_line))
    #print(spl_line)
    for i in range(len(spl_line)):
      Bins.append(float(spl_line[i]))
    Bins = np.array(Bins)
  #  print('Bins =',Bins)
  #  sys.exit('You SUCK!!!!!')
  #  year = int(SMPS_sheetname[12:15])
  #  month = int(SMPS_sheetname[15:17])
  #  day = int(SMPS_sheetname[17:19])
  #  print(year,month,day)
  #  year,month,day = 2016,5,14
    startT = dt.datetime(year,month,day,0,0) - dt.timedelta(hours=5)
    date = []
    for i in range(len(time3)):
  #    date.append(startT + dt.timedelta(seconds=time3[i]))
      date.append(dt.datetime.fromordinal(int(time3[i])) + dt.timedelta(days=time3[i]%1) - dt.timedelta(days = 366) - dt.timedelta(hours=5))
      #date.append(startT + dt.timedelta(seconds=int((time3[i] - time3[1])*24*60*60)))
  #    if flag[i]>0.:
  #      sizedist[i].fill('NaN')
  #  print(date[4]<SMPS_lower+dt.timedelta(minutes=60))
  #  print(SMPS_lower+dt.timedelta(minutes=60))
    
    SMPS_upper = SMPS_lower + dt.timedelta(hours=int(endtime[0]))
    init_date,init_sd= ([] for i in range(2))
    for i in range(len(date)):
      if date[i] >= SMPS_lower and date[i] <= SMPS_lower+dt.timedelta(minutes=20):
        init_date.append(date[i])
        #init_sd.append((np.log10(sizedist[i])))
        init_sd.append((sizedist[i]))
  init_date = np.array(init_date)
  init_sd = np.array(init_sd)
  #print('init_sd[0] =',init_sd[0])
  #print('init_sd =',init_sd)
  init_sd = np.nanmean(init_sd,axis=0)
#  init_sd = smoothdist(init_sd[65:])
#  init_sd = smoothdist(init_sd)
#  init_sd = smoothdist(init_sd)
#  init_sd = smoothdist(init_sd)
#  init_sd = smoothdist(init_sd)
#  init_sd = smoothdist(init_sd)
#  init_sd = init_sd[54:]
#  Bins = Bins[54:]
#  valid = ~np.isnan(init_sd)
#  final_date = np.array(final_date)
#  final_dsd = np.array(final_sd)


  stat=0
  guess1 = 100
  guess2 = 100
  guess3 = 100

#  while stat==0:
#    #print('Trying to fit:',i,guess1,guess2,guess3)
#    try:
#      #print('Trying')
#      p1,pcov = curve_fit(h, Bins, init_sd, p0=[guess2, guess1, guess3, 1.5, 1.5, 1.5, 4, 100, 400],maxfev=8000)
#      #print('p1 =',p1)
#      #stat = 1
#    except:
#      print('changing guess')
#      #sys.exit('Excepting')
#      guess1 = guess1+1000
#      if guess1>10000:
#        guess1=1000
#        guess2 = guess2+1000
#      if guess2>10000:
#        guess2=100
#        guess3 = guess3+1000
#      if guess3>10000:
#        pyplot.semilogx(Bins,init_sd)
#        pyplot.show( )
#        sys.exit('guess3 maxed out')
#        stat=1
#    else:
#      #print('p1 =',p1)
#      #pyplot.semilogx(fims_bins,temp,fims_bins,h(fims_bins,p1[0],p1[1],p1[2],p1[3],p1[4],p1[5],p1[6],p1[7],p1[8]))
#      #pyplot.show()
#      stat=1
  if day1 == '09112016' or day1 == '09172016' or day1 == '05012016':
    valid = ~np.isnan(init_sd)
    p1,pcov = curve_fit(h, Bins[valid], init_sd[valid] , p0=[100, 400, 400, 1.5, 1.5, 1.5, 40, 110, 400],maxfev=8000)
  else:
    init_sd = init_sd[54:]
    Bins = Bins[54:]
    valid = ~np.isnan(init_sd)
  #elif day1 == '05192016' or day1 == '05162016' or day1 == '05142016':
    while stat==0:
      #print('Trying to fit:',i,guess1,guess2,guess3)
      #try:
        #print('Trying')
      if day1 == '04272016':
        p1,pcov = curve_fit(h2, Bins[valid], init_sd[valid]*(np.pi*Bins[valid]**2), p0=[guess2, guess1, guess3, 1.5, 1.5, 1.1, 4, 100, 100],maxfev=8000, method = 'lm')
      if day1 == '04282016':
        p1,pcov = curve_fit(h2, Bins[valid], init_sd[valid]*(np.pi*Bins[valid]**2), p0=[guess2, guess1, guess3, 1.5, 1.5, 1.1, 4, 100, 100],maxfev=8000, method = 'lm')
      if day1 == '05112016':
        p1,pcov = curve_fit(h2, Bins[valid], init_sd[valid]*(np.pi*Bins[valid]**2), p0=[guess2, guess1, guess3, 1.5, 1.5, 2.0, 4, 100, 200],maxfev=8000, method = 'lm') 
      if day1 == '05142016':
        p1,pcov = curve_fit(h2, Bins[valid], init_sd[valid]*(np.pi*Bins[valid]**2), p0=[guess2, guess1, guess3, 1.5, 1.5, 1.5, 4, 100, 100], \
            maxfev=8000,bounds = [0,np.nanmax(init_sd)], method = 'trf') #5/14
      if day1 == '05162016':
        p1,pcov = curve_fit(h2, Bins[valid], init_sd[valid]*(np.pi*Bins[valid]**2), p0=[guess2, guess1, guess3, 1.5, 1.5, 2.0, 4, 100, 200],maxfev=8000, method = 'lm') 
      if day1 == '05192016':
        p1,pcov = curve_fit(h2, Bins[valid], init_sd[valid]*(np.pi*Bins[valid]**2), p0=[guess2, guess1, guess3, 1.5, 1.5, 2.0, 4, 100, 200],maxfev=8000, method = 'lm') 
      else: 
        p1,pcov = curve_fit(h2, Bins[valid], init_sd[valid]*(np.pi*Bins[valid]**2), p0=[guess2, guess1, guess3, 1.5, 1.5, 2.0, 4, 100, 200],maxfev=8000, method = 'lm') 

        #p1,pcov = curve_fit(h, Bins[valid], init_sd[valid], p0=[guess2, guess1, 1.5, 1.5, 100, 200],maxfev=8000, method = 'trf')
        #print('p1 =',p1)
      stat = 1
#      #except:
#        print('changing guess')
#        #sys.exit('Excepting')
#        guess1 = guess1+1000
#        if guess1>10000:
#          guess1=1000
#          guess2 = guess2+1000
#        if guess2>10000:
#          guess2=100
#          guess3 = guess3+1000
#        if guess3>10000:
#          pyplot.semilogx(Bins,init_sd)
#          pyplot.show( )
#          sys.exit('guess3 maxed out')
#          stat=1
  #else:
  #  p1,pcov = curve_fit(h, Bins[valid], init_sd[valid] , p0=[100, 400, 400, 1.5, 1.5, 1.5, 5, 40, 400],maxfev=20000)

  No_bg1 = ''
  Dpm_bg1 = ''
  sigma_bg1 = ''
  No_bg2 = ''
  Dpm_bg2 = ''
  sigma_bg2 = ''
  No_bg3 = ''
  Dpm_bg3 = ''
  sigma_bg3 = ''
  
#  print(year,month,day)
#  print('p1 =',p1)
#  print('max of init_sd = ',np.nanmax(init_sd))
#  sys.exit('You suck')

#  No1 = p1[0]
#  No2 = p1[1]
#  No3 = p1[2]
#  sigma1 = p1[3]
#  sigma2 = p1[4]
#  sigma3 = p1[5]
#  mu1 = p1[6]
#  mu2 = p1[7]
#  mu3 = p1[8]
#  No1 = p1[0]
#  No2 = p1[1]
#  sigma1 = p1[2]
#  sigma2 = p1[3]
#  mu1 = p1[4]
#  mu2 = p1[5]
  
  #N = h(Bins,No1,No2,No3,sigma1,sigma2,sigma3,mu1,mu2,mu3)
  #N = h(Bins,No1,No2,sigma1,sigma2,mu1,mu2)
  #pyplot.semilogx(Bins,N,Bins,init_sd)
  #pyplot.semilogx(Bins,N*((100**3)/((10**6)**3)*1600),Bins,init_sd*(np.pi/6*Bins**3)*((100**3)/((10**6)**3)*1600))
  #pyplot.semilogx(Bins,N*(np.pi*Bins**2)*((100**2)/((10**6)**2)),Bins,init_sd*(np.pi*Bins**2)*((100**2)/((10**6)**2)))
  #pyplot.show()
  #sys.exit('You suck')
#  print('day1 =',day1)
#  BG_SD_gc = gc.open('%s_FIMS_scale'%day1)
#  df_bgsd = np.array(BG_SD_gc.worksheet_by_title(title='Sheet1').get_as_df())
#  print('np.shape(df_bgsd) =',np.shape(df_bgsd))
#  lin_alt = np.array(df_bgsd[:,0])

#  if day1=='04272016' or day1=='09172016' or day1=='09112016' or day1=='05112016':
  for i in range(nlayers):
    No_bg1 = No_bg1 + '%15.7E'%g1(p1[0],pres_lay[i],temp_lay[i],stppres,stptemp)
    No_bg2 = No_bg2 + '%15.7E'%g1(p1[1],pres_lay[i],temp_lay[i],stppres,stptemp)
    No_bg3 = No_bg3 + '%15.7E'%g1(p1[2],pres_lay[i],temp_lay[i],stppres,stptemp)
    sigma_bg1 = sigma_bg1 + '%15.7E'%p1[3]
    sigma_bg2 = sigma_bg2 + '%15.7E'%p1[4]
    sigma_bg3 = sigma_bg3 + '%15.7E'%p1[5]
    Dpm_bg1 = Dpm_bg1 + '%15.7E'%g3(p1[6],i,1.0)
    Dpm_bg2 = Dpm_bg2 + '%15.7E'%g3(p1[7],i,1.0)
    Dpm_bg3 = Dpm_bg3 + '%15.7E'%g3(p1[8],i,1.0)
      
#  elif day1=='04282016': # or day1=='05112016': # or day1=='09112016':
#    BG_SD_gc = gc.open('%s_FIMS_scale'%day1)
#    df_bgsd = np.array(BG_SD_gc.worksheet_by_title(title='Sheet1').get_as_df())
#    print('np.shape(df_bgsd) =',np.shape(df_bgsd))
#    lin_alt = np.array(df_bgsd[:,0])
#    for i in range(nlayers):
#      indx = closest(lin_alt,Z_lay[i])
#      temp = np.nanmean(df_bgsd[indx-5:indx+5,1:],axis=0)
#      print('temp =',temp,'for layer:',i+1)
#      No_bg1 = No_bg1 + '%15.7E'%(p1[0]*temp[0])
#      No_bg2 = No_bg2 + '%15.7E'%(p1[1]*temp[1])
#      sigma_bg1 = sigma_bg1 + '%15.7E'%(p1[2]*temp[2])
#      sigma_bg2 = sigma_bg2 + '%15.7E'%(p1[3]*temp[3])
#      Dpm_bg1 = Dpm_bg1 + '%15.7E'%g3(p1[4],i,temp[4])
#      Dpm_bg2 = Dpm_bg2 + '%15.7E'%g3(p1[5],i,temp[5])
#    No_bg1 = No_bg1 + '%15.7E'%g1(p1[0],pres_lay[i],temp_lay[i],stppres,stptemp)
#    No_bg2 = No_bg2 + '%15.7E'%g1(p1[1],pres_lay[i],temp_lay[i],stppres,stptemp)
#    sigma_bg1 = sigma_bg1 + '%15.7E'%p1[2]
#    sigma_bg2 = sigma_bg2 + '%15.7E'%p1[3]
#    Dpm_bg1 = Dpm_bg1 + '%15.7E'%g3(p1[4],i,1.0)
#    Dpm_bg2 = Dpm_bg2 + '%15.7E'%g3(p1[5],i,1.0)
#  print(No_bg1,No_bg2,sigma_bg1,sigma_bg2,Dpm_bg1,Dpm_bg2)
#  sys.exit('Sucks')

  return No_bg1,Dpm_bg1,sigma_bg1,No_bg2,Dpm_bg2,sigma_bg2,No_bg3,Dpm_bg3,sigma_bg3

############################################################################################################

def SD_smps_fims(SMPS_sheetname,SMPS_lower,nlayers,endtime,gc,Z_lay,temp_lay,stppres,stptemp,day1,pres_lay,ML_top):
  if day1 == '0016':
    gsh_smps = gc.open('HiScaleSMPS')
    df_1 = gsh_smps.worksheet_by_title(title='%s'%SMPS_sheetname).get_as_df()
    
    temper = df_1.to_numpy()
    time3 = temper[:,0]
    Bins = temper[:110,-1]
    flag = temper[:,-2]
    sizedist = temper[:,1:-2].astype(float)
    
    year = int(SMPS_sheetname[17:-7])
    month = int(SMPS_sheetname[-7:-5])
    day = int(SMPS_sheetname[-5:-3])
    
    startT = dt.datetime(year,month,day,0,0) - dt.timedelta(hours=5)
    date = []
    
    for i in range(len(time3)):
      date.append(startT + dt.timedelta(seconds=time3[i]))
      if flag[i]>0.:
        sizedist[i].fill('NaN')
    SMPS_upper = SMPS_lower + dt.timedelta(hours=int(endtime[0]))
    
    init_date,init_sd= ([] for i in range(2))
    
    for i in range(len(date)):   
      if date[i] >= SMPS_lower and date[i] <= SMPS_lower+dt.timedelta(minutes=15): 
        init_date.append(date[i])
        init_sd.append((sizedist[i]))
  else:
    gsh_smps = gc.open('SMPS_merged')
    year = SMPS_lower.year
    month = SMPS_lower.month
    day = SMPS_lower.day
  
    df_1 = gsh_smps.worksheet_by_title(title='SMPS_merged_%s'%(str(year)+str(month).zfill(2)+str(day).zfill(2))).get_as_df()
    #df_1 = gsh_smps.worksheet_by_title(title='%s'%SMPS_sheetname).get_as_df()
    
    temper = df_1.to_numpy()
    time3 = temper[:,0]
  #  print('time3 =',time3)
  #  sys.exit('You suck')
  #  flag = temper[:,-2]
    sizedist = temper[:,1:].astype(float)
    
    Bins = []
    bin_data = 'SMPS_merged_dp.txt'
    bin_fid = open(bin_data, 'r')
    for line in bin_fid.readlines():
      #print(line)
      spl_line = line.split('	')
    #print(type(spl_line))
    #print(spl_line)
    for i in range(len(spl_line)):
      Bins.append(float(spl_line[i]))
    Bins = np.array(Bins)
  #  print('Bins =',Bins) 
  #  year = int(SMPS_sheetname[12:15])
  #  month = int(SMPS_sheetname[15:17])
  #  day = int(SMPS_sheetname[17:19])
  #  print(year,month,day)
  #  year,month,day = 2016,5,14
    startT = dt.datetime(year,month,day,0,0) - dt.timedelta(hours=5)
    date = []
    
    for i in range(len(time3)):
  #    date.append(startT + dt.timedelta(seconds=time3[i]))
      date.append(dt.datetime.fromordinal(int(time3[i])) + dt.timedelta(days=time3[i]%1) - dt.timedelta(days = 366) - dt.timedelta(hours=5))
      #date.append(startT + dt.timedelta(seconds=int((time3[i] - time3[1])*24*60*60)))
  #    if flag[i]>0.:
  #      sizedist[i].fill('NaN')
  #  print(date[4]<SMPS_lower+dt.timedelta(minutes=60))
  #  print(date[0:10])
  #  print(SMPS_lower+dt.timedelta(minutes=60))
    #sys.exit('You suck')
    SMPS_upper = SMPS_lower + dt.timedelta(hours=int(endtime[0]))
    #print('SMPS_lower =',SMPS_lower)    
    init_date,init_sd= ([] for i in range(2))
    for i in range(len(date)):   
      if date[i] >= SMPS_lower and date[i] <= SMPS_lower+dt.timedelta(minutes=20):
        init_date.append(date[i])
        init_sd.append((sizedist[i]))

  #  sys.exit('You really suck')
#    if date[i] >= SMPS_upper and date[i] <= SMPS_upper+dt.timedelta(minutes=15): 
#      final_date.append(date[i])
#      final_sd.append((sizedist[i]))
  
  #print('init_sd =',init_sd)
  init_date = np.array(init_date)
  init_sd = np.nanmean(np.array(init_sd),axis=0)
  #print('init_sd =',init_sd)
  valid = ~np.isnan(init_sd)
  #print('valid =',valid)
#  final_date = np.array(final_date)
#  final_dsd = np.array(final_sd)
  guess1 = 1000
  guess2 = 100
  guess3 = 100
  
  if day1 == '04272016':
    p1,pcov = curve_fit(h2, Bins[valid], init_sd[valid]*(np.pi*Bins[valid]**2), p0=[guess2, guess1, guess3, 1.5, 1.5, 1.1, 4, 100, 100],maxfev=8000, method = 'lm')
  if day1 == '04282016':
    p1,pcov = curve_fit(h2, Bins[valid], init_sd[valid]*(np.pi*Bins[valid]**2), p0=[guess2, guess1, guess3, 1.5, 1.5, 1.1, 4, 100, 100],maxfev=8000, method = 'lm')
  if day1 == '05112016':
    p1,pcov = curve_fit(h2, Bins[valid], init_sd[valid]*(np.pi*Bins[valid]**2), p0=[guess2, guess1, guess3, 1.5, 1.5, 2.0, 4, 100, 200],maxfev=8000, method = 'lm') 
  if day1 == '05142016':
    p1,pcov = curve_fit(h2, Bins[valid], init_sd[valid]*(np.pi*Bins[valid]**2), p0=[guess2, guess1, guess3, 1.5, 1.5, 1.5, 4, 100, 100], \
        maxfev=8000,bounds = [0,np.nanmax(init_sd)], method = 'trf') #5/14
  #p1,pcov = curve_fit(h, Bins[valid], init_sd[valid], p0=[100, 400, 400, 1.5, 1.5, 1.5, 5, 40, 110])
  else:
    p1,pcov = curve_fit(h2, Bins[valid], init_sd[valid]*(np.pi*Bins[valid]**2), p0=[guess2, guess1, guess3, 1.5, 1.5, 1.5, 4, 100, 100], \
        maxfev=8000,bounds = [0,np.nanmax(init_sd)], method = 'trf') #5/14
  
  No_bg1 = ''
  Dpm_bg1 = ''
  sigma_bg1 = ''
  No_bg2 = ''
  Dpm_bg2 = ''
  sigma_bg2 = ''
  No_bg3 = ''
  Dpm_bg3 = ''
  sigma_bg3 = ''
  
#  print('day1 =',day1)
#  BG_SD_gc = gc.open('%s_FIMS_scale'%day1)
#  df_bgsd = np.array(BG_SD_gc.worksheet_by_title(title='Sheet1').get_as_df())
#  print('np.shape(df_bgsd) =',np.shape(df_bgsd))
#  lin_alt = np.array(df_bgsd[:,0])
#  print('lin_alt =',lin_alt)
#  sys.exit('Sucker')
  BG_SD_gc = gc.open('%s_sortFIMS'%day1)
  df_bgsd = np.array(BG_SD_gc.worksheet_by_title(title='Sheet1').get_as_df())
  lin_alt = np.array(df_bgsd[:,1])

  #print('alt array:',lin_alt)
  #print('fims array:',df_bgsd[0:5,2:])
  #sys.exit('You suck')

  fims_bins = [ #10.0000, 
      #11.3565, 
      #12.8969, 
      #14.6464, 
      #16.6331, 
      #18.8894, 
      #21.4516, 
      #24.3615, 
      #27.6661, 
      #31.4189, 
      #35.6808, 
      #40.5208, 
      #46.0173, 
      #52.2594, 
      #59.3483, 
      #67.3987, 
      76.5412, 
      86.9238, 
      98.7148, 
      112.1051,
      127.3119,
      144.5814,
      164.1935,
      186.4658,
      211.7594,
      240.4840,
      273.1050,
      310.1509,
      352.2221,
      400.0000]

  for i in range(nlayers):
    if Z_lay[i] <= ML_top:
      No_bg1 = No_bg1 + '%15.7E'%g1(p1[0],pres_lay[i],temp_lay[i],stppres,stptemp)
      No_bg2 = No_bg2 + '%15.7E'%g1(p1[1],pres_lay[i],temp_lay[i],stppres,stptemp)
      No_bg3 = No_bg3 + '%15.7E'%g1(p1[2],pres_lay[i],temp_lay[i],stppres,stptemp)
      sigma_bg1 = sigma_bg1 + '%15.7E'%p1[3]
      sigma_bg2 = sigma_bg2 + '%15.7E'%p1[4]
      sigma_bg3 = sigma_bg3 + '%15.7E'%p1[5]
      Dpm_bg1 = Dpm_bg1 + '%15.7E'%g3(p1[6],i,1.0)
      Dpm_bg2 = Dpm_bg2 + '%15.7E'%g3(p1[7],i,1.0)
      Dpm_bg3 = Dpm_bg3 + '%15.7E'%g3(p1[8],i,1.0)
      
    elif Z_lay[i] > ML_top:
      indx = closest(lin_alt,Z_lay[i])
      #print('indx =',indx)
      #print('Z_lay[i] =',Z_lay[i],'lin_alt[indx] =',lin_alt[indx])
      if Z_lay[i]>60:
        #print('df_bgsd',df_bgsd[indx-10:indx+10,2:])
        temp = np.nanmean(df_bgsd[indx-10:indx+10,18:],axis=0)
      else:
        temp = np.nanmean(df_bgsd[0:indx+15,18:],axis=0)
      #print('shape of temp:',np.shape(temp))
      #valid = ~np.isnan(temp)
      #print('valid =',valid)
      #sys.exit('I hate you')
      #print('fims_bins =',fims_bins)
      #print('fims_bins[valid] =',fims_bins[valid])
      #try:
      stat=0
      guess1 = 1000
      guess2 = 100
      guess3 = 100
      valid = ~np.isnan(temp)
      fims_bins = np.array(fims_bins)
      #print('temp[valid] =',temp[valid])
      #print('len(temp) =',len(temp))
      while stat==0:
        #print('Trying to fit:',i,guess1,guess2,guess3)
        try:
          #print('Trying')
          p1,pcov = curve_fit(h2, fims_bins[valid], temp[valid]*(np.pi*fims_bins[valid]**2), p0=[guess1, guess2, guess3, 1.5, 1.5, 1.5, 4, 40, 200],maxfev=10000)
          #print('p1 =',p1)
          stat = 1
        except:
          print('changing guess',Z_lay[i])
          #sys.exit('Excepting')
          guess1 = guess1+1000
          if guess1>10000:
            guess1=1000
            guess2 = guess2+1000
          if guess2>10000:
            guess2=100
            guess3 = guess3+1000
          if guess3>10000:
            sys.exit('guess3 maxed out')
            stat=1
         #else:
          #print('p1 =',p1)
          #pyplot.semilogx(fims_bins,temp,fims_bins,h(fims_bins,p1[0],p1[1],p1[2],p1[3],p1[4],p1[5],p1[6],p1[7],p1[8]))
          #pyplot.show()
          #stat=1

      #p1,pcov = curve_fit(h, fims_bins, temp, p0=[100,400, 400, 1.5, 1.5, 1.5, 5, 10, 110],maxfev=5000)
      #except:
        #print('was not able to find fit parameters for the following:')
      #pyplot.semilogx(fims_bins,temp)
      #pyplot.show()
      #sys.exit('Also you suck')
      #print('p1 =',p1)
      #print('temp =',temp,'for layer:',i+1)
      No_bg1 = No_bg1 + '%15.7E'%(p1[0])
      No_bg2 = No_bg2 + '%15.7E'%(p1[1])
      No_bg3 = No_bg3 + '%15.7E'%(p1[2])
      sigma_bg1 = sigma_bg1 + '%15.7E'%(p1[3])
      sigma_bg2 = sigma_bg2 + '%15.7E'%(p1[4])
      sigma_bg3 = sigma_bg3 + '%15.7E'%(p1[5])
      Dpm_bg1 = Dpm_bg1 + '%15.7E'%g3(p1[6],i,1.0)
      Dpm_bg2 = Dpm_bg2 + '%15.7E'%g3(p1[7],i,1.0)
      Dpm_bg3 = Dpm_bg3 + '%15.7E'%g3(p1[8],i,1.0)
      #No_bg1 = No_bg1 + '%15.7E'%g1(p1[0],pres_lay[i],temp_lay[i],stppres,stptemp)
      #No_bg2 = No_bg2 + '%15.7E'%g1(p1[1],pres_lay[i],temp_lay[i],stppres,stptemp)
      #sigma_bg1 = sigma_bg1 + '%15.7E'%p1[2]
      #sigma_bg2 = sigma_bg2 + '%15.7E'%p1[3]
      #Dpm_bg1 = Dpm_bg1 + '%15.7E'%g3(p1[4],i,1.0)
      #Dpm_bg2 = Dpm_bg2 + '%15.7E'%g3(p1[5],i,1.0)
    #print(No_bg1,No_bg2,sigma_bg1,sigma_bg2,Dpm_bg1,Dpm_bg2)
  
  #sys.exit('Sucks for you, right?')
  return No_bg1,Dpm_bg1,sigma_bg1,No_bg2,Dpm_bg2,sigma_bg2,No_bg3,Dpm_bg3,sigma_bg3

############################################################################################################

def SD_nanosmps(SMPS_sheetname,SMPS_lower,nlayers,endtime,gc,Z_lay,temp_lay,stppres,stptemp,day1,pres_lay,ML_top):

  gsh_smps = gc.open('SMPS_merged_20160514')
  df_1 = gsh_smps.worksheet_by_title(title='%s'%SMPS_sheetname).get_as_df()
  
  temper = df_1.to_numpy()
  time3 = temper[:,0]
#  print('time3 =',time3)
#  sys.exit('You suck')
#  flag = temper[:,-2]
  sizedist = temper[:,1:].astype(float)
  
  Bins = []
  bin_data = 'SMPS_merged_dp.txt'
  bin_fid = open(bin_data, 'r')
  for line in bin_fid.readlines():
    #print(line)
    spl_line = line.split('	')
  #print(type(spl_line))
  #print(spl_line)
  for i in range(len(spl_line)):
    Bins.append(float(spl_line[i]))
  Bins = np.array(Bins)

#  year = int(SMPS_sheetname[12:15])
#  month = int(SMPS_sheetname[15:17])
#  day = int(SMPS_sheetname[17:19])
#  print(year,month,day)
  year,month,day = 2016,5,14
  startT = dt.datetime(year,month,day,0,0) - dt.timedelta(hours=5)
  date = []
  
  for i in range(len(time3)):
#    date.append(startT + dt.timedelta(seconds=time3[i]))
    date.append(startT + dt.timedelta(seconds=((time3[i] - time3[0])*24*60*60)))
#    if flag[i]>0.:
#      sizedist[i].fill('NaN')
  SMPS_upper = SMPS_lower + dt.timedelta(hours=int(endtime[0]))
  
  init_date,init_sd= ([] for i in range(2))
  
  for i in range(len(date)):   
    if date[i] >= SMPS_lower and date[i] <= SMPS_lower+dt.timedelta(minutes=15): 
      init_date.append(date[i])
      init_sd.append((sizedist[i]))
#    if date[i] >= SMPS_upper and date[i] <= SMPS_upper+dt.timedelta(minutes=15): 
#      final_date.append(date[i])
#      final_sd.append((sizedist[i]))
  
  init_date = np.array(init_date)
  init_sd = np.nanmean(np.array(init_sd),axis=0)
  valid = ~np.isnan(init_sd)
  print(init_sd)
#  final_date = np.array(final_date)
#  final_dsd = np.array(final_sd)
  
  p1,pcov = curve_fit(h, Bins[valid], init_sd[valid], p0=[400, 400, 1.5, 1.5, 40, 110])
  
  
  No_bg1 = ''
  Dpm_bg1 = ''
  sigma_bg1 = ''
  No_bg2 = ''
  Dpm_bg2 = ''
  sigma_bg2 = ''
  
#  print('day1 =',day1)
#  BG_SD_gc = gc.open('%s_FIMS_scale'%day1)
#  df_bgsd = np.array(BG_SD_gc.worksheet_by_title(title='Sheet1').get_as_df())
#  print('np.shape(df_bgsd) =',np.shape(df_bgsd))
#  lin_alt = np.array(df_bgsd[:,0])
#  print('lin_alt =',lin_alt)
#  sys.exit('Sucker')
  BG_SD_gc = gc.open('%s_SD_fit_params'%day1)
  df_bgsd = np.array(BG_SD_gc.worksheet_by_title(title='Sheet1').get_as_df())
  lin_alt = np.array(df_bgsd[:,0])

  for i in range(nlayers):
    if Z_lay[i] <= ML_top:
      No_bg1 = No_bg1 + '%15.7E'%g1(p1[0],pres_lay[i],temp_lay[i],stppres,stptemp)
      No_bg2 = No_bg2 + '%15.7E'%g1(p1[1],pres_lay[i],temp_lay[i],stppres,stptemp)
      sigma_bg1 = sigma_bg1 + '%15.7E'%p1[2]
      sigma_bg2 = sigma_bg2 + '%15.7E'%p1[3]
      Dpm_bg1 = Dpm_bg1 + '%15.7E'%g3(p1[4],i,1.0)
      Dpm_bg2 = Dpm_bg2 + '%15.7E'%g3(p1[5],i,1.0)
      
    elif Z_lay[i] > ML_top:
#      for j in range(nlayers):
        indx = closest(lin_alt,Z_lay[i])
        temp = np.nanmean(df_bgsd[indx-10:indx+10,2:],axis=0)
        #print('temp =',temp,'for layer:',i+1)
        No_bg1 = No_bg1 + '%15.7E'%(temp[0])
        No_bg2 = No_bg2 + '%15.7E'%(temp[1])
        sigma_bg1 = sigma_bg1 + '%15.7E'%(temp[2])
        sigma_bg2 = sigma_bg2 + '%15.7E'%(temp[3])
        Dpm_bg1 = Dpm_bg1 + '%15.7E'%g3(temp[4],i,1.0)
        Dpm_bg2 = Dpm_bg2 + '%15.7E'%g3(temp[5],i,1.0)
      #No_bg1 = No_bg1 + '%15.7E'%g1(p1[0],pres_lay[i],temp_lay[i],stppres,stptemp)
      #No_bg2 = No_bg2 + '%15.7E'%g1(p1[1],pres_lay[i],temp_lay[i],stppres,stptemp)
      #sigma_bg1 = sigma_bg1 + '%15.7E'%p1[2]
      #sigma_bg2 = sigma_bg2 + '%15.7E'%p1[3]
      #Dpm_bg1 = Dpm_bg1 + '%15.7E'%g3(p1[4],i,1.0)
      #Dpm_bg2 = Dpm_bg2 + '%15.7E'%g3(p1[5],i,1.0)
    #print(No_bg1,No_bg2,sigma_bg1,sigma_bg2,Dpm_bg1,Dpm_bg2)
  
  #sys.exit('Sucks for you, right?')
  return No_bg1,Dpm_bg1,sigma_bg1,No_bg2,Dpm_bg2,sigma_bg2
