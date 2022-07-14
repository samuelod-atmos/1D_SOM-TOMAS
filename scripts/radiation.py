import numpy as np
import pygsheets 
import datetime as dt
import pandas as pd 

#=====================================================================================================
# This requires two days worth of data since raw data are in UTC, so there can be carry-over from day to day. 
def rad(OH_sheetname1,OH_sheetname2,endtime,gc,OH_lower,delt):
  gsh_rad = gc.open('Radiation_data')
  df_1 = gsh_rad.worksheet_by_title(title='%s'%OH_sheetname1).get_as_df()
  df_2 = gsh_rad.worksheet_by_title(title='%s'%OH_sheetname2).get_as_df()
  
  index  = np.where(df_1.iloc[:,0]!='')[0].shape[0]
  sizedist1 = df_1['short_direct_normal'].iloc[:index].to_numpy()
  time1 = np.array(df_1['time_offset'].iloc[:index])
  
  index  = np.where(df_2.iloc[:,0]!='')[0].shape[0]
  sizedist2 = df_2['short_direct_normal'].iloc[:index].to_numpy()
  time2 = np.array(df_2['time_offset'].iloc[:index])
  
  sizedist1 = pd.to_numeric(sizedist1)
  sizedist2 = pd.to_numeric(sizedist2)
  
  year1 = int(OH_sheetname1[20:24])
  month1 = int(OH_sheetname1[24:26])
  day1 = int(OH_sheetname1[26:28])
  print('Running day: %s/%s/%s' %(year1,month1, day1))
  date1 = []
  
  for i in range(len(time1)):
    date1.append(dt.datetime(year1,month1,day1,int(time1[i][-8:-6]),int(time1[i][-5:-3]),int(time1[i][-2:-1])))
  date1 = np.array(date1)
  date1 = date1 - dt.timedelta(hours=5) # Convert from UTC to CDT 
  
  
  year2 = int(OH_sheetname2[20:24])
  month2 = int(OH_sheetname2[24:26])
  day2 = int(OH_sheetname2[26:28])
  #print(year2,month2,day2)
  date2 = []
  
  for i in range(len(time2)):
    date2.append(dt.datetime(year2,month2,day2,int(time2[i][-8:-6]),int(time2[i][-5:-3]),int(time2[i][-2:-1])))
  date2 = np.array(date2)
  date2 = date2 - dt.timedelta(hours=5) # Convert from UTC to CDT 
  
  date = np.append(date1,date2)
  sizedist = np.append(sizedist1,sizedist2)
  
  upper = OH_lower + dt.timedelta(hours = endtime[0])
  low = np.where(date>OH_lower)[0][0]
  up = np.where(date>upper)[0][0]
  SW_cut = sizedist[low:up]
  date_cut = date[low:up]
  
  x1 = np.linspace(0,len(SW_cut),int(endtime[0]*60*60/delt))
  x2 = np.arange(0,len(SW_cut))
  DWSW = np.interp(x1,x2,SW_cut)

  return DWSW
