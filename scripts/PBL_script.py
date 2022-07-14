# This program is intended to get 30-minute resolution PBL height 
# data from MERRA2 reanalysis data used in GEOS-Chem

#from pylab import * 
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt

def pbl(OH_lower,endtime,PBL_time,delt):

  pbldata = '../scripts/MERRA2_PBL_data'
  Date = []
  PBL_tmp = []
  
  fid1 = open(pbldata,'r')
  count = 0
  
  for line in fid1.readlines():
    spl_line = line.split(' ')
    Date.append(dt.datetime.strptime(str(spl_line[0]+'-'+spl_line[1]),'%Y-%m-%d-%H:%M:%S' )-dt.timedelta(hours=5))
    PBL_tmp.append(float(spl_line[2]))
  
  Date = np.array(Date)
  upper = OH_lower + dt.timedelta(hours = endtime[0])
  
  low = np.where(Date>OH_lower)[0][0]
  up = np.where(Date>upper)[0][0]
  
  PBL_cut = PBL_tmp[low:up]
  x1 = np.linspace(0,len(PBL_cut)-1,int(endtime[0]*60*60/delt))
  x2 = np.arange(0,len(PBL_cut))
  PBL = np.interp(x1,x2,PBL_cut)
  
  lower2 = dt.datetime(OH_lower.year,OH_lower.month,OH_lower.day) - dt.timedelta(days=1)
  upper2 = lower2 + dt.timedelta(days=1)

  low2 = np.where(Date>lower2)[0][0]
  up2 = np.where(Date>upper2)[0][0]
  RL_top = PBL_tmp[low2:up2]
  RL_top = np.nanmax(RL_top)

  low3 = np.where(Date>PBL_time)[0][0]
  ML_top = PBL_tmp[low3]

  return PBL,RL_top,ML_top

def pbl2(OH_lower,endtime,PBL_time,delt):
  
  
  data = 'sgppblhtsondeyr1mcfarlC1.s1.20160101.052700.csv'
  fid = open(data, 'r')

  time=[] #in UTC , seconds past midnight (2016-09-11)
  sizedist=[]#dN/dlogDm, shape time by 110 bins
  count = 0

  for line in fid.readlines():
    if count == 0:
      count=count + 1
    elif count>=1:
      spl_line=line.split(',')
      time.append(float(spl_line[3]))
      count=count+1
      temp=[]
      for i in range(10):
        if spl_line[i+10] == '':
          spl_line[i+10] = 'NaN'
        spl_line[i+10] = float(spl_line[i+10])
        temp.append(float(spl_line[i+10]))
      sizedist.append(temp)
  sizedist = np.transpose(sizedist)

  startT = dt.datetime(2016,1,1,0,0) - dt.timedelta(hours=5)
  date = []

  for i in range(len(time)):
    date.append(startT + dt.timedelta(seconds=time[i]))
  
  date = np.array(date)
  upper = OH_lower + dt.timedelta(hours = endtime[0])
  
  lower1 = dt.datetime(OH_lower.year,OH_lower.month,OH_lower.day)
  upper1 = lower1 + dt.timedelta(days=1)

  low = np.where(date>lower1)[0][0]
  up = np.where(date>upper1)[0][0]+1
#  print('low =',low,'up =',up) 
  pbl_liu_liang = sizedist[4]-318.
  
  PBL_cut = pbl_liu_liang[low:up]
#  print('PBL_cut =',PBL_cut)
  time_cut = time[low:up]
#  print('time_cut =',time_cut)

  x1 = np.linspace(time_cut[0],time_cut[-1],10000)
#  x1 = np.linspace(0,len(PBL_cut)-1,int(endtime[0]*60*60/delt))
#  x2 = np.arange(0,len(PBL_cut))
  PBL = np.interp(x1,time_cut,PBL_cut)
  #print('x1[0:50]=',x1[0:50])
  #print('x1[-20:]=',x1[-20:])
  #print('PBL[0:50] =',PBL[0:50])
  #print('PBL[-20:] =',PBL[-20:])
  date = []
  for i in range(len(x1)):
    date.append(startT + dt.timedelta(seconds=x1[i]))
  date = np.array(date)
  
#  print('date[0:10]',date[0:10])
#  print('date[-10:]',date[-10:])
  low = np.where(date>OH_lower)[0][0]
  up = np.where(date>upper)[0][0]
  
#  print('low =',low,'up =',up)

  PBL_cut = PBL[low:up]
  #print('PBL[0:50] =',PBL[0:50])
  #print('PBL[-20:] =',PBL[-20:])
  time_cut = x1[low:up]
  x1 = np.linspace(0,len(PBL_cut)-1,int(endtime[0]*60*60/delt)) 
  x2 = np.arange(0,len(PBL_cut))
  PBL = np.interp(x1,x2,PBL_cut)

#  lower2 = dt.datetime(OH_lower.year,OH_lower.month,OH_lower.day) - dt.timedelta(days=1)
#  upper2 = lower2 + dt.timedelta(days=1)
#
#  low2 = np.where(Date>lower2)[0][0]
#  up2 = np.where(Date>upper2)[0][0]
#  RL_top = PBL_tmp[low2:up2]
#  RL_top = np.nanmax(RL_top)
#
#  low3 = np.where(Date>PBL_time)[0][0]
#  ML_top = PBL_tmp[low3]

  return PBL
