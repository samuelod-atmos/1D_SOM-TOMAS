import numpy as np
import pygsheets 
import datetime as dt
import pandas as pd 
from functions import closest
import sys
from radiation import rad


def OH_profile(OH_lower,endtime,nlayers,aadt,eta_lay,day):
  ohfid = open('SGP_OH_Geos025x03125','r')
  etafid = open('47_layer_eta','r')
  
  Date,OH_full,eta_levs = [],[],[]
  for line in ohfid.readlines():
    spl_line = line.split(' ')
    OH_full.append(spl_line[2:-2:2])
    Date.append(dt.datetime.strptime(str(spl_line[0]+' '+spl_line[1]),'%Y-%m-%d %H:%M:%S')-dt.timedelta(hours=5))

  Date = np.array(Date)
  OH_full = np.array(OH_full)

  for line in etafid.readlines():
    spl_line = line.split(' ')
    eta_levs.append(float(spl_line[0]))
  
#  print('eta_levs =',eta_levs)
#  sys.exit('temporary')

  index = []
  for i in range(nlayers):
    index.append(closest(eta_levs,eta_lay[i]))
  
#  sys.exit('temporary')
  upper = OH_lower + dt.timedelta(hours=endtime[0])
  low = np.where(Date>OH_lower)[0][0]
  up = np.where(Date>upper)[0][0]
#  print(np.shape(OH_full[low:up]))
#  sys.exit('temporary')
  OHcut = np.empty([len(OH_full[low:up,0]),nlayers])
  for i in range(nlayers):
    OHcut[:,i] = OH_full[low:up,index[i]]

  x1 = np.linspace(0,len(OHcut[:,0])-1,int(endtime[0]*60*60/aadt)+1)
  x2 = np.arange(0,len(OHcut[:,0]))
  OH_interp = np.zeros([nlayers,len(x1)])
  for i in range(nlayers):
    OH_interp[i] = np.interp(x1,x2,OHcut[:,i])
   
  fil = open('../inputs/%s_OH'%day,'w')
  for i in range(len(OH_interp[0])):
    writer = ''
    for j in range(nlayers):
      writer = writer + '%6.2f '%OH_interp[j,i]
    fil.write('%s \n'%writer)
  fil.close()

  return 

