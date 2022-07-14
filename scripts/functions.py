import numpy as np


def h(x,N1,N2,N3,sig1,sig2,sig3,m1,m2,m3): return \
        (((2.303*N2)/(((2*np.pi)**.5)*np.log(sig2)))*np.exp(-((np.log(x)-np.log(m2))**2)/(2*np.log(sig2)**2)) +\
        ((2.303*N1)/(((2*np.pi)**.5)*np.log(sig1)))*np.exp(-((np.log(x)-np.log(m1))**2)/(2*np.log(sig1)**2))  +\
        ((2.303*N3)/(((2*np.pi)**.5)*np.log(sig3)))*np.exp(-((np.log(x)-np.log(m3))**2)/(2*np.log(sig3)**2))) #log10(Dh/Dl)
#def h(x,N1,N2,sig1,sig2,m1,m2): return \
#        (((2.303*N2)/(((2*np.pi)**.5)*np.log(sig2)))*np.exp(-((np.log(x)-np.log(m2))**2)/(2*np.log(sig2)**2)) +\
#        ((2.303*N1)/(((2*np.pi)**.5)*np.log(sig1)))*np.exp(-((np.log(x)-np.log(m1))**2)/(2*np.log(sig1)**2)))*np.pi/6*x**3
def h2(x,N1,N2,N3,sig1,sig2,sig3,m1,m2,m3): return \
        (((2.303*N2)/(((2*np.pi)**.5)*np.log(sig2)))*np.exp(-((np.log(x)-np.log(m2))**2)/(2*np.log(sig2)**2)) +\
        ((2.303*N1)/(((2*np.pi)**.5)*np.log(sig1)))*np.exp(-((np.log(x)-np.log(m1))**2)/(2*np.log(sig1)**2))  +\
        ((2.303*N3)/(((2*np.pi)**.5)*np.log(sig3)))*np.exp(-((np.log(x)-np.log(m3))**2)/(2*np.log(sig3)**2)))*np.pi*x**2 #log10(Dh/Dl)

# ====================================================================================================
def closest(lst, K): 
  lst = np.asarray(lst) 
  idx = (np.abs(lst - K)).argmin() 
  return idx
# ====================================================================================================

def g1(x,pre,tem,stpP,stpT): return x/(tem*stpP/(stpT*pre))   # SamO added this function for a vertical profile of background size distribution number concentration 
def g2(x,i1,i11): return np.abs(x*i1+i11)     # SamO added this function for a vertical profile of background size distribution sigma 
def g3(x,i1,i2): return i2*x*1e-3 # SamO added this function for a vertical profile of background size distribution median diametre
def ml_rh(ml_hgt,z_hght,sfcrh,toprh): return z_hght*(toprh-sfcrh)/(ml_hgt-0.0) + sfcrh

def rlft_rh(ml_hgt,top_hgt,z_hght,toprh,minrh,b):
  temp = z_hght*(minrh-toprh)/(top_hgt-ml_hgt)+b
  if temp<minrh:
    temp = minrh
  return temp
# ====================================================================================================

################################### Smoothing function ##################################################
def smoothdist(dist):
  sdist=np.empty(dist.shape)
  sdist[0]=dist[0]*0.75 + dist[1]*0.25
  sdist[1:-1]=dist[:-2]*0.25 + dist[1:-1]*0.5 + dist[2:]*0.25
  sdist[-1]=dist[-2]*0.25 + dist[-1]*0.75
  return sdist
#########################################################################################################
