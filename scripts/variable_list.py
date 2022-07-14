''' =================================================================================================
This is the variable list of inputs for SOM-TOMAS

Written by Ali Akherati MECH-E CSU Nov. 2018
================================================================================================= '''
#from google_som_parameters import rd_ggl_sht
import numpy as np
import pandas as pd
import pygsheets as pyg

# TOMAS inputs
# ===================================================================================================
fireID = 7
regime = 'highNOx' # [lowNOx, highNOx]
params = 'cappa2015'
ext    = 'Frn_PhnlT_UNSalk_oxy_OHprf_normal'
#mod    = ['MOD1']

aadt = 10.  # microphysics timestep [s]
nintern = 30 # frequency of printing values or number of internal microphysics steps
# adt = aadt*nintern # timestep for writing output [seconds]  

COAG = 1    # [0 or 1] the switch for On/Off coagulation
VWL  = [1]  # [0 or 1] the switch for On/Off vapor wall loss
PWL  = [0]  # [0 or 1] the switch for On/Off particle wall loss

ke = [0.0, 0.1, 0.156, 0.24]#, 2.5e-4, 1.4e-3] #, 1.0e-5, 1.0e-4, 2.5e-4, 5.0e-4, 1.0e-3] # loss rate constant
kw0     = 0.0        # loss rate constant

alpha   = 1.0        # accommodation coefficient 
storg   = 0.025      # [N/m] surface tension
Dbk     = [1.0E-10, 1.0E-15, 1.0E-18, 1.0E-20, 1.0E-22]  # particle-phase diffusion coefficient [m2/s]
kc      = [0.0E0]    # first-order loss rate of species in the particle phase [1/s]
nucrate = 0.0        # [cm-3 s-1] nucleation rate ***** We do not use it now
dilt    = 0.0        # dilution rate ***** Ali said it is not working
boxvol  = 10000000.0 # teflon [cm3] - [CalTech 24 m3, CSU 10 m3, CMU ?? m3]
endtime = [8.0]      # [hours] hours of run time (make even hours)

pres    = 88100.0    # [Pa] pressure
temp    = 298.0      # [K] temperature
rh      = 0.2        # relative humidity

# --- read input sheet -------------------------------------------------------------- 
gc = pyg.authorize(service_file='../inputs/SOM-TOMAS-905f18e05abd.json')
# open the google spreadsheet
gsh = gc.open('FIREX_fire_inputs')
# select the first sheet
wks_inp = gsh.worksheet_by_title(title='fire%03i'%(fireID))
df_inp = wks_inp.get_as_df()
inputs = np.array(df_inp.iloc[:23,:2])

No1     = inputs[17,1] # [# cm-3] background 1st number conc.
Dpm1    = inputs[18,1] # [microns] background 1st median diameter
sigma1  = inputs[19,1] # 1st background sigma 
No2     = inputs[20,1] # [# cm-3] background 2nd number conc.
Dpm2    = inputs[21,1] # [microns] background 2nd median diameter
sigma2  = inputs[22,1]    # 2nd background sigma 

# the following values are for OH concentration equation
# OH_conc = a_oh*exp(-1.*ax_oh*t) + b_oh*exp(-1.*bx_oh*t) - t in hour and unit is [Molecules/cm3]
lights_on = inputs[7,1]   # [s] of the day
a_oh      = inputs[10,1]
ax_oh     = inputs[11,1]
b_oh      = inputs[12,1]
bx_oh     = inputs[13,1]
OH0       = a_oh*np.exp(-1.*ax_oh*0.) + b_oh*np.exp(-1.*bx_oh*0.)
# 1.5e7*exp(-1.*t)+ 5e6*exp(-0.05*t) t in hour

nsomprec    = 9 # number of som precursors for parameterizations or SOM grids ***** 1 for now
somprecname = 'SVOSOMG UNSSOMG BNZSOMG AR1SOMG AR2SOMG FRNSOMG PHLSOMG ISPSOMG TRPSOMG' # AR1SOMG AR2SOMG' # som precursor parameterization's name
dlvp        = '   1.47    1.42    1.53    1.42    1.46   1.565   1.392  1.9302    1.91'

seed_dens = inputs[5,1]*1000.0
poa_1stname = 'SVO_C10'
poa_1st_lenname = len(poa_1stname)

iorg  = 796
ibins = 36   # number of bins
idiag = 2    # number of diagnostic species including sulfate and ammonia (2)
icomp = iorg + idiag + 1 # total number of gas phase species in TOMAS

# read the google sheet for emission
gc = pyg.authorize(service_file='../inputs/SOM-TOMAS-905f18e05abd.json')
# open the google spreadsheet
gsh = gc.open('FIREX_emission_chamber_Furan_PhenolTriple_UNSParo1')
# select the first sheet
wks = gsh.worksheet_by_title(title='fire%03i'%(fireID))
df_emiss = wks.get_as_df()

index  = np.where(df_emiss.iloc[:,0]!='')[0].shape[0] 
spname = df_emiss['Species Name'].iloc[:index]
g_ippm = df_emiss['Gas phase concentration - %6s [ppm]'%regime].iloc[:index]
p_frac = df_emiss['Initial Particle phase fraction - %6s'%regime].iloc[:index]
nspemiss = len(spname)
emiss_spname = ''
emiss_ippm   = ''
seed_frac    = ''
for i in range(len(spname)):
    emiss_spname = emiss_spname +   '%15s'%(spname.iloc[i])
    emiss_ippm   =   emiss_ippm + '%15.5E'%(g_ippm.iloc[i])
    seed_frac    =    seed_frac + '%15.5E'%(p_frac.iloc[i])
