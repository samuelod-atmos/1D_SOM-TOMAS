''' ==================================================================================
This code will extract the data from SOM_TOMAS (COSO) and calculate the yield and SOA
production

Written by Ali Akherati at CSU Oct 2018
================================================================================== '''

# Libraries
# ====================================================================================
import numpy as np
import pandas as pd
import pygsheets as pyg
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1 import make_axes_locatable

# Changing font
# ====================================================================================
#mpl.rcParams['font.family'] = 'Helvetica'
#params = {'mathtext.default': 'regular' }
#plt.rcParams.update(params)

from variable_list import *

# Google sheets
# ====================================================================================
gsh_name = 'comparison_coso'
sheet_number = 1
# authorization
gc = pyg.authorize(service_file='../inputs/SOM-TOMAS-905f18e05abd.json')
# open the google spreadsheet
gsh = gc.open(gsh_name)
# select the first sheet
wks_out = gsh[sheet_number-1]

# Data frame
df_out = pd.DataFrame()

df_out[0] = ['precursor name',
             'regime',
             'Number of Bins',
             'Number of Organics',
             'Coagulation',
             'Vapor Wall Losses',
             'Particle Wall Losses',
             'time [h]',
             'Initial total number concentration [#/cm3]',
             'Initial Median Diameter [um]',
             'sigma',
             'surface area [µm2/cm3]',
             'OH concentration [molecules/cm3]',
             'precursor concentration [ppm]',
             'ke',
             'kw0',
             'Alpha - Accom. Coef.',
             'Surface Tension [N/m]',
             'Db [m2/s]',
             'kc [1/s]',
             'Pressure [Pa]',
             'Temperature [K]',
             'RH',
             'Boxvol [cm3]',
             'SOM parameters',
             'mfrag',
             'dlvp',
             'PF1',
             'PF2',
             'PF3',
             'PF4',
             'kOH',
             'Result',
             'VOC_reacted',
             'SOM_TOMAS SOA Conc. [ug/m3]',
             'Yield',
             'O:C']

wks_out.set_dataframe(df_out,(1,1), copy_head=None)

# lognormal calculation for every diameter
# ====================================================================================
def lognorm_Dist(dp,Nt,dp_avg,sgma):
    N_dp = (Nt/(((2.*np.pi)**.5)*dp*np.log(sgma)))*\
           np.exp(-1.0*((np.log(dp)-np.log(dp_avg))**2.)\
                  /(2.*(np.log(sgma)**2.)))
    return N_dp
    
# inputs
# ====================================================================================
rho = 1400. # [kg/m3]
R = 8.314 # [j/mol/K]

# diameter
# ====================================================================================
print('number of bins = %s'%ibins)
Mo = 1.0e-21*2.e0**(-6) # [kg/particle]
dp = (Mo/rho*6/np.pi)**(1./3.) # [m]
#dp = dp*1.e9 # [nm]
print(dp, '[m]')

xk = []
for i in range(ibins+1):
    xk.append(Mo*(2**i))
xk=np.array(xk)

dpi = (xk/rho*6/np.pi)**(1./3.)*1.e9 # [nm]
print(dpi, '[nm]')
print('min dpi = %s'%dpi[0])
print('max dpi = %s'%dpi[-1])
tsurf=0.0

Dpv = (dpi[1:]*dpi[:-1])**.5

'''
# dN/dDp Calculation
# ====================================================================================
dNdDp = lognorm_Dist(dpi[:-1],No,Dpm,sigma)
dNdDp = np.array(dNdDp)
dNdlogDp = dNdlnDp*2.303

Dpv = (dpi[1:]*dpi[:-1])**.5
dMdlopDp = (np.pi*((Dpv*1e-3)**2)*dNdlogDp)
tsurf = np.nansum(dMdlopDp*np.log10(dpi[1]/dpi[0]))
dVdlogDp = (np.pi/6.*(Dpv**3)*dNdlogDp)
'''
# process
# ====================================================================================
src_path = '../src'
output_path = '../outputs'
ctr = 0
for vwl in VWL:
 for pwl in PWL:
  for ind_precn, precn in enumerate([precursor]):
   for ind_reg, reg in enumerate([regime]):
    for ind_tend, tend in enumerate(endtime):
     for ind_no, no in enumerate(No):
      for ind_dpm, dp_val in enumerate(Dpm):
       for ind_oh, oh in enumerate(OH_conc):
        for ind_db, db in enumerate(Dbk):
         for ind_kc, KC in enumerate(kc):
          for ind_ppm, ppm in enumerate(ippmprec):
            # dN/dDp Calculation
            # ====================================================================================
            dNdDp = lognorm_Dist(dpi[:-1],no,dp_val*1000.,sigma)
            dNdDp = np.array(dNdDp)
            dNdlnDp = dNdDp*dpi[:-1]
            dNdlogDp = dNdlnDp*2.303
            
            dMdlopDp = (np.pi*((Dpv*1e-3)**2)*dNdlogDp)
            tsurf = np.nansum(dMdlopDp*np.log10(dpi[1]/dpi[0]))
            dVdlogDp = (np.pi/6.*(Dpv**3)*dNdlogDp)

            ctr = ctr+1
            #rname = precn+'_'+reg+'_'+tend+'_'+no+'_'+dp_val+'_'+oh+'_'+ppm
            rname = '%s_%s_%04.1f_%08.1f_%05.3f_%7.1E_%5.3f'%(precn,reg,tend,no,dp_val,oh,ppm)
            if pwl==1:
                rname = 'PWL_'+rname
            if vwl==1:
                rname = 'VWL_'+rname
            print(rname, 'completed!')
          
            # Number Conc.
            # ------------------------------------------------------------------------------------
            df = np.array(pd.read_csv('%s/%s_noconc.dat'%(output_path,rname),
                                      header=None, delim_whitespace=True))
            Nk = df[:,1:] # [#/box]
            Ntot = np.sum(Nk, axis=1)
            Nk = Nk/boxvol # [#/cm3]
            # Need to conver to dN/dlogDp
            for i in range(Nk.shape[1]):
                Nk[:,i] = Nk[:,i]/np.log10(dpi[i+1]/dpi[i])
              
            #Ntot = []
            #for i in range(ibins):
            #    Ntot.append(np.log10(dpi[i+1]/dpi[i])*Nk[:,i])
            #Ntot = np.array(Ntot)

            # time
            # ------------------------------------------------------------------------------------
            time = df[:,0]
            #time = np.arange(0,time.shape[0])/3600.*10.

            # Aerosol Mass Conc.
            # ------------------------------------------------------------------------------------
            df = np.array(pd.read_csv('%s/%s_aemass.dat'%(output_path,rname),
                                      header=None, delim_whitespace=True))
            Mk = df[:,1:]
            Mk = np.reshape(Mk, (time.shape[0], int(Mk.shape[0]/time.shape[0]), Mk.shape[1]))
            Mkorg = Mk[:,1:iorg+1,:]
          
            # Gas Conc.
            # ------------------------------------------------------------------------------------
            df = np.array(pd.read_csv('%s/%s_gc.dat'%(output_path,rname),
                                      header=None, delim_whitespace=True))
            Gc = df[:,1:].astype(float)
            Gcorg = Gc[:,1:iorg+1]
            
            # SAPRC Gas Conc.
            # ------------------------------------------------------------------------------------
            df = np.array(pd.read_csv('%s/%s_saprcgc.dat'%(output_path,rname),
                                      header=None, delim_whitespace=True))
            saprcgc = df[:,:].astype(float)
            
            prec_name = 'GENVOC'
            saprc_doc = pd.read_csv('%s/saprc14_rev1.doc'%(src_path), header=None,
                                    skiprows=38, nrows=93+iorg, delim_whitespace=True)
            ind = np.where(np.array(saprc_doc)==prec_name)[0][0]
            prec_conc = saprcgc[:, ind]
            prec_mw = saprc_doc.loc[ind, 2]
            ppm_dvoc = (prec_conc[0] - prec_conc[-1])*1e3
            dvoc = (prec_conc[0] - prec_conc[-1])*1e-6* prec_mw*1e6*pres/R/temp # [ug/m3]

            # Specification (gas names, som names, carbon number, oxygen number)
            # ------------------------------------------------------------------------------------
            spec = pd.read_csv('%s/%s_spec.dat'%(output_path,rname),
                               header=None, delim_whitespace=True)
            act_sp = np.array(spec.loc[0,3:])
            CNO = np.array(spec.loc[2,3:iorg+2]).astype(int)
            ONO = np.array(spec.loc[3,3:iorg+2]).astype(int)
            mwt = np.array(spec.loc[4,3:iorg+2]).astype(float)
            mworg = np.array(spec.loc[5,3:iorg+2]).astype(float)
            mworg_hc = np.array(spec.loc[6,3:iorg+2]).astype(float)
            cstar = np.array(spec.loc[7,1:iorg]).astype(float)

            # PROCESSING
            # ====================================================================================
            Gcorgppm = np.zeros(Gcorg.shape)
            Mkorgppm = np.sum(Mkorg, axis=2)
            for i in range(len(time)):
                Gcorgppm[i,:] = Gcorg[i,:]*1.e9/pres*R*temp/mwt[:]*1.e3/boxvol*1.e6 # kg/bag -> ppm
                Mkorgppm[i,:] = Mkorgppm[i,:]*1.e3/pres*R*temp/mworg[:]/(boxvol*1.e-6)*1.e6 # kg/bag -> ppm
                Gcorg[i,:] = Gcorg[i,:]*1.e9/boxvol*1.e6 # kg/bag -> ug/m3
                for j in range(ibins):
                    Mkorg[i,:,j] = Mkorg[i,:,j]/boxvol*1.0e6*1.0e9 # kg/bag -> ug/m3
                
            mass_org = np.sum(Mkorg[1:,:,:], axis=2)
            mass_org = mass_org + Gcorg[1:]
            ppm_org = np.sum((Mkorgppm[:,:]+Gcorgppm[:,:]), axis=1)
            ppmc_org = np.sum((Mkorgppm[:,:]+Gcorgppm[:,:])*CNO[:80], axis=1)
            ppmc_prec = prec_conc*12

            #O2C = (np.sum((Mkorgppm*ONO + Gcorgppm*ONO)[1:,:75], axis=1))/\
            #     (np.sum((Mkorgppm*CNO + Gcorgppm*CNO)[1:,:75], axis=1))
            O2C = (np.sum((Mkorgppm*ONO)[1:,:], axis=1))/\
                  (np.sum((Mkorgppm*CNO)[1:,:], axis=1))

            # PRINTING THE RESULT
            # ====================================================================================
            soa = np.sum(Mkorg[1:,1:iorg+1,:], axis=(1,2))[-1]-\
                  np.sum(Mkorg[1:,1:iorg+1,:], axis=(1,2))[0]
            yld = soa/dvoc
            vcr = dvoc
            print('VOC Reacted =', dvoc, '[ug/m3] or', ppm_dvoc, '[ppb]')
            print('SOA Concentration =', soa, '[ug/m3]')
            print('Yield =', soa/dvoc)
            print('O:C =', O2C[-1])
            print('KOH = %s'%koh)
            
            df_out[ctr] = [precn,
                           reg,
                           ibins,
                           iorg,
                           COAG,
                           vwl,
                           pwl,
                           tend,
                           no,
                           dp_val,
                           sigma,
                           tsurf,
                           oh,
                           ppm,
                           ke,
                           kw0,
                           alpha,
                           storg,
                           db,
                           KC,
                           pres,
                           temp,
                           rh,
                           boxvol,
                           '',
                           cfrag,
                           dlvp,
                           pf1,
                           pf2,
                           pf3,
                           pf4,
                           koh,
                           ' ',
                           vcr,
                           soa,
                           yld,
                           O2C[-1]]
            
wks_out.set_dataframe(df_out,(1,1), copy_head=None)

#df_out.to_excel('%s_%s.xlsx'%(precname, regime), header=None, index=None)

