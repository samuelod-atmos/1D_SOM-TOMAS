''' ==================================================================================================
This code creates the SPARC11b.rxn file to add SOM precursor so that it runs prep.exe and creates
saprc14_rev1.mod and SAPRC14_rev1.f files

written by Ali Akherati October 2018 Colorado State University
================================================================================================== '''

import numpy as np
import pandas as pd
import pygsheets as pyg
#import subprocess as sp
import time
import re
import os
import sys

# Clean ipython
# ====================================================================================================
#from IPython import get_ipython
#get_ipython().magic('reset -sf') 

startTime = time.time()

# Run the python script to run the actual model
# ====================================================================================================
#njobs = 0 # what job you're on
#maxjobs = 100 # how many jobs you can run. For ozone, max is 16
queue = 'defaultfaculty.q@@students'
#queue = 'jathar.q'

# INPUTS
# ====================================================================================================
fireID = [3, 4, 7, 16, 22, 28, 37, 42, 54, 63, 67] #2, 4, 7, 16, 22, 28, 37, 42, 54, 63]

regime = 'highNOx' # [lowNOx, highNOx]
params = 'cappa2016'

extension = [
 'koss_highnox_bnzdiol_UNS-ALK',
 'koss_highnox_meffral_UNS-ALK',
 'koss_highnox_bnzdiol_UNS-AR1',
 'koss_highnox_meffral_UNS-AR1',
 #'koss_highnox_bnzdiol_UNS-AR2',
 #'koss_highnox_meffral_UNS-AR2',
]

OH_profile = ['butanol', 'terpene']
emission_mode = ['nonVolatile_IndividualFireProfile']

DLVP_v   = [
 #SVOC, ALK, UNS, BNZ, AR1, AR2, NPH, FRN, P01, P02, P03, P04, P05, P06, P07, P08, P09, P10,  ISP,  TRP
 '1.47 1.47 1.47 1.53 1.42 1.46 1.41 1.29 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.9302 1.91',
 '1.47 1.47 1.47 1.53 1.42 1.46 1.41 1.29 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.9302 1.91',
 '1.47 1.47 1.42 1.53 1.42 1.46 1.41 1.29 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.9302 1.91',
 '1.47 1.47 1.42 1.53 1.42 1.46 1.41 1.29 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.9302 1.91',
 '1.47 1.47 1.46 1.53 1.42 1.46 1.41 1.29 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.9302 1.91',
 '1.47 1.47 1.46 1.53 1.42 1.46 1.41 1.29 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.99 1.9302 1.91',
]

aadt = 10.  # microphysics timestep [s]
nintern = 6 # frequency of printing values or number of internal microphysics steps
# adt = aadt*nintern # timestep for writing output [seconds]  

COAG = 1    # [0 or 1] the switch for On/Off coagulation
VWL  = [1]  # [0 or 1] the switch for On/Off vapor wall loss
PWL  = [0]  # [0 or 1] the switch for On/Off particle wall loss

ke = [0.0, 0.00532, 0.192] #0.1, 0.156, 0.24]#, 2.5e-4, 1.4e-3] # loss rate constant
kw0     = 0.0        # loss rate constant

alpha   = 1.0        # accommodation coefficient 
storg   = 0.025      # [N/m] surface tension
Dbk     = [1.0E-10]#, 1.0E-18, 1.0E-19]#, 1.0E-18, 1.0E-20] #1.0E-15, 1.0E-18, 1.0E-20, 1.0E-22]  # particle-phase diffusion coefficient [m2/s]
kc      = [0.0E0]    # first-order loss rate of species in the particle phase [1/s]
nucrate = 0.0        # [cm-3 s-1] nucleation rate ***** We do not use it now
dilt    = 0.0        # dilution rate ***** Ali said it is not working
boxvol  = 10000000.0 # teflon [cm3] - [CalTech 24 m3, CSU 10 m3, CMU ?? m3]
endtime = [8.0]      # [hours] hours of run time (make even hours)

pres    = 90472.51   # [Pa] pressure
temp    = 298.0      # [K] temperature
rh      = 0.2        # relative humidity

# directory information
# ====================================================================================================
script_directory = os.popen('pwd').read()[:-1]
src_directory = script_directory[:-(len('scripts'))] + 'src'
run_directory = script_directory[:-(len('scripts'))] + 'runs'

gc = pyg.authorize(service_file='../inputs/SOM-TOMAS-905f18e05abd.json')
gsh_input = gc.open('FIREX_fire_inputs')
gsh_emission = gc.open('Emission_FIREX')

CTR = 0
for ctr_emiss, emiss_mod in enumerate(emission_mode):
 df_emiss = gsh_emission.worksheet_by_title(title='%s'%(emiss_mod)).get_as_df()
 for ctr_oh, oh_prf in enumerate(OH_profile):
  df_input = gsh_input.worksheet_by_title(title='OH_%s'%(oh_prf)).get_as_df()
  for ctr_ext, ext in enumerate(extension):
   for f_ctr, fID in enumerate(fireID):
       
    # --- read input sheet --------------------------------------------------------------  
    No1     = df_input.loc[df_input['variables']=='No_1', 'fire%03i'%fID].iloc[0]# [# cm-3] background 1st number conc.
    Dpm1    = df_input.loc[df_input['variables']=='Dpm_1', 'fire%03i'%fID].iloc[0]*1.0e-3 # [microns] background 1st median diameter
    sigma1  = df_input.loc[df_input['variables']=='sigma_1', 'fire%03i'%fID].iloc[0]# 1st background sigma 
    No2     = df_input.loc[df_input['variables']=='No_2', 'fire%03i'%fID].iloc[0]# [# cm-3] background 2nd number conc.
    Dpm2    = df_input.loc[df_input['variables']=='Dpm_2', 'fire%03i'%fID].iloc[0]*1.0e-3# [microns] background 2nd median diameter
    sigma2  = df_input.loc[df_input['variables']=='sigma_2', 'fire%03i'%fID].iloc[0]# 2nd background sigma 
    
    # the following values are for OH concentration equation
    # OH_conc = a_oh*exp(-1.*ax_oh*t) + b_oh*exp(-1.*bx_oh*t) - t in hour and unit is [Molecules/cm3]
    lights_on = df_input.loc[df_input['variables']=='lights on', 'fire%03i'%fID].iloc[0] # [s] of the day
    a_oh      = df_input.loc[df_input['variables']=='a_oh', 'fire%03i'%fID].iloc[0]
    ax_oh     = df_input.loc[df_input['variables']=='ax_oh', 'fire%03i'%fID].iloc[0]
    b_oh      = df_input.loc[df_input['variables']=='b_oh', 'fire%03i'%fID].iloc[0]
    bx_oh     = df_input.loc[df_input['variables']=='bx_oh', 'fire%03i'%fID].iloc[0]
    OH0       = a_oh*np.exp(-1.*ax_oh*0.) + b_oh*np.exp(-1.*bx_oh*0.)
    # 1.5e7*exp(-1.*t)+ 5e6*exp(-0.05*t) t in hour

    nsomprec    = 20 # number of som precursors for parameterizations or SOM grids ***** 1 for now
    somprecname = 'SVOSOMG ALKSOMG UNSSOMG BNZSOMG AR1SOMG AR2SOMG NPHSOMG FRNSOMG '+\
                  'P01SOMG P02SOMG P03SOMG P04SOMG P05SOMG P06SOMG P07SOMG P08SOMG P09SOMG P10SOMG '+\
                  'ISPSOMG TRPSOMG' # som precursor parameterization's name
    dlvp        = '%s'%DLVP_v[ctr_ext]

    seed_dens = df_input.loc[df_input['variables']=='POA density', 'fire%03i'%fID].iloc[0]*1.0e+3
    poa_1stname = 'SVO_C10'
    poa_1st_lenname = len(poa_1stname)

    iorg  = 1259
    ibins = 36   # number of bins
    idiag = 2    # number of diagnostic species including sulfate and ammonia (2)

    index  = np.where(df_emiss.iloc[:,0]!='')[0].shape[0] 
    spname = df_emiss['species'].iloc[:index]
    g_ippm = df_emiss['gas_fire%03i'%fID].iloc[:index]
    p_frac = df_emiss['pfrac_fire%03i'%fID].iloc[:index]
    nspemiss = len(spname)
    emiss_spname = ''
    emiss_ippm   = ''
    seed_frac    = ''
    for i in range(len(spname)):
        emiss_spname = emiss_spname +   '%15s'%(spname.iloc[i])
        emiss_ippm   =   emiss_ippm + '%15.5E'%(g_ippm.iloc[i])
        seed_frac    =    seed_frac + '%15.5E'%(p_frac.iloc[i])

    CTR = CTR+1
    # Running --------------------------------------
    print('%i - General information for the following simulations'%(CTR))
    print('======================================================')
    print('Fire%03i'%fID)
    print('exe_name = %s'%(ext))
    print('emission_mode = %s'%(emiss_mod))
    print('DLVP_v = %s'%(DLVP_v[ctr_ext]))
    print('OH0 =%s', OH0,'[molec/cm3]')
    print('inbins = %s'%ibins)
    print('COAG = %s'%COAG)
    print('boxvol = %s'%boxvol, '[cm3]')
    print('alpha = %s'%alpha)
    print('pressure = %s'%pres, '[Pa]')
    print('temperature = %s'%temp, '[K]')
    print('RH = %s'%rh)

    ctr = 0
    #for numb in range(len(mod)):
    for vwl in VWL:
     for pwl in PWL:
      for KE in ke:
       for tend in endtime:
        for db in Dbk:
         for KC in kc:
          #for ppm in ippmprec:
            ctr+=1

            rname = 'fire%03i_%s_OH%s_%s_%4.2e_%4.2e_%04.1f'%(fID,params,oh_prf,ext,KE,db,tend)
            if pwl==1:
                rname = 'PWL_'+rname
            #if vwl==1:
            #    rname = 'VWL_'+rname
                #rname = precursor+'_'+regime+'_'+tend+'_'+no+'_'+dp+'_'+oh+'_'+ppm

            print('Simulation #%s'%ctr)
            print('---------------')
            print('VWL = %s'%vwl)
            print('PWL = %s'%pwl)
            print('Db = %s [m2/s]'%db)
            print('Kc = %s [1/s]'%KC)
            print('ke = %s [1/s]'%KE)
            print('No1 = %s [#/cm3]'%No1)
            print('Dpm1 = %s [um]'%Dpm1)
            print('sigma1 = %s'%sigma1)
            print('No2 = %s [#/cm3]'%No2)
            print('Dpm2 = %s [um]'%Dpm2)
            print('sigma2 = %s'%sigma2)
            #print('ippm = %s [ppm]'%ppm)
            print('time = %s [h]'%tend)

            if os.path.exists('../runs/%s'%(rname)):
                os.system('rm -r ../runs/%s'%(rname))

            os.system('mkdir ../runs/%s'%(rname))
            os.system('cp ../src/box_%s.exe ../runs/%s/'%(ext, rname))
            os.system('cp ../src/saprc_files/saprc14_rev1.mod.%s ../runs/%s/saprc14_rev1.mod'%(ext, rname))

            f1 = open('../runs/%s/input'%(rname),'w')
            f1.write('%s\n'%rname)
            f1.write('%f\n'%aadt)
            f1.write('%i\n'%nintern)
            f1.write('%1i\n'%COAG)
            f1.write('%1i\n'%vwl)
            f1.write('%1i\n'%pwl)
            f1.write('%f\n'%KE)
            f1.write('%f\n'%kw0)
            f1.write('%f\n'%lights_on)
            f1.write('%e\n'%a_oh) # ***
            f1.write('%e\n'%ax_oh) # ***
            f1.write('%e\n'%b_oh) # ***
            f1.write('%e\n'%bx_oh) # ***
            f1.write('%5.3f\n'%nucrate)
            f1.write('%6.4f\n'%dilt)
            f1.write('%10.1f\n'%boxvol)
            f1.write('%5.2f\n'%tend) # ***
            f1.write('%8.5f\n'%alpha)
            f1.write('%4.2e\n'%db)
            f1.write('%f\n'%KC)
            f1.write('%8.5f\n'%storg)
            f1.write('%10.3f\n'%pres)
            f1.write('%6.2f\n'%temp)
            f1.write('%8.5f\n'%rh)
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
            f1.write('%s\n'%emiss_ippm)
            f1.write('%s\n'%seed_frac)
            f1.write('%s\n'%seed_dens)
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
            f2.write('ke          = %s [1/s]\n'%KE)
            f2.write('kw0         = %s [1/s]\n'%kw0)
            f2.write('LightsOn    = %s [s]\n'%lights_on)
            f2.write('OH eqn. = a_oh*exp(-1.*ax_oh*t) + b_oh*exp(-1.*bx_oh*t)') # ***
            f2.write('a_oh        = %s\n'%a_oh) # ***
            f2.write('ax_oh       = %s\n'%ax_oh) # ***
            f2.write('b_oh        = %s\n'%b_oh) # ***
            f2.write('bx_oh       = %s\n'%bx_oh) # ***

            f2.write('nucrate     = %s [cm-3 s-1]\n'%nucrate)
            f2.write('dilt        = %s\n'%dilt)
            f2.write('boxvol      = %s [cm3]\n'%boxvol)
            f2.write('endtime     = %s [hours]\n'%tend) # ***
            f2.write('alpha       = %s\n'%alpha)
            f2.write('Db          = %s [m2/s]\n'%db)
            f2.write('Kc          = %s [1/s]\n'%KC)
            f2.write('storg       = %s [N/m]\n'%storg)
            f2.write('pressure    = %s [Pa]\n'%pres)
            f2.write('temperature = %s [K]\n'%temp)
            f2.write('RH          = %s\n'%rh)
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
            f2.write('emiss_ippm  = %s\n'%emiss_ippm)
            f2.write('seed_frac   = %s\n'%seed_frac)
            f2.write('seed_dens   = %s\n'%seed_dens)

            f2.close()
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
            runme_line[ind3] = './box_%s.exe <input> /dev/null'%(ext)
            #runme_line[ind3] = './box_%s.exe <input> %s.out'%(ext, rname)

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



