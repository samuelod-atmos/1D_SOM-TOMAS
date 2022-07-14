''' ====================================================================================================
This code creates the SPARC11b.rxn file to add SOM precursor so that it runs prep.exe and creates
saprc14_rev1.mod and SAPRC14_rev1.f files

written by Ali Akherati October 2018 Colorado State University
==================================================================================================== '''

import pandas as pd
import numpy as np
import os
import time
import subprocess as sp

# Clean ipython
# ======================================================================================================
#from IPython import get_ipython
#get_ipython().magic('reset -sf') 

# Inputs
# ======================================================================================================
# precursor inputs
# ----------------------------------------------
conc = 0.0  # concentration which should be zero
mwt = 170.0 # molecular weight of the precursors
cno = 12.0  # carbon number of the species
nno = 0     # nitrogen number of the species
sno = 0     # sulfur number of the species
ono = 0     # oxygen number of the species
xno = 1     # 

# SOM grid inputs
# ----------------------------------------------
cstart = 1
cstop  = 30
ostart = 1
ostop  = 7
regime     = 'highnox' # [highnox or lownox]
precursor  = 'n-dodecane' # [n-dodecane, benzene, toluene, m-xylene, isoprene, alpha-pinene]
sheetn = 'Jathar_2015_GMD' # [Jathar_2015_GMD, Zhang_2014_PNAS]
filename1  = 'som_params' # excel file name which has contained SOM parameters

name = '%s_%s'%(precursor, regime) # name of the precursor and regime to find in the excel sheet
df1  = pd.read_excel('./%s.xlsx'%filename1, sheet_name=sheetn) # reading the excel file
df1  = np.array(df1) # converting to an array for searching
ind1 = np.where(df1==name)[0][0] # find the index for the specific precursor and regime

cfrag = df1[ind1][1] # fragmentation
dlvp  = df1[ind1][2] # decrease in volatility
pf1   = df1[ind1][3] # probability of adding 1 oxygen
pf2   = df1[ind1][4] # probability of adding 2 oxygen
pf3   = df1[ind1][5] # probability of adding 3 oxygen
pf4   = df1[ind1][6] # probability of adding 4 oxygen
koh   = df1[ind1][7] # rate constant with OH

# read saprc11b.rxn raw file
# ======================================================================================================
filename2 = 'raw_saprc11b.rxn'
df2 = pd.read_csv('./%s'%filename2, header=None) # reading the default saprc11b.rxn file
df2 = np.array(df2) # converting to an array to search and make changes

# Find the generic voc species and substitude the 
# ======================================================================================================
ind2 = np.where(df2=='GENVOC   0.0     0.0     0.0   0   0   0   1')[0][0]
# order        Name Conc(0) Mwt  #C   #N    #S   #O   #X
df2[ind2][0] = 'GENVOC   %3.1f   %6.2f   %4.1f   %1.0f   %1.0f   %1.0f   %1.0f'\
               %(conc, mwt, cno, nno, sno, ono, xno)
print('GENVOC   %3.1f   %6.2f   %4.1f   %1.0f   %1.0f   %1.0f   %1.0f'\
    %(conc, mwt, cno, nno, sno, ono, xno))

# Find the generic voc species and substitude the 
# ======================================================================================================
ind3 = np.where(df2=='GENSOMG     0      0     0      0     0.0001 0.0001 0.0001 0.0001  0.0001')[0][0]
df2[ind3][0] = 'GENSOMG    %2.0f      %2.0f   %2.0f     %2.0f     %6.4f %6.4f %6.4f %6.4f  %6.4f'\
    %(cstart, cstop, ostart, ostop, pf1, pf2, pf3, pf4, cfrag)
print('GENSOMG    %2.0f      %2.0f   %2.0f     %2.0f     %6.4f %6.4f %6.4f %6.4f  %6.4f'\
    %(cstart, cstop, ostart, ostop, pf1, pf2, pf3, pf4, cfrag))

# Finding the reaction and change the reaction rates based on the precursors
# ======================================================================================================
ind4 = np.where(df2=='BL20) 0.00e-00               ;GENVOC + OH = GENSOMG')[0][0]
df2[ind4][0] = 'BL20) %1.2E               ;GENVOC + OH = GENSOMG'%koh
print('BL20) %1.2E               ;GENVOC + OH = GENSOMG'%koh)

# Creating saprc11b.rxn file to run prep.exe
# ======================================================================================================
outfilename = 'saprc11b.rxn'
if os.path.exists('./%s'%outfilename):
    os.system('rm %s'%outfilename)
f1 = open('./%s'%outfilename, 'w')

for i in range(df2.shape[0]):
    f1.write(df2[i][0])
    f1.write('\n')

f1.close()

# Running SAPRC
# ======================================================================================================
saprc = 'saprc14_rev1'
os.system('cp %s ../v1_somtomas/'%outfilename) # Copying saprc11b.rxn to SAPRC folder
os.system('cd ../v1_somtomas/ ; echo %s | ./prep.exe'%(saprc)) # change the cwd and run prep.exe
os.system('mv ../v1_somtomas/%s.mod ../'%(saprc)) # move the .mod file to COSO 
os.system('mv ../v1_somtomas/%s.ftn ../v1_somtomas/%s.f'%(saprc, saprc)) # change .ftn to .f
os.system('mv ../v1_somtomas/%s.f ../'%(saprc)) # move the .f file to COSO

# Correcting the saprc14_rev1.f file
# ======================================================================================================
df3 = np.array(pd.read_table('../%s.f'%(saprc), header=None))
df3 = np.insert(df3, 113, df3[113][0][:40], 0)
df3[114] = '     &     ' + df3[114][0][40:]

output_file = '%s.f'%(saprc)
os.system('rm ../%s.f'%(saprc))
f2 = open('../%s.f'%(saprc), 'w')

for i in range(df3.shape[0]):
    f2.write(df3[i][0])
    f2.write('\n')

f2.close()


#f = pd.DataFrame(df3)
#f.to_csv('../saprc123.f', header=None, index=None)
