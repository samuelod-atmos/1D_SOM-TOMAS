

# Importing Libraries
# ====================================================================================================
import numpy as np
import pandas as pd
import pygsheets as pyg
import subprocess as sp
import time
import re
import os

# Clean ipython
# ====================================================================================================
from IPython import get_ipython
get_ipython().magic('reset -sf') 
'''
startTime = time.time()

from variable_list import *

# saprc
# ----------------------------------------------
saprc = 'saprc14_rev1'

# read saprc11b.rxn raw file
# ====================================================================================================
filename2 = 'raw_saprc11b.rxn'
df2 = pd.read_csv('../inputs/%s'%filename2, header=None) # reading the default saprc11b.rxn file
df2 = np.array(df2) # converting to an array to search and make changes

# Find the generic voc species and substitude the 
# ====================================================================================================
ind2 = np.where(df2=='GENVOC   0.0     0.0     0.0   0   0   0   1')[0][0]
# order        Name Conc(0) Mwt  #C   #N    #S   #O   #X
df2[ind2][0] = 'GENVOC   %3.1f   %6.2f   %4.1f   %1.0f   %1.0f   %1.0f   %1.0f'\
               %(conc, mwt, cno, nno, sno, ono, xno)
print('GENVOC   %3.1f   %6.2f   %4.1f   %1.0f   %1.0f   %1.0f   %1.0f'\
    %(conc, mwt, cno, nno, sno, ono, xno))

# Find the generic voc species and substitude the 
# ====================================================================================================
ind3 = np.where(df2=='GENSOMG     0      0     0      0     0.0001 0.0001 0.0001 0.0001  0.0001')[0][0]
df2[ind3][0] = 'GENSOMG    %2.0f      %2.0f   %2.0f     %2.0f     %6.4f %6.4f %6.4f %6.4f  %6.4f'\
    %(cstart, cstop, ostart, ostop, pf1, pf2, pf3, pf4, cfrag)
print('GENSOMG    %2.0f      %2.0f   %2.0f     %2.0f     %6.4f %6.4f %6.4f %6.4f  %6.4f'\
    %(cstart, cstop, ostart, ostop, pf1, pf2, pf3, pf4, cfrag))

# Finding the reaction and change the reaction rates based on the precursors
# ====================================================================================================
ind4 = np.where(df2=='BL20) 0.00e-00               ;GENVOC + OH = GENSOMG')[0][0]
df2[ind4][0] = 'BL20) %1.2E               ;GENVOC + OH = GENSOMG'%koh
print('BL20) %1.2E               ;GENVOC + OH = GENSOMG'%koh)
print('-----------------------------------------------------------------------------')

# Creating saprc11b.rxn file to run prep.exe
# ====================================================================================================
outfilename = 'saprc11b.rxn'
if os.path.exists('../inputs/%s'%outfilename):
    os.system('rm ../inputs/%s'%outfilename)
f1 = open('../inputs/%s'%outfilename, 'w')

for i in range(df2.shape[0]):
    f1.write(df2[i][0])
    f1.write('\n')

f1.close()
'''

sv_c_start = 5
sv_c_end = 30
for i in range(sv_c_start,sv_c_end+1):
    print('SVO_C%02i  0.0      0.0   %4.1f   0   0  0  1'%(i, i))

iv_c_start = 5
iv_c_end = 25
for i in range(iv_c_start,iv_c_end+1):
    print('IVO_C%02i  0.0      0.0   %4.1f   0   0  0  1'%(i, i))

alk_c_start = 5
alk_c_end = 20
for i in range(alk_c_start, alk_c_end+1):
    mw_alk = i*12.+(i*2+2)*1.
    print('ALK_C%02i  0.0    %6.2f  %4.1f   0   0  0  1'%(i, mw_alk, i))

'''
BENZENE  0.0    78.11   6.0   0   0  0  1
ARO1     0.0    92.14   7.0   0   0  0  1
ARO2     0.0   106.16   8.0   0   0  0  1
ISOPRENE 0.0     0.0    5.0   0   0  0  1
TERP     0.0     0.0   10.0   0   0  0  1
SESQ     0.0     0.0   15.0   0   0  0  1
SVOSOMG  0.0     0.0    0.0   0   0  0  1
IVOSOMG  0.0     0.0    0.0   0   0  0  1
ALKSOMG  0.0     0.0    0.0   0   0  0  1
BNZSOMG  0.0     0.0    0.0   0   0  0  1
AR1SOMG  0.0     0.0    0.0   0   0  0  1
AR2SOMG  0.0     0.0    0.0   0   0  0  1
ISPSOMG  0.0     0.0    0.0   0   0  0  1
TRPSOMG  0.0     0.0    0.0   0   0  0  1

!SPEC       CSTART CSTOP OSTART OSTOP PF1  PF2  PF3  PF4   CFRAG
SVOSOMG      1      30    0      7    0.9270 0.01010 0.0180 0.0445 0.0980
IVOSOMG      1      25    0      7    0.9270 0.01010 0.0180 0.0445 0.0980
ALKSOMG      1      20    0      7    0.9270 0.01010 0.0180 0.0445 0.0980
AR1SOMG      1       7    0      7    0.0029 0.00100 0.0010 1.0100 0.2220
AR2SOMG      1      15    0      7    0.9360 0.00100 0.0021 0.0609 0.0100
BNZSOMG      1       6    0      7    0.0792 0.00100 0.9190 0.0010 0.5350
ISPSOMG      1       5    0      7    0.9520 0.00109 0.0304 0.0163 0.0632
TRPSOMG      1      15    0      7    0.0680 0.63300 0.2750 0.0244 0.0353

'''
