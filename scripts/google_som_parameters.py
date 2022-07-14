# Libraries
# ===================================================================================================
import pygsheets as pyg
import os

# Function
# ===================================================================================================
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# Reading from a google sheet
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
def rd_ggl_sht(jsn, fln, shn, nm, hdn): # (json file, google file name, sheetname, name, header name)
    # authorization
    gc = pyg.authorize(service_file=jsn)
    # open the google spreadsheet
    gsh = gc.open(fln)
    # select the first sheet
    wks = gsh[shn] # [ 0-Jathar_2015_GMD, 1-Zhang_2014_PNAS ]
    df1 = wks.get_as_df()
    params = df1.where(df1[hdn]==nm) # find the voc and som parameters

    if len(params)==0:
        print('YO!!! pay attention! your either precursor name or regime name is wrong.'+
              'Take a look at the SOM parameters google sheets again.')
        os.stop()
    
    params = params.dropna() # delete NaNs

    cfrag = params['Fragmentation'].get_values()[0] # fragmentation
    dlvp  = params['dlvp'].get_values()[0] # decrease in volatility
    pf1   = params['P1'].get_values()[0] # probability of adding 1 oxygen
    pf2   = params['P2'].get_values()[0] # probability of adding 2 oxygen
    pf3   = params['P3'].get_values()[0] # probability of adding 3 oxygen
    pf4   = params['P4'].get_values()[0] # probability of adding 4 oxygen
    koh   = params['kOH'].get_values()[0] # rate constant with OH
    
    # precursor inputs
    # ----------------------------------------------
    conc = 0.0  # concentration which should be zero
    mwt = params['MW'].get_values()[0]  # molecular weight of the precursors
    cno = params['CNO'].get_values()[0]  # carbon number of the species
    nno = params['NNO'].get_values()[0] # nitrogen number of the species
    sno = params['SNO'].get_values()[0] # sulfur number of the species
    ono = params['ONO'].get_values()[0] # oxygen number of the species
    xno = params['XNO'].get_values()[0] # 

    return [cfrag, dlvp, pf1, pf2, pf3, pf4, koh, conc, mwt, cno, nno, sno, ono, xno]
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

