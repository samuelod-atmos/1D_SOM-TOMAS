# Libraries
# ====================================================================================
import numpy as np
import pandas as pd
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1 import make_axes_locatable

# Changing font
# ====================================================================================
#mpl.rcParams['font.family'] = 'Helvetica'
#params = {'mathtext.default': 'regular' }          
#plt.rcParams.update(params)

# inputs
# ====================================================================================
#name = 'test_3.0E06_4.0'
#name = 'c12_highnox_nowall_1.5E06_4.0'
#name = 'toluene_1.5E06_4.0'
name = 'test_n-dodecane_lownox-6_04.0_050000.0_0.100_1.5E+06_0.010'
ibins = 36
rho = 1400. # [kg/m3]
temp = 292.58 # [K]
pres = 101325.0*(1-2.25577e-5*977.0)**5.25588 # [Pa]
boxvol = 7000000.0 # [cm3]
#boxvol = 24000000.0 # [cm3]
R = 8.314 # [j/mol/K]

# diameter
# ------------------------------------------------------------------------------------
Mo = 1.0e-21*2.e0**(-6) # [kg/particle]
dp = (Mo/rho*6/np.pi)**(1./3.) # [m]
#dp = dp*1.e9 # [nm]
print(dp, '[m]')

xk = []
for i in range(37):
    xk.append(Mo*(2**i))
xk=np.array(xk)

dpi = (xk/rho*6/np.pi)**(1./3.)*1.e9 # [nm]
print(dpi, '[nm]')

# Number Conc.
# ------------------------------------------------------------------------------------
df = np.array(pd.read_csv('../outputs/%s_noconc.dat'%name, header=None, delim_whitespace=True))
#df = np.array(pd.read_csv('../tomas_bian_totallycleaned/burn37_160429_chemon2_0.943_2.493_1000000._5843.0_0.157_0._noconc.dat', header=None, delim_whitespace=True))
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
df = np.array(pd.read_csv('../outputs/%s_aemass.dat'%name, header=None, delim_whitespace=True))
Mk = df[:,1:]
Mk = np.reshape(Mk, (time.shape[0], int(Mk.shape[0]/time.shape[0]), Mk.shape[1]))
Mkorg = Mk[:,1:202,:]

# Gas Conc.
# ------------------------------------------------------------------------------------
df = np.array(pd.read_csv('../outputs/%s_gc.dat'%name, header=None, delim_whitespace=True))
Gc = df[:,1:].astype(float)
Gcorg = Gc[:,1:202]

# SAPRC Gas Conc.
# ------------------------------------------------------------------------------------
df = np.array(pd.read_csv('../outputs/%s_saprcgc.dat'%name, header=None, delim_whitespace=True))
saprcgc = df[:,:].astype(float)

precname = 'ALK5'
saprc_doc = pd.read_csv('../src/saprc14_rev1.doc', header=None, skiprows=38, nrows=303, delim_whitespace=True)
ind = np.where(np.array(saprc_doc)==precname)[0][0]
prec_conc = saprcgc[:, ind]
prec_mw = saprc_doc.loc[ind,2]

# Specification (gas names, som names, carbon number, oxygen number)
# ------------------------------------------------------------------------------------
spec = pd.read_csv('./%s_spec.dat'%name, header=None, delim_whitespace=True)
act_sp = np.array(spec.loc[0,3:])
CNO = np.array(spec.loc[2,3:203]).astype(int)
ONO = np.array(spec.loc[3,3:203]).astype(int)
mwt = np.array(spec.loc[4,3:203]).astype(float)
mworg = np.array(spec.loc[5,3:203]).astype(float)
cstar = np.array(spec.loc[6,1:201]).astype(float)

# PROCESSING
# ====================================================================================
Gcorgppm = np.zeros(Gcorg.shape)
for i in range(len(time)):
    Gcorgppm[i,:] = Gcorg[i,:]*1.e9/pres*R*temp/mwt[:]*1.e3/boxvol*1.e6 # kg/bag -> ppm
    Gcorg[i,:] = Gcorg[i,:]*1.e9/boxvol*1.e6 # kg/bag -> ug/m3
    for j in range(ibins):
        Mkorg[i,:,j] = Mkorg[i,:,j]/boxvol*1.0e6*1.0e9 # kg/bag -> ug/m3

mass_org = np.sum(Mkorg[1:,:,:], axis=2)
mass_org = mass_org + Gcorg[1:]
ppm_org = mass_org/mwt/pres*R*temp
ppmc_org = ppm_org*CNO
c_org = np.sum(ppm_org, axis=1)

print('SOA Concentration =', np.sum(Mkorg[1:,1:202,:], axis=(1,2))[-1]-\
                             np.sum(Mkorg[1:,1:202,:], axis=(1,2))[0],\
                             '[ug/m3]')

# PLOTS
# ====================================================================================
'''
# Banana plot
# ------------------------------------------------------------------------------------
T,D = np.meshgrid(time,dpi[:-1])

plt.figure()
p = plt.pcolormesh(T,D,np.transpose(Nk))
plt.yscale('log')
plt.colorbar(p)
plt.show()

# Size distribution plot
# ------------------------------------------------------------------------------------
plt.figure()

plt.plot(dpi[:-1], Nk[0,:], lw=1.6, label='initial')
plt.plot(dpi[:-1], Nk[-1,:], lw=1.6, label='final')

plt.xlabel('$D_p\ [nm]$')
plt.ylabel('$dN\ dlogD_p\ [cm^{-3}]$')
plt.xscale('log')

plt.legend()
plt.tight_layout()
plt.show()

# Gas concentration
# ------------------------------------------------------------------------------------
plt.figure()

plt.plot(time[1:], np.sum(Gcorg[1:,1:202], axis=1), '.')

plt.xlabel('time [hr]')
plt.ylabel('total gas concentration [$\mu g\ m^{-3}$]')

plt.tight_layout()
plt.show()

# SOA concentration
# ------------------------------------------------------------------------------------
plt.figure()

plt.plot(time[1:], np.sum(Mkorg[1:,1:202,:], axis=(1,2)), '.')

plt.xlabel('time [hr]')
plt.ylabel('total particle concentration [$\mu g\ m^{-3}$]')

plt.tight_layout()
plt.show()
print('SOA Concentration =', np.sum(Mkorg[1:,1:202,:], axis=(1,2))[-1]-\
                             np.sum(Mkorg[1:,1:202,:], axis=(1,2))[0],\
                             '[ug/m3]')
# total gas and particle concentration
# ------------------------------------------------------------------------------------
plt.figure()

plt.plot(time[1:], np.sum(Mkorg[1:,1:202,:], axis=(1,2))+np.sum(Gcorg[1:,1:202], axis=1), '.')

plt.xlabel('time [hr]')
plt.ylabel('total aerosol concentration [$\mu g\ m^{-3}$]')

plt.tight_layout()
plt.show()

# CSTAR plot
# ------------------------------------------------------------------------------------
fig, ax = plt.subplots(1, 1, facecolor='w', figsize=(11,3))

ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["left"].set_visible(False)

p = ax.scatter(CNO, ONO, c=np.log(cstar), marker='s', s=420, cmap='jet')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.05)
cb = fig.colorbar(p, cax=cax)

cb.outline.set_edgecolor(None)
plt.ylabel('$log\ C^*$', fontsize=14)
plt.tick_params(labelsize=12)

plt.axes(ax)
plt.xlabel('Carbon Number', fontsize=14)
plt.ylabel('Oxygen Number', fontsize=14)
plt.tick_params(labelsize=12)
plt.ylim([.5,7.5])
plt.xlim([.5,30.5])

plt.tight_layout()#h_pad=.01)
#plt.savefig('./figures/Cstar_COspace.pdf')
plt.show()
'''
'''
# Gc in Carbon/Oxygen Space plot
# ------------------------------------------------------------------------------------
for i in range(Gcorg.shape[0]/60):
    fig, ax = plt.subplots(1, 1, facecolor='w', figsize=(11,3))
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    
    p = ax.scatter(CNO, ONO, c=Gcorg[i*60,:], marker='s', s=420, cmap='jet')#,
                   #vmax=0.002, vmin=0.)
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    cb = fig.colorbar(p, cax=cax)
    
    cb.outline.set_edgecolor(None)
    plt.ylabel('Gas Concentration [$\mu g\ m^{-3}$]', fontsize=14)
    plt.tick_params(labelsize=12)
    
    plt.axes(ax)
    plt.xlabel('Carbon Number', fontsize=14)
    plt.ylabel('Oxygen Number', fontsize=14)
    plt.tick_params(labelsize=12)
    plt.ylim([.5,7.5])
    plt.xlim([.5,30.5])
    
    plt.tight_layout()#h_pad=.01)
#    plt.savefig('./figures/Gcorg_time_C_O_minute_%s.jpg'%(i*60))
    plt.show()

'''
'''
# Aerosol in Carbon/Oxygen Space plot
# ------------------------------------------------------------------------------------
for i in range(Gcorg.shape[0]):
    fig, ax = plt.subplots(1, 1, facecolor='w', figsize=(11,3))
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    
    #p = ax.scatter(CNO, ONO, c=np.sum(Mkorg[i*60,:,:], axis=1)+Gcorg[i*60,:], marker='s', s=420, cmap='jet')
    p = ax.scatter(CNO, ONO, c=np.sum(Mkorg[i,:,:], axis=1), marker='s', s=420, cmap='jet')    
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    cb = fig.colorbar(p, cax=cax)
    
    cb.outline.set_edgecolor(None)
    plt.ylabel('Aerosol Concentration [$\mu g\ m^{-3}$]', fontsize=14)
    plt.tick_params(labelsize=12)
    
    plt.axes(ax)
    plt.xlabel('Carbon Number', fontsize=14)
    plt.ylabel('Oxygen Number', fontsize=14)
    plt.tick_params(labelsize=12)
    plt.ylim([.5,7.5])
    plt.xlim([.5,30.5])
    
    plt.tight_layout()#h_pad=.01)
#    plt.savefig('./figures/Gcorg_time_C_O_minute_%s.jpg'%(i*60))
    plt.show()
'''
'''
# Gc in Carbon/Oxygen Space plot
# ------------------------------------------------------------------------------------
for i in range(Gcorg.shape[0]/60):
    fig, ax = plt.subplots(1, 1, facecolor='w', figsize=(11,3))
    
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    
    p = ax.scatter(CNO, ONO, c=saprcgc[i*60,105:], marker='s', s=420, cmap='jet',
                   vmax=0.005, vmin=0.)
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    cb = fig.colorbar(p, cax=cax)
    
    cb.outline.set_edgecolor(None)
    plt.ylabel('Gas Concentration [ppm]', fontsize=14)
    plt.tick_params(labelsize=12)
    
    plt.axes(ax)
    plt.xlabel('Carbon Number', fontsize=14)
    plt.ylabel('Oxygen Number', fontsize=14)
    plt.tick_params(labelsize=12)
    plt.ylim([.5,7.5])
    plt.xlim([.5,30.5])
    
    plt.tight_layout()#h_pad=.01)
#    plt.savefig('./figures/saprcgc_time_C_O_minute_%s.jpg'%(i*60))
plt.show()
'''
#plt.close('all')
