# -*- coding: utf-8 -*-
"""
Created on Fri May  9 09:19:15 2025

@author: User
"""

#imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
import os



"""
THINGS TO PLOT:

Transfer function (Theta_l) and integrand of Cl (Theta_l^2/k) for 3 values of l

CMB power spectrum (maybe include the separate contributions)

EE and TE spectra as well

Matter power spectrum with equality scale k_eq. marked in
"""

#Units
m           = 1.0;                         #Length (in meters)
s           = 1.0                       #Time (in seconds)
kg          = 1.0                       #Kilo (in kilos)
K           = 1.0                       #Temperature (in Kelvins)

#Derived units
km          = 1e3 * m                 #Kilometers
N           = kg*m/(s*s)                #Newton
J           = N*m                       #Joule
W           = J/s                       #Watt
Mpc         = 3.08567758e22 * m         #Megaparsec
eV          = 1.60217653e-19 * J        #Electronvolt
yr = 365 * 24 * 60 * 60 * s #year
Gyr = 1e9 * yr  #Gigayear

# Physical constants    
k_b         = 1.38064852e-23 * J/K      #Bolzmanns constant
m_e         = 9.10938356e-31 * kg       #Mass of electron
m_H         = 1.6735575e-27 * kg        #Mass of hydrogen atom
c           = 2.99792458e8 * m/s        #Speed of light
G           = 6.67430e-11 * N*m*m/(kg*kg)#Gravitational constant
hbar        = 1.054571817e-34 * J*s;        # Reduced Plancks constant
sigma_T     = 6.6524587158e-29 * m*m;       # Thomas scattering cross-section
lambda_2s1s = 8.227 / s;                    # Transition time between 2s and 1s in Hydrogen
H0_over_h   = 100 * km/s/Mpc;               # H0 / h
epsilon_0   = 13.605693122994 * eV;         # Ionization energy for the ground state of hydrogen
xhi0        = 24.587387 * eV;               # Ionization energy for neutral Helium
xhi1        = 4.0 * epsilon_0;              # Ionization energy for singly ionized Helium

# Small h (needed for plotting the matter power spectrum)
h           = 0.67

#opening the recombination.txt file
with open("cells.txt", 'r') as file:
    #makes a list of the lines   
    lines = file.readlines()
    #make a dictionary where the keys are the variables and the values are lists of their corresponding values
    valuesdic1 = {'ell': [float(line.split(' ')[0]) for line in lines], 'CellTT': [float(line.split(' ')[1]) for line in lines],'CellEE': [float(line.split(' ')[2]) for line in lines], 'CellTE': [float(line.split(' ')[3]) for line in lines]}
    
    
with open("pofk.txt", 'r') as file:
    #makes a list of the lines   
    lines = file.readlines()
    #make a dictionary where the keys are the variables and the values are lists of their corresponding values
    valuesdic2 = {'k': [float(line.split(' ')[0]) for line in lines[:-2]], 'pofk': [float(line.split(' ')[1]) for line in lines[:-2]], 'Theta10': [float(line.split(' ')[9]) for line in lines[:-2]], 'Theta100': [float(line.split(' ')[21]) for line in lines[:-2]], 'Theta300': [float(line.split(' ')[30]) for line in lines[:-2]]}
    
    k_eq = float(lines[-1])

print(f"k at matter radiation equality: {k_eq*Mpc}")
with open("contributions.txt", "r") as file:
    lines = file.readlines()
    valuesdic3 = {'ell': [float(line.split(' ')[0]) for line in lines], 'CellSW': [float(line.split(' ')[1]) for line in lines], 'CellISW': [float(line.split(' ')[2]) for line in lines], 'CellDoppler': [float(line.split(' ')[3]) for line in lines], 'CellQuad': [float(line.split(' ')[4]) for line in lines]}
"""
Define the directories for plots. If they don't yet exist, create them.
"""

PROJECT_ROOT_DIR = "Project_plots"
if not os.path.exists(PROJECT_ROOT_DIR):
    os.makedirs(PROJECT_ROOT_DIR)

SPECTRUM_DIR = os.path.join(PROJECT_ROOT_DIR, "PowerSpectra")
if not os.path.exists(SPECTRUM_DIR):
    os.makedirs(SPECTRUM_DIR)

max_ell_tt = valuesdic1['CellTT'].index(max(valuesdic1['CellTT']))
first_trough_TT = 0
second_peak_TT = 0
trough_counter = 1
peak_counter = 1

for i,j in zip(valuesdic1['ell'][max_ell_tt+1:-1], valuesdic1['CellTT'][max_ell_tt+1:-1]):
    index = valuesdic1['ell'].index(i)

    if (j < valuesdic1['CellTT'][index+1] and j < valuesdic1['CellTT'][index-1] and trough_counter):
        trough_counter = 0
        first_trough_TT = i
    if (j > valuesdic1['CellTT'][index+1] and j > valuesdic1['CellTT'][index-1] and peak_counter):
        peak_counter = 0
        second_peak_TT = i

print(max_ell_tt, first_trough_TT, second_peak_TT)

with open("data/planck_cell_low.txt", 'r') as file:
    next(file)
    lines = file.readlines()

    low_ell_errors = {'ell': [float(line.lstrip(' ').split(' ')[0]) for line in lines], 'Cell':[float(line.lstrip(' ').split(' ')[7]) for line in lines], 'ErrorUp':[float(line.lstrip(' ').split(' ')[14]) for line in lines],'ErrorDown':[float(line.lstrip(' ').split(' ')[21].rstrip('\n')) for line in lines]}

with open("data/planck_cell_high.txt", 'r') as file:
    next(file)
    lines = file.readlines()
    
    high_ell_errors = {'ell': [float(line.lstrip(' ').split(' ')[0]) for line in lines], 'Cell':[float(line.lstrip(' ').split(' ')[3]) for line in lines], 'ErrorDown':[float(line.lstrip(' ').split(' ')[6]) for line in lines],'ErrorUp':[float(line.lstrip(' ').split(' ')[9].rstrip('\n')) for line in lines]}

with open("data/planck_cell_EE.txt", 'r') as file:
    next(file)
    lines = [line.split(' ') for line in file.readlines()]
    for line in lines:
        while '' in line:
            line.remove('')
    EE_ell_errors = {'ell': [float(line[0]) for line in lines], 'Cell':[float(line[1]) for line in lines], 'ErrorDown':[float(line[2]) for line in lines],'ErrorUp':[float(line[3].rstrip('\n')) for line in lines]}
    
with open("data/planck_cell_TE.txt", 'r') as file:
    next(file)
    lines = [line.split(' ') for line in file.readlines()]
    for line in lines:
        for i in line:
            while '' in line:
                line.remove('')
    TE_ell_errors = {'ell': [float(line[0]) for line in lines], 'Cell':[float(line[1]) for line in lines], 'ErrorDown':[float(line[2]) for line in lines],'ErrorUp':[float(line[3].rstrip('\n')) for line in lines]}
    
with open("data/WMAP_ACT.txt", 'r') as file:
    next(file)
    lines = [line.split(' ') for line in file.readlines()]
    for line in lines:
        for i in line:
            while '' in line:
                line.remove('')
    WMAP_ACT_errors = {'k': [float(line[0]) for line in lines], 'p_of_k':[float(line[1]) for line in lines], 'Error':[float(line[2])-float(line[1]) for line in lines]}

with open("data/SDSS_DR7_LRG.txt", 'r') as file:
    next(file)
    lines = [line.split(' ') for line in file.readlines()]
    for line in lines:
        for i in line:
            while '' in line:
                line.remove('')
    SDSS_DR7_LRG_errors = {'k': [float(line[0]) for line in lines], 'p_of_k':[float(line[1]) for line in lines], 'Error':[float(line[2]) for line in lines]}



#%%
fig, axes = plt.subplots(nrows=2, figsize=(9.5, 13.5), sharex=False)
xticks = [2, 10, 50, 500, 1000, 2000]
error_array_low = np.array([low_ell_errors['ErrorDown'], low_ell_errors['ErrorUp']])
error_array_high = np.array([high_ell_errors['ErrorDown'], high_ell_errors['ErrorUp']])

for i in axes:
    i.set_xlabel(r'$\ell$', fontsize = 18)
    i.tick_params(axis='both', labelsize=12)
    i.set_xlim([2,2000])
    
    i.yaxis.labelpad = 20
    i.set_xscale('log')
    i.set_xticks(xticks, labels=[str(x) for x in xticks])
    i.tick_params(axis='both', labelsize=15)
axes[0].set_title('CMB Temperature power spectrum', fontsize = 25, pad = 10)
axes[0].set_ylabel(r'$\ell(\ell+1)C_\ell^{TT}/\,2\pi\ (\mu K^2)$', fontsize = 23)
axes[0].plot(valuesdic1['ell'], valuesdic1['CellTT'], label = 'Numerical prediction')
axes[0].errorbar(low_ell_errors['ell'], low_ell_errors['Cell'], yerr = error_array_low, fmt = ".", color= 'orange', label ='Planck 2018')
axes[0].errorbar(high_ell_errors['ell'], high_ell_errors['Cell'], yerr = error_array_high, fmt = ".", color = 'orange')

axes[1].set_ylabel(r'Individual terms $(\mu K^2)$', fontsize = 23)
axes[1].plot(valuesdic3['ell'], valuesdic3['CellSW'], label = 'Sachs Wolfe effect')
axes[1].plot(valuesdic3['ell'], valuesdic3['CellISW'], label = 'ISW effect')
axes[1].plot(valuesdic3['ell'], valuesdic3['CellDoppler'], label = 'Doppler effect')
axes[1].plot(valuesdic3['ell'], valuesdic3['CellQuad'], label = 'Quadrupolar correction')

axes[0].legend(fontsize = 20)
axes[1].legend(fontsize = 20)
plt.tight_layout()

plt.savefig(SPECTRUM_DIR + "/TT_Spectrum.png", bbox_inches = "tight")
#%%
fig, ax = plt.subplots(figsize=(10, 6.5))

# Main plot

error_array_EE = np.array([EE_ell_errors['ErrorDown'], EE_ell_errors['ErrorUp']])

ax.set_title('CMB Polarization power spectrum', fontsize = 25, pad = 10)
ax.set_xlabel(r'$\ell$', fontsize = 25)
ax.set_ylabel(r'$\ell(\ell+1)C_\ell^{EE}/\,2\pi\ (\mu K^2)$', fontsize = 23)
ax.plot(valuesdic1['ell'], [i*(j+2)*(j+1)*j*(j-1) for i,j in zip(valuesdic1['CellEE'], valuesdic1['ell'])], label = 'Numerical prediction')
ax.errorbar(EE_ell_errors['ell'], EE_ell_errors['Cell'], yerr = error_array_EE, fmt = ".", color= 'orange', label ='Planck 2018')
ax.set_xlim([2,2000])
xticks = [500, 1000, 1500, 2000]
ax.set_xticks(xticks, labels=[str(x) for x in xticks])
ax.tick_params(axis='both', labelsize=17)

axins = inset_axes(ax, width="30%", height="35%", loc="lower center")

x = valuesdic1['ell']
y = [i*(j+2)*(j+1)*j*(j-1) for i,j in zip(valuesdic1['CellEE'], valuesdic1['ell'])]
# Plot the same data on the inset
axins.plot(x, y)
# Set limits for zoom region
axins.set_xlim(0, 40)
axins.set_ylim(min(y[:42]), max(y[:42])+0.1)

#axins.set_xticks([])
axins.tick_params(labelsize=8)
axins.xaxis.set_ticks_position('top')
axins.tick_params(axis='x', labelbottom=False, labeltop=True)

ax.legend(fontsize = 16, loc = 'upper left')
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
#plt.tight_layout()

plt.savefig(SPECTRUM_DIR + "/EE_Spectrum.png", bbox_inches = "tight")
#%%
plt.figure(9, figsize = (12, 6))

plt.title('CMB temperature-polarization cross spectrum', fontsize = 25, pad = 10)
plt.xlabel(r'$\ell$', fontsize = 25)
plt.ylabel(r'$\ell(\ell+1)C_\ell^{TE}/\,2\pi\ (\mu K^2)$', fontsize = 23)
error_array_TE = np.array([TE_ell_errors['ErrorDown'], TE_ell_errors['ErrorUp']])


plt.plot(valuesdic1['ell'], [i*((j+2)*(j+1)*j*(j-1))**(1/2) for i,j in zip(valuesdic1['CellTE'],valuesdic1['ell'])])
plt.errorbar(TE_ell_errors['ell'], TE_ell_errors['Cell'], yerr = error_array_TE, fmt = ".")
plt.xlim([2,2000])
xticks = [500, 1000, 1500, 2000]
plt.xticks(xticks, labels=[str(x) for x in xticks])
plt.tick_params(axis='both', labelsize=17)

plt.tight_layout()

plt.savefig(SPECTRUM_DIR + "/TE_Spectrum.png", bbox_inches = "tight")
#%%

#### PLOT OF MATTER POWER SPECTRUM ####

plt.figure(3, figsize = (9, 6.5))

plt.title('Matter Power Spectrum', fontsize = 25, pad = 10)
plt.xlabel(r'$k\,(h/\rm{Mpc})$', fontsize = 20)
plt.ylabel(r'$P(k)\,(\rm{Mpc}/h)^3$', fontsize = 20)

#
plt.yscale('log')
plt.xscale('log')
plt.xlim([1e-3, max(valuesdic2['k']) * Mpc/h])
plt.ylim([2.5e2, 4e4])
plt.plot([i * Mpc/h for i in valuesdic2['k']], [i * (h/Mpc)**3 for i in valuesdic2['pofk']], label = 'Numerical prediction')

#Plot error bars
plt.errorbar(WMAP_ACT_errors['k'], WMAP_ACT_errors['p_of_k'], yerr = WMAP_ACT_errors['Error'], fmt = ".", color= 'orange', label ='WMAP + ACT')
plt.errorbar(SDSS_DR7_LRG_errors['k'], SDSS_DR7_LRG_errors['p_of_k'], yerr = SDSS_DR7_LRG_errors['Error'], fmt = ".", color= 'black', label ='SDSS survey (DR7 LRG)')
#plt.errorbar(Lyman_alpha_errors['k'], Lyman_alpha_errors['p_of_k'], yerr = Lyman_alpha_errors['Error'], fmt = ".", color= 'green', label =r'$\text{Lyman}-\alpha forest')


#Plot vertical line of equality scale
plt.plot([k_eq * Mpc/h for i in range(100)], np.linspace(100, max([i * (h/Mpc)**3 for i in valuesdic2['pofk']]), 100), linestyle='--', color='red')

plt.tick_params(axis='both', labelsize=15)
plt.legend(fontsize = 17)
plt.tight_layout()
plt.savefig(SPECTRUM_DIR + "/Matter_power_Spectrum.png", bbox_inches = "tight")
#%%
plt.figure(4)

eta0 = 14.19163177767912 * 1000 * Mpc

fig, axes = plt.subplots(nrows=2, figsize=(8.5, 13.5), sharex=False)

for i in axes:
    i.set_xlabel(r'$k\eta_0$', fontsize = 18)
    i.tick_params(axis='both', labelsize=12)
    i.set_xlim([0,200])

axes[0].set_ylabel(r'$\Theta_{10}(k)$', fontsize = 18)
axes[0].plot([i * eta0 for i in valuesdic2['k']], valuesdic2['Theta10'], label = r'$\ell=10$')
axes[0].yaxis.labelpad = 20

axes[1].set_ylabel(r'$\Theta_{100}(k)$', fontsize = 18)
axes[1].plot([i * eta0 for i in valuesdic2['k']], valuesdic2['Theta100'], label = r'$\ell=100$')
axes[1].yaxis.labelpad = 20

plt.savefig(SPECTRUM_DIR + "/Transfer_functions.png", bbox_inches = "tight")
#%%
eta0 = 14.19163177767912 * 1000 * Mpc

fig, axes = plt.subplots(nrows=2, figsize=(8.5, 13.5), sharex=False)

for i in axes:
    i.set_xlabel(r'$k\eta_0$', fontsize = 18)
    i.tick_params(axis='both', labelsize=12)
    i.set_xlim([0,200])
 
    
axes[0].set_ylabel(r'$|\Theta_{10}|^2(k)/k\,(10^{-7} \eta_0)$', fontsize = 18)
axes[0].plot([i * eta0 for i in valuesdic2['k']], [i**2*10**7/(j*eta0) for i,j in zip(valuesdic2['Theta10'], valuesdic2['k'])], label = r'$\ell=100$')
axes[0].yaxis.labelpad = 20

axes[1].set_ylabel(r'$|\Theta_{100}|^2(k)/k\,(10^{-7} \eta_0)$', fontsize = 18)
axes[1].yaxis.labelpad = 20
axes[1].plot([i * eta0 for i in valuesdic2['k']], [i**2*10**7/(j*eta0) for i,j in zip(valuesdic2['Theta100'], valuesdic2['k'])], label = r'$\ell=100$')
#plt.yscale('log')

#plt.plot([i * eta0 for i in valuesdic2['k']], [i**2*h*H0_over_h/j for i,j in zip(valuesdic2['Theta10'], valuesdic2['k'])], label = r'$\ell=10$')

#plt.plot([i * eta0 for i in valuesdic2['k']], [i**2/j for i,j in zip(valuesdic2['Theta300'], valuesdic2['k'])], label = r'$\ell=300$')

plt.savefig(SPECTRUM_DIR + "/C_ell_integrands.png", bbox_inches = "tight")

