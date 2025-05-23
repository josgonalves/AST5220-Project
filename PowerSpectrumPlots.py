# -*- coding: utf-8 -*-
"""
Created on Fri May  9 09:19:15 2025

@author: User
"""

#imports
import numpy as np
import matplotlib.pyplot as plt
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
    valuesdic2 = {'k': [float(line.split(' ')[0]) for line in lines], 'pofk': [float(line.split(' ')[1]) for line in lines], 'Theta10': [float(line.split(' ')[9]) for line in lines], 'Theta100': [float(line.split(' ')[21]) for line in lines], 'Theta300': [float(line.split(' ')[30]) for line in lines]}
"""
Define the directories for plots. If they don't yet exist, create them.
"""

PROJECT_ROOT_DIR = "Project_plots"
if not os.path.exists(PROJECT_ROOT_DIR):
    os.makedirs(PROJECT_ROOT_DIR)

SPECTRUM_DIR = os.path.join(PROJECT_ROOT_DIR, "PowerSpectra")
if not os.path.exists(SPECTRUM_DIR):
    os.makedirs(SPECTRUM_DIR)


#%%
plt.figure(0, figsize = (11, 6))

plt.xlabel(r'$\ell$', fontsize = 25)
plt.title(r'$\frac{\ell(\ell+1)}{2\pi}C_\ell^{TT} (\mu K^2)$', fontsize = 25, pad = 10)

plt.plot(valuesdic1['ell'], valuesdic1['CellTT'])
plt.xscale('log')

plt.xlim([2,2000])
xticks = [2, 10, 50, 500, 1000, 2000]
plt.xticks(xticks, labels=[str(x) for x in xticks])
plt.tick_params(axis='both', labelsize=15)

plt.tight_layout()

#%%
plt.figure(1, figsize = (11, 7))

plt.xlabel(r'$\ell$', fontsize = 25)
plt.title(r'$\frac{\ell(\ell+1)}{2\pi}C_\ell^{EE} (\mu K^2)$', fontsize = 25, pad = 10)

plt.plot(valuesdic1['ell'], [i*(j+2)*(j+1)*j*(j-1) for i,j in zip(valuesdic1['CellEE'], valuesdic1['ell'])])
plt.xlim([2,2000])
xticks = [500, 1000, 1500, 2000]
plt.xticks(xticks, labels=[str(x) for x in xticks])
plt.tick_params(axis='both', labelsize=17)

plt.tight_layout()

#%%
plt.figure(2, figsize = (12, 6))

plt.xlabel(r'$\ell$', fontsize = 25)
plt.title(r'$\frac{\ell(\ell+1)}{2\pi}C_\ell^{TE} (\mu K^2)$', fontsize = 25, pad = 10)

plt.plot(valuesdic1['ell'], [i*((j+2)*(j+1)*j*(j-1))**(1/2) for i,j in zip(valuesdic1['CellTE'],valuesdic1['ell'])])
plt.xlim([2,2000])
xticks = [500, 1000, 1500, 2000]
plt.xticks(xticks, labels=[str(x) for x in xticks])
plt.tick_params(axis='both', labelsize=17)

plt.tight_layout()
#%%
plt.figure(3)

plt.xlabel(r'$k\,(h/\rm{Mpc})$', fontsize = 12)
plt.title(r'$P(k)\,(\rm{Mpc}/h)^3$', fontsize = 15)

plt.yscale('log')
plt.xscale('log')
plt.xlim([min(valuesdic2['k']) * Mpc/h, max(valuesdic2['k']) * Mpc/h])
plt.plot([i * Mpc/h for i in valuesdic2['k']], [i * (h/Mpc)**3 for i in valuesdic2['pofk']])
#%%
plt.figure(4)

eta0 = 14.19163177767912 * 1000 * Mpc

plt.xlabel(r'$k\eta_0$')
plt.title(r'$\Theta_\ell(k)$', fontsize = 15)
plt.xlim([0,500])

plt.plot([i * eta0 for i in valuesdic2['k']], valuesdic2['Theta10'], label = r'$\ell=10$')
plt.plot([i * eta0 for i in valuesdic2['k']], valuesdic2['Theta100'], label = r'$\ell=100$')
plt.plot([i * eta0 for i in valuesdic2['k']], valuesdic2['Theta300'], label = r'$\ell=300$')
plt.legend()

#%%
plt.figure(4)

eta0 = 14.19163177767912 * 1000 * Mpc

plt.xlabel(r'$k\eta_0$')
plt.title(r'$\Theta_\ell(k)$', fontsize = 15)

#plt.yscale('log')

plt.plot([i * eta0 for i in valuesdic2['k']], [i**2/j for i,j in zip(valuesdic2['Theta10'], valuesdic2['k'])], label = r'$\ell=10$')
plt.plot([i * eta0 for i in valuesdic2['k']], [i**2/j for i,j in zip(valuesdic2['Theta100'], valuesdic2['k'])], label = r'$\ell=100$')
plt.plot([i * eta0 for i in valuesdic2['k']], [i**2/j for i,j in zip(valuesdic2['Theta300'], valuesdic2['k'])], label = r'$\ell=300$')