# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 12:15:03 2025

@author: User
"""

#imports
import numpy as np
import matplotlib.pyplot as plt
import os



"""
THINGS TO PLOT:

Xe(x), tau(x), -dtaudx(x), ddtauddx(x) 


g(x), dgdx, ddgddx (scaled to the same plot)
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


#opening the cosmology.txt file
with open("recombination.txt", 'r') as file:
    #skips the first two lines of the file
    for line in range(2):
        next(file)

    #makes a list of the remaining lines    
    lines = file.readlines()
    #make a dictionary where the keys are the variables and the values are lists of their corresponding values
    valuesdic = {'x': [float(line.split(' ')[0]) for line in lines[:-6]], 'Xe': [float(line.split(' ')[1]) for line in lines[:-6]],'ne': [float(line.split(' ')[2]) for line in lines[:-6]],
                 'tau': [float(line.split(' ')[3]) for line in lines[:-6]],'dtaudx': [float(line.split(' ')[4]) for line in lines[:-6]],'ddtauddx': [float(line.split(' ')[5]) for line in lines[:-6]],'g_tilde': [float(line.split(' ')[6]) for line in lines[:-6]],'dgdx': [float(line.split(' ')[7]) for line in lines[:-6]],'ddgddx': [float(line.split(' ')[8]) for line in lines[:-6]],
                 'Xe_saha': [float(line.split(' ')[9]) for line in lines[:-6]], 's':[float(line.split(' ')[10]) for line in lines[:-6]], 'Xe_reion': [float(line.split(' ')[11]) for line in lines[:-6]],
                 'tau_reion': [float(line.split(' ')[12]) for line in lines[:-6]],'dtaudx_reion': [float(line.split(' ')[13]) for line in lines[:-6]],'ddtauddx_reion': [float(line.split(' ')[14]) for line in lines[:-6]],'g_tilde_reion': [float(line.split(' ')[15]) for line in lines[:-6]],'dgdx_reion': [float(line.split(' ')[16]) for line in lines[:-6]],'ddgddx_reion': [float(line.split(' ')[17]) for line in lines[:-6]]}
    
    # Fetch relevant values (last scattering, recombination, etc.)
    last_scattering_xtz = [float(i) for i in lines[-5].split(':')[1].split(' ')]
    recombination_xtzT =[float(i) for i in lines[-4].split(':')[1].split(' ')]
    last_scattering_Saha_xtz = [float(i) for i in lines[-3].split(':')[1].split(' ')]
    recombination_Saha_xtzT =[float(i) for i in lines[-2].split(':')[1].split(' ')]
    sound_horizon = float(lines[-1].split(':')[1])
    freeze_out = valuesdic['Xe'][-1]
    
    # Obtain the plateau of tau_reion 
    plateau_index = valuesdic['x'].index(-6.0012) # Use the index of some x value within the plateau
    tau_plateau = valuesdic['tau_reion'][plateau_index]

print("Optical depth plateau: ", tau_plateau)
print("Freeze-out free electron fraction: ", freeze_out)
print("Sound horizon at decoupling (Mpc): ", sound_horizon/Mpc)    
print("Last scattering surface (x, t, z): ", last_scattering_xtz)
print("Recombination (x, t, z, T): ", recombination_xtzT)
print("Last scattering surface Saha (x, t, z): ", last_scattering_Saha_xtz)
print("Recombination Saha (x, t, z): ", recombination_Saha_xtzT)

#Get cosmic time in Kyr and the recombination CMB temperature in eV
print("Previously determined values of t converted to Kyr:")
print(last_scattering_xtz[1]/(1000*yr))
print(recombination_xtzT[1]/(1000*yr))
print(recombination_xtzT[3]*k_b/eV)

print(last_scattering_Saha_xtz[1]/(1000*yr))
print(recombination_Saha_xtzT[1]/(1000*yr))

print(np.exp(2.25)-1)

PROJECT_ROOT_DIR = "Project_plots"
if not os.path.exists(PROJECT_ROOT_DIR):
    os.makedirs(PROJECT_ROOT_DIR)

RECOMBINATION_DIR = os.path.join(PROJECT_ROOT_DIR, "Recombination_history")
if not os.path.exists(RECOMBINATION_DIR):
    os.makedirs(RECOMBINATION_DIR)
#print(valuesdic['Xe_saha'])
#%%

plt.figure(0, figsize = (7,4))
plt.grid()
plt.xlabel('x', fontsize = 14)
plt.title(r'$X_e$', fontsize = 17, pad = 10)
plt.tick_params(axis='y', labelsize=11.5)

plt.xlim([-15,0])
plt.yscale('log')
plt.ylim([min(valuesdic['Xe'])*0.9, max(valuesdic['Xe']) * 1.5])

plt.plot(valuesdic['x'], valuesdic['Xe'], label = "Saha + Peebles", lw = 2)
plt.plot(valuesdic['x'], valuesdic['Xe_reion'], '--', label = "Reionization", alpha = 0.9)
plt.plot(valuesdic['x'], valuesdic['Xe_saha'], '--', color = 'purple', label = "Saha only", alpha = 0.9)
plt.scatter(recombination_xtzT[0], 0.1, color='green', alpha = 0.9, s = 35)

plt.legend(fontsize=11)
plt.savefig(RECOMBINATION_DIR + "/Xe_of_x.png")

#%%

fig, axes = plt.subplots(nrows=2, figsize=(6,9), sharex=False)
for i in axes:
    i.set_xlabel('x', fontsize = 16)
    i.set_xlim([-15,0])
    i.grid()
    i.tick_params(axis='both', labelsize=11)
    i.set_ylim([10**-8,max(valuesdic['tau'])])
    i.set_yscale('log')


axes[0].set_title("No reionization", fontsize = 16, pad = 10)
axes[0].axvline(last_scattering_xtz[0],linestyle = '--', alpha = 0.6, color='red')
axes[0].plot(valuesdic['x'], valuesdic['tau'], label = r"$\tau$")
axes[0].plot(valuesdic['x'], [-i for i in valuesdic['dtaudx']], label = r"-$\tau'$")
axes[0].plot(valuesdic['x'], valuesdic['ddtauddx'], label = r"$\tau''$")

axes[1].set_title("Reionization", fontsize = 16, pad = 10)
axes[1].plot(valuesdic['x'], valuesdic['tau_reion'], label = r"$\tau$")
axes[1].plot(valuesdic['x'], [-i for i in valuesdic['dtaudx_reion']], label = r"-$\tau'$")
axes[1].plot(valuesdic['x'], valuesdic['ddtauddx_reion'], label = r"$\tau''$")
axes[1].legend(fontsize=14)

plt.subplots_adjust(hspace=0.3)
plt.savefig(RECOMBINATION_DIR + "/tau_and_derivs_of_x.png", bbox_inches ='tight')

#%%

plt.figure(3)
plt.grid()

plt.xlabel('x', fontsize = 14)
plt.title(r'$s(x)$ (Mpc)', fontsize = 15, pad = 10)

plt.xlim([-15,0])
plt.yscale('log')

plt.plot(valuesdic['x'], [i/(Mpc) for i in valuesdic['s']])
plt.savefig(RECOMBINATION_DIR + "/s_of_x.png")

#%%

fig, axes = plt.subplots(nrows=3, figsize=(7, 14), sharex=False)


for i in axes:
    i.set_xlabel('x', fontsize = 14)
    i.set_xlim([-15,0])
    i.grid()
    i.tick_params(axis='both', labelsize=12)
    
axes[0].set_title(r"$\tilde{g}$", fontsize = 17, pad = 10)
axes[1].set_title(r"$\tilde{g}'$", fontsize = 17, pad = 10) 
axes[2].set_title(r"$\tilde{g}''$", fontsize = 17, pad = 10) 

axes[0].plot(valuesdic['x'], valuesdic['g_tilde'], color = 'orange', label = "No reionization")
axes[0].plot(valuesdic['x'], valuesdic['g_tilde_reion'], '--', color = 'blue', alpha=0.6, label = "Reionization")
axes[0].axvline(last_scattering_xtz[0], linestyle='-', color='red')

axes[1].plot(valuesdic['x'], valuesdic['dgdx'], color = 'orange', label = "No reionization")
axes[1].plot(valuesdic['x'], valuesdic['dgdx_reion'],'--', color = 'blue', alpha=0.6, label = "Reionization")

axes[2].plot(valuesdic['x'], valuesdic['ddgddx'], color = 'orange', label = "No reionization")
axes[2].plot(valuesdic['x'], valuesdic['ddgddx_reion'],'--', color = 'blue', alpha=0.6, label = "Reionization")
axes[2].legend(fontsize = 14)

plt.subplots_adjust(hspace=0.4)
plt.savefig(RECOMBINATION_DIR + "/g_and_derivs_of_x.png", bbox_inches = "tight")

