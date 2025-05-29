# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 10:14:52 2025

@author: User
"""

#imports
import numpy as np
import matplotlib.pyplot as plt
import os



"""
THINGS TO PLOT:

v_gamma(x), delta_gamma(x), v_b(x), delta_b(x), v_CDM(x), delta_CDM(x)

Theta_2(x), Phi(x) , Phi+Psi(x)


for 3-4 values of k (large, medium and small scales)
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


#opening the perturbations files
with open("perturbations_k0.1.txt", 'r') as file:
   lines = file.readlines()
   
   valuesdic1 = {'x': [float(line.split(' ')[0]) for line in lines[:-4]], 'Theta0': [float(line.split(' ')[1]) for line in lines[:-4]],'Theta1': [float(line.split(' ')[2]) for line in lines[:-4]],
                'Theta2': [float(line.split(' ')[3]) for line in lines[:-4]], 'Phi': [float(line.split(' ')[4]) for line in lines[:-4]],'Psi': [float(line.split(' ')[5]) for line in lines[:-4]],'Pi': [float(line.split(' ')[6]) for line in lines[:-4]],
                'Thetap0': [float(line.split(' ')[7]) for line in lines[:-4]], 'Thetap1': [float(line.split(' ')[8]) for line in lines[:-4]], 'Thetap2': [float(line.split(' ')[9]) for line in lines[:-4]], 'Nu0': [float(line.split(' ')[10]) for line in lines[:-4]], 'Nu1': [float(line.split(' ')[11]) for line in lines[:-4]], 'Nu2': [float(line.split(' ')[12]) for line in lines[:-4]],
                'Source_T_10': [float(line.split(' ')[14]) for line in lines[:-4]], 'Source_T_100': [float(line.split(' ')[15]) for line in lines[:-4]], 'delta_b': [float(line.split(' ')[17]) for line in lines[:-4]], 'delta_CDM': [float(line.split(' ')[18]) for line in lines[:-4]], 'v_b': [float(line.split(' ')[19]) for line in lines[:-4]], 'v_CDM': [float(line.split(' ')[20]) for line in lines[:-4]], 'Source_E': [float(line.split(' ')[-2]) for line in lines[:-4]]}
   matradequalxtz = [float(i) for i in lines[-3].split(':')[1].split(' ')]
   matlambdaequalxtz = [float(i) for i in lines[-2].split(':')[1].split(' ')]
   last_scattering_xtz = [float(i) for i in lines[-1].split(':')[1].split(' ')]
   
with open("perturbations_k0.01.txt", 'r') as file:
   lines = file.readlines()
   
   valuesdic2 = {'x': [float(line.split(' ')[0]) for line in lines[:-4]], 'Theta0': [float(line.split(' ')[1]) for line in lines[:-4]],'Theta1': [float(line.split(' ')[2]) for line in lines[:-4]],
                'Theta2': [float(line.split(' ')[3]) for line in lines[:-4]], 'Phi': [float(line.split(' ')[4]) for line in lines[:-4]],'Psi': [float(line.split(' ')[5]) for line in lines[:-4]],'Pi': [float(line.split(' ')[6]) for line in lines[:-4]],
                'Thetap0': [float(line.split(' ')[7]) for line in lines[:-4]], 'Thetap1': [float(line.split(' ')[8]) for line in lines[:-4]], 'Thetap2': [float(line.split(' ')[9]) for line in lines[:-4]], 'Nu0': [float(line.split(' ')[10]) for line in lines[:-4]], 'Nu1': [float(line.split(' ')[11]) for line in lines[:-4]], 'Nu2': [float(line.split(' ')[12]) for line in lines[:-4]],
                'Source_T_10': [float(line.split(' ')[14]) for line in lines[:-4]], 'Source_T_100': [float(line.split(' ')[15]) for line in lines[:-4]], 'delta_b': [float(line.split(' ')[17]) for line in lines[:-4]], 'delta_CDM': [float(line.split(' ')[18]) for line in lines[:-4]], 'v_b': [float(line.split(' ')[19]) for line in lines[:-4]], 'v_CDM': [float(line.split(' ')[20]) for line in lines[:-4]], 'Source_E': [float(line.split(' ')[-2]) for line in lines[:-4]]}
with open("perturbations_k0.001.txt", 'r') as file:
   lines = file.readlines()
   
   valuesdic3 = {'x': [float(line.split(' ')[0]) for line in lines[:-4]], 'Theta0': [float(line.split(' ')[1]) for line in lines[:-4]],'Theta1': [float(line.split(' ')[2]) for line in lines[:-4]],
                'Theta2': [float(line.split(' ')[3]) for line in lines[:-4]], 'Phi': [float(line.split(' ')[4]) for line in lines[:-4]],'Psi': [float(line.split(' ')[5]) for line in lines[:-4]],'Pi': [float(line.split(' ')[6]) for line in lines[:-4]],
                'Thetap0': [float(line.split(' ')[7]) for line in lines[:-4]], 'Thetap1': [float(line.split(' ')[8]) for line in lines[:-4]], 'Thetap2': [float(line.split(' ')[9]) for line in lines[:-4]], 'Nu0': [float(line.split(' ')[10]) for line in lines[:-4]], 'Nu1': [float(line.split(' ')[11]) for line in lines[:-4]], 'Nu2': [float(line.split(' ')[12]) for line in lines[:-4]],
                'Source_T_10': [float(line.split(' ')[14]) for line in lines[:-4]], 'Source_T_100': [float(line.split(' ')[15]) for line in lines[:-4]], 'delta_b': [float(line.split(' ')[17]) for line in lines[:-4]], 'delta_CDM': [float(line.split(' ')[18]) for line in lines[:-4]], 'v_b': [float(line.split(' ')[19]) for line in lines[:-4]], 'v_CDM': [float(line.split(' ')[20]) for line in lines[:-4]], 'Source_E': [float(line.split(' ')[-2]) for line in lines[:-4]]}
   
"""
Define the directories for plots. If they don't yet exist, create them.
"""
print(valuesdic1['x'].index(-6.00174))
PROJECT_ROOT_DIR = "Project_plots"
if not os.path.exists(PROJECT_ROOT_DIR):
    os.makedirs(PROJECT_ROOT_DIR)

PERTURBATIONS_DIR = os.path.join(PROJECT_ROOT_DIR, "Perturbations")
if not os.path.exists(PERTURBATIONS_DIR):
    os.makedirs(PERTURBATIONS_DIR)
#%%

#plot density perturbations for photons and neutrinos

fig, axes = plt.subplots(nrows=2, figsize=(8.5, 14), sharex=False)

for i in axes:
    i.set_xlabel('x', fontsize = 14)
    i.set_xlim([-13,0])
    i.grid()
    i.tick_params(axis='both', labelsize=12)


axes[0].set_title(r'$ \delta_\gamma= 4\Theta_0$', fontsize = 20, pad = 10)
axes[0].plot(valuesdic1['x'], [4*i for i in valuesdic1['Theta0']], color = 'green', label = 'k = 0.1/Mpc')
axes[0].plot(valuesdic2['x'], [4*i for i in valuesdic2['Theta0']], color = 'orange', label = 'k = 0.01/Mpc')
axes[0].plot(valuesdic3['x'], [4*i for i in valuesdic3['Theta0']], color = 'blue', label = 'k = 0.001/Mpc')
axes[0].axvline(last_scattering_xtz[0], linestyle='--', color='black')

axes[1].set_title(r'$ \delta_\nu= 4\mathcal{N}_0$', fontsize = 20, pad = 10)
axes[1].plot(valuesdic1['x'], [4*i for i in valuesdic1['Nu0']], color = 'green', label = 'k = 0.1/Mpc')
axes[1].plot(valuesdic2['x'], [4*i for i in valuesdic2['Nu0']], color = 'orange', label = 'k = 0.01/Mpc')
axes[1].plot(valuesdic3['x'], [4*i for i in valuesdic3['Nu0']], color = 'blue', label = 'k = 0.001/Mpc')
axes[1].axvline(last_scattering_xtz[0], linestyle='--', color='black')

axes[1].legend(fontsize = 14)

plt.subplots_adjust(hspace=0.4)

plt.savefig(PERTURBATIONS_DIR + "/photon_neutrino_density.png", bbox_inches = "tight")
#%%

#plot density perturbations for baryons and CDM

plt.figure(2, figsize = (7,4))
plt.grid()
plt.xlabel('x', fontsize = 14)
plt.title(r'$ \delta_{\rm CDM}, \delta_{b}$', fontsize = 17, pad = 10)
plt.tick_params(axis='y', labelsize=11.5)

plt.xlim([-15,0])
plt.yscale('log')
plt.ylim([0.1, 1e5])

plt.axvline(last_scattering_xtz[0], linestyle='--', color='black')
plt.plot(valuesdic1['x'], valuesdic1['delta_b'], linestyle = '--', color='green')
plt.plot(valuesdic1['x'], valuesdic1['delta_CDM'], color='green', label = 'k = 0.1/Mpc')
plt.plot(valuesdic2['x'], valuesdic2['delta_b'], linestyle = '--', color='orange')
plt.plot(valuesdic2['x'], valuesdic2['delta_CDM'], color='orange', label = 'k = 0.01/Mpc')
plt.plot(valuesdic3['x'], valuesdic3['delta_b'], linestyle = '--', color='blue')
plt.plot(valuesdic3['x'], valuesdic3['delta_CDM'], color='blue', label = 'k = 0.001/Mpc')
plt.legend()

plt.savefig(PERTURBATIONS_DIR + "/baryon_CDM_density.png", bbox_inches = "tight")
#%%

#plot velocity profiles for photons and neutrinos

fig, axes = plt.subplots(nrows=2, figsize=(8.5, 14), sharex=False)

for i in axes:
    i.set_xlabel('x', fontsize = 14)
    i.set_xlim([-13,0])
    i.grid()
    i.tick_params(axis='both', labelsize=12)


axes[0].set_title(r'$v_\gamma= -3\Theta_1$', fontsize = 20, pad = 10)
axes[0].plot(valuesdic1['x'], [-3*i for i in valuesdic1['Theta1']], color = 'green', label = 'k = 0.1/Mpc')
axes[0].plot(valuesdic2['x'], [-3*i for i in valuesdic2['Theta1']], color = 'orange', label = 'k = 0.01/Mpc')
axes[0].plot(valuesdic3['x'], [-3*i for i in valuesdic3['Theta1']], color = 'blue', label = 'k = 0.001/Mpc')
axes[0].axvline(last_scattering_xtz[0], linestyle='--', color='black')


axes[1].set_title(r'$ v_\nu= -3\mathcal{N}_1$', fontsize = 20, pad = 10)
axes[1].plot(valuesdic1['x'], [-3*i for i in valuesdic1['Nu1']], color = 'green', label = 'k = 0.1/Mpc')
axes[1].plot(valuesdic2['x'], [-3*i for i in valuesdic2['Nu1']], color = 'orange', label = 'k = 0.01/Mpc')
axes[1].plot(valuesdic3['x'], [-3*i for i in valuesdic3['Nu1']], color = 'blue', label = 'k = 0.001/Mpc')
axes[1].axvline(last_scattering_xtz[0], linestyle='--', color='black')

axes[1].legend(fontsize = 14)

plt.subplots_adjust(hspace=0.4)

plt.savefig(PERTURBATIONS_DIR + "/photon_neutrino_velocity.png", bbox_inches = "tight")
#%%

#plot velocity profiles for baryons and CDM

plt.figure(5, figsize = (7,4))
plt.grid()
plt.xlabel('x', fontsize = 14)
plt.title(r'$ v_{\rm CDM}, v_{b}$', fontsize = 17, pad = 10)
plt.tick_params(axis='y', labelsize=11.5)

plt.xlim([-15,0])
plt.yscale('log')

plt.axvline(last_scattering_xtz[0], linestyle='--', color='black')


plt.plot(valuesdic1['x'], valuesdic1['v_b'], linestyle = '--', color='green')
plt.plot(valuesdic1['x'], valuesdic1['v_CDM'], color='green', label = 'k = 0.1/Mpc')
plt.plot(valuesdic2['x'], valuesdic2['v_b'], linestyle = '--', color='orange')
plt.plot(valuesdic2['x'], valuesdic2['v_CDM'], color='orange', label = 'k = 0.01/Mpc')
plt.plot(valuesdic3['x'], valuesdic3['v_b'], linestyle = '--', color='blue')
plt.plot(valuesdic3['x'], valuesdic3['v_CDM'], color='blue', label = 'k = 0.001/Mpc')
plt.legend()

plt.savefig(PERTURBATIONS_DIR + "/baryon_CDM_velocity.png", bbox_inches = "tight")
#%%

# Plot Phi
plt.figure(6, figsize = (7,4))
plt.grid()
plt.xlabel('x', fontsize = 14)
plt.title(r'$ \Phi$', fontsize = 17, pad = 10)
plt.tick_params(axis='y', labelsize=11.5)

plt.xlim([-15,0])

plt.axvline(matradequalxtz[0], linestyle='--', color='black')
plt.axvline(matlambdaequalxtz[0], linestyle='--', color='red')
plt.plot(valuesdic1['x'], valuesdic1['Phi'], color = 'green', label = 'k = 0.1/Mpc')
plt.plot(valuesdic2['x'], valuesdic2['Phi'], color = 'orange', label = 'k = 0.01/Mpc')
plt.plot(valuesdic3['x'], valuesdic3['Phi'], color = 'blue', label = 'k = 0.001/Mpc')
plt.legend()

plt.savefig(PERTURBATIONS_DIR + "/Phi.png", bbox_inches = "tight")
#%%

# Plot Phi + Psi
plt.figure(7, figsize = (7,4))
plt.grid()
plt.xlabel('x', fontsize = 14)
plt.title(r'$ \Phi + \Psi$', fontsize = 17, pad = 10)
plt.tick_params(axis='y', labelsize=11.5)

plt.xlim([-15,0])

plt.axvline(matradequalxtz[0], linestyle='--', color='black')
plt.axvline(matlambdaequalxtz[0], linestyle='--', color='red')
plt.plot(valuesdic1['x'], [i+j for i,j in zip (valuesdic1['Phi'], valuesdic1['Psi'])], color = 'green', label = 'k = 0.1/Mpc')
plt.plot(valuesdic2['x'], [i+j for i,j in zip (valuesdic2['Phi'], valuesdic2['Psi'])], color = 'orange', label = 'k = 0.01/Mpc')
plt.plot(valuesdic3['x'], [i+j for i,j in zip (valuesdic3['Phi'], valuesdic3['Psi'])], color = 'blue', label = 'k = 0.001/Mpc')
plt.legend()

plt.savefig(PERTURBATIONS_DIR + "/Anisotropic_stress.png", bbox_inches = "tight")
#%%

# Plot polarization multipoles

fig, axes = plt.subplots(nrows=3, figsize=(8.5, 21), sharex=False)


for i in axes:
    i.set_xlabel('x', fontsize = 14)
    i.set_xlim([-12,0])
    i.grid()
    i.tick_params(axis='both', labelsize=12)
    
axes[0].set_title(r"$\Theta_{P0}$", fontsize = 20, pad = 10)
axes[1].set_title(r"$\Theta_{P1}$", fontsize = 20, pad = 10) 
axes[2].set_title(r"$\Theta_{P2}$", fontsize = 20, pad = 10) 

axes[0].plot(valuesdic1['x'], valuesdic1['Thetap0'], color = 'green', label = 'k = 0.1/Mpc', lw = 2)
axes[0].plot(valuesdic2['x'], valuesdic2['Thetap0'], color = 'orange', label = 'k = 0.01/Mpc', lw = 2)
axes[0].plot(valuesdic3['x'], valuesdic3['Thetap0'], color = 'blue', label = 'k = 0.001/Mpc', lw = 2)
axes[0].axvline(last_scattering_xtz[0], linestyle='--', color='black')

axes[1].plot(valuesdic1['x'], valuesdic1['Thetap1'], color = 'green', label = 'k = 0.1/Mpc', lw = 2)
axes[1].plot(valuesdic2['x'], valuesdic2['Thetap1'], color = 'orange', label = 'k = 0.01/Mpc', lw = 2)
axes[1].plot(valuesdic3['x'], valuesdic3['Thetap1'], color = 'blue', label = 'k = 0.001/Mpc', lw = 2)
axes[1].axvline(last_scattering_xtz[0], linestyle='--', color='black')


axes[2].plot(valuesdic1['x'], valuesdic1['Thetap2'], color = 'green', label = 'k = 0.1/Mpc', lw = 2)
axes[2].plot(valuesdic2['x'], valuesdic2['Thetap2'], color = 'orange', label = 'k = 0.01/Mpc', lw = 2)
axes[2].plot(valuesdic3['x'], valuesdic3['Thetap2'], color = 'blue', label = 'k = 0.001/Mpc', lw = 2)
axes[2].axvline(last_scattering_xtz[0], linestyle='--', color='black')

axes[2].legend(fontsize = 14)

plt.subplots_adjust(hspace=0.4)

plt.savefig(PERTURBATIONS_DIR + "/Polarization_multipoles.png", bbox_inches = "tight")
#%%

#plot quadrupole for photons and neutrinos

fig, axes = plt.subplots(nrows=2, figsize=(8.5, 14), sharex=False)

for i in axes:
    i.set_xlabel('x', fontsize = 14)
    i.set_xlim([-13,0])
    i.grid()
    i.tick_params(axis='both', labelsize=12)


axes[0].set_title(r'$\Theta_2$', fontsize = 20, pad = 10)
axes[0].plot(valuesdic1['x'], valuesdic1['Theta2'], color = 'green', label = 'k = 0.1/Mpc')
axes[0].plot(valuesdic2['x'], valuesdic2['Theta2'], color = 'orange', label = 'k = 0.01/Mpc')
axes[0].plot(valuesdic3['x'], valuesdic3['Theta2'], color = 'blue', label = 'k = 0.001/Mpc')
axes[0].axvline(last_scattering_xtz[0], linestyle='--', color='black')


axes[1].set_title(r'$\mathcal{N}_2$', fontsize = 20, pad = 10)
axes[1].plot(valuesdic1['x'], valuesdic1['Nu2'], color = 'green', label = 'k = 0.1/Mpc')
axes[1].plot(valuesdic2['x'], valuesdic2['Nu2'], color = 'orange', label = 'k = 0.01/Mpc')
axes[1].plot(valuesdic3['x'], valuesdic3['Nu2'], color = 'blue', label = 'k = 0.001/Mpc')
axes[1].axvline(last_scattering_xtz[0], linestyle='--', color='black')

axes[1].legend(fontsize = 14)

plt.subplots_adjust(hspace=0.4)

plt.savefig(PERTURBATIONS_DIR + "/photon_neutrino_quadrupole.png", bbox_inches = "tight")

#%%

#plot overdensity and velocity at k=0.1 for photons and baryons

fig, axes = plt.subplots(nrows=2, figsize=(8.5, 14), sharex=False)

for i in axes:
    i.set_xlabel('x', fontsize = 14)
    i.set_xlim([-13,last_scattering_xtz[0]+1])
    i.grid()
    i.tick_params(axis='both', labelsize=12)


axes[0].set_title(r'$ \delta_{\gamma}, \delta_{b}$', fontsize = 25, pad = 10)
axes[0].plot(valuesdic1['x'][:3697], [4*i for i in valuesdic1['Theta0']][:3697], color = 'green', label = r'$\delta_{\gamma}$')
axes[0].plot(valuesdic1['x'][:3697], valuesdic1['delta_b'][:3697], color = 'orange', label = r'$\delta_{b}$')
axes[0].axvline(last_scattering_xtz[0], linestyle='--', color='black')


axes[1].set_title(r'$v_{\gamma}, v_{b}$', fontsize = 25, pad = 10)
axes[1].plot(valuesdic1['x'][:3697], [-3*i for i in valuesdic1['Theta1']][:3697], color = 'green', label = r'$v_{\gamma}$')
axes[1].plot(valuesdic1['x'][:3697], valuesdic1['v_b'][:3697], color = 'orange', label = r'$v_{b}$')
axes[1].axvline(last_scattering_xtz[0], linestyle='--', color='black')

axes[0].legend(fontsize = 20)
axes[1].legend(fontsize = 20)

plt.subplots_adjust(hspace=0.4)

plt.savefig(PERTURBATIONS_DIR + "/photon_baryon_k01.png", bbox_inches = "tight")

#%%

# Plot Pi
plt.figure(18, figsize = (7,4))
plt.grid()
plt.xlabel('x', fontsize = 14)
plt.title(r'$ \Pi$', fontsize = 17, pad = 10)
plt.tick_params(axis='y', labelsize=11.5)

plt.xlim([-15,0])

plt.axvline(matradequalxtz[0], linestyle='--', color='black')
plt.axvline(matlambdaequalxtz[0], linestyle='--', color='red')
plt.plot(valuesdic1['x'], [i+j for i,j in zip(valuesdic1['Theta0'],valuesdic1['Psi'])], color = 'green', label = 'k = 0.1/Mpc')
plt.plot(valuesdic2['x'], [i+j for i,j in zip(valuesdic2['Theta0'],valuesdic2['Psi'])], color = 'orange', label = 'k = 0.01/Mpc')
plt.plot(valuesdic3['x'], [i+j for i,j in zip(valuesdic3['Theta0'],valuesdic3['Psi'])], color = 'blue', label = 'k = 0.001/Mpc')
plt.legend()
#%%

# Plot Pi
plt.figure(19, figsize = (7,4))
plt.grid()
plt.xlabel('x', fontsize = 14)
plt.title(r'$ \Pi$', fontsize = 17, pad = 10)
plt.tick_params(axis='y', labelsize=11.5)

plt.xlim([-15,-0.1])

plt.axvline(matradequalxtz[0], linestyle='--', color='black')
plt.axvline(matlambdaequalxtz[0], linestyle='--', color='red')
#plt.plot(valuesdic1['x'], valuesdic1['Source_T_10'], color = 'green', label = 'k = 0.1/Mpc')
plt.plot(valuesdic2['x'], valuesdic2['Source_T_100'], color = 'orange', label = 'k = 0.01/Mpc')
#plt.plot(valuesdic3['x'], valuesdic3['Source_T_10'], color = 'blue', label = 'k = 0.001/Mpc')
plt.legend()

