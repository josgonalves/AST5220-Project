# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 13:54:51 2025

@author: User
"""
#imports

import numpy as np
import matplotlib.pyplot as plt
import os

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





with open('results_supernovafitting.txt', 'r') as file:
    for line in range(200):
        next(file)
    lines = [i.split(' ') for i in file.readlines()]
    """
    directly extract the values from splitting each line with spaces and removing
    those spaces
    put them in np arrays to apply argmin
    """
    def filterspaces(x):
        if x=='':
            return False
        return True
    filteredlines = [list(filter(filterspaces, i)) for i in lines]
    
    chi2_array = np.array([float(i[0]) for i in filteredlines])
    h_array = np.array([float(i[1]) for i in filteredlines])
    omegaM_array = np.array([float(i[2]) for i in filteredlines])
    omegaK_array = np.array([float(i[3]) for i in filteredlines])
    
    #compute omega lambda values from those of omega M and omega K
    omegaL_array = np.full(omegaK_array.shape, 1) - omegaM_array - omegaK_array

best_chi2 = chi2_array[np.argmin(chi2_array)]
best_h = h_array[np.argmin(chi2_array)]
best_omegaM = omegaM_array[np.argmin(chi2_array)]
best_omegaL = omegaL_array[np.argmin(chi2_array)]
print(best_h, best_omegaM, best_omegaL)



accepted_omegaM = omegaM_array[chi2_array < best_chi2 + 3.53]
accepted_omegaL = omegaL_array[chi2_array < best_chi2 + 3.53]
accepted_h = h_array[chi2_array < best_chi2 + 3.53]



#means, standard deviations (for posteriors)
avg_h = np.average(h_array)
avg_omegaM = np.average(omegaM_array)
avg_omegaK = np.average(omegaK_array)
avg_omegaL = np.average(omegaL_array)


d_h = np.std(h_array)
d_omegaM = np.std(omegaM_array)
d_omegaK = np.std(omegaK_array)
d_omegaL = np.std(omegaL_array)

#Planck best fit values

h           = 0.67
OmegaB      = 0.05
OmegaCDM    = 0.267
OmegaK      = 0.0
Neff        = 3.046
TCMB        = 2.7255



H0 = h * H0_over_h
OmegaR = (2 * pow(np.pi,3) * pow(k_b*TCMB,4) * 8 * G)/(90 * pow(hbar,3) * pow(c,5) * pow(H0,2))
OmegaNu = Neff * (7/8) * pow(4/11, 4/3) * OmegaR
OmegaLambda = 1 - (OmegaB + OmegaCDM + OmegaK + OmegaR + OmegaNu)

"""
Gaussian function for fitting the posteriors
"""

def gauss(x, dev, avg):
    return (1/((2*np.pi)**(1/2)*dev))*np.exp(-(1/2)*((x-avg)**2/dev**2))



    
"""
Define the directories for plots. If they don't yet exist, create them.
"""

PROJECT_ROOT_DIR = "Project_plots"
if not os.path.exists(PROJECT_ROOT_DIR):
    os.makedirs(PROJECT_ROOT_DIR)

SUPERNOVA_DIR = os.path.join(PROJECT_ROOT_DIR, "Supernova_fit")
if not os.path.exists(SUPERNOVA_DIR):
    os.makedirs(SUPERNOVA_DIR)

#%%
#1-sigma confidence OmegaM vs. OmegaLambda plot
        
plt.figure(0)
plt.title(r'Constraints on the $\Omega_{\rm M0} \times \Omega_{\Lambda0}$ plane', fontsize = 17)
plt.ylabel(r'$\Omega_{\Lambda0}$', fontsize = 14)
plt.tick_params(axis='both', labelsize=11.5)
plt.xlim([0,1])
plt.ylim([0,1])
plt.xlabel(r'$\Omega_{\rm M0}$', fontsize = 14)
plt.scatter(accepted_omegaM, accepted_omegaL, label = '1$\sigma$ region')
plt.plot(np.linspace(0,1,100),[1-i for i in np.linspace(0,1,100)], '--', color='black', label = 'Flat universe')
plt.scatter(OmegaB+OmegaCDM, OmegaLambda, color='orange', label = 'Fiducial values')

plt.legend(bbox_to_anchor=(1.45, 1), fontsize = 12)
plt.savefig(SUPERNOVA_DIR + "/Omega_M_Lambda_Plane.png", bbox_inches='tight')

#%%
#OmegaM posterior plot
fig, axes = plt.subplots(nrows=3, figsize=(8, 12), sharex=False)

#h posterior plot

axes[0].set_title(r'Posterior for h', fontsize = 17)
counts, bins, patches = axes[0].hist(h_array, bins = 30, density = True, alpha = 0.6)
axes[0].tick_params(axis='both', labelsize=13, pad = 5)
axes[0].set_xlim([0.665,0.725])
#Define the 1-sigma region 
lower_bound = min(accepted_h)
upper_bound = max(accepted_h)

#Color bars within the 1-sigma region
for count, patch, bin_left, bin_right in zip(counts, patches, bins[:-1], bins[1:]):
    #Color the 1-sigma region differently
    #print(bin_left, bin_right)
    if bin_right >= lower_bound >= bin_left:
        axes[0].fill_between(y1 = 0, y2 = count, x = [lower_bound, bin_right], color = 'blue', label='1$\sigma$ region')
        axes[0].fill_between(y1 = 0, y2 = count, x = [bin_left, lower_bound], color = 'gray')
    elif bin_right >= upper_bound >= bin_left:
        axes[0].fill_between(y1 = 0, y2 = count, x = [bin_left, upper_bound], color = 'blue')
        axes[0].fill_between(y1 = 0, y2 = count, x = [upper_bound, bin_right], color = 'gray')
    elif lower_bound <= bin_left <= upper_bound:
        axes[0].fill_between(y1 = 0, y2 = count, x = [bin_left, bin_right], color = 'blue')
    else:
        axes[0].fill_between(y1 = 0, y2 = count, x = [bin_left, bin_right], color = 'gray')  # Different color outside

# Vertical reference lines (planck expectation and 1 sigma confidence)
axes[0].axvline(h, linestyle='--', color='black', label='Planck best fit value')
axes[0].axvline(lower_bound, linestyle='--', color='red')
axes[0].axvline(upper_bound, linestyle='--', color='red')

# Gaussian fit
axes[0].plot(np.sort(h_array), gauss(np.sort(h_array), d_h, avg_h))

# Labels and legend
axes[0].set_xlabel('h', fontsize = 17)



axes[1].set_title(r'Posterior for $\Omega_{\rm M0}$', fontsize = 17)
counts, bins, patches = axes[1].hist(omegaM_array, bins=30, density=True, alpha=0.6)
axes[1].tick_params(axis='both', labelsize=13, pad = 5)
axes[1].set_xlim([0,0.6])

#Define the 1-sigma region 
lower_bound = min(accepted_omegaM)
upper_bound = max(accepted_omegaM)

#Color bars within the 1-sigma region
for count, patch, bin_left, bin_right in zip(counts, patches, bins[:-1], bins[1:]):
    #Color the 1-sigma region differently
    #print(bin_left, bin_right)
    if bin_right >= lower_bound >= bin_left:
        axes[1].fill_between(y1 = 0, y2 = count, x = [lower_bound, bin_right], color = 'blue', label='1$\sigma$ region')
        axes[1].fill_between(y1 = 0, y2 = count, x = [bin_left, lower_bound], color = 'gray')
    elif bin_right >= upper_bound >= bin_left:
        axes[1].fill_between(y1 = 0, y2 = count, x = [bin_left, upper_bound], color = 'blue')
        axes[1].fill_between(y1 = 0, y2 = count, x = [upper_bound, bin_right], color = 'gray')
    elif lower_bound <= bin_left <= upper_bound:
        axes[1].fill_between(y1 = 0, y2 = count, x = [bin_left, bin_right], color = 'blue')
    else:
        axes[1].fill_between(y1 = 0, y2 = count, x = [bin_left, bin_right], color = 'gray')  # Different color outside

# Vertical reference lines (planck expectation and 1 sigma confidence)
axes[1].axvline(OmegaB + OmegaCDM, linestyle='--', color='black', label='Planck best fit value')
axes[1].axvline(lower_bound, linestyle='--', color='red')
axes[1].axvline(upper_bound, linestyle='--', color='red')

# Gaussian fit
axes[1].plot(np.sort(omegaM_array), gauss(np.sort(omegaM_array), d_omegaM, avg_omegaM), label = 'Gaussian fit')

# Labels and legend
axes[1].set_xlabel(r'$\Omega_{\rm M0}$', fontsize = 17)
axes[1].legend(bbox_to_anchor=(1, 1), fontsize = 15)






#OmegaLambda posterior plot

axes[2].set_title(r'Posterior for $\Omega_{\Lambda0}$', fontsize = 17)
counts, bins, patches = axes[2].hist(omegaL_array, bins = 30, density = True, alpha = 0.6)
axes[2].set_xlim([0.1, 1.2])
axes[2].tick_params(axis='both', labelsize=13, pad = 5)

#Define the 1-sigma region 
lower_bound = min(accepted_omegaL)
upper_bound = max(accepted_omegaL)

#Color bars within the 1-sigma region
for count, patch, bin_left, bin_right in zip(counts, patches, bins[:-1], bins[1:]):
    #Color the 1-sigma region differently
    #print(bin_left, bin_right)
    if bin_right >= lower_bound >= bin_left:
        axes[2].fill_between(y1 = 0, y2 = count, x = [lower_bound, bin_right], color = 'blue', label='1$\sigma$ region')
        axes[2].fill_between(y1 = 0, y2 = count, x = [bin_left, lower_bound], color = 'gray')
    elif bin_right >= upper_bound >= bin_left:
        axes[2].fill_between(y1 = 0, y2 = count, x = [bin_left, upper_bound], color = 'blue')
        axes[2].fill_between(y1 = 0, y2 = count, x = [upper_bound, bin_right], color = 'gray')
    elif lower_bound <= bin_left <= upper_bound:
        axes[2].fill_between(y1 = 0, y2 = count, x = [bin_left, bin_right], color = 'blue')
    else:
        axes[2].fill_between(y1 = 0, y2 = count, x = [bin_left, bin_right], color = 'gray')  # Different color outside

# Vertical reference lines (planck expectation and 1 sigma confidence)
axes[2].axvline(OmegaLambda, linestyle='--', color='black', label='Planck best fit value')
axes[2].axvline(lower_bound, linestyle='--', color='red')
axes[2].axvline(upper_bound, linestyle='--', color='red')

# Gaussian fit
axes[2].plot(np.sort(omegaL_array), gauss(np.sort(omegaL_array), d_omegaL, avg_omegaL))

# Labels and legend
axes[2].set_xlabel(r'$\Omega_{\Lambda0}$', fontsize = 17)

plt.subplots_adjust(hspace=0.3)
plt.tight_layout()
plt.savefig(SUPERNOVA_DIR + "/Posterior_distributions.png", bbox_inches='tight')
