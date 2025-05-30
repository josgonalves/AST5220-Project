
#imports
import numpy as np
import matplotlib.pyplot as plt
import os



"""
THINGS TO PLOT:

1/Hp dHp/dx, 1/Hp ddHp/ddx, eta*Hp/c

Hp(x)

All of the omegas together

t(x), eta(x)/c

Express time in Gyr
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



#opening the cosmology.txt file
with open("cosmology.txt", 'r') as file:
    #skips the first two lines of the file
    for line in range(2):
        next(file)

    #makes a list of the remaining lines    
    lines = file.readlines()
    #make a dictionary where the keys are the variables and the values are lists of their corresponding values
    valuesdic = {'x': [float(line.split(' ')[0]) for line in lines[:-5]], 'eta': [float(line.split(' ')[1]) for line in lines[:-5]], 
    't': [float(line.split(' ')[2]) for line in lines[:-5]],'Hp': [float(line.split(' ')[3]) for line in lines[:-5]], 'dHpdx':[float(line.split(' ')[4]) for line in lines[:-5]], 
    'ddHpddx': [float(line.split(' ')[5]) for line in lines[:-5]], 'dL': [float(line.split(' ')[6]) for line in lines[:-5]],
    'OmegaB': [float(line.split(' ')[7]) for line in lines[:-5]], 'OmegaCDM': [float(line.split(' ')[8]) for line in lines[:-5]],
    'OmegaLambda': [float(line.split(' ')[9]) for line in lines[:-5]], 'OmegaR': [float(line.split(' ')[10]) for line in lines[:-5]],
    'OmegaNu': [float(line.split(' ')[11]) for line in lines[:-5]], 'OmegaK': [float(line.split(' ')[12]) for line in lines[:-5]],
    'etaprime': [float(line.split(' ')[13]) for line in lines[:-5]]}
    
    """
    Current conformal time in gigaparsec, current age of the Universe in gigayears
    """
    currenthorizon = float(lines[-5].split(':')[1].split(' ')[0])/(Mpc*1000)
    currentage = float(lines[-5].split(':')[1].split(' ')[1])/Gyr
    
    """
    Lists that contain the x,t,z values for matter radiation eq.,
    matter dark energy eq. and start of acceleration
    """
    matradequalxtz = [float(i) for i in lines[-3].split(':')[1].split(' ')]
    matlambdaequalxtz = [float(i) for i in lines[-2].split(':')[1].split(' ')]
    acceleration_startxtz = [float(i) for i in lines[-1].split(':')[1].split(' ')]

"""
define the redshift for all values of x as z = 1/a - 1
"""
valuesdic['z'] = [np.exp(-i) - 1 for i in valuesdic['x']]

print(currenthorizon)
print(currentage)

with open("data/supernovadata.txt", 'r') as file:
    #skips the first two lines of the file
    for line in range(1):
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
    supernovadic = {'z': [float(i[0]) for i in filteredlines], 'dL':[float(i[1]) for i in filteredlines], 'sigma':[float(i[2]) for i in filteredlines]}
    supernovadic['x'] =[np.log(1/(1+i)) for i in supernovadic['z']]

"""
Sum of all density parameters, 
OmegaM = OmegaB + OmegaCDM,
OmegaRel = OmegaR + OmegaNu
"""
    
Totalomega = [i+j+k+l+m+n for i,j,k,l,m,n in zip(valuesdic['OmegaB'], valuesdic['OmegaCDM'], valuesdic['OmegaR'], valuesdic['OmegaNu'], valuesdic['OmegaK'], valuesdic['OmegaLambda'])]
OmegaM = [i+j for i,j in zip(valuesdic['OmegaB'], valuesdic['OmegaCDM'])]
OmegaRel = [i+j for i,j in zip(valuesdic['OmegaR'], valuesdic['OmegaNu'])]

"""
Scale factors a = exp(x)
Redshift 1/a = 1+z then z = 1/a -1
"""
a = [np.exp(i) for i in valuesdic['x']] 
z = [1/i - 1 for i in a]


"""
Define the directories for plots. If they don't yet exist, create them.
"""

PROJECT_ROOT_DIR = "Project_plots"
if not os.path.exists(PROJECT_ROOT_DIR):
    os.makedirs(PROJECT_ROOT_DIR)

BACKGROUND_DIR = os.path.join(PROJECT_ROOT_DIR, "Background_cosmology")
if not os.path.exists(BACKGROUND_DIR):
    os.makedirs(BACKGROUND_DIR)

SUPERNOVA_DIR = os.path.join(PROJECT_ROOT_DIR, "Supernova_fit")
if not os.path.exists(SUPERNOVA_DIR):
    os.makedirs(SUPERNOVA_DIR)


"""
A small note:
    
    For the plots of the evolution of the cosmological parameters, Hp'/Hp and Hp''/Hp, to better visualize
    the domination of Lambda, x is made to run all the way to 5. However, all
    other plots only run up to x=0, which is closest to the 1126th element of each of the
    lists I made.
"""

print(matradequalxtz)
print(matlambdaequalxtz)
print(acceleration_startxtz)

print(matradequalxtz[1]/Gyr)

#%%
#Plot of all of the cosmological params summed over x
plt.figure(0)
plt.grid()
plt.xlabel('x', fontsize = 14)
plt.title(r'$\sum\, \Omega_i$', fontsize = 15, pad = 10)

plt.xlim([-15,0])
plt.ylim([0,2])

plt.plot(valuesdic['x'][:1126], Totalomega[:1126])
plt.savefig(BACKGROUND_DIR + "/Sum_of_cosmo_params.png")


#%%

#Plot of Hp(x) in 100 Km/s/Mpc

plt.figure(2)
plt.title(r'$\mathcal{H}\,\left(\frac{100\,km/s}{Mpc}\right)$', fontsize = 15, pad = 10)
plt.grid()
plt.xlabel("x", fontsize = 14)
plt.tick_params(axis='y', labelsize=11.5)
plt.xlim([-15,0])

rad_domination = np.linspace(-15, matradequalxtz[0], 1000)
mat_domination = np.linspace(matradequalxtz[0], matlambdaequalxtz[0], 1000)



plt.yscale('log')
plt.plot(valuesdic['x'][:1126], [i*Mpc/(100*km) for i in valuesdic['Hp']][:1126])
plt.plot(rad_domination, [np.exp(-i)*H0*(OmegaR+OmegaNu)**(1/2)*Mpc/(100*km) for i in rad_domination], linestyle ='--',alpha = 0.9, label = 'Radiation domination expectation')
plt.plot(mat_domination, [np.exp(-i/2)*H0*(OmegaB+OmegaCDM)**(1/2)*Mpc/(100*km) for i in mat_domination], linestyle ='--', alpha = 0.9, label= 'Matter domination expectation')
plt.legend()
plt.savefig(BACKGROUND_DIR + "/Hp_of_x.png")



#%%

#Plots of eta(x)/c in gigaparsec and t(x) in gigayears
fig1, axes1 = plt.subplots(nrows=2, figsize=(8, 12), sharex=False)
axes1[0].grid()
axes1[0].set_xlabel('x', fontsize = 17)
axes1[0].set_title(r'$\eta$ (Gpc)', fontsize = 18, pad = 5)
axes1[0].tick_params(axis='both', labelsize=13, pad = 10)

axes1[0].set_yscale('log')
axes1[0].set_xlim([-15,0])

axes1[0].plot(valuesdic['x'][:1126], [i/(1000*Mpc) for i in valuesdic['eta']][:1126])

axes1[0].set_ylim([min([i/(1000*Mpc) for i in valuesdic['eta'][:1126]]), 1.5*max([i/(1000*Mpc) for i in valuesdic['eta'][:1126]])])

axes1[1].set_title('t (Gyr)', fontsize = 18, pad = 10)
axes1[1].grid()
axes1[1].set_xlabel("x", fontsize = 17)
axes1[1].tick_params(axis='both', labelsize=13, pad = 10)

axes1[1].set_yscale('log')
axes1[1].set_xlim([-15,0])
axes1[1].set_ylim([min([i/Gyr for i in valuesdic['t'][:1126]]), 1.5*max([i/Gyr for i in valuesdic['t'][:1126]])])


axes1[1].plot(valuesdic['x'][:1126], [i/Gyr for i in valuesdic['t'][:1126]])

plt.subplots_adjust(hspace=0.3)
plt.savefig(BACKGROUND_DIR + "/t_and_eta_of_x.png")



#%%

#plot H'p/Hp

fig, axes = plt.subplots(nrows=2, figsize=(8, 12), sharex=False)

axes[0].set_title(r"$\mathcal{H}\,'/\mathcal{H}$", fontsize = 18, pad = 10)
axes[0].grid()
axes[0].set_xlabel('x', fontsize = 16)
axes[0].tick_params(axis='y', labelsize=13)

axes[0].set_xlim([-15,5])

axes[0].plot(valuesdic['x'],[i/j for i,j in zip(valuesdic['dHpdx'], valuesdic['Hp'])])




#plot H''p/Hp


axes[1].set_title(r"$\mathcal{H}\,''/\mathcal{H}$", fontsize = 18, pad = 10)
axes[1].grid()
axes[1].set_xlabel('x', fontsize = 16)
axes[1].tick_params(axis='y', labelsize=13)
axes[1].set_xlim([-15,5])

axes[1].plot(valuesdic['x'],[i/j for i,j in zip(valuesdic['ddHpddx'], valuesdic['Hp'])])

plt.subplots_adjust(hspace=0.3)
plt.savefig(BACKGROUND_DIR + "/Hp_numerical_tests.png")


#%%


plt.figure(6)
plt.title('dL/z (Gpc)', fontsize = 17, pad = 10)
plt.grid()
plt.xlabel('x', fontsize = 17)
plt.tick_params(axis='both', labelsize=13, pad = 5)

plt.xlim([-1,0])

plt.plot(valuesdic['x'][1049:1126], [i/(1000*Mpc*j) for i,j in zip(valuesdic['dL'],valuesdic['z'])][1049:1126], label = 'Fiducial prediction')
plt.scatter(supernovadic['x'], [i/j for i,j in zip(supernovadic['dL'], supernovadic['z'])], s = 15, color = 'orange', label = 'Observations')
plt.errorbar(supernovadic['x'], [i/j for i,j in zip(supernovadic['dL'], supernovadic['z'])], fmt = 'none' ,yerr = [i/j for i,j in zip(supernovadic['sigma'], supernovadic['z'])], color = 'orange')


plt.legend(loc='upper left', fontsize=11, bbox_to_anchor=(1, 1))
plt.savefig(SUPERNOVA_DIR + "/Luminosity_distance_comparison.png", bbox_inches='tight')

#%%


plt.figure(7)
plt.title("$\Omega_i$", fontsize = 17)
plt.grid()
plt.xlabel('z', fontsize = 17)
plt.tick_params(axis='both', labelsize=13, pad = 5)
#plt.xlim([-15,5])
plt.xscale('log')
plt.xlim([z[0], z[-376]])
plt.xticks(np.logspace(-2, 6, 5))


plt.plot(z, OmegaM, label = '$\Omega_{M} = \Omega_b + \Omega_{CDM}$')
plt.plot(z, OmegaRel, label = r'$\Omega_{r} = \Omega_{\gamma} + \Omega_{\nu}$')
plt.plot(z, valuesdic['OmegaLambda'], label = '$\Omega_{\Lambda}$')

plt.axvline(matradequalxtz[2], linestyle='--', color='#BA8E23')
plt.axvline(matlambdaequalxtz[2], linestyle='--', color='red')
plt.axvline(acceleration_startxtz[2], linestyle='--', color='black')

plt.legend(loc='upper left', fontsize=13, bbox_to_anchor=(1, 1))
plt.savefig(BACKGROUND_DIR + "/Evolution_of_cosmo_params.png", bbox_inches='tight')  


#%%

#test etaHp/c

plt.figure(8)
plt.grid()
plt.tick_params(axis='y', labelsize=11.5)
plt.xlim([-15,0])

plt.title(r"$\eta\mathcal{H}$", fontsize = 16, pad = 10)
plt.xlabel('x', fontsize = 14)

plt.plot(valuesdic['x'][:1126], [i*j/c for i,j in zip(valuesdic['eta'],valuesdic['Hp'])][:1126])

plt.savefig(BACKGROUND_DIR + "/Eta_times_Hp.png", bbox_inches='tight')  

#%%

#test etaprimeHp/c

plt.figure(9)
plt.grid()
plt.ylim([0,2])
plt.xlim([-15,0])

plt.tick_params(axis='y', labelsize=11.5)


plt.xlabel('x', fontsize = 14)
plt.title(r"$\frac{\eta\,'\mathcal{H}}{c}$", fontsize = 16, pad = 10)


plt.plot(valuesdic['x'][:1126], [(i*j)/c for i,j in zip(valuesdic['etaprime'],valuesdic['Hp'])][:1126])

plt.savefig(BACKGROUND_DIR + "/Etaprime_times_Hp.png", bbox_inches='tight')  
