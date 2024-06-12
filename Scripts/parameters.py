"""
Parameters of the model
"""
import numpy as np
from dvm_functions import length_weight_func, swimming_speed_func, micron_dict, natural_mortality_L22, mortality_allometricF, speedF
from pathlib import Path
path_to_plots = Path('~/simul_dvm2/Plots').expanduser()

### SIMULATION PARAMETERS

# Taxonomy ID
bm_idL = ["bm_fish", "bm_crust", "bm_ceph"]
bm_id = bm_idL[0]
# Mean standard length of a micronekton community in mm
L = 35
# Weight
Wg, Wc_ind, Wd = length_weight_func(L, bm_id) #Wg and Wd in mg and Wc in mgC
#Gut evacuation rate, conversion efficiency and mortality rate
d = micron_dict[bm_id][7]
e = micron_dict[bm_id][13]
f = micron_dict[bm_id][15]

taxo_ratio = micron_dict[bm_id][14]

if bm_id == "bm_fish":
    mu = natural_mortality_L22(L*0.1)
elif bm_id == "bm_crust":
    mu = mortality_allometricF(Wd*1e-3)
else:
    mu = 9e-5

#Swimming speed of micronekton
vmax = swimming_speed_func(L*0.1, bm_id) #in m/h

# AMR = 4*RMR
w_ref = vmax/2

dt_sup = 0.05
#dt_save = 0.1
dt_save = 0.8
tmax = 24 * 334 #For season: 24*334

dz = 0.2
Zmax = 500
watercolumn = np.arange(0, Zmax, dz)
zz = len(watercolumn)  # number of depth layers
CFL = 0.99  # Courant–Friedrichs–Lewy condition

temps = np.arange(0, tmax, dt_save)
month = np.arange(0, 24*31, dt_save)

mm = len(month)
tt = len(temps)
dt_inf = CFL * dz / vmax  # plus petit pas de temps possible avec la vitesse max
time12 = np.arange(0, 12, dt_inf)
Nt12 = len(time12)

# The relative gradient of light is computed
V = speedF(Nt12, dt_inf, vmax)

coef_alpha = 3

#beta = 1
#beta = 0.001
beta = 0.01
rho = 0.0062 # SEE Jager et Ravagnan 2016
Frac_PP_Z = 0.32

Csda = 0.175

# Initial conditions !! Concentration of consumer and their resource NEEDS TO BE IN CARBON WEIGHT (mg C m-3 h-1)
NPPi = 10 #34
Ki = NPPi * Frac_PP_Z
Ri = NPPi * Frac_PP_Z
Ci = (Ri*0.1) * taxo_ratio

cen = 40
sigma = 10
sigma2 = 10

#List of the parameters
alphaL = np.arange(0.05, 0.8, 0.06)
betaL = np.arange(8e-3, 1, 0.1)
sizeL = np.arange(20, 85, 3)
#sizeL_fish = np.arange(20, 85, 5)

coef_fish_alphaL = np.arange(2, 5.5, 0.2) # for fish 20-80mm
coef_crust_alphaL = np.arange(0.5, 3, 0.08) #for crust 10-50mm
coef_ceph_alphaL = np.arange(0.5, 3.5, 0.08) #for cephalopods 20-80mm 600M

#Respiration coefficients for respi_func
a0 = micron_dict[bm_id][8]
a1 = micron_dict[bm_id][9]
a2 = micron_dict[bm_id][10]
a3 = micron_dict[bm_id][11]
RQ = micron_dict[bm_id][12]

## Dictionary of the parameters with their range for the sensitivity analysis in sensitivity_analysis.py : (min, max)
problem_crust = {    "num_vars": 11,
                    "names": ["d", "e", "Mort_coef", "Csda", "fac_swim", "a_swim", 'RQ', 'a0', 'a1', 'a2', 'a3'],
                    "bounds": [
                               [0.5, 0.84],
                               [0.7, 0.9],
                               [3.95, 6.6],
                               [0.13, 0.22],
                               [1, 3],
                               [0.825, 1.375],
                               [0.61, 1.62],
                               [22.109, 24.049],
                               [0.8, 0.826],
                               [-6.528, -5.968],
                               [-0.325, -0.147]]
                    }


problem_ceph = {    "num_vars": 11,
                    "names": ["d", "e", "mu", "Csda", "fac_swim", "a_swim", 'RQ', 'a0', 'a1', 'a2', 'a3'],
                    "bounds": [
                               [0.375, 0.625],
                               [0.8, 0.96],
                               [5.7e-5, 1e-4],
                               [0.13, 0.22],
                               [1, 4],
                               [0.94, 1.56],
                               [0.64, 1.06],
                               [18.6, 30.4],
                               [0.814, 0.922],
                               [-8.074, -4.774],
                               [-0.325, -0.197]]
                    }


problem_fish = {    "num_vars": 12,
                    "names": ["d", "e", "L_inf", "K", "Csda", "fac_swim", "a_swim", 'RQ', 'a0', 'a1', 'a2', 'a3'],
                    "bounds": [
                               [0.15, 0.25],
                               [0.7, 0.9],
                               [7.68, 12.8],
                               [0.75, 1.25],
                               [0.13, 0.22],
                               [1, 4],
                               [1.35, 2.25],
                               [0.675, 1.125],
                               [28, 33.2],
                               [0.85, 0.89],
                               [-9.252, -7.778],
                               [-0.119, -0.057]]
                    }