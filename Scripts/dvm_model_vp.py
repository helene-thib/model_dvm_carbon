"""
Script simulating the DVM of micronekton population (as a macro-organism) and carbon production
              ----------Establishment of the water column----------
              ---------Computation of the swimming speed----------
              ------Computation of the visual capture rate-------
Importing functions containing numerical schemes from num_scheme_model_vp.py
Plotting the results with Resultss.py : states variables and carbon detritus

-- From Hélène Thibault - October 2023 / April 2024 --
"""

### IMPORTS
import matplotlib.pyplot as plt
import os
import pandas as pd

#Parameters and functions needed
from num_scheme_model_vp import AN_RGCD_VP, AN_RGC_VP2
from parameters import *
from plotting_funcs import contour_levels_func
from dvm_functions import I_richards, beer_lambert, windowed_mean, gaussian, normalize
import cmocean

path_to_plot = Path('~/simul_dvm2/Plots').expanduser()
path_to_results = Path('~/simul_dvm2/Results').expanduser()
# Import data of temperature in NA from Copernicus
path_to_data = Path('~/simul_dvm2/Data').expanduser()
os.chdir(path_to_data)
data_filename_temp = "thetao.csv"
os.chdir("..")
path_to_file_temp = path_to_data / data_filename_temp

df_temp = pd.read_csv(path_to_file_temp, header=7)
df_temp = df_temp.loc[0:18] #select until 500m depth

### --- PRELIMINARY WORK FOR THE MODEL --- ###

# Interpolation of the temperature on the watercolumn grid
Tc = np.interp(watercolumn, df_temp["depth"], df_temp["thetao"])
Tc = windowed_mean(Tc, 500)
TkL = np.asarray(Tc) + 273.15 # temperature in Kelvin

# Surface irradiance as a function of time of the day
time12 = np.arange(0, 12, dt_inf)
time24 = np.arange(0, 24 + dt_inf, dt_inf)
Nt12 = len(time12)

I0 = np.zeros(Nt12)
for t in range(Nt12):
    I0[t] = I_richards(t * dt_inf + dt_inf)
I0 = np.hstack((I0, np.flip(I0)))

# Create a 2d array (temps*depth) to represent the irradiance along depth
I_depth = np.zeros((zz, Nt12*2))
for z in range(zz):
    I_depth[z] = beer_lambert(I0, z)

# Normalized the irradiance between 0 and 1 to set a visual predation parameter by micronekton
I_norm = normalize(I_depth)
path_Inorm = path_to_results / "I_norm.txt"
np.savetxt(path_Inorm, I_norm)

# Implementation of the initial conditions
KL = np.zeros(zz)
KL[0:500] = gaussian(watercolumn[0:500], Ki, cen, sigma)

#If the system is at equilibrium, set eq to 1
eq = 0
if eq == 1:
    path_Req = path_to_results / 'Req_fish.txt'
    path_Ceq = path_to_results / 'Ceq_fish.txt'
    R0 = np.loadtxt(path_Req)
    C0 = np.loadtxt(path_Ceq)
else:
    C0 = gaussian(watercolumn, Ci, cen, sigma)
    R0 = gaussian(watercolumn, Ri, cen, sigma)

### --- EXECUTION of the model --- ###
det = 1
if det == 1:
    results = AN_RGCD_VP(R0, C0, V, TkL, d, e, mu, w_ref, beta, dt_inf, KL, I_norm, coef_alpha, RQ, a0, a1, a2, a3, Wc_ind, Csda)
    output_names = ['results_R.txt', 'results_G.txt', 'results_C.txt', 'results_Dg.txt', 'results_M.txt',
                    'results_SDA.txt', 'results_AMR.txt']
else:
    results = AN_RGC_VP2(R0, C0, Wg, V, TkL, d, e, mu, w_ref, Csda, beta, dt_inf, KL, I_norm, coef_alpha, RQ, a0, a1, a2, a3)
    output_names = ['results_R.txt', 'results_G.txt', 'results_C.txt']

### --- Save the results : States variables (R, G, C) and carbon detritus --- ###
for i, output in enumerate(results):
    output_path = path_to_results / output_names[i]
    np.savetxt(output_path, output)
