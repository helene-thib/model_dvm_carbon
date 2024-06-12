"""
Script simulating the DVM of fish population and seasonal carbon production : respiration, excretion, dead falls and fecal pellets

Seasonal variability v1 :
- Carrying capacity (K) of the resource being primary production (mgC m-3) monthly mean from Marine Copernicus at PAP Station
- Constant temperature over time (monthly mean of June at the PAP station, represented along depth)
Establishment of the watercolumn - Computation of the swimming speed

Hélène Thibault
"""

### IMPORTS ###
from plotting_funcs import nonlinear_colormap
cm = nonlinear_colormap()

# Parameters and functions needed
from parameters import *
from dvm_functions import gaussian, I_richards, normalize, beer_lambert
from read_netCDF import int_annualNPP

from num_scheme_season import AN_RGCD_K2
import os
import pandas as pd
import numpy as np

path_to_plot = Path('~/model_dvm_carbon/Plots').expanduser()
path_to_results = Path('~/model_dvm_carbon/Results/env_K').expanduser()
# Import data of temperature in NA
path_to_data = Path('~/model_dvm_carbon/Data').expanduser()
os.chdir(path_to_data)
data_filename_temp = "thetao.csv"

os.chdir("..")

path_to_file_temp = path_to_data / data_filename_temp
df_temp = pd.read_csv(path_to_file_temp, header=7)
df_temp = df_temp.loc[0:18] #select until 1000m depth

# --- PRELIMINARY WORK FOR THE MODEL
dt_inf = CFL * dz / vmax  # Courant–Friedrichs–Lewy condition to ensure stability
time12 = np.arange(0, 12, dt_inf)
Nt12 = len(time12)

# --- Computation of the irradiance
I0 = np.zeros(Nt12)
for t in range(Nt12):
    I0[t] = I_richards(t * dt_inf + dt_inf)
I0 = np.hstack((I0, np.flip(I0)))
# Create a 2d array (time*depth) to represent the irradiance along depth
I_depth = np.zeros((zz, Nt12*2))
for z in range(zz):
    I_depth[z] = beer_lambert(I0, z)

# Normalized the irradiance between 0 and 1 to set a visual predation parameter by micronekton
I_norm = normalize(I_depth)
path_Inorm = path_to_results / "I_norm.txt"
np.savetxt(path_Inorm, I_norm)

alphaA = I_norm * coef_alpha

# --- SET ENVIRONMENTAL PARAMETERS

#Interpolation of the temperature on the watercolumn grid
Tc = np.interp(watercolumn, df_temp["depth"], df_temp["thetao"])
TkL = Tc + 273.15 # temperature in Kelvin

## INITIAL CONDITIONS
C0 = np.zeros(zz)
R0 = np.zeros(zz)
Zchl = np.zeros((zz, tmax))
for t in range(24*334):
    Kt = int_annualNPP[t]
    Zchl[:, t] = np.zeros(zz)
    Zchl[:, t][0:500] = gaussian(watercolumn[0:500], Kt, cen, sigma)

NPPi = int_annualNPP[0]
Ri = NPPi * Frac_PP_Z
Ci = Ri*0.1

# Set eq=1 if at equilibrium
eq = 0
if eq == 1:
    path_Req = path_to_results / 'Req_fishK.txt'
    path_Ceq = path_to_results / 'Ceq_fishK.txt'
    R0 = np.loadtxt(path_Req)
    C0 = np.loadtxt(path_Ceq)
else:
    C0 = gaussian(watercolumn, Ci, cen, sigma)
    R0 = gaussian(watercolumn, Ri, cen, sigma)

# EXECUTION
results = AN_RGCD_K2(R0, C0, Wc_ind, V, TkL, d, coef_alpha, I_norm, beta, dt_inf, Zchl)

# Results files at the given respective paths; the results can be read and used as numpy arrays
output_names = ['results_R.txt', 'results_G.txt', 'results_C.txt', 'results_Dg.txt', 'results_M.txt', 'results_K.txt']
for i, output in enumerate(results):
    output_path = path_to_results / output_names[i]
    np.savetxt(output_path, output)