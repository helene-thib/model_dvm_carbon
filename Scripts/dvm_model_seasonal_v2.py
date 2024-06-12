"""
Script simulating the DVM of fish population and seasonal carbon production : respiration, excretion, dead falls and fecal pellets

Seasonal variability v2:
- MODIFICATION from v1 : Variation of temperature over time (monthly mean of June at the PAP station, represented along depth)
+ Growth rate of the resource varying with temperature and remineralization of FP

Hélène Thibault
"""

### IMPORTS
import os
import matplotlib.pyplot as plt
# Parameters and functions needed
from parameters import *
from dvm_functions import I_richards, gaussian, vant_hoff, normalize, beer_lambert
from num_scheme_season import AN_RGCD_T2
from read_netCDF import annualT

## PATHS
path_to_plot = Path('~/simul_dvm2/Plots').expanduser()
path_to_results = Path('~/simul_dvm2/Results').expanduser()
path_to_results_K = Path('~/simul_dvm2/Results/env_K').expanduser()
path_to_results_T = Path('~/simul_dvm2/Results/env_T').expanduser()
# Import data of temperature in NA
path_to_data = Path('~/simul_dvm2/Data').expanduser()
os.chdir(path_to_data)
os.chdir("..")

## LOAD DATA
# Convert the temperature from Celsius to Kelvin
annual_Tk = annualT + 273.15

# --- PRELIMINARY WORK FOR THE MODEL ---#
dt_inf = CFL * dz / vmax  # plus petit pas de temps possible avec la vitesse max
time12 = np.arange(0, 12, dt_inf)
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
alphaA = I_norm * coef_alpha

### Implementation of the initial conditions ###
KL = np.zeros(zz)
Ki = 150
KL[0:500] = gaussian(watercolumn[0:500], Ki, cen, sigma)

Ri = Ki * Frac_PP_Z
Ci = Ri*0.1

# Set eq=1 if at equilibrium
eq = 0
if eq == 1:
    path_Req = path_to_results_T / 'Req_fishT.txt'
    path_Ceq = path_to_results_T / 'Ceq_fishT.txt'
    R0 = np.loadtxt(path_Req)
    C0 = np.loadtxt(path_Ceq)
else:
    C0 = np.zeros(zz)
    R0 = np.zeros(zz)
    # Concentration de depart des 100m surface
    C0[0:500] = gaussian(watercolumn[0:500], Ci, cen, sigma)
    R0[0:500] = gaussian(watercolumn[0:500], Ri, cen, sigma)

fig = plt.figure()
plt.plot(KL, -watercolumn, label="K")
plt.plot(C0, -watercolumn, label="C")
plt.plot(R0, -watercolumn, label="R")
plt.title("Initial concentrations")
plt.legend()
plot_RCKi = path_to_plot / 'RCKi.pdf'
plt.savefig(plot_RCKi)
plt.close()
#Growth rate of zooplankton dependant of temperature
rhoT = np.vectorize(vant_hoff)(annualT)

# EXECUTION
results = AN_RGCD_T2(R0, C0, V, annual_Tk, d, coef_alpha, beta, I_norm, dt_inf, KL, rhoT)
## Save the results
output_names = ['results_R.txt', 'results_G.txt', 'results_C.txt', 'results_Dg.txt', 'results_M.txt']
for i, output in enumerate(results):
    output_path = path_to_results_T / output_names[i]
    np.savetxt(output_path, output)


