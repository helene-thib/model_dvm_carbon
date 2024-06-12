"""
Script simulating the DVM of fish population and seasonal carbon production : respiration, excretion, dead falls and fecal pellets

Seasonal variability: All environmental parameters varying together
- Carrying capacity (K) of the resource being phytoplankton concentrations (mgC m-3) monthly mean from Marine Copernicus at PAP Station
- Varying temperature over one year (monthly mean of June at the PAP station, represented along depth)
- Computation of the variation of irradiance along depth and over a year
Plot the results of the numerical scheme with Results_envKTI.py

By Helene Thibault
"""

### IMPORTS
from plotting_funcs import nonlinear_colormap
import cmocean
import matplotlib.pyplot as plt
import pvlib
from pvlib.location import Location
import pandas as pd
import os
cm = nonlinear_colormap()

#Parameters and functions needed
from parameters import *
from dvm_functions import gaussian, I_sinus_annualF, I_sinusF, normalize, vant_hoff
from num_scheme_season import AN_RGCD_KTI2
from read_netCDF import annualT, int_annualNPP, sum_annualNPP

path_to_plot = Path('~/simul_dvm2/Plots').expanduser()
path_to_results = Path('~/simul_dvm2/Results').expanduser()
path_to_results_K = Path('~/simul_dvm2/Results/env_K').expanduser()
path_to_results_KTI = Path('~/simul_dvm2/Results/env_KTI').expanduser()
# Import data of temperature in NA
path_to_data = Path('~/simul_dvm2/Data').expanduser()
os.chdir(path_to_data)
os.chdir("..")

path_to_file_chl = path_to_data / "chl.csv"
df_chl = pd.read_csv(path_to_file_chl, header=7)

## Calculation of the Kd PAR at 490nm
df_chl['time'] = df_chl['time'].str[5:-17]
df_chl = df_chl.groupby(['time'])['chl'].max().reset_index(name='monthly_chl')
df_chl["Kd490"] = 0.0166 + 0.072*df_chl["monthly_chl"]**0.69
df_chl['time'] = pd.to_numeric(df_chl['time'])

#Interpolation for each day of a year
year = np.arange(0, 365, 1)
year2 = np.arange(0, 24*334, 1)
Kd_annual = np.interp(year, df_chl['time']*30.5, df_chl["Kd490"])
sum_NPP = np.interp(temps, year2, sum_annualNPP)
# Convert the temperature from Celsius to Kelvin
annual_Tk = annualT + 273.15
# Growth rate of zooplankton temperature-dependant
rhoT = np.vectorize(vant_hoff)(annualT)

### PRELIMINARY WORK FOR THE MODEL
dt_inf = CFL * dz / vmax  # plus petit pas de temps possible avec la vitesse max
time12 = np.arange(0, 12, dt_inf)
Nt12 = len(time12)
Nt24 = Nt12 * 2
time24 = np.arange(0, 24, dt_save)

# Definition of Location object. Coordinates and elevation of Madrid Ciemat Headquarters (Spain)
site = Location(48.5, 16.5, 'Etc/GMT+0', 0) # latitude, longitude, time_zone, altitude, name
# Definition of a time range of simulation
times = pd.date_range('2022-01-01 00:00:00', '2022-12-31 23:59:00', freq='D', tz=site.tz)

# Estimate Solar Position with the 'Location' object
solpos = pvlib.solarposition.get_solarposition(times, site.latitude, site.longitude, site.altitude)
annual_ZA = np.array(np.cos(np.radians(solpos['zenith'])))

### --- Annual variation of surface irradiance --- ###
path_speed = path_to_results / 'V_annual.txt'
V_annual365 = np.loadtxt(path_speed)
V_annual = V_annual365[15:349]

### --- Compute visual coefficient alpha at the surface --- ###
shape = (len(year), Nt24)
I0 = np.zeros(shape)
for j in year:
    n = int(I_sinus_annualF(j))
    I_day = np.zeros(Nt24)
    for t in range(Nt24):
        ZA = annual_ZA[j-15]
        I_day[t] = I_sinusF(t * dt_inf + dt_inf, n=n, a=200)*ZA
        #I_depth = np.zeros((zz, Nt12 * 2))
    I0[j-15] = I_day

# Only at the surface !!
I_norm = normalize(I0[15:349])
im = plt.imshow(I_norm, aspect='auto', extent=[0, 24, 0, 364], cmap=cmocean.cm.solar)
plt.close()

## INITIAL CONDITIONS : NPP, R and C
Zchl = np.zeros((zz, tmax))
for t in range(24*334):
    Kt = int_annualNPP[t] * Frac_PP_Z
    Zchl[:, t] = np.zeros(zz)
    Zchl[:, t][0:500] = gaussian(watercolumn[0:500], Kt, cen, sigma2)

Ri = NPPi * Frac_PP_Z
Ci = Ri*0.1

#If the system is at equilibrium, set eq to 1
eq = 0
if eq == 1:
    path_Req = path_to_results_KTI / 'Req_fish.txt'
    path_Ceq = path_to_results_KTI / 'Ceq_fish.txt'
    R0 = np.loadtxt(path_Req)
    C0 = np.loadtxt(path_Ceq)
else:
    C0 = gaussian(watercolumn, Ci, cen, sigma)
    R0 = gaussian(watercolumn, Ri, cen, sigma2)

### --- EXECUTION of the model --- ###
results = AN_RGCD_KTI2(R0, C0, V_annual, annual_Tk, Kd_annual, I_norm, d, coef_alpha, beta, dt_inf, rhoT, Zchl)

# Calculation of pe_ratio below 200m
id_200 = int(200 / dz)
Dg = results[3]
Dg_int_200 = np.sum(results[3][:, id_200:zz], axis=1)
Mort_int_200 = np.sum(results[2][:, id_200:zz] * mu * dt_save, axis=1)

pe_ratioL = (Dg_int_200 + Mort_int_200) / sum_NPP
pe_ratioL = np.bincount(np.arange(len(pe_ratioL))//30, pe_ratioL)
output_path_pr = path_to_results_KTI / "pe_ratio.txt"
np.savetxt(output_path_pr, pe_ratioL)

# Results files at the given respective paths; the results can be read and used as numpy arrays
output_names = ['results_R.txt', 'results_G.txt', 'results_C.txt', 'results_Dg.txt', 'results_M.txt', "alphaL.txt"]
### --- Save the results : States variables (R, G, C) and carbon detritus --- ###
for i, output in enumerate(results):
    output_path = path_to_results_KTI / output_names[i]
    np.savetxt(output_path, output)



