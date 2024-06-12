"""
Script simulating the DVM of fish population and seasonal carbon production : respiration, excretion, dead falls and fecal pellets

Seasonal variability v3 :
- Carrying capacity (K) of the resource being phytoplankton concentrations (mgC m-3) monthly mean from Marine Copernicus at PAP Station
- Varying temperature over one year (monthly mean of June at the PAP station, represented along depth)
- New : Computation of the variation of irradiance along depth and over a year

Commentaries:
-Set Zmax=700 in parameters.py

Hélène Thibault
"""

### IMPORTS
from plotting_funcs import nonlinear_colormap
from pylab import colorbar
import cmocean
import matplotlib.pyplot as plt
import pvlib
from pvlib.location import Location
import pandas as pd

cm = nonlinear_colormap()

#Parameters and functions needed

from parameters import *
from dvm_functions import gaussian, I_sinus_annualF, I_sinusF, I_richards, normalize, SDA_func

from num_scheme_season import AN_RGCD_I2

import os

path_to_plot = Path('~/simul_dvm2/Plots').expanduser()
path_to_results = Path('~/simul_dvm2/Results').expanduser()
path_to_results_I = Path('~/simul_dvm2/Results/env_I').expanduser()
data_filename_temp = "thetao.csv"
# Import data of temperature in NA
path_to_data = Path('~/simul_dvm2/Data').expanduser()
os.chdir(path_to_data)
os.chdir("..")

path_to_file_chl = path_to_data / "chl.csv"
df_chl = pd.read_csv(path_to_file_chl, header=7)

## Calculation of the Kd PAR at 490nm
df_chl['time'] = df_chl['time'].str[5:-17]
df_chl = df_chl.groupby(['time'])['chl'].mean().reset_index(name='monthly_chl')
df_chl["Kd490"] = 0.0166 + 0.072*df_chl["monthly_chl"]**0.69
df_chl['time'] = pd.to_numeric(df_chl['time'])

#Interpolation for each day of a year
year = np.arange(0, 365, 1)
Kd_annual = np.interp(year, df_chl['time']*30.5, df_chl["Kd490"])

path_to_file_temp = path_to_data / data_filename_temp
df_temp = pd.read_csv(path_to_file_temp, header=7)
df_temp = df_temp.loc[0:18] #select until 1000m depth
#Interpolation of the temperature on the watercolumn grid
Tc = np.interp(watercolumn, df_temp["depth"], df_temp["thetao"])
TkL = Tc + 273.15 # temperature in Kelvin

### PRELIMINARY WORK FOR THE MODEL
dt_inf = CFL * dz / vmax  # plus petit pas de temps possible avec la vitesse max
time12 = np.arange(0, 12, dt_inf)
Nt12 = len(time12)
Nt24 = Nt12 * 2
time24 = np.arange(0, 24, 0.2)

# Definition of Location oject. Coordinates and elevation of Madrid Ciemat Headquarters (Spain)
site = Location(48.5, 16.5, 'Etc/GMT+0', 0) # latitude, longitude, time_zone, altitude, name
# Definition of a time range of simulation
times = pd.date_range('2022-01-01 00:00:00', '2022-12-31 23:59:00', freq='D', tz=site.tz)

# Estimate Solar Position with the 'Location' object
solpos = pvlib.solarposition.get_solarposition(times, site.latitude, site.longitude, site.altitude)
annual_ZA = np.array(np.cos(np.radians(solpos['elevation'])))

### --- Annual variation of surface irradiance
compute_irr = 1
count = 0
if compute_irr == 1:
    I_annual = np.zeros((len(year), len(time24)))
    for j in year:
        n = I_sinus_annualF(j)
        ZA = annual_ZA[j]
        I_annual[j] = (I_sinusF(time24, n)*ZA)
    fig = plt.figure(figsize=(4, 4), layout="constrained")
    im = plt.imshow(I_annual.T, aspect='auto', extent=[1, 12, 0, 24], cmap=cmocean.cm.solar)
    plt.xticks([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
               ["Jan", "Fev", "Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec"], size=10)
    plt.tick_params(axis='x', labelrotation=45)
    plt.ylabel('Hours of a day')
    plt.xlabel('Months')
    colorbar(im, label='Relative gradient of surface irradiance', location='top')
    name_irran = 'Annual_irr.pdf'
    plot_irran = path_to_plot / name_irran
    plt.savefig(plot_irran)
    plt.close()

    # The relative gradient of light is computed for the first 12h of the day, symmetric is taken and concatenated to it
    dI_annual = np.zeros((len(year), Nt24))
    dI_r_maxL = np.zeros(len(year))

    for j in year:
        #Solar zenith angle
        ZA = annual_ZA[j]
        # Parameter that change the length of day for computing the daily surface irradiance
        n = int(I_sinus_annualF(j))

        dI_r_max = 0
        dI_r = np.zeros(Nt24)
        It_1 = np.zeros(Nt24)
        It = np.zeros(Nt24)
        for t in range(Nt24):
            It_1[t] = I_sinusF(t * dt_inf + dt_inf, n)*ZA
            It[t] = I_sinusF(t * dt_inf, n)*ZA
            # Computing the rate of change of irradiance
            dI_r[t] = (1 / dt_inf) * (It_1[t] / It[t] - 1)
            if np.abs(dI_r[t]) > dI_r_max:
                dI_r_max = np.abs(dI_r[t])
        dI_r_maxL[j] = dI_r_max
        dI_annual[j] = dI_r

    # Computation of the speed for each day of one year
    V_annual = np.zeros((len(year), Nt24))
    plot_time24 = []
    for j in year:
        dI_r_max = dI_r_maxL[j]
        V_day = np.zeros(Nt24)
        for t in range(Nt24):
            V_day[t] = vmax * dI_annual[j,t] / dI_r_max
        V_annual[j] = V_day
    np.savetxt(path_to_results / "V_annual.txt", V_annual)

    fig = plt.figure()
    im2 = plt.imshow(V_annual, aspect='auto', extent=[0, 24, 0, 364], cmap=cmocean.cm.balance)
    plt.xlabel('Hours of a day')
    plt.ylabel('Days of the year')
    cbar = fig.colorbar(im2, label='Relative daily speed [$m$ $h^{-1}$]')
    name_speed = 'Annual_speed.pdf'
    plot_speed = path_to_plot / name_speed
    plt.savefig(plot_speed)
    plt.close()
else:
    path_speed = path_to_results / 'V_annual.txt'
    V_annual = np.loadtxt(path_speed)
V_annual = V_annual[15:349]

### --- Compute visual coefficient alpha
shape = (len(year), Nt24)
I0 = np.zeros(shape)
for j in year:
    n = int(I_sinus_annualF(j))
    I_day = np.zeros(Nt24)
    Kd = Kd_annual[j]
    for t in range(Nt24):
        ZA = annual_ZA[j]
        I_day[t] = I_sinusF(t * dt_inf + dt_inf, n=n, a=200)*ZA
        #I_depth = np.zeros((zz, Nt12 * 2))
    I0[j] = I_day

## Normalized the irradiance between 0 and 1 to set a visual predation parameter by micronekton
# Only at the surface !!
I_norm = normalize(I0[15:349])
im = plt.imshow(I_norm, aspect='auto', extent=[0, 24, 0, 364], cmap=cmocean.cm.solar)

# Implementation of the initial conditions
KL = np.zeros(zz)
Ki = 30
KL[0:500] = gaussian(watercolumn[0:500], Ki, cen, sigma)
Ri = Ki * Frac_PP_Z
Ci = Ri*0.1

#eq = int(input('Enter 1 at equilibrium else 0 :'))
eq = 0
if eq == 1:
    path_Req = path_to_results_I / 'Req_fishI.txt'
    path_Ceq = path_to_results_I / 'Ceq_fishI.txt'
    Req = np.loadtxt(path_Req)
    Ceq = np.loadtxt(path_Ceq)
    C0 = gaussian(watercolumn, np.max(Ceq), cen, sigma)
    R0 = gaussian(watercolumn, np.max(Req), cen, sigma)
else:
    C0 = gaussian(watercolumn, Ci, cen, sigma)
    R0 = gaussian(watercolumn, Ri, cen, sigma)

Kd_annual = Kd_annual[15:349]
# EXECUTION
results = AN_RGCD_I2(R0, C0, V_annual, TkL, Kd_annual, I_norm, d, coef_alpha, beta, dt_inf, KL)
# Results files at the given respective paths; the results can be read and used as numpy arrays
output_names = ['results_R.txt', 'results_G.txt', 'results_C.txt', 'results_Dg.txt', 'results_M.txt', "alphaL.txt"]
for i, output in enumerate(results):
    output_path = path_to_results_I / output_names[i]
    np.savetxt(output_path, output)