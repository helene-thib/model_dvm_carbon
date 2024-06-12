"""
Script to perform a sensitivity analysis of bioenergetic parameters for a community of crustaceans (3cm)
- Generate parameters values with a quasi-random MonteCarlo method
- Generate the outputs : biomass, detritus and pe-ratio200m
> Then use the script sensitivity_analysis.py to plot the results

Outputs: biomass of C and R, detritus and pe-ratio200
"""

import pandas as pd
import numpy as np

from SALib.sample import sobol
from num_scheme_model_vp import AN_RGC_VP2, AN_RGCD_VP
from parameters import *
from dvm_functions import gaussian, windowed_mean, I_richards, beer_lambert, normalize
import os
path_to_results = Path('~/simul_dvm2/Results/Sobol').expanduser()

path_to_plot = Path('~/simul_dvm2/Plots').expanduser()

# Import data of temperature in NA
path_to_data = Path('~/simul_dvm2/Data').expanduser()
data_filename_temp = "thetao.csv"
os.chdir("..")
path_to_file_temp = path_to_data / data_filename_temp
df_temp = pd.read_csv(path_to_file_temp, header=7)
df_temp = df_temp.loc[0:18] #select until 500m depth

#Interpolation of the temperature on the watercolumn grid
Tc = np.interp(watercolumn, df_temp["depth"], df_temp["thetao"])
Tc = windowed_mean(Tc, 500)
TkL = np.asarray(Tc) + 273.15 # temperature in Kelvin

# Generate Samples
param_values = sobol.sample(problem_crust, 2**6) #N must be power of 2
Ctot = np.zeros([param_values.shape[0]])
Rtot = np.zeros([param_values.shape[0]])
Cmax = np.zeros([param_values.shape[0]])
Rmax = np.zeros([param_values.shape[0]])

# Implementation of the initial conditions
KL = np.zeros(zz)
KL[0:500] = gaussian(watercolumn[0:500], Ki, cen, sigma)
C0 = gaussian(watercolumn, Ci, cen, sigma)
R0 = gaussian(watercolumn, Ri, cen, sigma)

#Run the model for crust
coef_alpha = 3
detritus = 0
if detritus == 0:
    for i, X in enumerate(param_values):
        bm_id = "bm_crust"
        d, e, Mort_coef, Csda, fac_swim, a_swim, RQ, a0, a1, a2, a3 = X

        mu = mortality_allometricF(Wd * 1e-3, Mort_coef)

        vmax = (a_swim * L*0.1)/100 * 3600

        dt_inf = CFL * dz / vmax  # plus petit pas de temps possible avec la vitesse max
        time12 = np.arange(0, 12, dt_inf)
        time24 = np.arange(0, 24 + dt_inf, dt_inf)
        Nt12 = len(time12)

        I0 = np.zeros(Nt12)
        for t in range(Nt12):
            I0[t] = I_richards(t * dt_inf + dt_inf)
        I0 = np.hstack((I0, np.flip(I0)))

        # Create a 2d array (temps*depth) to represent the irradiance along depth
        I_depth = np.zeros((zz, Nt12 * 2))
        for z in range(zz):
            I_depth[z] = beer_lambert(I0, z)

        # Normalized the irradiance between 0 and 1 to set a visual predation parameter by micronekton
        I_norm = normalize(I_depth)

        Vm = speedF(Nt12, dt_inf, vmax)
        w_ref = vmax / fac_swim

        Ctot[i], Rtot[i], Cmax[i], Rmax[i] = AN_RGC_VP2(R0, C0, Wc_ind, Vm, TkL, d, e, mu, w_ref, Csda, beta, dt_inf, KL, I_norm, coef_alpha, RQ, a0, a1, a2, a3)

    ### SAVE THE OUTPUTS
    output_names = ["sobol_crust_Ctot.txt", "sobol_crust_Rtot.txt", "sobol_crust_Cmax.txt", "sobol_crust_Rmax.txt", "sobol_param_crust.txt"]
    outputs = [Ctot, Rtot, Cmax, Rmax, param_values]

    for i, output in enumerate(outputs):
        output_path = path_to_results / output_names[i]
        np.savetxt(output_path, output)

else:
    ##Load data
    path_param = path_to_results / "sobol_param_crust.txt"
    path_Rmax = path_to_results / "sobol_crust_Rmax.txt"
    path_Cmax = path_to_results / "sobol_crust_Cmax.txt"
    param_values = np.loadtxt(path_param)
    RmaxL = np.loadtxt(path_Rmax)
    CmaxL = np.loadtxt(path_Cmax)

    ##Create 2d zero array to save the results
    Resp = np.zeros(param_values.shape[0])
    POC = np.zeros(param_values.shape[0])
    pe_ratio_200 = np.zeros(param_values.shape[0])
    NPP = np.sum(KL*dz)/Frac_PP_Z

    for i, X in enumerate(param_values):
        bm_id = "bm_crust"
        d, e, Mort_coef, Csda, fac_swim, a_swim, RQ, a0, a1, a2, a3 = X
        mu = mortality_allometricF(Wd * 1e-3, Mort_coef)
        vmax = (a_swim * L * 0.1) / 100 * 3600

        dt_inf = CFL * dz / vmax  # plus petit pas de temps possible avec la vitesse max
        time12 = np.arange(0, 12, dt_inf)
        time24 = np.arange(0, 24 + dt_inf, dt_inf)
        Nt12 = len(time12)

        I0 = np.zeros(Nt12)
        for t in range(Nt12):
            I0[t] = I_richards(t * dt_inf + dt_inf)
        I0 = np.hstack((I0, np.flip(I0)))

        # Create a 2d array (temps*depth) to represent the irradiance along depth
        I_depth = np.zeros((zz, Nt12 * 2))
        for z in range(zz):
            I_depth[z] = beer_lambert(I0, z)

        # Normalized the irradiance between 0 and 1 to set a visual predation parameter by micronekton
        I_norm = normalize(I_depth)

        Vm = speedF(Nt12, dt_inf, vmax)
        w_ref = vmax / fac_swim

        R0 = np.zeros(zz)
        C0 = np.zeros(zz)
        R0[0:500] = gaussian(watercolumn[0:500], RmaxL[i], cen, sigma)
        C0[0:500] = gaussian(watercolumn[0:500], CmaxL[i], cen, sigma)

        results = AN_RGCD_VP(R0, C0, Vm, TkL, d, e, mu, w_ref, beta, dt_inf, KL, I_norm, coef_alpha, RQ, a0, a1, a2, a3, Wc_ind, Csda)
        POC[i] = np.sum(results[3] + results[2] * mu * dt_save)
        Resp[i] = np.sum(results[4] + results[5] + results[6])

        # Calculate the pe-ratio : efficiency of carbon transport below 200m
        id_200 = int(200 / dz)
        Dg_int_200 = np.sum(results[3][:, id_200:zz])
        Mort_int_200 = np.sum(results[2][:, id_200:zz] * mu * dt_save)
        pe_ratio_200[i] = (Dg_int_200 + Mort_int_200) / (NPP * 24)

    ### SAVE THE OUTPUTS
    output_names = ["sobol_crust_Resp.txt", "sobol_crust_POC.txt", "sobol_crust_pe_ratio200.txt"]
    outputs = [Resp, POC, pe_ratio_200]
    for i, output in enumerate(outputs):
        output_path = path_to_results / output_names[i]
        np.savetxt(output_path, output)