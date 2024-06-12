"""
Script to perform a sensitivity analysis of bioenergetic parameters for a community of fish (3cm)
- Generate parameters values with a quasi-random MonteCarlo method
- Generate the outputs : biomass, detritus and pe-ratio200m
> Then use the script sensitivity_analysis.py to plot the results

Outputs: biomass of C and R, detritus and pe-ratio200
"""

import pandas as pd
import matplotlib.pyplot as plt

from SALib.sample import sobol
from num_scheme_model_vp import AN_RGC_VP2, AN_RGCD_VP
from parameters import *
from dvm_functions import gaussian, windowed_mean, natural_mortality_L22, I_richards, beer_lambert, normalize
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

#Swimming speed of different taxonomic groups
swim = 0
if swim == 1:
    swim_dict = {    "bm_ceph":    [0.5, 2],
                     "bm_crust":    [0.2, 2],
                     "bm_fish_inf":    [2, 5]}

    taxoL = swim_dict.keys()
    fish_lengthL = []
    crust_lengthL = []
    ceph_lengthL = []
    fish_vL = []
    crust_vL = []
    ceph_vL = []

    for bm_id in taxoL:
        for i in range(1000000):
            length = np.random.uniform(1, 5)
            a = np.random.uniform(swim_dict[bm_id][0], swim_dict[bm_id][1])
            V = a * length
            if bm_id == "bm_fish_inf":
                fish_lengthL.append(length)
                fish_vL.append(V)
            elif bm_id == "bm_crust":
                crust_lengthL.append(length)
                crust_vL.append(V)
            else:
                ceph_lengthL.append(length)
                ceph_vL.append(V)


    zipped = list(zip(fish_lengthL, fish_vL))
    df_fish = pd.DataFrame(zipped, columns=['BL', 'SS'])
    zipped = list(zip(crust_lengthL, crust_vL))
    df_crust = pd.DataFrame(zipped, columns=['BL', 'SS'])
    zipped = list(zip(ceph_lengthL, ceph_vL))
    df_ceph = pd.DataFrame(zipped, columns=['BL', 'SS'])

    df_fish['BL'] = np.round(df_fish['BL'], 1)
    df_ceph['BL'] = np.round(df_ceph['BL'], 1)
    df_crust['BL'] = np.round(df_crust['BL'], 1)
    df_fish_mean = df_fish.groupby(['BL'])['SS'].mean().reset_index(name='mean_SS')
    df_fish_std = df_fish.groupby(['BL'])['SS'].std().reset_index(name='std_SS')
    df_ceph_mean = df_ceph.groupby(['BL'])['SS'].mean().reset_index(name='mean_SS')
    df_ceph_std = df_ceph.groupby(['BL'])['SS'].std().reset_index(name='std_SS')
    df_crust_mean = df_crust.groupby(['BL'])['SS'].mean().reset_index(name='mean_SS')
    df_crust_std = df_crust.groupby(['BL'])['SS'].std().reset_index(name='std_SS')

    #plt.scatter(fish_lengthL, fish_vL, label='Fish', alpha=.01)
    plt.plot(df_fish_mean['BL'], df_fish_mean['mean_SS'], label='Fish')
    plt.fill_between(df_fish_mean['BL'], df_fish_mean['mean_SS']-df_fish_std['std_SS'],
                     df_fish_mean['mean_SS']+df_fish_std['std_SS'],alpha=.1)
    #plt.errorbar(df_fish_mean['BL'], df_fish_mean['mean_SS'], df_fish_std['std_SS'], label='Fish')

    #plt.scatter(crust_lengthL, crust_vL, label='Crustacea', alpha=.01)
    plt.plot(df_crust_mean['BL'], df_crust_mean['mean_SS'], label='Crustacea')
    plt.fill_between(df_crust_mean['BL'], df_crust_mean['mean_SS']-df_crust_std['std_SS'],
                     df_crust_mean['mean_SS']+df_crust_std['std_SS'],alpha=.1)

    #plt.scatter(ceph_lengthL, ceph_vL, label='Cephalopods', alpha=.01)
    plt.plot(df_ceph_mean['BL'], df_ceph_mean['mean_SS'], label='Cephalopods')
    plt.fill_between(df_ceph_mean['BL'], df_ceph_mean['mean_SS']-df_ceph_std['std_SS'],
                     df_ceph_mean['mean_SS']+df_ceph_std['std_SS'],alpha=.1)
    plt.legend(loc='upper left')
    plt.xlabel("Body length [cm]")
    plt.ylabel("Migration speed [cm/s]")
    name_speed = "BL_speed_mc.pdf"
    plot_name_speed = path_to_plot / name_speed
    plt.savefig(plot_name_speed)


## Dictionary of the parameters with their range : (min, value, max)

sens_dict_fish = { "rho":       [0.00496, 0.0062, 0.00744], # ww in g to cw (conversion from g to mg)
                   "d":         [0.16, 0.2, 0.24], # dw in mg to cw for Euphausiids (R2=0.98)
                   "e":         [0.14, 0.17, 0.2],
                   "L_inf":     [3.3, 10.24, 26.5],
                   "K":         [0.17, 0.53, 5.62],
                   "phi":       [1.04, 1.75, 2.52],
                   "fac_feed":  [1.6, 2, 2.4],
                   "fac_swim":  [3.2, 4, 4.8],
                   "RQ":        [0.72, 0.9, 1.08],
                   "a0":        [28.316, 30.767, 33.218],
                   "a1":        [0.85, 0.87, 0.89],
                   "a2":        [-7.778, -8.515, -9.252],
                   "a3":        [0.057, 0.088, 0.119]}

problem_fish2 = {    "num_vars": 12,
                    "names": ["d", "e", "L_inf", "K", "Csda", "fac_swim", "a_swim", 'RQ', 'a0', 'a1', 'a2', 'a3'],
                    "bounds": [
                               [0.15, 0.25],
                               [0.13, 0.21],
                               [7.68, 12.8],
                               [0.75, 1.25],
                               [0.12, 0.16],
                               [1, 4],
                               [1.35, 2.25],
                               [0.675, 1.125],
                               [28, 33.2],
                               [0.85, 0.89],
                               [7.778, 9.252],
                               [0.057, 0.119]
]
                    }
# Generate Samples
param_values = sobol.sample(problem_fish, 2**6) #N must be power of 2

Ctot = np.zeros([param_values.shape[0]])
Rtot = np.zeros([param_values.shape[0]])
Cmax = np.zeros([param_values.shape[0]])
Rmax = np.zeros([param_values.shape[0]])

# Implementation of the initial conditions
KL = np.zeros(zz)
KL[0:500] = gaussian(watercolumn[0:500], Ki, cen, sigma)
C0 = gaussian(watercolumn, Ci, cen, sigma)
R0 = gaussian(watercolumn, Ri, cen, sigma)

#Run the model for fish
coef_alpha = 4
detritus = 1
if detritus == 0:
    for i, X in enumerate(param_values):
        bm_id = "bm_fish"
        #d, e, L_inf, K, Csda, fac_swim, a_swim, RQ, a0, a1, a2, a3 = X
        d, e, L_inf, K, Csda, fac_swim, a_swim = X
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
        mu = natural_mortality_L22(L*0.1, L_inf, K)

        w_ref = vmax / fac_swim

        Ctot[i], Rtot[i], Cmax[i], Rmax[i] = AN_RGC_VP2(R0, C0, Wc_ind, Vm, TkL, d, e, mu, w_ref, Csda, beta, dt_inf, KL, I_norm, coef_alpha, RQ, a0, a1, a2, a3)

    ### SAVE THE OUTPUTS
    output_names = ["sobol_fish_Ctot.txt", "sobol_fish_Rtot.txt", "sobol_fish_Cmax.txt", "sobol_fish_Rmax.txt",
                    "sobol_param_fish.txt"]
    outputs = [Ctot, Rtot, Cmax, Rmax, param_values]

    for i, output in enumerate(outputs):
        output_path = path_to_results / output_names[i]
        np.savetxt(output_path, output)

else:
    ##Load data
    path_param = path_to_results / "sobol_param_fish.txt"
    path_Rmax = path_to_results / "sobol_fish_Rmax.txt"
    path_Cmax = path_to_results / "sobol_fish_Cmax.txt"
    param_values = np.loadtxt(path_param)
    RmaxL = np.loadtxt(path_Rmax)
    CmaxL = np.loadtxt(path_Cmax)

    ##Create 2d zero array to save the results
    Resp = np.zeros(param_values.shape[0])
    POC = np.zeros(param_values.shape[0])
    pe_ratio_200 = np.zeros(param_values.shape[0])
    NPP = np.sum(KL*dz)/Frac_PP_Z

    for i, X in enumerate(param_values):
        bm_id = "bm_fish"
        #d, e, L_inf, K, Csda, fac_swim, a_swim, RQ, a0, a1, a2, a3 = X
        d, e, L_inf, K, Csda, fac_swim, a_swim = X
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
        mu = natural_mortality_L22(L * 0.1, L_inf, K)

        w_ref = vmax / fac_swim

        R0 = np.zeros(zz)
        C0 = np.zeros(zz)
        R0[0:500] = gaussian(watercolumn[0:500], RmaxL[i], cen, sigma)
        C0[0:500] = gaussian(watercolumn[0:500], CmaxL[i], cen, sigma)

        results = AN_RGCD_VP(R0, C0, Vm, TkL, d, e, mu, w_ref, beta, dt_inf, KL, I_norm, coef_alpha, RQ, a0, a1, a2, a3, Wc_ind, Csda)
        POC[i] = np.sum(results[3] + results[2] * mu*dt_save)
        Resp[i] = np.sum(results[4] + results[5] + results[6])

        # Calculate the pe-ratio : efficiency of carbon transport below 200m
        id_200 = int(200 / dz)
        Dg_int_200 = np.sum(results[3][:, id_200:zz])
        Mort_int_200 = np.sum(results[2][:, id_200:zz] * mu*dt_save)
        pe_ratio_200[i] = (Dg_int_200 + Mort_int_200) / (NPP * 24)

    ### SAVE THE OUTPUTS
    output_names = ["sobol_fish_Resp.txt", "sobol_fish_POC.txt", "sobol_fish_pe_ratio200.txt"]
    outputs = [Resp, POC, pe_ratio_200]
    for i, output in enumerate(outputs):
        output_path = path_to_results / output_names[i]
        np.savetxt(output_path, output)