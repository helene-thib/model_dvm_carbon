"""
Global descriptors of detritus to test the influence of the size and taxonomy on the model
Stable values generated with global_descriptors.py

-Set Zmax=700, tmax=24 in parameters.py
"""
import os

import pandas as pd
from matplotlib import pyplot as plt
from sklearn import linear_model
import matplotlib as mpl
import seaborn as sns
from matplotlib import cm

from parameters import *
from dvm_functions import length_weight_func, I_richards, gaussian, beer_lambert, normalize
from num_scheme_model_vp import AN_RGCD_VP

path_to_result = Path('~/simul_dvm2/Results/GB').expanduser()
path_to_temp = Path('~/simul_dvm2/Results').expanduser()
path_to_plots = Path('~/simul_dvm2/Plots').expanduser()
path_to_data = Path('~/simul_dvm2/Data').expanduser()

# Import data of temperature in NA
data_filename_temp = "temp_profil.txt"
os.chdir("..")
path_to_file_temp = path_to_temp / data_filename_temp

### LOAD FILES
# Load and concatenate data
ab_fish = pd.read_csv(path_to_result /"alpha_size_fish.txt")
ab_crust = pd.read_csv(path_to_result /"alpha_size_crust.txt")
ab_ceph = pd.read_csv(path_to_result /"alpha_size_ceph.txt")

### PRELIMINARY WORK FOR THE MODEL
# Temperature profile
TkL = np.loadtxt(path_to_file_temp)
kn = np.arange(20, 80)
kn_crust = np.arange(10, 50)

## EXECUTION : testing the influence of size and taxonomic groups on the production of detritus

# Select R and C initial concentrations based en ratio C/R
ab_fish["C/R"] = ab_fish["C"] / ab_fish["R"]
ab_crust["C/R"] = ab_crust["C"] / ab_crust["R"]
ab_ceph["C/R"] = ab_ceph["C"] / ab_ceph["R"]

ab_fish = ab_fish.loc[(ab_fish['C/R'] >= 0.07) & (ab_fish['C/R'] <= 0.13)]
ab_fish = ab_fish.loc[ab_fish.groupby('size')['C/R'].idxmax()]

ab_crust = ab_crust.loc[(ab_crust['C/R'] >= 0.07) & (ab_crust['C/R'] <= 0.13)]
ab_crust = ab_crust.loc[ab_crust.groupby('size')['C/R'].idxmin()]

ab_ceph = ab_ceph.loc[(ab_ceph['C/R'] >= 0.07) & (ab_ceph['C/R'] <= 0.13)]
ab_ceph = ab_ceph.loc[ab_ceph.groupby('size')['C/R'].idxmin()]

run_sizeL = ab_fish["size"].unique()
run_size_crustL = ab_crust["size"].unique()

length = len(ab_fish['alpha'])
sizeL = ab_fish['size'].unique()

## Generating coef alpha coefficient for fish
alphaL = ab_fish.alpha.values
# Linear regression
lim_alpha_fish = alphaL.reshape(length, 1)
lim_size = sizeL.reshape(length, 1)
lm = linear_model.LinearRegression()
lm.fit(lim_size, lim_alpha_fish)

pred_fish_sizeL = np.array(run_sizeL).reshape(-1, 1)
coef_alpha_fish = lm.predict(pred_fish_sizeL)
a_fish = lm.intercept_
b_fish = lm.coef_
equation_fish = 'Fish: y='+str(np.round(a_fish[0], 2))+"x"+"+"+str(np.round(b_fish[0][0], 2))

## Generating coef alpha coefficient for crust
length = len(ab_crust['alpha'])
sizeL = ab_crust['size'].unique()
alphaL = ab_crust.alpha.values
pred_crust_sizeL = np.array(run_size_crustL).reshape(-1, 1)
# Linear regression
lim_alpha_crust = alphaL.reshape(length, 1)
lim_size = sizeL.reshape(length, 1)
lm = linear_model.LinearRegression()
lm.fit(lim_size, lim_alpha_crust)

coef_alpha_crust = lm.predict(pred_crust_sizeL)
a_crust = lm.intercept_
b_crust = lm.coef_
equation_crust = 'Crustacean: y='+str(np.round(a_crust[0], 2))+"x"+"+"+str(np.round(b_crust[0][0], 2))

## Generating coef alpha coefficient for ceph
pred_ceph_sizeL = np.array(run_sizeL).reshape(-1, 1)
length = len(ab_ceph['alpha'])
alphaL = ab_ceph.alpha.values
sizeL = ab_ceph['size'].unique()
# Linear regression
lim_alpha_ceph = alphaL.reshape(length, 1)
lim_size = sizeL.reshape(length, 1)
lm = linear_model.LinearRegression()
lm.fit(lim_size, lim_alpha_ceph)

coef_alpha_ceph = lm.predict(pred_ceph_sizeL)
a_ceph = lm.intercept_
b_ceph = lm.coef_
equation_ceph = 'Squids: y='+str(np.round(a_ceph[0], 2))+"x"+"+"+str(np.round(b_ceph[0][0], 2))

# Implementation of the initial conditions : carrying capacity
KL = np.zeros(zz)
KL[0:500] = gaussian(watercolumn[0:500], Ki, cen, sigma)

##Create empty dataframe for
ix_fishL = [0, 2, 6, 7, 9, 11]
num_fish = len(ix_fishL)
fish_pe_ratioL = np.zeros(num_fish)
crust_pe_ratioL = np.zeros(len(run_size_crustL[::4]))
ceph_pe_ratioL = np.zeros(num_fish)

fish_bio_ratioL = np.zeros(num_fish)
crust_bio_ratioL = np.zeros(len(run_size_crustL[::4]))
ceph_bio_ratioL = np.zeros(num_fish)

## Create empty dict for diff pools of detritus : FP, resp, mort
fish_detD = {}
crust_detD = {}
ceph_detD = {}
detL = ["FP", "Resp", "Mort"]
for det in detL:
    fish_detD[det] = np.zeros((num_fish, zz))
    crust_detD[det] = np.zeros((len(run_size_crustL[::4]), zz))
    ceph_detD[det] = np.zeros((num_fish, zz))

## Select only 6 sizes
ab_fish = ab_fish[ab_fish['size'].isin(list(run_sizeL[ix_fishL]))]
ab_crust = ab_crust[ab_crust['size'].isin(run_size_crustL[::4])]
ab_ceph = ab_ceph[ab_ceph['size'].isin(list(run_sizeL[ix_fishL]))]

crust_sizeL = pred_crust_sizeL.ravel().tolist()[::4]
fish_sizeL_ = pred_fish_sizeL.ravel().tolist()
fish_sizeL = list(np.array(fish_sizeL_)[ix_fishL])
ceph_sizeL_ = pred_ceph_sizeL.ravel().tolist()
ceph_sizeL = list(np.array(ceph_sizeL_)[ix_fishL])

for taxo in ["bm_fish", "bm_crust", "bm_ceph"]: #bm_idL
    # Respiration coefficient
    a0 = micron_dict[taxo][8]
    a1 = micron_dict[taxo][9]
    a2 = micron_dict[taxo][10]
    a3 = micron_dict[taxo][11]
    RQ = micron_dict[taxo][12]

    # Gut evacuation rate
    d = micron_dict[taxo][7]
    e = micron_dict[taxo][13]

    if taxo == "bm_crust":
        coef_alpha_pred = coef_alpha_crust[::4]
        pred_sizeL = crust_sizeL
        RiniL = ab_crust['R']
        CiniL = ab_crust['C']

    elif taxo == "bm_fish":
        coef_alpha_pred = list(np.array(coef_alpha_fish)[ix_fishL])
        pred_sizeL = fish_sizeL
        RiniL = ab_fish['R']
        CiniL = ab_fish['C']
    else:
        coef_alpha_pred = list(np.array(coef_alpha_ceph)[ix_fishL])
        pred_sizeL = ceph_sizeL
        RiniL = ab_ceph['R']
        CiniL = ab_ceph['C']

    for count, size in enumerate(pred_sizeL):
        print(count, size)
        # Swimming speed of micronekton
        Vmax = swimming_speed_func(size*0.1, taxo)  # in m/h

        dt_inf = CFL * dz / Vmax  # plus petit pas de temps possible avec la vitesse max
        time12 = np.arange(0, 12, dt_inf)
        Nt12 = len(time12)

        # The relative gradient of light is computed for the first 12h of the day, symmetric is taken and concatenated to it
        V = speedF(Nt12, dt_inf, Vmax)
        # Normalized the irradiance between 0 and 1 to set a visual predation parameter by micronekton
        I0 = np.zeros(Nt12)
        for t in range(Nt12):
            I0[t] = I_richards(t * dt_inf + dt_inf)
        I0 = np.hstack((I0, np.flip(I0)))

        # Create a 2d array (temps*depth) to represent the irradiance along depth
        I_depth = np.zeros((zz, Nt12 * 2))
        for z in range(zz):
            I_depth[z] = beer_lambert(I0, z)
        I_norm = normalize(I_depth)

        # Weight
        Wg, Wc, Wd = length_weight_func(size, taxo)  # Wg in mg and Wc in mgC

        if bm_id == "bm_fish":
            mu = natural_mortality_L22(size * 0.1)
        elif bm_id == "bm_crust":
            mu = mortality_allometricF(Wd * 1e-3)
        else:
            mu = 9e-5

        # Implementation of the initial conditions
        C0 = gaussian(watercolumn, 0.3, cen, sigma)
        R0 = gaussian(watercolumn, 3, cen, sigma)
        KL = gaussian(watercolumn, Ki, cen, sigma)
        NPP = np.sum(KL)

        coef_alpha = coef_alpha_pred[count]

        # EXECUTION
        results = AN_RGCD_VP(R0, C0, V, TkL, d, e, mu, w_ref, beta, dt_inf, KL, I_norm, coef_alpha, RQ, a0, a1, a2, a3, Wc_ind, Csda)

        id_200 = int(200 / dz)
        Ctot = np.sum(C0)

        mort = mu * results[2] * dt_save
        Det = results[3]
        Resp = results[4] + results[5] + results[6]

        # Calculation of pe_ratio below 200m
        Dg_int_200 = np.sum(results[3][:, id_200:zz])*7
        Mort_int_200 = np.sum(results[2][:, id_200:zz] * mu * dt_save)*7
        Resp_int_200 = np.sum(Resp[:, id_200:zz])

        if taxo == "bm_fish":
            FPsum = Det.sum(axis=0) * dt_save * 0.7 * 365
            Rsum = Resp.sum(axis=0) * dt_save * 0.7 * 365
            Mortsum = mort.sum(axis=0) * dt_save * 0.7 * 365
            axs[0, 0].plot(FPsum, -watercolumn, color=cmap(float(size)/kn.max()))
            axs[1, 0].plot(Rsum, -watercolumn, color=cmap(float(size) / kn.max()))
            axs[2, 0].plot(Mortsum, -watercolumn, color=cmap(float(size) / kn.max()))
            fish_detD["Resp"][count] = Rsum
            fish_detD["FP"][count] = FPsum
            fish_detD["Mort"][count] = Mortsum

            fish_pe_ratioL[count] = (Dg_int_200 + Mort_int_200) / NPP
            fish_bio_ratioL[count] = (Dg_int_200 + Mort_int_200 + Resp_int_200)/Ctot

        elif taxo == "bm_crust":
            FPsum = Det.sum(axis=0) * dt_save * 0.1 * 365
            Rsum = Resp.sum(axis=0) * dt_save * 0.1 * 365
            Mortsum = mort.sum(axis=0) * dt_save * 0.1 * 365
            axs[0, 1].plot(FPsum, -watercolumn, color=cmap(float(size)/kn_crust.max()))
            map = axs[1, 1].plot(Rsum, -watercolumn, color=cmap(float(size) / kn_crust.max()))
            axs[2, 1].plot(Mortsum, -watercolumn, color=cmap(float(size) / kn_crust.max()))
            crust_detD["Resp"][count] = Rsum
            crust_detD["FP"][count] = FPsum
            crust_detD["Mort"][count] = Mortsum

            crust_pe_ratioL[count] = (Dg_int_200 + Mort_int_200) / NPP
            crust_bio_ratioL[count] = (Dg_int_200 + Mort_int_200 + Resp_int_200) / Ctot

        else:
            FPsum = Det.sum(axis=0) * dt_save * 0.2 * 365
            Rsum = Resp.sum(axis=0) * dt_save * 0.2 * 365
            Mortsum = mort.sum(axis=0) * dt_save * 0.2 * 365
            axs[0, 2].plot(FPsum, -watercolumn, color=cmap(float(size)/kn.max()))
            axs[1, 2].plot(Rsum, -watercolumn, color=cmap(float(size) / kn.max()))
            axs[2, 2].plot(Mortsum, -watercolumn, color=cmap(float(size) / kn.max()))
            ceph_detD["Resp"][count] = Rsum
            ceph_detD["FP"][count] = FPsum
            ceph_detD["Mort"][count] = Mortsum

            ceph_pe_ratioL[count] = (Dg_int_200 + Mort_int_200) / NPP
            ceph_bio_ratioL[count] = (Dg_int_200 + Mort_int_200 + Resp_int_200) / Ctot

axs[0, 0].set_title("Fish")
axs[0, 1].set_title("Crustacean")
axs[0, 2].set_title("Cephalopod")
#fig.suptitle('Detritus production of different taxonomic group')
axs[0, 1].set_xlabel('Fecal pellets [$mgC$ $m^{-3}$ $y^{-1}$]')
axs[1, 1].set_xlabel('Respiration [$mgC$ $m^{-3}$ $y^{-1}$]')
axs[2, 1].set_xlabel('Dead bodies [$mgC$ $m^{-3}$ $y^{-1}$]')

axs[0, 0].set_ylabel('Depth [$m$]')
axs[1, 0].set_ylabel('Depth [$m$]')
axs[2, 0].set_ylabel('Depth [$m$]')

for ax in axs.flat:
    ax.set_ylim(-700, 0)
#xlim FP
axs[0, 0].set_xlim(0.0001, 0.95)
axs[0, 1].set_xlim(0.0001, 0.37)
axs[0, 2].set_xlim(0.0001, 0.28)
#xlim respiration
axs[1, 0].set_xlim(0.0001, 0.95)
axs[1, 1].set_xlim(0.0001, 0.37)
axs[1, 2].set_xlim(0.0001, 0.28)
#xlim DB
axs[2, 0].set_xlim(0.0001, 0.95)
axs[2, 1].set_xlim(0.0001, 0.37)
axs[2, 2].set_xlim(0.0001, 0.28)

norm = mpl.colors.Normalize(vmin=20, vmax=77)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

norm_crust = mpl.colors.Normalize(vmin=10, vmax=50)
sm_crust = plt.cm.ScalarMappable(cmap=cmap, norm=norm_crust)
sm_crust.set_array([])
# Colorbars
fig.colorbar(sm, ticks=np.linspace(20, 80, 6), ax=[axs[2, 0]], label="Size[mm]", orientation="horizontal")
fig.colorbar(sm_crust, ticks=np.linspace(10, 50, 6), ax=[axs[2, 1]], label="Size[mm]", orientation="horizontal")
fig.colorbar(sm, ticks=np.linspace(20, 80, 6), ax=[axs[2, 2]], label="Size[mm]", orientation="horizontal")

name = "Det_taxo_size.pdf"
file_name_det = path_to_plots / name
plt.savefig(file_name_det)
plt.close()

## Plot the detritus production as function of size, sum all det
fish_detA = fish_detD["Resp"] + fish_detD["FP"] + fish_detD["Mort"]
crust_detA = crust_detD["Resp"] + crust_detD["FP"] + crust_detD["Mort"]
ceph_detA = ceph_detD["Resp"] + ceph_detD["FP"] + ceph_detD["Mort"]

# Integrate along the watercolumn for each size : len(sum_det_fishL) = len(sizeL)
sum_det_fishL = np.sum(fish_detA, axis=1)
sum_det_crustL = np.sum(crust_detA, axis=1)
sum_det_cephL = np.sum(ceph_detA, axis=1)

for key, value in fish_detD.items():
    fish_detD[key] = np.sum(value, axis=1)
for key, value in crust_detD.items():
    crust_detD[key] = np.sum(value, axis=1)
for key, value in ceph_detD.items():
    ceph_detD[key] = np.sum(value, axis=1)

dict_det_fish = {
          'R': fish_detD["Resp"],
           'FP': fish_detD["FP"],
           'DB': fish_detD["Mort"]}
dict_det_crust = {
          'R': crust_detD["Resp"],
           'FP': crust_detD["FP"],
           'DB': crust_detD["Mort"]}
dict_det_ceph = {
          'R': ceph_detD["Resp"],
           'FP': ceph_detD["FP"],
           'DB': ceph_detD["Mort"]}

#Transforme dictionary in df
df_det_fish = pd.DataFrame(data=dict_det_fish)
df_det_crust = pd.DataFrame(data=dict_det_crust)
df_det_ceph = pd.DataFrame(data=dict_det_ceph)

row_sum_fish = df_det_fish.sum(axis=1)
row_sum_crust = df_det_crust.sum(axis=1)
row_sum_ceph = df_det_ceph.sum(axis=1)
# divide each column by row_sum and multiply by 100 to get percentage
df_perc_fish = df_det_fish.div(row_sum_fish, axis=0).multiply(100)
df_perc_crust = df_det_crust.div(row_sum_crust, axis=0).multiply(100)
df_perc_ceph = df_det_ceph.div(row_sum_ceph, axis=0).multiply(100)

print(df_perc_fish.describe())
print(df_perc_crust.describe())
print(df_perc_ceph.describe())

fig, axes = plt.subplots(1, 3, figsize=(10, 3), layout="constrained", sharex=True)
color_map = ["papayawhip", "rosybrown", "dimgrey"]
sns.boxplot(x="variable", y="value", data=pd.melt(df_perc_fish), ax=axes[0], palette=color_map)
sns.boxplot(x="variable", y="value", data=pd.melt(df_perc_crust), ax=axes[1], palette=color_map)
sns.boxplot(x="variable", y="value", data=pd.melt(df_perc_ceph), ax=axes[2], palette=color_map)
data = pd.melt(df_det_fish)

axes[0].set_xlabel('')
axes[0].set_ylabel('Relative carbon production')
axes[1].set_xlabel('Detritus pools')
axes[1].set_ylabel('')
axes[2].set_xlabel('')
axes[2].set_ylabel('')
axes[2].legend()

axes[0].set_title("Fish")
axes[1].set_title("Crustacean")
axes[2].set_title("Cephalopod")

name = "boxplot_det_taxo.pdf"
file_name_det = path_to_plots / name
plt.savefig(file_name_det)
plt.close()

# Create Plot
fig, axes = plt.subplots(1, 3, figsize=(10, 3), layout="constrained", sharey=True)
axes0 = axes[0].twinx()
axes1 = axes[1].twinx()
axes2 = axes[2].twinx()
#For fish
axes[0].plot(fish_sizeL, fish_detD["Resp"], color=color_map[0], label="Respiration", marker='o')
axes[0].plot(fish_sizeL, fish_detD["FP"], color=color_map[1], label="Fecal pellets", marker='o')
axes[0].plot(fish_sizeL, fish_detD["Mort"], color=color_map[2], label="Dead bodies", marker='o')
axes0.plot(fish_sizeL, fish_pe_ratioL, color="black", label="$pe-ratio_{200}$", linestyle='--', alpha=0.4)

axes[1].set_xlabel("Size [mm]")
axes[0].set_ylabel("$mgC$ $m^{-2}$ $y^{-1}$")

#For crustaceans
axes[1].plot(crust_sizeL, crust_detD["Resp"], color=color_map[0], label="Respiration", marker='o')
axes[1].plot(crust_sizeL, crust_detD["FP"], color=color_map[1], label="Fecal pellets", marker='o')
axes[1].plot(crust_sizeL, crust_detD["Mort"], color=color_map[2], label="Dead bodies", marker='o')
axes1.plot(crust_sizeL, crust_pe_ratioL, color="black", label="$pe_ratio_{200}$", linestyle='--', alpha=0.4)

#For cephalopods
axes[2].plot(ceph_sizeL, ceph_detD["Resp"], color=color_map[0], label="Respiration", marker='o')
axes[2].plot(ceph_sizeL, ceph_detD["FP"], color=color_map[1], label="Fecal pellets", marker='o')
axes[2].plot(ceph_sizeL, ceph_detD["Mort"], color=color_map[2], label="Dead bodies", marker='o')
axes2.plot(ceph_sizeL, ceph_pe_ratioL, color="black", label="$pe-ratio_{200}$", linestyle='--', alpha=0.4)
axes[2].legend()
axes0.legend()
axes2.set_ylabel("$pe-ratio_{200}$")
axes2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

axes0.set_ylim(0, 0.009)
axes1.set_ylim(0, 0.009)
axes2.set_ylim(0, 0.009)

axes[0].set_xlim(20, 77)
axes[1].set_xlim(10, 50)
axes[2].set_xlim(20, 77)

axes0.get_yaxis().set_visible(False)
axes1.get_yaxis().set_visible(False)

name = "Det_taxo_size_func.pdf"
file_name_det = path_to_plots / name
plt.savefig(file_name_det)
plt.close('all')

fig, axes = plt.subplots()
plt.plot(fish_sizeL, fish_bio_ratioL, color="black", label="$bio-ratio_{200}$", linestyle='--', alpha=0.4)
plt.plot(crust_sizeL, crust_bio_ratioL, color="black", label="$bio-ratio_{200}$", linestyle='--', alpha=0.4)
plt.plot(ceph_sizeL, ceph_bio_ratioL, color="black", label="$bio-ratio_{200}$", linestyle='--', alpha=0.4)
plt.show()




