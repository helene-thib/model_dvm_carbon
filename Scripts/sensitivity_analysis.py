"""
Sensivity analysis of metabolic parameters
Plot of Sobol's indices : heatmap and radar-plot
"""

import sys
sys.path.append("../..")

from SALib.analyze import sobol
import pandas as pd
from math import pi
import seaborn as sns

import matplotlib.pyplot as plt
from plotting_funcs import get_continuous_cmap
import numpy as np
from pathlib import Path
path_to_plot = Path('~/simul_dvm2/Plots').expanduser()
path_to_results = Path('~/simul_dvm2/Results/Sobol').expanduser()

from parameters import problem_crust, problem_fish, problem_ceph

# --- LOAD FILES
file_names_fish = ["sobol_fish_Ctot.txt", "sobol_fish_Rtot.txt", "sobol_fish_Cmax.txt", "sobol_fish_Rmax.txt", "sobol_param_fish.txt", "sobol_fish_POC.txt", "sobol_fish_Resp.txt", "sobol_fish_pe_ratio200.txt"]
file_names_crust = ["sobol_crust_Ctot.txt", "sobol_crust_Rtot.txt", "sobol_crust_Cmax.txt", "sobol_crust_Rmax.txt", "sobol_param_crust.txt", "sobol_crust_POC.txt", "sobol_crust_Resp.txt", "sobol_crust_pe_ratio200.txt"]
file_names_ceph = ["sobol_ceph_Ctot.txt", "sobol_ceph_Rtot.txt", "sobol_ceph_Cmax.txt", "sobol_ceph_Rmax.txt", "sobol_param_ceph.txt", "sobol_ceph_POC.txt", "sobol_ceph_Resp.txt",  "sobol_ceph_pe_ratio200.txt"]

outputs_fish = []
outputs_crust = []
outputs_ceph = []

# --- Crustaceans
for file_name in file_names_crust:
    file_path = path_to_results / file_name
    loaded_data = np.loadtxt(file_path)
    outputs_crust.append(loaded_data)

POC_tot_crust = outputs_crust[5]
Resp_tot_crust = outputs_crust[6]
pe_ratio_crust200 = outputs_crust[7]

## Generate sobol indices
Si_crust_bio = sobol.analyze(problem_crust, outputs_crust[2])
Si_crust_POC = sobol.analyze(problem_crust, POC_tot_crust)
Si_crust_Resp = sobol.analyze(problem_crust, Resp_tot_crust)
Si_crust_ratio200 = sobol.analyze(problem_crust, pe_ratio_crust200)

data_crust_bio = [Si_crust_bio['S1'], Si_crust_bio['ST']]
data_crust_ratio200 = [Si_crust_ratio200['S1'], Si_crust_ratio200['ST']]
data_crust_POC = [Si_crust_POC['S1'], Si_crust_POC['ST']]
print("biomass_crust ST: ", Si_crust_bio['ST'])

# --- Cephalopods
for file_name in file_names_ceph:
    file_path = path_to_results / file_name
    loaded_data = np.loadtxt(file_path)
    outputs_ceph.append(loaded_data)


outputs_ceph[5][np.isnan(outputs_ceph[5])] = 0
POC_tot_ceph = outputs_ceph[5]
outputs_ceph[6][np.isnan(outputs_ceph[6])] = 0
Resp_tot_ceph = outputs_ceph[6]
outputs_ceph[7][np.isnan(outputs_ceph[7])] = 0
pe_ratio_ceph200 = outputs_ceph[7]
#outputs_ceph[8][np.isnan(outputs_ceph[8])] = 0
#outputs_ceph[2][np.isnan(outputs_ceph[2])] = 0
bio_ceph = outputs_ceph[2]

## Generate sobol indices
Si_ceph_bio = sobol.analyze(problem_ceph, bio_ceph)
Si_ceph_POC = sobol.analyze(problem_ceph, POC_tot_ceph)
Si_ceph_Resp = sobol.analyze(problem_ceph, Resp_tot_ceph)
Si_ceph_ratio200 = sobol.analyze(problem_ceph, pe_ratio_ceph200)


data_ceph = [Si_ceph_bio['S1'],
              Si_ceph_bio['ST']]

# --- Fish
for file_name in file_names_fish:
    file_path = path_to_results / file_name
    loaded_data = np.loadtxt(file_path)
    outputs_fish.append(loaded_data)

POC_tot_fish = outputs_fish[5]
Resp_tot_fish = outputs_fish[6]
pe_ratio_fish200 = outputs_fish[7]

Si_fish_bio = sobol.analyze(problem_fish, outputs_fish[2])
Si_fish_POC = sobol.analyze(problem_fish, POC_tot_fish)
Si_fish_Resp = sobol.analyze(problem_fish, Resp_tot_fish)
Si_fish_ratio200 = sobol.analyze(problem_fish, pe_ratio_fish200)

data_fish_ratio200 = [Si_fish_ratio200['S1'], Si_fish_ratio200['ST']]

data_fish_bio = [Si_fish_bio['S1'],
              Si_fish_bio['ST']]

print("biomass_fish ST: ", Si_fish_bio['ST'])

"""Si_crust_ratio200.plot()
plt.show()"""

# -- First heatmap
viridis = sns.cubehelix_palette(start=.5, rot=-.5, as_cmap=True)
plt.figure(figsize=(10, 5))
ax1 = sns.heatmap(data_ceph, yticklabels=['S1', 'ST'], xticklabels=problem_ceph["names"], cmap='coolwarm', annot=True)
#ax2 = sns.heatmap(data_ceph, yticklabels=['S1', 'ST'], xticklabels =problem_ceph["names"], norm=LogNorm())
#plt.show()
plt.close()

# Set data : ST indice
df_biomass = pd.DataFrame({
    'group': ['Fish', 'Crust', 'Squid'],
    '$d$':        [Si_fish_bio['ST'][0], Si_crust_bio['ST'][0], Si_ceph_bio['ST'][0]],
    '$e$':        [Si_fish_bio['ST'][1], Si_crust_bio['ST'][1], Si_ceph_bio['ST'][1]],
    '$L_{inf}$':  [Si_fish_bio['ST'][2], 0, 0],
    '$k$':        [Si_fish_bio['ST'][3], 0, 0],
    '$\mu$':      [0, Si_crust_bio['ST'][2], Si_ceph_bio['ST'][2]],
    '$C_{sda}$':  [Si_fish_bio['ST'][4], Si_crust_bio['ST'][3], Si_ceph_bio['ST'][3]],
    '$A_f$':      [Si_fish_bio['ST'][5], Si_crust_bio['ST'][4], Si_ceph_bio['ST'][4]],
    '$a_{swim}$': [Si_fish_bio['ST'][6], Si_crust_bio['ST'][5], Si_ceph_bio['ST'][5]],
    '$RQ$':       [Si_fish_bio['ST'][7], Si_crust_bio['ST'][6], Si_ceph_bio['ST'][6]],
    '$a_0$':      [Si_fish_bio['ST'][8], Si_crust_bio['ST'][7], Si_ceph_bio['ST'][7]],
    '$a_1$':      [Si_fish_bio['ST'][9], Si_crust_bio['ST'][8], Si_ceph_bio['ST'][8]],
    '$a_2$':      [Si_fish_bio['ST'][10], Si_crust_bio['ST'][9], Si_ceph_bio['ST'][9]],
    '$a_3$':      [Si_fish_bio['ST'][11], Si_crust_bio['ST'][10], Si_ceph_bio['ST'][10]]
})



df_det = pd.DataFrame({
    'group': ['Fish', 'Crust', 'Squid'],
    '$d$':        [Si_fish_POC['ST'][0], Si_crust_POC['ST'][0], Si_ceph_POC['ST'][0]],
    '$e$':        [Si_fish_POC['ST'][1], Si_crust_POC['ST'][1], Si_ceph_POC['ST'][1]],
    '$L_{inf}$':  [Si_fish_POC['ST'][2], 0, 0],
    '$k$':        [Si_fish_POC['ST'][3], 0, 0],
    '$\mu$':      [0, Si_crust_POC['ST'][2], Si_ceph_POC['ST'][2]],
    '$C_{sda}$':  [Si_fish_POC['ST'][4], Si_crust_POC['ST'][3], Si_ceph_POC['ST'][3]],
    '$A_f$':      [Si_fish_POC['ST'][5], Si_crust_POC['ST'][4], Si_ceph_POC['ST'][4]],
    '$a_{swim}$': [Si_fish_POC['ST'][6], Si_crust_POC['ST'][5], Si_ceph_POC['ST'][5]],
    '$RQ$':       [Si_fish_POC['ST'][7], Si_crust_POC['ST'][6], Si_ceph_POC['ST'][6]],
    '$a_0$':      [Si_fish_POC['ST'][8], Si_crust_POC['ST'][7], Si_ceph_POC['ST'][7]],
    '$a_1$':      [Si_fish_POC['ST'][9], Si_crust_POC['ST'][8], Si_ceph_POC['ST'][8]],
    '$a_2$':      [Si_fish_POC['ST'][10], Si_crust_POC['ST'][9], Si_ceph_POC['ST'][9]],
    '$a_3$':      [Si_fish_POC['ST'][11], Si_crust_POC['ST'][10], Si_ceph_POC['ST'][10]]
})


df_ratio200 = pd.DataFrame({
    'group': ['Fish', 'Crust', 'Squid'],
    '$d$':        [Si_fish_ratio200['ST'][0], Si_crust_ratio200['ST'][0], Si_ceph_ratio200['ST'][0]],
    '$e$':        [Si_fish_ratio200['ST'][1], Si_crust_ratio200['ST'][1], Si_ceph_ratio200['ST'][1]],
    '$L_{inf}$':  [Si_fish_ratio200['ST'][2], 0, 0],
    '$k$':        [Si_fish_ratio200['ST'][3], 0, 0],
    '$\mu$':      [0, Si_crust_ratio200['ST'][2], Si_ceph_ratio200['ST'][2]],
    '$C_{sda}$':  [Si_fish_ratio200['ST'][4], Si_crust_ratio200['ST'][3], Si_ceph_ratio200['ST'][3]],
    '$A_f$':      [Si_fish_ratio200['ST'][5], Si_crust_ratio200['ST'][4], Si_ceph_ratio200['ST'][4]],
    '$a_{swim}$': [Si_fish_ratio200['ST'][6], Si_crust_ratio200['ST'][5], Si_ceph_ratio200['ST'][5]],
    '$RQ$':       [Si_fish_ratio200['ST'][7], Si_crust_ratio200['ST'][6], Si_ceph_ratio200['ST'][6]],
    '$a_0$':      [Si_fish_ratio200['ST'][8], Si_crust_ratio200['ST'][7], Si_ceph_ratio200['ST'][7]],
    '$a_1$':      [Si_fish_ratio200['ST'][9], Si_crust_ratio200['ST'][8], Si_ceph_ratio200['ST'][8]],
    '$a_2$':      [Si_fish_ratio200['ST'][10], Si_crust_ratio200['ST'][9], Si_ceph_ratio200['ST'][9]],
    '$a_3$':      [Si_fish_ratio200['ST'][11], Si_crust_ratio200['ST'][10], Si_ceph_ratio200['ST'][10]]
})


# --- HEATMAP --- #
#Fish biomass
sobol_bio_fishL = df_biomass.loc[0].drop('group').drop('$\mu$').tolist()
sobol_bio_fishL += sobol_bio_fishL[:0]
#Fish POC
sobol_det_fishL = df_det.loc[0].drop('group').drop('$\mu$').tolist()
sobol_det_fishL += sobol_det_fishL[:0]
#Fish pe-ratio 200
sobol_ratio200_fishL = df_ratio200.loc[0].drop('group').drop('$\mu$').tolist()
sobol_ratio200_fishL += sobol_ratio200_fishL[:0]

#Crustaceans biomass
sobol_bio_crustL = df_biomass.loc[1].drop('group').drop('$L_{inf}$').drop('$k$').tolist()
sobol_bio_crustL += sobol_bio_fishL[:0]
#Crustaceans POC
sobol_det_crustL = df_det.loc[1].drop('group').drop('$L_{inf}$').drop('$k$').tolist()
sobol_det_crustL += sobol_det_crustL[:0]
#Crustaceans pe-ratio 200
sobol_ratio200_crustL = df_ratio200.loc[1].drop('group').drop('$L_{inf}$').drop('$k$').tolist()
sobol_ratio200_crustL += sobol_ratio200_crustL[:0]

#Squids biomass
sobol_bio_cephL = df_biomass.loc[2].drop('group').drop('$L_{inf}$').drop('$k$').tolist()
sobol_bio_cephL += sobol_bio_cephL[:0]
#Squids POC
sobol_det_cephL = df_det.loc[2].drop('group').drop('$L_{inf}$').drop('$k$').tolist()
sobol_det_cephL += sobol_det_cephL[:0]
#Squids pe-ratio200
sobol_ratio200_cephL = df_ratio200.loc[2].drop('group').drop('$L_{inf}$').drop('$k$').tolist()
sobol_ratio200_cephL += sobol_ratio200_cephL[:0]

# --- HEATMAP
colors = ["#ffffff", "#a5c6e2", "#80afd6", "#5b97ca", "#3b7fb9", "#2f6694", "#234c6f"]
cmap = get_continuous_cmap(colors)

x_names_fish = ['$d$', '$e$', '$L_{inf}$', '$k$', '$C_{sda}$', '$A_f$', '$a_{swim}$', '$RQ$', '$a_0$', '$a_1$', '$a_2$', '$a_3$']
x_names_ceph = ['$d$', '$e$', '$\mu$', '$C_{sda}$', '$A_f$', '$a_{swim}$', '$RQ$', '$a_0$', '$a_1$', '$a_2$', '$a_3$']
x_names_crust = ['$d$', '$e$', '$\lambda$', '$C_{sda}$', '$A_f$', '$a_{swim}$', '$RQ$', '$a_0$', '$a_1$', '$a_2$', '$a_3$']

fig, (ax0, ax1, ax2) = plt.subplots(3, 1, layout='constrained')
heat_fish = [sobol_bio_fishL, sobol_det_fishL, sobol_ratio200_fishL]
heat_crust = [sobol_bio_crustL, sobol_det_crustL, sobol_ratio200_crustL]
heat_ceph = [sobol_bio_cephL, sobol_det_cephL, sobol_ratio200_cephL]
viridis = sns.cubehelix_palette(start=.5, rot=-.5, as_cmap=True)
x = [0.5, 1.5, 2.25, 2.75, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5]
y = [0.5, 1.5, 2.5]
im0 = ax0.pcolor(x, y, heat_fish, cmap=cmap, edgecolors='dimgrey', linewidths=1, vmin=0, vmax=0.65)
ax0.set_title('Fish')
ax0.set_xticks(x, x_names_fish, color="dimgrey", size=8)
ax0.set_yticks([0.5, 1.5, 2.5], ['Biomass', 'POC', 'pe-ratio$_{200}$'], color="dimgrey", size=8)

im1 = ax1.pcolor(heat_crust, cmap=cmap, edgecolors='dimgrey', linewidths=1, vmin=0, vmax=0.65)
ax1.set_title('Crustacean')
ax1.set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5], x_names_crust, color="dimgrey", size=8)
ax1.set_yticks([0.5, 1.5, 2.5], ['Biomass', 'POC', 'pe-ratio$_{200}$'], color="dimgrey", size=8)

im2 = ax2.pcolor(heat_ceph, cmap=cmap, edgecolors='dimgrey', linewidths=1, vmin=0, vmax=0.65)
ax2.set_title('Squid')
ax2.set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5], x_names_ceph, color="dimgrey", size=8)
ax2.set_yticks([0.5, 1.5, 2.5], ['Biomass', 'POC', 'pe-ratio$_{200}$'], color="dimgrey", size=8)

fig.colorbar(im1, ax=ax1)

# Add legend
plot_heat = path_to_plot / "sobol_heatmap.pdf"
plt.savefig(plot_heat)

# --- PART 1: Create background

# number of variable
categories = list(df_biomass)[1:]
N = len(categories)

# What will be the angle of each axis in the plot (we divide the plot / number of variable)
angles = [n / float(N) * 2 * pi for n in range(N)]
angles += angles[:1]

# Initialise the spider plot
fig = plt.figure(figsize=(11, 5), layout='constrained')
ax1 = fig.add_subplot(131, projection='polar')
ax2 = fig.add_subplot(132, projection='polar')
ax3 = fig.add_subplot(133, projection='polar')

# First axis on top:
ax1.set_theta_offset(pi / 2)
ax1.set_theta_direction(-1)
ax2.set_theta_offset(pi / 2)
ax2.set_theta_direction(-1)
ax3.set_theta_offset(pi / 2)
ax3.set_theta_direction(-1)

# Draw one axe per variable + add labels
ax1.set_xticks(angles[:-1], categories)
ax2.set_xticks(angles[:-1], categories)
ax3.set_xticks(angles[:-1], categories)

# Draw ylabels
ax1.set_rlabel_position(0)
ax2.set_rlabel_position(0)
ax3.set_rlabel_position(0)

ax1.set_yticks([0, 0.25, 0.5, 0.75, 1], ["0", "0.25", "0.5", "0.75", "1"], color="dimgrey", size=6)
ax2.set_yticks([0, 0.25, 0.5, 0.75, 1], ["0", "0.25", "0.5", "0.75", "1"], color="dimgrey", size=6)
ax3.set_yticks([0, 0.25, 0.5, 0.75, 1], ["0", "0.25", "0.5", "0.75", "1"], color="dimgrey", size=6)

ax1.set_ylim(0, 0.8)
ax2.set_ylim(0, 0.8)
ax3.set_ylim(0, 0.8)

ax1.set_title('C biomass')
ax2.set_title('Total POC')
ax3.set_title('$pe-ratio_{200m}$')

# ------- PART 2: Add plots
# Plot each individual = each line of the data

# Ind3
sobol_bio = df_biomass.loc[2].drop('group').values.flatten().tolist()
sobol_bio += sobol_bio[:1]
#ax1.plot(angles, sobol_bio, linewidth=1.5, linestyle='solid', label="Squids", color="gold")
#ax1.fill(angles, sobol_bio, 'r', alpha=0.1)
ax1.bar(angles, sobol_bio, width=0.6, alpha=0.2, label="Squid", color="gold", edgecolor="darkorange")
# Ind3 - Det
sobol_det = df_det.loc[2].drop('group').values.flatten().tolist()
sobol_det += sobol_det[:1]
#ax2.plot(angles, sobol_det, linewidth=1.5, linestyle='solid', label="Squids", color="gold")
#ax2.fill(angles, sobol_det, 'r', alpha=0.1)
ax2.bar(angles, sobol_det, width=0.6, alpha=0.2, label="Squid", color="gold", edgecolor="darkorange")
# Ind3 - pe-ratio200
sobol_ratio = df_ratio200.loc[2].drop('group').values.flatten().tolist()
sobol_ratio += sobol_ratio[:1]
#ax3.plot(angles, sobol_ratio, linewidth=1.5, linestyle='solid', label="Squids", color="gold")
#ax3.fill(angles, sobol_ratio, 'r', alpha=0.1)
ax3.bar(angles, sobol_ratio, width=0.6, alpha=0.2, label="Squid", color="gold", edgecolor="darkorange")

# Ind1 - Biomass
sobol_bio = df_biomass.loc[0].drop('group').values.flatten().tolist()
sobol_bio += sobol_bio[:1]
#ax1.plot(angles, sobol_bio, linewidth=1.5, linestyle='solid', label="Fish", color="lightblue")
ax1.bar(angles, sobol_bio, width=0.6, color='lightblue', alpha=0.7, label="Fish", edgecolor="steelblue")
#ax1.fill(angles, sobol_bio, 'b', alpha=0.1)
# Ind1 - Det
sobol_det = df_det.loc[0].drop('group').values.flatten().tolist()
sobol_det += sobol_det[:1]
#ax2.plot(angles, sobol_det, linewidth=1.5, linestyle='solid', label="Fish", color="lightblue")
#ax2.fill(angles, sobol_det, 'b', alpha=0.1)
ax2.bar(angles, sobol_det, width=0.6, color='lightblue', alpha=0.7, label="Fish", edgecolor="steelblue")
# Ind1 - pe-ratio200
sobol_ratio = df_ratio200.loc[0].drop('group').values.flatten().tolist()
sobol_ratio += sobol_ratio[:1]
#ax3.plot(angles, sobol_ratio, linewidth=1.5, linestyle='solid', label="Fish", color="lightblue")
#ax3.fill(angles, sobol_ratio, 'b', alpha=0.1)
ax3.bar(angles, sobol_ratio, width=0.6, color='lightblue', alpha=0.7, label="Fish", edgecolor="steelblue")

# Ind2 - Biomass
sobol_bio = df_biomass.loc[1].drop('group').values.flatten().tolist()
sobol_bio += sobol_bio[:1]
#ax1.plot(angles, sobol_bio, linewidth=1.5, linestyle='solid', label="Crustaceans", color="salmon")
#ax1.fill(angles, sobol_bio, 'r', alpha=0.1)
ax1.bar(angles, sobol_bio, width=0.6, alpha=0.4, label="Crustacean", color="salmon", edgecolor="crimson")
# Ind2 - Det
sobol_det = df_det.loc[1].drop('group').values.flatten().tolist()
sobol_det += sobol_det[:1]
#ax2.plot(angles, sobol_det, linewidth=1.5, linestyle='solid', label="Crustaceans", color="salmon")
#ax2.fill(angles, sobol_det, 'r', alpha=0.1)
ax2.bar(angles, sobol_det, width=0.6, alpha=0.4, label="Crustacean", color="salmon", edgecolor="crimson")
# Ind2 - pe-ratio200
sobol_ratio = df_ratio200.loc[1].drop('group').values.flatten().tolist()
sobol_ratio += sobol_ratio[:1]
#ax3.plot(angles, sobol_ratio, linewidth=1.5, linestyle='solid', label="Crustaceans", color="salmon")
#ax3.fill(angles, sobol_ratio, 'r', alpha=0.1)
ax3.bar(angles, sobol_ratio, width=0.6, alpha=0.4, label="Crustacean", color="salmon", edgecolor="crimson")

# Add legend
plt.suptitle("Total Sobol indices of the model outputs")
plt.legend(loc='upper left', bbox_to_anchor=(0.85, 1.3))
plot_Dg_d = path_to_plot / "sobol_radar.pdf"
plt.savefig(plot_Dg_d)