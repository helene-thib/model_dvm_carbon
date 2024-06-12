"""
Plot the results of the seasonal model from dvm_model_seasonal_v4.py
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from pylab import meshgrid, imshow, colorbar
from parameters import *
from plotting_funcs import nonlinear_colormap, contour_levels_func, get_continuous_cmap


import cmocean
cm = nonlinear_colormap()

path_to_results_KTI = Path('~/simul_dvm2/Results/env_KTI').expanduser()
path_to_results_KTI2 = Path('~/simul_dvm2/Results/env_KTI2').expanduser()
path_to_results_KTI3 = Path('~/simul_dvm2/Results/env_KTI3').expanduser()
path_to_results_I = Path('~/simul_dvm2/Results/env_I').expanduser()
path_to_plot = Path('~/simul_dvm2/Plots').expanduser()
path_to_data = Path('~/simul_dvm2/Data').expanduser()

# Load results
path_R = path_to_results_KTI / 'results_R.txt'
path_G = path_to_results_KTI / 'results_G.txt'
path_C = path_to_results_KTI / 'results_C.txt'
path_RMR = path_to_results_KTI / 'results_M.txt'
path_Dg = path_to_results_KTI / 'results_Dg.txt'

path_resp_KTI2 = path_to_results_KTI2 / 'Maint_time_KTI2.txt'
path_Dg_KTI2 = path_to_results_KTI2 / 'Dg_time_KTI2.txt'
path_Du_KTI2 = path_to_results_KTI2 / 'Maint_time_KTI2.txt'

path_resp_KTI3 = path_to_results_KTI3 / 'Maint_time_KTI3.txt'
path_Dg_KTI3 = path_to_results_KTI3 / 'Dg_time_KTI3.txt'
path_Du_KTI3 = path_to_results_KTI3 / 'Maint_time_KTI3.txt'

path_AL_KTI = path_to_results_KTI / 'alphaL.txt'
path_AL_I = path_to_results_I / 'alphaL.txt'
path_PR = path_to_results_KTI / 'pe_ratio.txt'

R = np.loadtxt(path_R)
Req = R[-1, :]
C = np.loadtxt(path_C)
Ceq = C[-1, :]
G = np.loadtxt(path_G)

print("Req=", np.max(Req))
print("Ceq=", np.max(Ceq))
print("C/R=", np.max(Ceq)/np.max(Req))

if bm_id == "bm_fish":
    path_Ceq = path_to_results_KTI / 'Ceq_fish.txt'
    path_Req = path_to_results_KTI / 'Req_fish.txt'
elif bm_id == "bm_crust":
    path_Ceq = path_to_results_KTI / 'Ceq_crust.txt'
    path_Req = path_to_results_KTI / 'Req_crust.txt'
else:
    path_Ceq = path_to_results_KTI / 'Ceq_ceph.txt'
    path_Req = path_to_results_KTI / 'Req_ceph.txt'
np.savetxt(path_Req, Req)
np.savetxt(path_Ceq, Ceq)

# Integration of the state variables and detritus along depth
Ctot = C.sum(axis=1)
Rtot = R.sum(axis=1)
Gtot = G.sum(axis=1)

Ctot *= dz
Rtot *= dz
Gtot *= dz

Ctot_end = Ctot[-1]
Rtot_end = Rtot[-1]

# Last 24h values of R and C of the simulation
plt.ion()
fig = plt.figure()
for i in range(-int((24/dt_save)), -1, 2):
    plt.title("R and C at equilibrium")
    plt.plot(R[i, :], -watercolumn, label='Resource')
    plt.plot(C[i, :], -watercolumn, label='Consumers')
    plt.legend()

    # to flush the GUI events
    plt.pause(0.001)
    plot_name_RCeq = path_to_plot / "RCeq.pdf"
    plt.savefig(plot_name_RCeq)
    plt.clf()
plt.close('all')

# Plot of the integrated biomass over time of the 3 state variables
month = temps[:len(Ctot)]/24/30.5
fig, ax = plt.subplots()
ax.plot(month, Rtot, label='Resource')
ax.plot(month, Gtot, label='Gut')
ax.plot(month, Ctot, label='Consumer')
#ax.plot(month, zoo, label='Zooplankton')
ax.set_xlabel('Months')
ax.set_ylabel('Carbon mass ($mgC$ $m^{-2}$)')
textstr = '\n'.join((
    r'$\mathrm{Vmax (m/h)}=%.1f$' % (int(vmax), ),
    r'$\mathrm{µ (h^{-1})}=%.4f$' % (mu, ),
    r'$\mathrm{\alpha (h^{-1})}=%.3f$' % (coef_alpha, ),
    r'$\mathrm{\beta (h R^{-1})}=%.2f$' % (beta, ),
    r'$\mathrm{e}=%.2f$' % (e, )))
# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', edgecolor='none', facecolor='blue', alpha=0.05)
# place a text box in upper left in axes coords
ax.text(0.05, 0.55, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
plt.xlim(0, 11)
ax.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], ["Jan", "Fev", "Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec"], color="dimgrey", size=9)
plt.legend()
name_RGCtot = 'RGCtot_KTI.pdf'
plot_RGCtot = path_to_plot / name_RGCtot
plt.savefig(plot_RGCtot)
plt.close()

# Plot of the integrated biomass over time of the 3 state variables
fig, ax = plt.subplots(layout='constrained')
month = month[ : : 30]
ax.plot(month, Rtot[ : : 30], label='Resource')
ax.plot(month, Ctot[ : : 30], label='Consumer')

ax.set_xlabel('Months')
ax.set_ylabel('Carbon mass ($mgC$ $m^{-2}$)')
ax.set_xlim(0, 11)
ax2 = ax.twinx()
ax2.plot(month, Ctot[ : : 30]/Rtot[ : : 30], label='ratio C/R', color='black', alpha=0.4, linestyle='--')
ax.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], ["Jan", "Fev", "Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec"], size=10)
ax.tick_params(axis='x', labelrotation=45)
ax.legend(loc='upper left')
ax2.legend()
name_RGCtot = 'RC_ratio_KTI.pdf'
plot_RGCtot = path_to_plot / name_RGCtot
plt.savefig(plot_RGCtot)
plt.close()

alphaL_KTI = np.loadtxt(path_AL_KTI)
alphaL_I = np.loadtxt(path_AL_I)
## Plot alpha at 40m, 00h over year
fig, ax = plt.subplots(layout='constrained', sharex=True, figsize=(4.5, 3.5))
x = np.linspace(1, 12, len(alphaL_KTI))
x2 = np.linspace(1, 12, len(alphaL_I))
ax.plot(x2, alphaL_I, label="normal")
ax.plot(x, alphaL_KTI, label="with Chl $\it{a}$")
ax.set_xlabel('Months')
ax.set_ylabel('α')
ax.legend(loc="upper right")
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
name_AI = 'alpha_KTI.pdf'
plot_AI = path_to_plot / name_AI
plt.savefig(plot_AI)
plt.close()

### DETRITUS ###
Dg = np.loadtxt(path_Dg)
Dresp = np.loadtxt(path_RMR)
#Dead bodies of consumers and their gut
Du = (C+G) * mu*dt_save

# Integration along depth
Dg_tot = Dg.sum(axis=0)*dz
Maint_tot = Dresp.sum(axis=0)*dz
Mort_tot = Du.sum(axis=0)*dz
# Plot of the integrated biomass along depth of the detritus
det_dict = {'Respiration': Maint_tot.ravel().tolist(),
            'Fecal pellets': Dg_tot.ravel().tolist(),
            'Dead bodies': Mort_tot.ravel().tolist()}
detA = np.array((det_dict['Respiration'], det_dict['Fecal pellets'], det_dict['Dead bodies']))
sum_det = np.sum(detA, axis=0)
##Figure of detritus along depth
fig, ax = plt.subplots()
color_map = ["papayawhip", "rosybrown", "dimgrey"]
ax.fill_betweenx(-watercolumn, 0, det_dict['Respiration'], label='Respiration', color="papayawhip")
ax.fill_betweenx(-watercolumn, det_dict['Respiration'], np.add(det_dict['Respiration'], det_dict['Fecal pellets']), label='Fecal pellets', color="rosybrown")
ax.fill_betweenx(-watercolumn, np.add(det_dict['Respiration'], det_dict['Fecal pellets']), sum_det, label='Dead bodies', color="dimgrey")
#ax.stackplot([det_dict['Respiration'], det_dict['Fecal pellets'],det_dict['Dead bodies']], [watercolumn, watercolumn, watercolumn], labels=det_dict.keys(), alpha=0.8, colors=color_map)
plt.ylabel('Depth [$m$]')
plt.xlabel('$mgC$ $m^{-2}$')
ax.set_xlim(0)
plt.legend(loc="lower right")
fig.tight_layout()
name_det_tot = 'Det_tot_KTI.pdf'
plot_det_tot = path_to_plot / name_det_tot
plt.savefig(plot_det_tot)
plt.close()

# Nombre de jour
nj = tmax/24
# Integration over time
#FP
Dg_time = Dg.sum(axis=1)*dt_save
Dg_time = np.add.reduceat(Dg_time, np.arange(0, Dg_time.size, int(24/dt_save))) # Sum values every day to only see the annual variation
#Dead bodies
Du_time = Du.sum(axis=1)*dt_save
Du_time = np.add.reduceat(Du_time, np.arange(0, Du_time.size, int(24/dt_save)))
#RMR
Maint_time = Dresp.sum(axis=1)*dt_save
Maint_time = np.add.reduceat(Maint_time, np.arange(0, Maint_time.size, int(24/dt_save)))

year = np.linspace(0.5, 11.5, 334)
pe_ratio = np.loadtxt(path_PR)
year2 = np.linspace(0.5, 11.5, len(pe_ratio))
# Plot of the integrated biomass over time of the detritus

det_dict_time = {'Respiration': Maint_time.ravel().tolist(),
            'Fecal pellets': Dg_time.ravel().tolist(),
            #'Excretion': Exc_tot.ravel().tolist(),
            'Dead bodies': Du_time.ravel().tolist()}
detT_dict = pd.DataFrame({
                           'Respiration': Maint_time.ravel().tolist(),
                           'Fecal pellets': Dg_time.ravel().tolist(),
                           'Dead bodies': Du_time.ravel().tolist(), })
det_perc = detT_dict.divide(detT_dict.sum(axis=1), axis=0)*100

fig, axs = plt.subplots(2, 1, layout='constrained', sharex=True, figsize=(5, 6))
axs[0].stackplot(year, det_dict_time.values(),
             labels=det_dict_time.keys(), alpha=0.8, colors=color_map)
axs2 = axs[1].twinx()
axs2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
axs2.plot(year2, pe_ratio, color="black", label="$pe-ratio_{200}$", linestyle='--', alpha=0.7)
axs[0].set_xlim(0.5, 11.5)
#axs[0].set_ylim(0, 900)
axs[1].set_ylim(0, 100)
axs[1].set_xlabel('Months')
axs[0].set_ylabel('$mgC$ $m^{-2}$')
axs[1].set_ylabel('Relative carbon contribution [%]')
axs2.set_ylabel('$pe-ratio_{200}$')
axs[0].set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5], ["Jan", "Fev", "Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec"])
axs[0].tick_params(axis='x', labelrotation=45)
axs[1].tick_params(axis='x', labelrotation=45)
#plt.ylabel('Detritus concentrations [$mgC$ $m^{-2}$ $d^{-1}$]')
axs[1].stackplot(year,  det_perc['Respiration'],  det_perc['Fecal pellets'],  det_perc['Dead bodies'], colors=color_map)
axs[0].legend(loc='upper left', prop={'size': 9})
axs2.legend(loc='upper left', prop={'size': 9})
name_det_tot = 'Det_tot_timeKTI.pdf'
plot_det_tot = path_to_plot / name_det_tot
plt.savefig(plot_det_tot)
plt.close()

## All detritus
all_det_KTI = (Dg_time + Du_time + Maint_time)*1e-3*dt_save/365

Maint_time_KTI2 = np.loadtxt(path_resp_KTI2)
Dg_time_KTI2 = np.loadtxt(path_Dg_KTI2)
Du_time_KTI2 = np.loadtxt(path_Du_KTI2)
all_det_KTI2 = (Maint_time_KTI2 + Dg_time_KTI2 + Du_time_KTI2)*1e-3*dt_save/365

Maint_time_KTI3 = np.loadtxt(path_resp_KTI3)
Dg_time_KTI3 = np.loadtxt(path_Dg_KTI3)
Du_time_KTI3 = np.loadtxt(path_Du_KTI3)
all_det_KTI3 = (Maint_time_KTI3 + Dg_time_KTI3 + Du_time_KTI3)*1e-3*dt_save/365

##Sum of the detritus for present and 2100
sum_resp_KTI = np.sum(Maint_time)
sum_FP_KTI = np.sum(Dg_time)
sum_mort_KTI = np.sum(Du_time)

sum_resp_KTI2 = np.sum(Maint_time_KTI2)
sum_FP_KTI2 = np.sum(Dg_time_KTI2)
sum_mort_KTI2 = np.sum(Du_time)

sum_resp_KTI3 = np.sum(Maint_time_KTI3)
sum_FP_KTI3 = np.sum(Dg_time_KTI3)
sum_mort_KTI3 = np.sum(Du_time_KTI3)

legend = ['Respiration', 'Fecal pellets', 'Dead bodies']
x = ["Present", "2100 NPP", "2100 NPP+T"]
resp = [sum_resp_KTI, sum_resp_KTI2, sum_resp_KTI3]
FP = [sum_FP_KTI, sum_FP_KTI2, sum_FP_KTI3]
Mort = [sum_mort_KTI, sum_mort_KTI2, sum_mort_KTI3]

fig, axs = plt.subplots(1, 2, layout='constrained', figsize=(10, 4))
axs[0].plot(year, all_det_KTI, label="Present")
axs[0].plot(year, all_det_KTI2, label="2100 NPP")
axs[0].plot(year, all_det_KTI3, label="2100 NPP+T")
axs[0].set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5], ["Jan", "Fev", "Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec"], size=10)
axs[0].set_ylabel('Detritus concentrations [$gC$ $m^{-2}$]')
axs[0].set_xlabel('Months')
axs[0].set_xlim(0.5, 11.5)

axs[1].bar(x, resp, color=color_map[0])
axs[1].bar(x, FP, bottom=resp, color=color_map[1])
axs[1].bar(x, Mort, bottom=np.add(resp, FP), color=color_map[2])

axs[0].legend()
axs[1].legend(legend)
#fig.tight_layout()
name_det_tot = 'All_det_timeKTI.pdf'
plot_det_tot = path_to_plot / name_det_tot
plt.savefig(plot_det_tot)

### Seasonal production of detritus along depth
winter_idx = int((24*1)/dt_save)
spring_idx = int((24*90)/dt_save)
summer_idx = int((24*190)/dt_save)
fall_idx = int((24*250)/dt_save)
nb_days = 30
one_month = int(24*nb_days/dt_save)
## Select the day
Dg_winter = Dg[winter_idx:winter_idx+one_month].sum(axis=0)/nb_days
Dg_spring = Dg[spring_idx:spring_idx+one_month].sum(axis=0)/nb_days
Dg_summer = Dg[summer_idx:summer_idx+one_month].sum(axis=0)/nb_days
Dg_fall = Dg[fall_idx:fall_idx+one_month].sum(axis=0)/nb_days

Dr_winter = Dresp[winter_idx:winter_idx+one_month].sum(axis=0)/nb_days
Dr_spring = Dresp[spring_idx:spring_idx+one_month].sum(axis=0)/nb_days
Dr_summer = Dresp[summer_idx:summer_idx+one_month].sum(axis=0)/nb_days
Dr_fall = Dresp[fall_idx:fall_idx+one_month].sum(axis=0)/nb_days

Du_winter = Du[winter_idx:winter_idx+one_month].sum(axis=0)/nb_days
Du_spring = Du[spring_idx:spring_idx+one_month].sum(axis=0)/nb_days
Du_summer = Du[summer_idx:summer_idx+one_month].sum(axis=0)/nb_days
Du_fall = Du[fall_idx:fall_idx+one_month].sum(axis=0)/nb_days

#Sum det winter
det_winter_A = np.array((Dr_winter, Dg_winter, Du_winter))
sum_winter_det = np.sum(det_winter_A, axis=0)
#Sum det spring
det_spring_A = np.array((Dr_spring, Dg_spring, Du_spring))
sum_spring_det = np.sum(det_spring_A, axis=0)
#Sum det summer
det_summer_A = np.array((Dr_summer, Dg_summer, Du_summer))
sum_summer_det = np.sum(det_summer_A, axis=0)
#Sum det Fall
det_fall_A = np.array((Dr_fall, Dg_fall, Du_fall))
sum_fall_det = np.sum(det_fall_A, axis=0)

fig, axs = plt.subplots(1, 4, layout='constrained', sharey=True, sharex=True, figsize=(10, 4))
#Winter
axs[0].fill_betweenx(-watercolumn, 0, Dr_winter, label='Respiration', color="papayawhip")
axs[0].fill_betweenx(-watercolumn, Dr_winter, np.add(Dr_winter, Dg_winter), label='Fecal pellets', color="rosybrown")
axs[0].fill_betweenx(-watercolumn, np.add(Dr_winter, Dg_winter), sum_winter_det, label='Dead bodies', color="dimgrey")
#Spring
axs[1].fill_betweenx(-watercolumn, 0, Dr_spring, label='Respiration', color="papayawhip")
axs[1].fill_betweenx(-watercolumn, Dr_spring, np.add(Dr_spring, Dg_spring), label='Fecal pellets', color="rosybrown")
axs[1].fill_betweenx(-watercolumn, np.add(Dr_spring, Dg_spring), sum_spring_det, label='Dead bodies', color="dimgrey")
#Summer
axs[2].fill_betweenx(-watercolumn, 0, Dr_summer, label='Respiration', color="papayawhip")
axs[2].fill_betweenx(-watercolumn, Dr_summer, np.add(Dr_summer, Dg_summer), label='Fecal pellets', color="rosybrown")
axs[2].fill_betweenx(-watercolumn, np.add(Dr_summer, Dg_summer), sum_summer_det, label='Dead bodies', color="dimgrey")
#Fall
axs[3].fill_betweenx(-watercolumn, 0, Dr_fall, label='Respiration', color="papayawhip")
axs[3].fill_betweenx(-watercolumn, Dr_fall, np.add(Dr_fall, Dg_fall), label='Fecal pellets', color="rosybrown")
axs[3].fill_betweenx(-watercolumn, np.add(Dr_fall, Dg_fall), sum_fall_det, label='Dead bodies', color="dimgrey")

axs[0].set_ylabel('Depth [$m$]')
axs[0].set_ylim(-450, 0)
axs[0].set_xlim(0)
axs[1].set_xlabel('$mgC$ $m^{-3}$ $d^{-1}$')
ax.set_xlim(0)
axs[0].legend(loc="lower right")

axs[0].set_title("Winter")
axs[1].set_title("Spring")
axs[2].set_title("Summer")
axs[3].set_title("Autumn")

name_det_season = 'Det_seasons_KTI.pdf'
plot_det_season = path_to_plot / name_det_season
plt.savefig(plot_det_season)