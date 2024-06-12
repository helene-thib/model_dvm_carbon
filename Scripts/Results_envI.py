"""
Plot the results of the seasonal model from dvm_model_seasonal_v3.py
"""
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from pylab import meshgrid, imshow, colorbar
from parameters import *
from plotting_funcs import nonlinear_colormap, contour_levels_func, get_continuous_cmap


import cmocean
cm = nonlinear_colormap()

path_to_results_I = Path('~/simul_dvm2/Results/env_I').expanduser()
path_to_plot = Path('~/simul_dvm2/Plots').expanduser()
path_to_data = Path('~/simul_dvm2/Data').expanduser()

# Load results
#path_zoo = path_to_data / 'zoo_annual.txt'
path_R = path_to_results_I / 'results_R.txt'
path_G = path_to_results_I / 'results_G.txt'
path_C = path_to_results_I / 'results_C.txt'
path_RMR = path_to_results_I / 'results_M.txt'
path_Dg = path_to_results_I / 'results_Dg.txt'
path_AL = path_to_results_I / 'alphaL.txt'

# Read files at the given respective paths
output_names = ['results_R.txt', 'results_G.txt', 'results_C.txt', 'results_Dg.txt', 'results_M.txt']
results_K = []
results_T = []
results_I = []

#zoo = np.loadtxt(path_zoo)
R = np.loadtxt(path_R)
Req = R[-1, :]
C = np.loadtxt(path_C)
Ceq = C[-1, :]
G = np.loadtxt(path_G)
alphaL = np.loadtxt(path_AL)

## Plot alpha at 40m, 00h over year
fig, ax = plt.subplots(layout="constrained")
x = np.linspace(1, 12, len(alphaL))
ax.plot(x, alphaL)
ax.set_xlabel('Months')
ax.set_ylabel('α')
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
name_AI = 'alpha_I.pdf'
plot_AI = path_to_plot / name_AI
plt.savefig(plot_AI)

print("Req=", np.max(Req))
print("Ceq=", np.max(Ceq))
print("C/R=", np.max(Ceq)/np.max(Req))

if bm_id == "bm_fish":
    path_Ceq = path_to_results_I / 'Ceq_fishI.txt'
    path_Req = path_to_results_I / 'Req_fishI.txt'
elif bm_id == "bm_crust":
    path_Ceq = path_to_results_I / 'Ceq_crust.txt'
    path_Req = path_to_results_I / 'Req_crust.txt'
else:
    path_Ceq = path_to_results_I / 'Ceq_ceph.txt'
    path_Req = path_to_results_I / 'Req_ceph.txt'
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
ax.set_xlabel('Month ($h$)')
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
plt.xlim(0.5, 11.5)
ax.set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5], ["Jan", "Fev", "Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec"], color="dimgrey", size=9)
plt.legend()
name_RGCtot = 'RGCtot_I.pdf'
plot_RGCtot = path_to_plot / name_RGCtot
plt.savefig(plot_RGCtot)

fig, ax = plt.subplots(1, 1, figsize=(5, 3), layout="constrained")
# Pattern of consumers migration over a day
colors = ["#1f0954", "#1e488f", "#7bc8f6", "#6fc276", "#fcb001", "#fff917"]
cmap_dvm = get_continuous_cmap(colors)
dt = dt_save
props = dict(boxstyle='round', edgecolor='black', facecolor='white', alpha=0.9)
ims = []
tmax = len(Ctot)*dt
for i in range(int(tmax/24)):
    im = ax.imshow(np.transpose(C[int(i*24/dt):int((i*24/dt)+24/dt)]), aspect='auto', cmap=cmap_dvm, extent=[-12, 12, Zmax, 0], animated=True)
    t = ax.text(-10, 150, "Day " + str(i), bbox=props)
    plt.xticks([-12, -6, 0, 6, 12])
    plt.xlabel('Time relative to local noon [$h$]')
    plt.ylabel('Depth [$m$]')
    """plt.pause(0.01)
    plot_dvmC = path_to_plot / 'colormapDVMC.pdf'
    plt.savefig(plot_dvmC)
    plt.clf()"""
    if i == 0:
        #cbar = fig.colorbar(im, boundaries=np.linspace(0,0.05,6))
        cbar = fig.colorbar(im)
        cbar.set_label("$mgC$ $m^{-2}$")
    ims.append([im, t])

ani = animation.ArtistAnimation(fig, ims, interval=100, blit=True, repeat=False)
plt.show()
plot_dvmC = path_to_plot / 'colormapDVMC_I.gif'
ani.save(plot_dvmC, writer="pillow")
plt.close('all')

### DETRITUS ###
Dg = np.loadtxt(path_Dg)
Dresp = np.loadtxt(path_RMR)
#Dead bodies of consumers and their gut
Du = (C+G) * mu*dt_save

# Integration along depth
Dg_tot = Dg.sum(axis=0)[:2000]*dz
Maint_tot = Dresp.sum(axis=0)[:2000]*dz
Mort_tot = Du.sum(axis=0)[:2000]*dz

# Plot of the integrated biomass along depth of the detritus
det_dict = {'Respiration': Maint_tot.ravel().tolist(),
            'Fecal pellets': Dg_tot.ravel().tolist(),
            #'Excretion': Exc_tot.ravel().tolist(),
            'Dead bodies': Mort_tot.ravel().tolist()}
##Figure of detritus along depth
fig, ax = plt.subplots()
color_map = ["papayawhip", "rosybrown", "dimgrey"]
ax.stackplot(watercolumn[:2000], det_dict.values(),
             labels=det_dict.keys(), alpha=0.8, colors=color_map)
plt.ylabel('Detritus concentrations [$mgC$ $m^{-2}$ $h^{-1}$]')
plt.xlabel('Depth [$m$]')
plt.legend()
fig.tight_layout()
name_det_tot = 'Det_tot_I.pdf'
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

fig, axs = plt.subplots(2, 1, layout='constrained', sharex=True, figsize=(4.5, 6))
axs[0].stackplot(year, det_dict_time.values(),
             labels=det_dict_time.keys(), alpha=0.8, colors=color_map)

axs[0].set_xlim(0.5, 11.5)
axs[0].set_ylim(0, 140)
axs[1].set_ylim(0, 100)
axs[1].set_xlabel('Months')
axs[1].tick_params(axis='x', labelrotation=45)
axs[0].set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5], ["Jan", "Fev", "Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec"], color="dimgrey", size=9)
#plt.ylabel('Detritus concentrations [$mgC$ $m^{-2}$ $d^{-1}$]')
axs[1].stackplot(year,  det_perc['Respiration'],  det_perc['Fecal pellets'],  det_perc['Dead bodies'], colors=color_map)

name_det_tot = 'Det_tot_timeI.pdf'
plot_det_tot = path_to_plot / name_det_tot
plt.savefig(plot_det_tot)
plt.close()